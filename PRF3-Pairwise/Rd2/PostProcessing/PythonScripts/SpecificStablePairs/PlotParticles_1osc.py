###MODULES###

import numpy as np
import pandas as pd
import os, sys
import time as t
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.ticker import MaxNLocator
import pathlib
from matplotlib.colors import Normalize
from scipy import interpolate
norm = Normalize()
from resource import getrusage, RUSAGE_SELF
import random
import scipy.ndimage as ndimage

mpl.rcParams['axes.linewidth'] = 1.5 #set the value globally
mpl.rcParams['contour.negative_linestyle'] = 'solid'

#CONSTANTS
cwd_PYTHON = os.getcwd() + '/'
RHO = 1000.0
NX = 512
PERIOD = 0.1
RADIUSLARGE = 0.002
RADIUSSMALL = 0.5*RADIUSLARGE
maxR = 0.025/RADIUSLARGE

csfont = {'fontname':'Times New Roman'}

#System Arguments
config = sys.argv[1]
Re = sys.argv[2]#"2"
perNumber = int(sys.argv[3])#5
local = int(sys.argv[4])

#PARTICLE STREAMLINE CONSTANTS
nFrames = 240
minVal, maxVal = -6.0,6.0
rows, cols = 100, 100
nPart, nTrail = rows*cols, 80
timestep = 0.4/60.0
dX = 2.0*maxR/(1.0*NX)
seed = random.seed(11235)

if local:
    cwd_FIGS = cwd_PYTHON+"../../Figures/VorticityDetection/{0}/".format(config)
    pathlib.Path(cwd_FIGS).mkdir(parents=True, exist_ok=True)
    cwd_Re = cwd_PYTHON+'../../FieldData/TestField/'
    cwd_POS = cwd_Re
else:
    cwd_FIGS = cwd_PYTHON+'../Figures/Bifurcation/{0}/Particles/'.format(config)
    pathlib.Path(cwd_FIGS).mkdir(parents=True, exist_ok=True)
    cwd_Re = cwd_PYTHON+'../{0}/Re{1}/VTK/AVG/'.format(config,Re)
    cwd_POS = cwd_PYTHON+'../PosData/{0}/Re{1}/'.format(config,Re)

# constructs a filepath for the pos data of Re = $Re
def pname(cwd,idx):
    return cwd+"pd_rot_%04d.csv"%idx

def GetPosData(cwd,idx):
    data = pd.read_csv(pname(cwd,idx),delimiter=' ')
    return data

def GetAvgFieldData(cwd,idx):
    #Load position data
    #Columns
    #mx.flat my.flat avgW.flat avgP.flat avgUx.flat avgUy.flat
    fieldData = pd.read_csv(cwd+'AVGRot_%04d.csv'%idx,delimiter=' ')
    print(fieldData.head())
    #All field values to a list
    mxList = fieldData['mx'].values.tolist()
    myList = fieldData['my'].values.tolist()
    WList  = fieldData['avgW'].values.tolist()
    UxList = fieldData['avgUx'].values.tolist()
    UyList = fieldData['avgUy'].values.tolist()
    #Convert lists to numpy arrays
    #Reshape them to be Nx x Ny
    Nx, Ny = 512, 512
    mxArr = np.array(mxList).reshape((Nx,Ny))
    myArr = np.array(myList).reshape((Nx,Ny))
    WArr  = np.array(WList).reshape((Nx,Ny))
    UxArr = np.array(UxList).reshape((Nx,Ny))
    UyArr = np.array(UyList).reshape((Nx,Ny))
    return (mxArr, myArr, WArr, UxArr, UyArr)

def AddDiscsToPlot(ax,pos):
    #Add Discs
    circle1 = Circle((pos.loc[0,'aXU_rot'], pos.loc[0,'aYU_rot']), 1.0,         facecolor=(0.25,)*3,
                     linewidth=1,alpha=1.0,zorder=6)
    ax.add_patch(circle1)
    circle2 = Circle((pos.loc[0,'aXL_rot'], pos.loc[0,'aYL_rot']), 0.5, facecolor=(0.25,)*3,
                     linewidth=1,alpha=1.0,zorder=6)
    ax.add_patch(circle2)
    circle3 = Circle((pos.loc[0,'bXU_rot'], pos.loc[0,'bYU_rot']), 1.0, facecolor=(0.75,)*3,
                     linewidth=1,alpha=1.0,zorder=6)
    ax.add_patch(circle3)
    circle4 = Circle((pos.loc[0,'bXL_rot'], pos.loc[0,'bYL_rot']), 0.5, facecolor=(0.75,)*3,
                     linewidth=1,alpha=1.0,zorder=6)
    ax.add_patch(circle4)
    #Add Swimmer "springs"
    ax.plot([pos.loc[0,'aXU_rot'],pos.loc[0,'aXL_rot']],
            [pos.loc[0,'aYU_rot'],pos.loc[0,'aYL_rot']],
            color=(0.25,)*3,linewidth=3,zorder=6)
    ax.plot([pos.loc[0,'bXU_rot'],pos.loc[0,'bXL_rot']],
            [pos.loc[0,'bYU_rot'],pos.loc[0,'bYL_rot']],
            color=(0.75,)*3,linewidth=3,zorder=6)
    return

def set_size(w,h, ax=None):
    """ w, h: width, height in inches """
    if not ax: ax=plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w)/(r-l)
    figh = float(h)/(t-b)
    ax.figure.set_size_inches(figw, figh)
    return ax

def CalcPsi2D(fx,fy,NX,DX):
    #From here, we are going to calculate the stream function
    psi = np.zeros((NX,NX))
    for idx in range(1,NX):
        psi[idx,0] = psi[idx-1,0] - fy[idx,0]*DX
    for idy in range(1,NX):
        psi[:,idy] = psi[:,idy-1] + fx[:,idy]*DX
    return psi

#Plot New mesh and interpolated velocity field Ux and Uy
def PlotAvgW(cwd,mx,my,W,Ux,Uy,pos,space,scale):
    global FIGNUM, PERIOD,minVal,maxVal, Re, perNumber
    #Here, we will visualize the velocity field on the new coordinate system
    nRows, nCols = 1, 1
    fig, ax = plt.subplots(nrows=nRows, ncols=nCols, num=0,figsize=(6,6),dpi=200)
    #ax.set_title(r'Average Velocity Field',fontsize=12)
    #Plot Streamlines
    #Use two grids and combine them
    UxT, UyT, WT = Ux.T, Uy.T, W.T
    psi = CalcPsi2D(Ux,Uy,NX,dX)
    print('psi.min() = ',psi.min())
    print('psi.max() = ',psi.max())
    sys.stdout.flush()
    #Psi Contour
    psi2 = ndimage.gaussian_filter(psi, sigma=5.0, order=0)
    levels = MaxNLocator(nbins=61).tick_values(-1.0*max(abs(psi2.min()),psi2.max()), max(abs(psi2.min()),psi2.max()))
    #levels = MaxNLocator(nbins=21).tick_values(-1.0*max(abs(psi2.min()),psi2.max()), max(abs(psi2.min()),psi2.max()))
    ax.contour(mx,my,psi2,colors='k',extend='both',levels=levels)
    #PlotVorticity with imshow (interpolate to smooth)
    ax.imshow(W.T,cmap='bwr',extent=(-1.0*maxR-0.5*dX,maxR+0.5*dX,
                                   -1.0*maxR-0.5*dX,maxR+0.5*dX),
            origin='lower',vmin=-1.0,vmax=1.0,interpolation='bilinear')
    #Add swimmer
    AddDiscsToPlot(ax,pos)
    xmin = min(pos.loc[0,'aXU_rot'],pos.loc[0,'aXL_rot'],
               pos.loc[0,'bXU_rot'],pos.loc[0,'bXL_rot'])
    xmax = max(pos.loc[0,'aXU_rot'],pos.loc[0,'aXL_rot'],
               pos.loc[0,'bXU_rot'],pos.loc[0,'bXL_rot'])
    ymin = min(pos.loc[0,'aYU_rot'],pos.loc[0,'aYL_rot'],
               pos.loc[0,'bYU_rot'],pos.loc[0,'bYL_rot'])
    ymax = max(pos.loc[0,'aYU_rot'],pos.loc[0,'aYL_rot'],
               pos.loc[0,'bYU_rot'],pos.loc[0,'bYL_rot'])
    #ax.axis([xmin-0.25,xmax+0.25,ymin-2.0,ymax+2.0])
    ax.axis([-5.0,5.0,-5.0,5.0])
    fig.tight_layout()
    fig.savefig(cwd+'W_{0}_Re{1}_per{2}_.png'.format(config,Re,perNumber))
    fig.clf()
    plt.close()
    return

#Plot New mesh and interpolated velocity field Ux and Uy
def PlotAvgU(cwd,mx,my,Ux,Uy,pos,space,config,Re,per):
    global FIGNUM, PERIOD,minVal,maxVal
    #Here, we will visualize the velocity field on the new coordinate system
    nRows, nCols = 1, 1
    fig, ax = plt.subplots(nrows=nRows, ncols=nCols, num=0,figsize=(6,6),dpi=200)
    ax.set_title(r'Average Velocity Field',fontsize=12)
    normUx,normUy = Ux/np.hypot(Ux,Uy),Uy/np.hypot(Ux,Uy)
    magU = np.hypot(Ux,Uy)
    #Plot Vector field with quiver
    ax.quiver(mx[::space,::space],my[::space,::space],
              Ux[::space,::space],Uy[::space,::space],
              color='white',pivot='mid',angles='xy',scale_units='xy', scale=10,zorder=5)
    #Plot magnitude with contourf
    ax.contourf(mx,my,magU,cmap='viridis')

    AddDiscsToPlot(ax,pos)
    #print('RSMALL = ',RSMALL)
    ax.axis([minVal,maxVal,minVal,maxVal])
    fig.tight_layout()
    #plt.show()
    fig.savefig(cwd+'avgU_{0}_Re{1}_per{2}_.png'.format(config,Re,per))
    fig.clf()
    plt.close()
    return

#Plot New mesh and interpolated velocity field Ux and Uy
def PlotParticles(mx,my,pos,particles,frame,cwd,config,Re,per):
    global FIGNUM, PERIOD,minVal,maxVal
    #Here, we will visualize the velocity field on the new coordinate system
    nRows, nCols = 1, 1
    fig, ax = plt.subplots(nrows=nRows, ncols=nCols, num=1,figsize=(6,6),dpi=200)
    
    alpha = np.linspace(0.2,1.0,particles.nTime)
    for idTime in range(particles.nTime):
        pointColor = (1.0-(1.0*idTime/(1.0*particles.nTime)),)*3
        alphaValue = alpha[idTime]
        markerSize = particles.size[idTime,0]
        ax.plot(particles.x[idTime].flatten(),particles.y[idTime].flatten(),
                marker='o',ms=markerSize,color=pointColor,zorder=5,alpha=1,linewidth=0)
    
    #Plot Swimmer
    AddDiscsToPlot(ax,pos)
    #print('RSMALL = ',RSMALL)
    #ax.axis([minVal,maxVal,minVal,maxVal])
    ax.axis([-5,5,-5,5])
    fig.tight_layout()
    fig.savefig(cwd+'Part_{0}_Re{1}_per{2}_.png'.format(config,Re,per))
    #plt.show()
    fig.clf()
    plt.close()
    return

class Particles:
    def __init__(self,rows,cols,minVal,maxVal,nPart,nTime):
        self.nPart = nPart
        print('nTime = ',nTime)
        self.nTime = nTime
        print('self.nTime = ',self.nTime)
        xvals = np.linspace(minVal,maxVal,rows)
        yvals = np.linspace(minVal,maxVal,cols)
        mx, my = np.meshgrid(xvals,yvals)
        self.x  = np.array([mx.flatten()]*nTime).reshape((nTime,nPart))
        self.y  = np.array([my.flatten()]*nTime).reshape((nTime,nPart))
        self.xinit = self.x.copy()
        self.yinit = self.y.copy()
        self.vx  = np.zeros((nTime,nPart))
        self.vy  = np.zeros((nTime,nPart))
        self.vxinit = self.vx.copy()
        self.vyinit = self.vy.copy()
        self.idx = np.array([[int((self.x[a,b] + maxR)/dX) for b in range(self.nPart)] for a in range(self.nTime)])
        self.idy = np.array([[int((self.y[a,b] + maxR)/dX) for b in range(self.nPart)] for a in range(self.nTime)])
        self.age  = np.array([random.randrange(10,240) for b in range(nPart)]*nTime).reshape((nTime,nPart))
        self.life = self.age.copy()
        self.size = np.array([np.linspace(0.001,.1,self.nTime) for a in range(self.nPart)]).T.reshape((nTime,nPart))
        self.curr_age  = np.zeros((nTime,nPart))
    def CalcMeshIndex(self,idTime):
        self.idx[idTime] = [int((self.x[idTime,b] + maxR)/dX) for b in range(self.nPart)]
        self.idy[idTime] = [int((self.y[idTime,b] + maxR)/dX) for b in range(self.nPart)]
    def AssignVelocity(self,idTime,nX,avgUx,avgUy):
        indices = nX*self.idx[idTime]+self.idy[idTime]
        self.vx[idTime] = avgUx[indices]
        self.vy[idTime] = avgUy[indices]

if __name__ == '__main__':
    #Get AvgVel Field and Rotate Frame
    #Save Vel Field as AvgUx and AvgUy
    #READ ALL AVG FILES IN A SIMULATION DIRECTORY
    #EXTRACT AVERAGE FIELD DATA INTO NUMPY ARRAYS
    #PLOT AVERAGED FIELD DATA
    #Simulation Parameters
    #Extract Position Data
    #Paths to data and plots
    cwd_DATA = cwd_Re
    countPer = perNumber
    
    AVGPlot = pathlib.Path(cwd_DATA+'AVGRot_%04d.csv'%countPer)
    if AVGPlot.exists ():
        start = t.clock()
        #Get Avg Field Data
        mx,my,avgW,avgUx,avgUy = GetAvgFieldData(cwd_DATA,countPer)
        #Extract Position and Time Data
        posData = GetPosData(cwd_POS,countPer)
        #Plot Averaged Field Data
        #Vorticity And Streamlines
        stend = t.clock()
        diff = stend - start
        print('Time to run for 1 period = %.5fs'%diff)
        sys.stdout.flush()

    #Plot Flow Field Visual
    PlotAvgU(cwd_FIGS,mx,my,avgUx,avgUy,posData,6,config,Re,perNumber)
    #PlotAvgW(cwd_FIGS,mx,my,avgW,avgUx,avgUy,posData,4,5)

    #Now that we have the abg velocity field, we can calculate particle trajectories
    #Let's start with 10
    # Initialize 10 random coordinates
    # Each coordinate will have a lifetime from 60 -> 240 frames
    # For each, calculate the new position based on the velocity field (do this for twenty timesteps)
    # The velocity field does not evolve
    # Create array of points (use interpolate.griddata to find velocities)
    # pos += dt*velocity_interp
    # Plot scatter which decreases in opacity and point size for timesteps going backward
    # Advance time for each new frame

    #Initialize Uniform Distribution of points. The structure should be a 2D ndarray
    #Choose points in range (-7.5,7.5) for both x and y

    #Flatten avg Velocity field
    avgU = np.array([avgUx,avgUy])
    print('shape = ',avgU.shape)
    print('avgU = ',avgU[1,281,248])
    avgUx = avgUx.flatten()
    avgUy = avgUy.flatten()
    print('avgUy = ',avgUy[512*248 + 281])
    print('avgUy = ',avgUy[512*281 + 248])
    magU = np.hypot(avgUx,avgUy)
    while np.amax(magU)*timestep > 0.75*dX:
        timestep *= 0.95
    print('timestep = ',timestep)
    #print('dX = ',dX)
    #print('max dX = ',timestep*np.amax(magU))
    assert 0.75*dX >= np.amax(magU)*timestep
    particles = Particles(rows,cols,minVal,maxVal,nPart,nTrail)

    #Initialize Particles
    #Find velocity by index value (no interpolation)
    for idTime in range(nTrail):
        #Calculate Mesh Index for idTime
        particles.CalcMeshIndex(idTime)
        #Assign velocity by index (no interp)
        particles.AssignVelocity(idTime,NX,avgUx,avgUy)
        #Update position idTime+1
        if idTime < nTrail - 1:
            changeX = timestep*particles.vx[idTime]
            changeY = timestep*particles.vy[idTime]
            particles.x[idTime+1:nTrail] += changeX
            particles.y[idTime+1:nTrail] += changeY
            #Increase age of particles
            particles.curr_age[:idTime+1] -= 1
    #Save Initial particle stream pos and vel
    particles.xinit = particles.x[0,:].copy()#particles.x.copy()
    particles.yinit = particles.y[0,:].copy()#particles.y.copy()
    particles.vxinit = particles.vx[0,:].copy()#particles.vx.copy()
    particles.vyinit = particles.vy[0,:].copy()#particles.vy.copy()

    #Loop over # of frames
    for idxFrame in range(nFrames+1):
        
        #Check Age
        particles.x = np.where(particles.curr_age >= particles.age, particles.xinit, particles.x)
        particles.y = np.where(particles.curr_age >= particles.age, particles.yinit, particles.y)
        particles.vx = np.where(particles.curr_age >= particles.age, particles.vxinit, particles.vx)
        particles.vy = np.where(particles.curr_age >= particles.age, particles.vyinit, particles.vy)
        particles.age = np.where(particles.curr_age >= particles.age, nFrames, particles.age)
        particles.curr_age = np.where(particles.curr_age >= particles.life, 0, particles.curr_age)
        particles.life = particles.age.copy()

        #Plot Particle Stream
        if idxFrame >= nFrames:
            PlotParticles(mx,my,posData,particles,idxFrame-nFrames,cwd_FIGS,config,Re,perNumber)
        #Increase particle age
        particles.curr_age += 1
        #Roll positions and velocities back 1. Will need to calculate last position each particle
        particles.x = np.roll(particles.x,-1,axis=0)
        particles.y = np.roll(particles.y,-1,axis=0)
        particles.vx = np.roll(particles.vx,-1,axis=0)
        particles.vy = np.roll(particles.vy,-1,axis=0)
        #Change Position of last value by using velocity of 2nd to last
        changeX = timestep*particles.vx[-2]
        #print('vy[-2,1] = ',particles.vy[-2,1])
        changeY = timestep*particles.vy[-2]
        #print('changeY = ',changeY[1])
        particles.x[-1] = particles.x[-2] + changeX
        particles.y[-1] = particles.y[-2] + changeY
        #Update Mesh Index
        particles.CalcMeshIndex(-1)
        #Update Velocity of first particle in trail
        particles.AssignVelocity(-1,NX,avgUx,avgUy)

        if idxFrame %60 == 0:
            print('Frame {0} is complete'.format(idxFrame))
            sys.stdout.flush()

