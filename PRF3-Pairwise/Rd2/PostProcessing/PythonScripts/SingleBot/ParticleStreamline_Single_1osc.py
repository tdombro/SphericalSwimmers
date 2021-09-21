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

mpl.rcParams['axes.linewidth'] = 1.5 #set the value globally

#CONSTANTS
cwd_PYTHON = os.getcwd() + '/'
RHO = 1000.0
DX = 0.025/256.0
NX = 512
PERIOD = 0.1
RADIUSLARGE = 0.002
RADIUSSMALL = 0.5*RADIUSLARGE
maxR = 0.025/RADIUSLARGE

csfont = {'fontname':'Times New Roman'}

#System Arguments
Re = sys.argv[1]#"2"
perNumber = int(sys.argv[2])#5

#PARTICLE STREAMLINE CONSTANTS
nFrames = 240
minVal, maxVal = -6.0,6.0
rows, cols = 100, 100
nPart, nTrail = rows*cols, 80
timestep = 0.2/60.0
dX = 2.0*maxR/(1.0*NX)
seed = random.seed(11235)
cwd_FIGS = cwd_PYTHON+'../Figures/ParticleStreamline/Re{0}/per{1}/'.format(Re,perNumber)
pathlib.Path(cwd_FIGS).mkdir(parents=True, exist_ok=True)
cwd_Re = cwd_PYTHON+'../Fig4_VisitFiles_Single/Re{0}/VTK/AVG/'.format(Re)
cwd_POS = cwd_PYTHON+'../PosData/'

# constructs a filepath for the pos data of Re = $Re
def pname(cwd,Re):
    #return cwd+"/pd.txt"
    #cwd = cwd_PYTHON
    return cwd+"/pd_Re{0}.txt".format(Re)

def GetPosData(cwd,time,Re):
    global RADIUSLARGE
    data = pd.read_csv(pname(cwd,Re),delimiter=' ')
    topData = data[data['idx'] == 6].copy()
    botData = data[data['idx'] == 19].copy()
    topData = topData.sort_values(by=['time'])
    botData = botData.sort_values(by=['time'])
    topData = topData.reset_index(drop=True)
    botData = botData.reset_index(drop=True)
    dictPos = {'xU':topData['x'],'yU':topData['y'],'xL':botData['x'],'yL':botData['y'],'time':topData['time']}
    pos = pd.DataFrame(data=dictPos)
    pos = pos[pos['time'] == time]
    pos = pos.reset_index(drop=True)
    pos['xU'] /= RADIUSLARGE
    pos['xL'] /= RADIUSLARGE
    pos['yU'] /= RADIUSLARGE
    pos['yL'] /= RADIUSLARGE
    return pos

def GetPosDataLength(cwd,Re):
    data = pd.read_csv(pname(cwd,Re),delimiter=' ')
    topData = data[data['idx'] == 6].copy()
    botData = data[data['idx'] == 19].copy()
    topData = topData.sort_values(by=['time'])
    botData = botData.sort_values(by=['time'])
    topData = topData.reset_index(drop=True)
    botData = botData.reset_index(drop=True)
    assert len(topData['time']) == len(botData['time'])
    return len(topData['time'])

def GetAvgFieldData(cwd,idx):
    global RADIUSLARGE
    #Load position data
    #Columns
    #mx.flat my.flat avgW.flat avgP.flat avgUx.flat avgUy.flat
    fieldData = pd.read_csv(cwd+'AVG_%04d.csv'%idx,delimiter=' ')
    print(fieldData.head())
    #All field values to a list
    mxList = fieldData['mx'].values.tolist()
    myList = fieldData['my'].values.tolist()
    WList  = fieldData['avgW'].values.tolist()
    PList  = fieldData['avgP'].values.tolist()
    UxList = fieldData['avgUx'].values.tolist()
    UyList = fieldData['avgUy'].values.tolist()
    #Convert lists to numpy arrays
    #Reshape them to be Nx x Ny
    Nx, Ny = 1024, 1024
    mxArr = np.array(mxList).reshape((Nx,Ny))/RADIUSLARGE
    myArr = np.array(myList).reshape((Nx,Ny))/RADIUSLARGE
    WArr  = np.array(WList).reshape((Nx,Ny))
    PArr  = np.array(PList).reshape((Nx,Ny))
    UxArr = np.array(UxList).reshape((Nx,Ny))/RADIUSLARGE
    UyArr = np.array(UyList).reshape((Nx,Ny))/RADIUSLARGE
    return (mxArr, myArr, WArr, PArr, UxArr, UyArr)

def AddDiscsToPlot(ax,pos):
    #Add Discs
    circle1 = Circle((pos.loc[0,'aXU_rot'], pos.loc[0,'aYU_rot']), 1.0, facecolor=(0.0,)*3,
                     linewidth=1,alpha=1.0,zorder=6)
    ax.add_patch(circle1)
    circle2 = Circle((pos.loc[0,'aXL_rot'], pos.loc[0,'aYL_rot']), 0.5, facecolor=(0.0,)*3,
                  linewidth=1,alpha=1.0,zorder=6)
    ax.add_patch(circle2)
    #Add Swimmer "springs"
    ax.plot([pos.loc[0,'aXU_rot'],pos.loc[0,'aXL_rot']],
         [pos.loc[0,'aYU_rot'],pos.loc[0,'aYL_rot']],
         color='black',linewidth=3,zorder=6)
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

def Rotate(xy, theta):
    # https://en.wikipedia.org/wiki/Rotation_matrix#In_two_dimensions
    #First Rotate based on Theta
    #Allocate Arrays
    rotationMatrix = np.zeros((2,2))
    #Calculate rotation matrix
    rotationMatrix[0,0] = np.cos(theta)
    rotationMatrix[0,1] = -1.0*np.sin(theta)
    rotationMatrix[1,0] = np.sin(theta)
    rotationMatrix[1,1] = np.cos(theta)
    return rotationMatrix.dot(xy)

def CalcLabAngle(pos):
    #Find swimming axis (normal y-axis)
    xU, xL = pos.loc[0,'xU'], pos.loc[0,'xL']
    yU, yL = pos.loc[0,'yU'], pos.loc[0,'yL']
    labX   = xU - xL
    labY   = yU - yL
    length = np.hypot(labX,labY)
    normX = labX/length
    normY = labY/length
    #2) Calculate Theta
    if(normX <= 0.0):
        theta = np.arccos(normY)
    else:
        theta = -1.0*np.arccos(normY)+2.0*np.pi
    print('theta = ',theta*180.0/np.pi)
    return 2.0*np.pi - theta

def InterpolateToNewCoordinateSystem(x,y,mx,my,arrayUx,arrayUy):
    #Create a uniform mesh for the interpolated velocity vectors!
    mx_new, my_new = np.meshgrid(x,y)
    
    #Interpolate Ux and Uy from original cartesian coordainates to new ones
    #Griddata
    print('About to inteprolate field data')
    print('peak memory = ',getrusage(RUSAGE_SELF).ru_maxrss)
    sys.stdout.flush()
    arrayUx_new=interpolate.griddata((mx.flatten(),my.flatten()),arrayUx.flatten() , (mx_new,my_new),method='linear')
    print('X transformation complete')
    print('peak memory = ',getrusage(RUSAGE_SELF).ru_maxrss)
    sys.stdout.flush()
    arrayUy_new=interpolate.griddata((mx.flatten(),my.flatten()),arrayUy.flatten() , (mx_new,my_new),method='linear')
    print('Coordinate Transformation Complete!')
    print('peak memory = ',getrusage(RUSAGE_SELF).ru_maxrss)
    sys.stdout.flush()
    return (arrayUx_new,arrayUy_new)

def RotateSimulation(cwd,time,mx,my,Ux,Uy,pos):
    global RADIUSLARGE
    #Shift x and y by the CM location
    xCM = 0.8*pos.loc[0,'xU'] + 0.2*pos.loc[0,'xL']
    yCM = 0.8*pos.loc[0,'yU'] + 0.2*pos.loc[0,'yL']
    #Do the same for mx and my
    mx -= xCM
    my -= yCM
    #Shift pos data by xCM and yCM
    pos['xU'] -= xCM
    pos['xL'] -= xCM
    pos['yU'] -= yCM
    pos['yL'] -= yCM
    #Rotate Reference frame by swimmer 1's axis
    #Calculate Theta (Rotate by -Theta)
    theta_rotate = CalcLabAngle(pos)
    print('theta_rotate = ',theta_rotate*180.0/np.pi)
    mxy = np.array([mx.flatten(),my.flatten()])
    mxy_rot = np.zeros((2,1024*1024))
    #Do the same for the U field
    Uxy = np.array([Ux.flatten(),Uy.flatten()])
    Uxy_rot = np.zeros((2,1024*1024))
    for jdx in range(1024*1024):
        mxy_rot[:,jdx] = Rotate(mxy[:,jdx],theta_rotate)
        Uxy_rot[:,jdx] = Rotate(Uxy[:,jdx],theta_rotate)
    mx_rot = mxy_rot[0,:].reshape((1024,1024))
    my_rot = mxy_rot[1,:].reshape((1024,1024))
    Ux_rot = Uxy_rot[0,:].reshape((1024,1024))
    Uy_rot = Uxy_rot[1,:].reshape((1024,1024))

    aU_pos = np.array([pos.loc[0,'xU'],pos.loc[0,'yU']])
    aL_pos = np.array([pos.loc[0,'xL'],pos.loc[0,'yL']])
    aU_rot = Rotate(aU_pos,theta_rotate)
    print('aU = ',aU_pos)
    print('aU_rot = ',aU_rot)
    aL_rot = Rotate(aL_pos,theta_rotate)
    pos['aXU_rot'], pos['aYU_rot'] = aU_rot[0], aU_rot[1]
    pos['aXL_rot'], pos['aYL_rot'] = aL_rot[0], aL_rot[1]
    #Interpolate onto a new coordinate system
    x = np.linspace(-0.025/RADIUSLARGE,0.025/RADIUSLARGE,512)
    y = np.linspace(-0.025/RADIUSLARGE,0.025/RADIUSLARGE,512)
    mx_stream, my_stream = np.meshgrid(x,y)
    interpUx, interpUy = InterpolateToNewCoordinateSystem(x,y,mx_rot,my_rot,Ux_rot,Uy_rot)
    
    return (mx_stream, my_stream, interpUx, interpUy, pos)

#Plot New mesh and interpolated velocity field Ux and Uy
def PlotAvgU(cwd,mx,my,Ux,Uy,pos,space):
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
    fig.savefig(cwd+'avgU.png')
    fig.clf()
    plt.close()
    return

#Plot New mesh and interpolated velocity field Ux and Uy
def PlotParticles(mx,my,U,pos,particles,frame,cwd):
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
    ax.axis([minVal,maxVal,minVal,maxVal])
    fig.tight_layout()
    fig.savefig(cwd+'PartStream_{0}_.png'.format(frame))
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
        indices = nX*self.idy[idTime]+self.idx[idTime]
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
    #Calculate # Periods
    DUMP_INT = 20.0
    nTime = GetPosDataLength(cwd_POS,Re)
    nPer = int(np.trunc(1.0*nTime/DUMP_INT))
    #nPer = 2
    #Paths to data and plots
    cwd_DATA = cwd_Re
    countPer = 0
    for countPer in range(nPer):
        if(countPer == perNumber):
            AVGPlot = pathlib.Path(cwd_DATA+'AVG_%04d.csv'%countPer)
            if AVGPlot.exists ():
                start = t.clock()
                #Get Avg Field Data
                mx,my,avgW,avgP,avgUx,avgUy = GetAvgFieldData(cwd_DATA,countPer)
                #Extract Position and Time Data
                time = np.round(0.05 + countPer*PERIOD,2)
                posData = GetPosData(cwd_POS,time,Re)
                #Plot Averaged Field Data
                #Vorticity And Streamlines
                mx,my,avgUx,avgUy,posData = RotateSimulation(cwd_PYTHON,time,mx,my,avgUx,avgUy,posData)
                stend = t.clock()
                diff = stend - start
                print('Time to run for 1 period = %.5fs'%diff)
                sys.stdout.flush()

    #Visual Check of vel field's data (any idx)
    PlotAvgU(cwd_FIGS,mx,my,avgUx,avgUy,posData,6)

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
    #print(avgU[1,281,248])
    avgU = np.array([[avgUx],[avgUy]])
    avgUx = avgUx.flatten()
    avgUy = avgUy.flatten()
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
            particles.x[idTime+1:nTime] += changeX
            particles.y[idTime+1:nTime] += changeY
            #Increase age of particles
            particles.curr_age[:idTime+1] -= 1
    #Save Initial particle stream pos and vel
    particles.xinit = particles.x[0,:].copy()#particles.x.copy()
    particles.yinit = particles.y[0,:].copy()#particles.y.copy()
    particles.vxinit = particles.vx[0,:].copy()#particles.vx.copy()
    particles.vyinit = particles.vy[0,:].copy()#particles.vy.copy()

    '''
    print('B4')
    print('life[100] = ',particles.age[:,100])
    print('curr_age[100] = ',particles.curr_age[:,100])
    print('x[100] = ',particles.x[:,100])
    print('y[100] = ',particles.y[:,100])
    print('vx[100] = ',particles.vx[:,100])
    print('vy[100] = ',particles.vy[:,100])
    '''
    #Loop over # of frames
    for idxFrame in range(2*nFrames):
        
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
            PlotParticles(mx,my,avgU,posData,particles,idxFrame-nFrames,cwd_FIGS)
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

    os.chdir(cwd_FIGS)
    os.system("ffmpeg -r 60 -i PartStream_%d_.png -vcodec libx264 -pix_fmt yuv420p -y PartStream.mp4")
    os.system("rm -rf PartStream_*")
    os.chdir(cwd_PYTHON)
    print('Movie {0}: Re{1} is complete'.format(config,Re))
