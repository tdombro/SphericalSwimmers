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
Theta = sys.argv[1]
Hx = sys.argv[2]
Hy = sys.argv[3]
perNumber = int(sys.argv[4])#5

#PARTICLE STREAMLINE CONSTANTS
nFrames = 240
minVal, maxVal = -6.0,6.0
rows, cols = 100, 100
nPart, nTrail = rows*cols, 80
timestep = 0.2/60.0
dX = 2.0*maxR/(1.0*NX)
seed = random.seed(11235)
cwd_FIGS = cwd_PYTHON+'../Figures/ParticleStreamline/L/per{0}/'.format(perNumber)#cwd_PYTHON+"../../Figures/ParticleStreamline/TestField2/"
pathlib.Path(cwd_FIGS).mkdir(parents=True, exist_ok=True)
cwd_Re = cwd_PYTHON+'Theta{0}/Hx{1}/Hy{2}/VTK/AVG/'.format(Theta,Hx,Hy)#cwd_PYTHON+'../../FieldData/TestField/'
cwd_POS = cwd_PYTHON+'Theta{0}/Hx{1}/Hy{2}/'.format(Theta,Hx,Hy)#cwd_POS = cwd_PYTHON+'../../FieldData/TestField/'

# constructs a filepath for the pos data of Re = $Re
def pname(cwd):
    #return cwd+"/pd.txt"
    #cwd = cwd_PYTHON
    return cwd+"pd.txt"

def GetPosData(cwd,time,parTheta,parHx,parHy):
    global RADIUSLARGE
    
    #Load position data
    pdData = pd.read_csv(pname(cwd),delimiter=' ')
    #Split up individual sphere data by given index
    UAdata = pdData[pdData['idx'] == 6].copy()
    LAdata = pdData[pdData['idx'] == 19].copy()
    UBdata = pdData[pdData['idx'] == 32].copy()
    LBdata = pdData[pdData['idx'] == 45].copy()
    #Sort data by time and reset indices
    UAdata = UAdata.sort_values(by=['time'])
    LAdata = LAdata.sort_values(by=['time'])
    UBdata = UBdata.sort_values(by=['time'])
    LBdata = LBdata.sort_values(by=['time'])
    UAdata = UAdata.reset_index(drop=True)
    LAdata = LAdata.reset_index(drop=True)
    UBdata = UBdata.reset_index(drop=True)
    LBdata = LBdata.reset_index(drop=True)
    #Rename columns to previous data frames
    UAdata = UAdata.rename(columns={"x":"aXU", "y":"aYU"})
    LAdata = LAdata.rename(columns={"x":"aXL", "y":"aYL"})
    UBdata = UBdata.rename(columns={"x":"bXU", "y":"bYU"})
    LBdata = LBdata.rename(columns={"x":"bXL", "y":"bYL"})
    #Combine separate dataframes to create previous dataframe used
    splitDict = {'aXU':UAdata['aXU'],'aYU':UAdata['aYU'],'aXL':LAdata['aXL'],'aYL':LAdata['aYL'],
        'bXU':UBdata['bXU'],'bYU':UBdata['bYU'],'bXL':LBdata['bXL'],'bYL':LBdata['bYL'],'time':UAdata['time']}
    posData = pd.DataFrame(data=splitDict)
    pos = posData[posData['time'] == time]
    pos = pos.reset_index(drop=True)

    #Renormalize
    pos['aXU'] /= RADIUSLARGE
    pos['aXL'] /= RADIUSLARGE
    pos['aYU'] /= RADIUSLARGE
    pos['aYL'] /= RADIUSLARGE
    pos['bXU'] /= RADIUSLARGE
    pos['bXL'] /= RADIUSLARGE
    pos['bYU'] /= RADIUSLARGE
    pos['bYL'] /= RADIUSLARGE
    return pos

def GetPosDataLength(cwd):
    data = pd.read_csv(pname(cwd),delimiter=' ')
    return len(data['time'])

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
    circle3 = Circle((pos.loc[0,'bXU_rot'], pos.loc[0,'bYU_rot']), 1.0, facecolor=(0.5,)*3,
                  linewidth=1,alpha=1.0,zorder=6)
    ax.add_patch(circle3)
    circle4 = Circle((pos.loc[0,'bXL_rot'], pos.loc[0,'bYL_rot']), 0.5, facecolor=(0.5,)*3,
                  linewidth=1,alpha=1.0,zorder=6)
    ax.add_patch(circle4)
    #Add Swimmer "springs"
    ax.plot([pos.loc[0,'aXU_rot'],pos.loc[0,'aXL_rot']],
         [pos.loc[0,'aYU_rot'],pos.loc[0,'aYL_rot']],
         color='black',linewidth=3,zorder=6)
    ax.plot([pos.loc[0,'bXU_rot'],pos.loc[0,'bXL_rot']],
         [pos.loc[0,'bYU_rot'],pos.loc[0,'bYL_rot']],
         color=(0.5,)*3,linewidth=3,zorder=6)
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
    xU, xL = pos.loc[0,'aXU'], pos.loc[0,'aXL']
    yU, yL = pos.loc[0,'aYU'], pos.loc[0,'aYL']
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
    xCM = 0.25*(pos.loc[0,'aXU'] + pos.loc[0,'bXU'] + pos.loc[0,'aXL'] + pos.loc[0,'bXL'])
    yCM = 0.25*(pos.loc[0,'aYU'] + pos.loc[0,'bYU'] + pos.loc[0,'aYL'] + pos.loc[0,'bYL'])
    #Do the same for mx and my
    mx -= xCM
    my -= yCM
    #Shift pos data by xCM and yCM
    pos['aXU'] -= xCM
    pos['aXL'] -= xCM
    pos['bXU'] -= xCM
    pos['bXL'] -= xCM
    pos['aYU'] -= yCM
    pos['aYL'] -= yCM
    pos['bYU'] -= yCM
    pos['bYL'] -= yCM
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

    aU_pos = np.array([pos.loc[0,'aXU'],pos.loc[0,'aYU']])
    aL_pos = np.array([pos.loc[0,'aXL'],pos.loc[0,'aYL']])
    bU_pos = np.array([pos.loc[0,'bXU'],pos.loc[0,'bYU']])
    bL_pos = np.array([pos.loc[0,'bXL'],pos.loc[0,'bYL']])
    aU_rot = Rotate(aU_pos,theta_rotate)
    print('aU = ',aU_pos)
    print('aU_rot = ',aU_rot)
    aL_rot = Rotate(aL_pos,theta_rotate)
    bU_rot = Rotate(bU_pos,theta_rotate)
    bL_rot = Rotate(bL_pos,theta_rotate)
    pos['aXU_rot'], pos['aYU_rot'] = aU_rot[0], aU_rot[1]
    pos['aXL_rot'], pos['aYL_rot'] = aL_rot[0], aL_rot[1]
    pos['bXU_rot'], pos['bYU_rot'] = bU_rot[0], bU_rot[1]
    pos['bXL_rot'], pos['bYL_rot'] = bL_rot[0], bL_rot[1]
    #Interpolate onto a new coordinate system
    x = np.linspace(-0.025/RADIUSLARGE,0.025/RADIUSLARGE,512)
    y = np.linspace(-0.025/RADIUSLARGE,0.025/RADIUSLARGE,512)
    mx_stream, my_stream = np.meshgrid(x,y)
    interpUx, interpUy = InterpolateToNewCoordinateSystem(x,y,mx_rot,my_rot,Ux_rot,Uy_rot)
    
    return (mx_stream, my_stream, interpUx, interpUy, pos)

#Plot New mesh and interpolated velocity field Ux and Uy
def PlotAvgU(cwd,mx,my,Ux,Uy,pos,space):
    global FIGNUM, PERIOD,minVal,maxVal,Theta,Hx,Hy
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
    fig.savefig(cwd+'avgU_T{0}_Hx{1}_Hy{2}_.png'.format(Theta,np.round(float(Hx)/2.0,2),np.round(float(Hy)/2.0,1)))
    fig.clf()
    plt.close()
    return

#Plot New mesh and interpolated velocity field Ux and Uy
def PlotParticles(mx,my,U,pos,particles,frame,cwd):
    global FIGNUM, PERIOD,minVal,maxVal,Theta,Hx,Hy
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
    fig.savefig(cwd+'PartStream_T{0}_Hx{1}_Hy{2}_{3}_.png'.format(Theta,np.round(float(Hx)/2.0,2),np.round(float(Hy)/2.0,1),frame))
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
    nTime = GetPosDataLength(cwd_POS)
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
                posData = GetPosData(cwd_POS,time,float(Theta),float(Hx),float(Hy))
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
    strMovie = "ffmpeg -r 60 -i PartStream_T{0}_Hx{1}_Hy{2}_%d_.png -vcodec libx264 -pix_fmt yuv420p -y PartMov_T{0}_Hx{1}_Hy{2}_.mp4".format(Theta,np.round(float(Hx)/2.0,2),np.round(float(Hy)/2.0,1))
    os.system(strMovie)
    os.system("rm -rf PartStream_T{0}_Hx{1}_Hy{2}_*".format(Theta,np.round(float(Hx)/2.0,2),np.round(float(Hy)/2.0,1)))
    os.chdir(cwd_PYTHON)
    print('Movie T{0}: Hx{1}: Hy{2} is complete'.format(Theta,Hx,Hy))
