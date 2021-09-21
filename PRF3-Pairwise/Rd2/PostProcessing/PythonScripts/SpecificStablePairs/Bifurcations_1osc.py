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

minVal, maxVal = -6.0,6.0
dX = 2.0*maxR/(1.0*NX)
if local:
    cwd_FIGS = cwd_PYTHON+"../../Figures/VorticityDetection/{0}/".format(config)
    pathlib.Path(cwd_FIGS).mkdir(parents=True, exist_ok=True)
    cwd_Re = cwd_PYTHON+'../../FieldData/TestField/'
    cwd_POS = cwd_Re
else:
    cwd_FIGS = cwd_PYTHON+'../Figures/Bifurcation/{0}/'.format(config)
    pathlib.Path(cwd_FIGS).mkdir(parents=True, exist_ok=True)
    cwd_Re = cwd_PYTHON+'../{0}/Re{1}/VTK/AVG/'.format(config,Re)
    cwd_POS = cwd_PYTHON+'../PosData/{0}/Re{1}/'.format(config,Re)

# constructs a filepath for the pos data of Re = $Re
def pname(cwd):
    return cwd+"pd.txt"

def GetPosData(cwd,time,config):
    global RADIUSLARGE
    data = pd.read_csv(pname(cwd),delimiter=' ')
    if(config == 'V' or config == 'O'):
        pos = data[data['time'] == time*2.0]
    else:
        pos = data[data['time'] == time]
    #pos = data[data['time'] == time]
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

def InterpolateToNewCoordinateSystem(x,y,mx,my,arrayUx,arrayUy,arrayW):
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
    arrayW_new=interpolate.griddata((mx.flatten(),my.flatten()),arrayW.flatten() , (mx_new,my_new),method='linear')
    print('Vorticity Transformation Complete!')
    print('peak memory = ',getrusage(RUSAGE_SELF).ru_maxrss)
    sys.stdout.flush()
    return (arrayUx_new,arrayUy_new, arrayW_new)

def RotateSimulation(cwd,time,mx,my,Ux,Uy,W,pos):
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
    interpUx, interpUy, interpW = InterpolateToNewCoordinateSystem(x,y,mx_rot,my_rot,Ux_rot,Uy_rot, W)
    
    return (mx_stream.T, my_stream.T, interpW.T, interpUx.T, interpUy.T, pos)

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
    levels = MaxNLocator(nbins=21).tick_values(-1.0*max(abs(psi2.min()),psi2.max()), max(abs(psi2.min()),psi2.max()))
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
    ax.axis([xmin-1.0,xmax+1.0,ymin-1.0,ymax+1.0])
    fig.tight_layout()
    fig.savefig(cwd+'W_{0}_Re{1}_per{2}_.png'.format(config,Re,perNumber))
    fig.clf()
    plt.close()
    return

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
                posData = GetPosData(cwd_POS,time,config)
                #Plot Averaged Field Data
                #Vorticity And Streamlines
                mx,my,avgW,avgUx,avgUy,posData = RotateSimulation(cwd_PYTHON,time,mx,my,avgUx,avgUy,avgW,posData)
                rotatedDict = {'mx':mx.flatten(),'my':my.flatten(),
                               'avgUx':avgUx.flatten(),'avgUy':avgUy.flatten(),
                               'avgW':avgW.flatten()
                              }
                rotatedData = pd.DataFrame(data=rotatedDict)
                rotatedData.to_csv(cwd_DATA+'AVGRot_%04d.csv'%countPer,index=False,sep=' ',float_format='%.5e')
                posData.to_csv(cwd_POS+'pd_rot_%04d.csv'%countPer,index=False,sep=' ',float_format='%.5e')
                stend = t.clock()
                diff = stend - start
                print('Time to run for 1 period = %.5fs'%diff)
                sys.stdout.flush()

    #Plot Flow Field Visual
    PlotAvgW(cwd_FIGS,mx,my,avgW,avgUx,avgUy,posData,4,5)

