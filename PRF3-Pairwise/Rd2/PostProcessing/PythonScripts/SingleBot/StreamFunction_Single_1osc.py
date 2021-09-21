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
Re = sys.argv[1]#"2"
perNumber = int(sys.argv[2])#5
local = int(sys.argv[3])

minVal, maxVal = -6.0,6.0
dX = 2.0*maxR/(1.0*NX)
if local:
    cwd_FIGS = cwd_PYTHON+'../../Figures/VorticityDetection/Test/SingleBot/'
    pathlib.Path(cwd_FIGS).mkdir(parents=True, exist_ok=True)
    cwd_Re = cwd_PYTHON+'../../FieldData/SingleBot/'#
    cwd_POS = cwd_Re
else:
    cwd_FIGS = cwd_PYTHON+'../Figures/VorticityDetection/'
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
    #Add Saddles at edge of spheres
    #Large Sphere
    xSaddle = [pos.loc[0,'aXU_rot']-1.0,pos.loc[0,'aXU_rot'],pos.loc[0,'aXU_rot'],pos.loc[0,'aXU_rot']+1.0]
    ySaddle = [pos.loc[0,'aYU_rot'],pos.loc[0,'aYU_rot']-1.0,pos.loc[0,'aYU_rot']+1.0,pos.loc[0,'aYU_rot']]
    ax.scatter(xSaddle,ySaddle,color='purple',s=20,zorder=10,marker='D')
    xSaddle = [pos.loc[0,'aXL_rot']-0.5,pos.loc[0,'aXL_rot'],pos.loc[0,'aXL_rot'],pos.loc[0,'aXL_rot']+0.5]
    ySaddle = [pos.loc[0,'aYL_rot'],pos.loc[0,'aYL_rot']-0.5,pos.loc[0,'aYL_rot']+0.5,pos.loc[0,'aYL_rot']]
    ax.scatter(xSaddle,ySaddle,color='purple',s=20,zorder=10,marker='D')
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

def InterpolateToNewCoordinateSystem(x,y,mx,my,arrayW,arrayUx,arrayUy):
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
    arrayW_new=interpolate.griddata((mx.flatten(),my.flatten()),arrayW.flatten() , (mx_new,my_new),method='linear')
    print('Coordinate Transformation Complete!')
    print('peak memory = ',getrusage(RUSAGE_SELF).ru_maxrss)
    sys.stdout.flush()
    return (arrayW_new,arrayUx_new,arrayUy_new)

def RotateSimulation(cwd,time,mx,my,W,Ux,Uy,pos):
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
    interpW,interpUx, interpUy = InterpolateToNewCoordinateSystem(x,y,mx_rot,my_rot,W,Ux_rot,Uy_rot)
    
    return (mx_stream.T, my_stream.T, interpW.T, interpUx.T, interpUy.T, pos)

#Plot New mesh and interpolated velocity field Ux and Uy
def PlotStreamAndPsi(cwd,mx,my,W,Ux,Uy,psi,pos):
    global FIGNUM, PERIOD,minVal,maxVal, Re, perNumber,maxR
    #Here, we will visualize the velocity field on the new coordinate system
    nRows, nCols = 1, 2
    fig, ax = plt.subplots(nrows=nRows, ncols=nCols, num=0,figsize=(12,6),dpi=200)
    ax[0].set_title(r'U,V Streamplot',fontsize=12)
    ax[1].set_title(r'Psi Contour',fontsize=12)
    x = np.linspace(-6.125,6.125,256)
    y = np.linspace(-6.125,6.125,256)
    mx_stream, my_stream = np.meshgrid(x,y)
    UxT, UyT, WT = Ux.T, Uy.T, W.T
    Umag = np.hypot(UxT,UyT)
    lw = 2 * Umag/Umag.max()
    lw = np.where(lw < 1.0, 1.0, lw)
    #Streamplot
    ax[0].streamplot(x,y,UxT[128:384,128:384],UyT[128:384,128:384],
                  color = 'k',
                  linewidth=1.5,minlength=0.3,
                  arrowsize=0.01,density=2.0)
    #Plot magnitude with contourf
    levels = MaxNLocator(nbins=21).tick_values(-1.0, 1.0)
    #ax[0].contourf(mx,my,W,cmap='bwr',levels=levels,extend='both')
    ax[0].imshow(W.T,cmap='bwr',extent=(-1.0*maxR-0.5*dX,maxR+0.5*dX,
                                        -1.0*maxR-0.5*dX,maxR+0.5*dX),
                 origin='lower',vmin=-1.0,vmax=1.0,interpolation='bilinear')
    #levels = [-2.5,-2.0,-1.5,1.5,2.0,2.5]
    #ax[0].contour(mx,my,psi,colors='k',extend='both',levels=levels,linewidth=1.0)
    #Add swimmer
    AddDiscsToPlot(ax[0],pos)
    
    #Psi Contour
    psi2 = ndimage.gaussian_filter(psi, sigma=5.0, order=0)
    levels = MaxNLocator(nbins=11).tick_values(0.0, max(abs(psi2.min()),psi2.max()))
    #levels = [-2.5,-2.3,-2.1,-1.9,-1.7,-1.5,1.5,1.7,1.9,2.1,2.3,2.5]
    ax[1].contour(mx,my,abs(psi2),colors='k',extend='both',levels=levels)
    #Plot magnitude with contourf
    #levels = MaxNLocator(nbins=21).tick_values(-psi2.max(), psi2.max())
    #ax[1].contourf(mx,my,W,cmap='bwr',levels=levels,extend='both')
    ax[1].imshow(W.T,cmap='bwr',extent=(-1.0*maxR-0.5*dX,maxR+0.5*dX,
                                        -1.0*maxR-0.5*dX,maxR+0.5*dX),
                 origin='lower',vmin=-1.0,vmax=1.0,interpolation='bilinear')
    #Add swimmer
    AddDiscsToPlot(ax[1],pos)
    ax[0].axis([minVal,maxVal,minVal,maxVal])
    ax[1].axis([minVal,maxVal,minVal,maxVal])
    fig.tight_layout()
    pathlib.Path(cwd+'Zoom/').mkdir(parents=True, exist_ok=True)
    fig.savefig(cwd+'Zoom/Psi_Re{0}_per{1}_.png'.format(Re,perNumber))
    fig.clf()
    plt.close()
    return

def CalcPsi2D(fx,fy,NX,DX):
    #From here, we are going to calculate the stream function
    psi = np.zeros((NX,NX))
    for idx in range(1,NX):
        psi[idx,0] = psi[idx-1,0] - fy[idx,0]*DX
    for idy in range(1,NX):
        psi[:,idy] = psi[:,idy-1] + fx[:,idy]*DX

    return psi

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
                mx,my,avgW,avgUx,avgUy,posData = RotateSimulation(cwd_PYTHON,time,mx,my,avgW,avgUx,avgUy,posData)
                stend = t.clock()
                diff = stend - start
                print('Time to run for 1 period = %.5fs'%diff)
                sys.stdout.flush()

    psi = CalcPsi2D(avgUx,avgUy,NX,dX)
    print('min psi = ',psi.min())
    print('max psi = ',psi.max())
    PlotStreamAndPsi(cwd_FIGS,mx,my,avgW,avgUx,avgUy,psi,posData)




















