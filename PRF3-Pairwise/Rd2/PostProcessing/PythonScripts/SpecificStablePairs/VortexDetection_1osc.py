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
    cwd_FIGS = cwd_PYTHON+'../Figures/VorticityDetection/{0}/'.format(config)
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
            
    #Add Saddles at edge of spheres
    #Swimmer1
    xSaddle = [pos.loc[0,'aXU_rot']-1.0,pos.loc[0,'aXU_rot'],pos.loc[0,'aXU_rot'],pos.loc[0,'aXU_rot']+1.0]
    ySaddle = [pos.loc[0,'aYU_rot'],pos.loc[0,'aYU_rot']-1.0,pos.loc[0,'aYU_rot']+1.0,pos.loc[0,'aYU_rot']]
    ax.scatter(xSaddle,ySaddle,color='lime',s=30,zorder=10,marker='o')
    xSaddle = [pos.loc[0,'aXL_rot']-0.5,pos.loc[0,'aXL_rot'],pos.loc[0,'aXL_rot'],pos.loc[0,'aXL_rot']+0.5]
    ySaddle = [pos.loc[0,'aYL_rot'],pos.loc[0,'aYL_rot']-0.5,pos.loc[0,'aYL_rot']+0.5,pos.loc[0,'aYL_rot']]
    ax.scatter(xSaddle,ySaddle,color='lime',s=30,zorder=10,marker='o')
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

def MaximumVorticityMethod(avgW):
    #Let's first try identifying vortex centers with the maximum vorticity method
    #Here is the pseudocode for it
    #1: compute ω at all grid nodes
    #2: for all cell faces do
    #3: examine its 4×4 surrounding nodes
    #4: if ∃ maximum |ω| in central nodes then
    #5: mark grid face as candidate face
    #6: end if
    #7: end for
    #8: for all candidate faces do
    #9: compute ∇|ω| using central difference at nodes
    #10: compute solution points where ∇|ω| = 0
    #11: if points within face and are local maxima then
    #12: mark them as vortex core points
    #13: end if
    #14: end for
    
    #Calculate the vorticity of the rotated velocity field
    #gradavgUx_x, gradavgUx_y = grad(avgUx,h=dX)
    #gradavgUy_x, gradavgUy_y = grad(avgUy,h=dX)
    #avgW = gradavgUy_x - gradavgUx_y
    
    centerCandidates = []
    #Loop over all cells:
    NXW = NX
    for idx in range(1,NXW-1):
        for idy in range(1,NXW-1):
            #Create a small nodegroup (mesh) of vorticities
            #Check to see if (idx,idy) is a maximum
            #If it is, then save the index in a candidateList
            nodeGroup = avgW[idx-1:idx+2,idy-1:idy+2].copy()
            maxValue = 0.0
            #Make sure vorticity is not close to zero
            if abs(avgW[idx,idy]) >= 0.1:
                #Check if center cell is maximum value
                maxValue = np.amax(abs(nodeGroup))
                if maxValue == abs(avgW[idx,idy]):
                    centerCandidates.append([idx,idy])
    #print('centerCandidates = ',centerCandidates)

    savedCenters = []
    #Now that all candidates have been chosen, calculate del(|w|) using
    #central difference at nodes
    for candidates in centerCandidates:
        idx, idy = candidates[0], candidates[1]
        delWx = (abs(avgW[idx+1,idy]) - abs(avgW[idx-1,idy]))/(2.0*dX)
        delWy = (abs(avgW[idx,idy+1]) - abs(avgW[idx,idy-1]))/(2.0*dX)
        delW = np.hypot(delWx,delWy)
        #print('delW = ',delW)
        tolerance = 1.0e-2
        if delW <= tolerance:
            #Save as vorticity center
            savedCenters.append([idx,idy])
    savedCenters = np.array(savedCenters)
    print('saved centers = ',savedCenters)
    return savedCenters

def CombinatorialMethod(avgUx,avgUy,nLabels):
    global NX, dX
    #Combinatorial method
    #1: for all grid cells do
    #2: compute swirl plane normal n at cell center
    #3: project v from surrounding nodes
    #4: for all vp in swirl plane do
    #5: compute its angle α from local x-axis
    #6: label direction range for α
    #7: end for
    #8: if all direction ranges are labeled then
    #9: mark grid cell as vortex core
    #10: end if
    #11: end for
    
    #Swirl plane is just (Ux,Uy)
    #For each cell, calculate velocity angle
    #label direction range: 3=[0-120,120-240,240-360]
    #4 = [0-90,90-180,180-270,270-360]
    Umag = np.hypot(avgUx,avgUy)
    nUx = avgUx/Umag
    nUy = avgUy/Umag
    thetaArr = np.zeros((NX,NX))
    labelArr = np.zeros((NX,NX))
    
    for idx in range(NX):
        for idy in range(NX):
            if nUy[idx,idy] >= 0.0:
                thetaArr[idx,idy] = np.arccos(nUx[idx,idy])
            else:
                thetaArr[idx,idy] = -1.0*np.arccos(nUx[idx,idy])+2.0*np.pi
    #label 4 direction spans
    if nLabels == 4:
        labelArr[(thetaArr <= np.pi/2.0)] = 0
        labelArr[(thetaArr > np.pi/2.0) & (thetaArr <= np.pi)] = 1
        labelArr[(thetaArr > np.pi) & (thetaArr <= 1.5*np.pi)] = 2
        labelArr[(thetaArr > 1.5*np.pi)] = 3
    #label 3 direction spans
    elif nLabels == 3:
        labelArr[(thetaArr <= 2.0*np.pi/3.0)] = 0
        labelArr[(thetaArr > 2.0*np.pi/3.0) & (thetaArr <= 4.0*np.pi/3.0)] = 1
        labelArr[(thetaArr > 4.0*np.pi/3.0)] = 2
    
    print('thetaArr = ',thetaArr)
    print('labelArr = ',labelArr)
    sys.stdout.flush()
    PlotCheckLabels(cwd_FIGS,mx,my,nUx,nUy,labelArr,posData,8)
    fixedPoints = []
    mask = np.ones(9,dtype=bool)
    for idx in range(1,NX-1):
        for idy in range(1,NX-1):
            #Create Neighbor group
            neighborLabel = labelArr[idx-1:idx+2,idy-1:idy+2].reshape(9)
            mask[5] = False
            neighborLabel = neighborLabel[mask]
            neighborSet = list(np.unique(neighborLabel))
            neighborSet = [int(a) for a in neighborSet]
            #Check to see if direction spanning is complete
            
            if (nLabels == 4 and neighborSet == [0,1,2,3]):
                #We have found a vortex candidate
                print('Fixed point ({0},{1}) has been found'.format(idx,idy))
                fixedPoints.append([idx,idy])
                sys.stdout.flush()
            elif (nLabels == 3 and neighborSet == [0,1,2]):
                #We have found a vortex candidate
                print('Fixed point ({0},{1}) has been found'.format(idx,idy))
                fixedPoints.append([idx,idy])
                sys.stdout.flush()
    return fixedPoints

def CalcPsi2D(fx,fy,NX,DX):
    #From here, we are going to calculate the stream function
    psi = np.zeros((NX,NX))
    for idx in range(1,NX):
        psi[idx,0] = psi[idx-1,0] - fy[idx,0]*DX
    for idy in range(1,NX):
        psi[:,idy] = psi[:,idy-1] + fx[:,idy]*DX
    return psi

#Plot New mesh and interpolated velocity field Ux and Uy
def PlotAvgW(cwd,mx,my,W,Ux,Uy,pos,space,scale,cores,saddles,zoom):
    global FIGNUM, PERIOD,minVal,maxVal, Re, perNumber
    #Here, we will visualize the velocity field on the new coordinate system
    nRows, nCols = 1, 1
    fig, ax = plt.subplots(nrows=nRows, ncols=nCols, num=0,figsize=(6,6),dpi=200)
    ax.set_title(r'Average Velocity Field',fontsize=12)
    #Plot Streamlines
    #Use two grids and combine them
    UxT, UyT, WT = Ux.T, Uy.T, W.T
    if zoom:
        psi = CalcPsi2D(Ux,Uy,NX,dX)
        print('psi.min() = ',psi.min())
        print('psi.max() = ',psi.max())
        sys.stdout.flush()
        #Psi Contour
        psi2 = ndimage.gaussian_filter(psi, sigma=5.0, order=0)
        levels = MaxNLocator(nbins=21).tick_values(-1.0*max(abs(psi2.min()),psi2.max()), max(abs(psi2.min()),psi2.max()))
        ax.contour(mx,my,psi2,colors='k',extend='both',levels=levels)
    else:
        x = np.linspace(-12.5,12.5,512)
        y = np.linspace(-12.5,12.5,512)
        #Apply variable density to
        ax.streamplot(x,y,UxT,UyT,
                      #color=WT,cmap='bwr',
                      color = 'k',
                      #norm=Normalize(vmax=0.01, vmin=-0.01),
                      linewidth=1.0,arrowsize=0.01,density=2.0)
    #PlotVorticity with imshow (interpolate to smooth)
    ax.imshow(W.T,cmap='bwr',extent=(-1.0*maxR-0.5*dX,maxR+0.5*dX,
                                   -1.0*maxR-0.5*dX,maxR+0.5*dX),
            origin='lower',vmin=-1.0,vmax=1.0,interpolation='bilinear')
    #Add swimmer
    AddDiscsToPlot(ax,pos)
    #Add vorticity centers
    if list(cores):
        ax.scatter(cores[:,0],cores[:,1],color='purple',s=30,zorder=10,marker='D')
    if list(saddles):
        ax.scatter(saddles[:,0],saddles[:,1],color='lime',s=30,zorder=10,marker='o')
    #print('RSMALL = ',RSMALL)
    if zoom:
        ax.axis([minVal,maxVal,minVal,maxVal])
        fig.tight_layout()
        pathlib.Path(cwd+'Zoom/').mkdir(parents=True, exist_ok=True)
        fig.savefig(cwd+'Zoom/W_{0}_Re{1}_per{2}_.png'.format(config,Re,perNumber))
    else:
        ax.axis([-12.5,12.5,-12.5,12.5])
        fig.tight_layout()
        pathlib.Path(cwd+'WholeSim/').mkdir(parents=True, exist_ok=True)
        fig.savefig(cwd+'WholeSim/W_{0}_Re{1}_per{2}_.png'.format(config,Re,perNumber))
    fig.clf()
    plt.close()
    return

#Plot New mesh and interpolated velocity field Ux and Uy
def PlotCheckLabels(cwd,mx,my,Ux,Uy,labels,pos,space):
    global NX, dX, Re, perNumber, minVal, maxVal, config
    #Here, we will visualize the velocity field on the new coordinate system
    nRows, nCols = 1, 1
    fig, ax = plt.subplots(nrows=nRows, ncols=nCols, num=0,figsize=(6,6),dpi=200)
    ax.set_title(r'CheckAngleLabels',fontsize=12)
    #Plot Vector field with quiver
    #ax.scatter(mx[::space,::space],my[::space,::space],color='k',s=12)
    ax.quiver(mx[::space,::space],my[::space,::space],
              Ux[::space,::space],Uy[::space,::space],
              color='k',pivot='mid',angles='xy',scale_units='xy', scale=5,zorder=5)
    image = np.ones((NX,NX,3))
    for idx in range(NX):
        for jdx in range(NX):
            labelValue = labels[idx,jdx]
            if labelValue == 0:
                image[jdx,idx,:] = [1.0,0.0,0.0]
            elif labelValue == 1:
                image[jdx,idx,:] = [0.0,1.0,0.0]
            elif labelValue == 2:
                image[jdx,idx,:] = [0.0,0.0,1.0]
            elif labelValue == 3:
                image[jdx,idx,:] = (0.5,)*3
            else:
                image[jdx,idx] = (1.0,)*3
    ax.imshow(image,cmap=None,extent=(-1.0*maxR-0.5*dX,maxR+0.5*dX,
                                      -1.0*maxR-0.5*dX,maxR+0.5*dX),
              alpha=0.75,origin='lower',zorder=2)

    #Add swimmer
    AddDiscsToPlot(ax,pos)
    axis = [minVal,maxVal,minVal,maxVal]
    ax.axis(axis)
    fig.tight_layout()
    #plt.show()
    pathlib.Path(cwd+'CheckLabels/').mkdir(parents=True, exist_ok=True)
    fig.savefig(cwd+'CheckLabels/CheckLabels_{0}_Re{1}_per{2}_.png'.format(config,Re,perNumber))
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
                stend = t.clock()
                diff = stend - start
                print('Time to run for 1 period = %.5fs'%diff)
                sys.stdout.flush()

    #From here, we need to detect vorticity centers and saddle points
    #We will start with vortex centers

    #Maximum Vorticity Method
    #savedCenters = MaximumVorticityMethod(avgW)
    #Combinatorial Method
    fixedPoints = CombinatorialMethod(avgUx,avgUy,3)
    print('len(fixedPoints) = ',len(fixedPoints))

    #Average fixed points
    #Most likely, multiple fixed points have been identified around a specific region
    #We need to average those points to get 1 fixed point from the region
    #We will consider fixed points to be identical if their distance apart < threshold
    #Count up and average those in the same region
    #Remove fixed points that were averaged from the list
    #Repeat with next fixed point entry
    cleanedFixedPoints = []
    threshold = 5.0
    while fixedPoints:
        trackedIndices = []
        count = 1
        fP = fixedPoints[0]
        for idx in range(1,len(fixedPoints)):
            fP2 = fixedPoints[idx]
            #Calculate Distance apart from one another
            dist = np.hypot(fP[0] - fP2[0],fP[1] - fP2[1])
            #Check if dist is within threshold
            if dist <= threshold:
                trackedIndices.append(idx)
                count += 1
        #Check if trackedIndices exist
        sumFixed = np.array([mx[fP[0],fP[1]],my[fP[0],fP[1]]])
        if trackedIndices:
            for index in sorted(trackedIndices, reverse=True):
                fP2 = fixedPoints[index]
                sumFixed += np.array([mx[fP2[0],fP2[1]],my[fP2[0],fP2[1]]])
                fixedPoints.pop(index)
            avgFixed = sumFixed/(1.0*count)
            cleanedFixedPoints.append(avgFixed)
        else:
            #No other fixed points near
            cleanedFixedPoints.append([mx[fP[0],fP[1]],my[fP[0],fP[1]]])
        fixedPoints.pop(0)
                
    print('len(cleaned) = ',len(cleanedFixedPoints))
    print(cleanedFixedPoints)
    print('cleaned[0,:] = ',cleanedFixedPoints[0])
    print('cleaned[:,0] = ',cleanedFixedPoints[:][0])

    #Now, we need to determine if a fixedPoint is considered a vortex
    #We will test this using the velocity gradient tensor
    #If 2D and incompressible, then U_y*V_x - U_x*V_y < 0 indicates vortex
    #First, we must identify closest mesh point to identified fixedPoint
    fixedPointsIndex = []
    for fP in cleanedFixedPoints:
        print('fP = ',fP)
        print('fP[0] = ',fP[0])
        idx = int((fP[0] + maxR)/dX)
        idy = int((fP[1] + maxR)/dX)
        #meshX, meshY = mx[idx,idy], my[idx,idy]
        fixedPointsIndex.append([idx,idy])
    #Now we calculate gradients (using central difference) around fixed point

    savedVortexCores = []
    savedSaddles = []
    for fP in fixedPointsIndex:
        gradUx_y = (avgUx[fP[0],fP[1]+1] - avgUx[fP[0],fP[1]-1])/(2.0*dX)
        gradUx_x = (avgUx[fP[0]+1,fP[1]] - avgUx[fP[0]-1,fP[1]])/(2.0*dX)
        gradUy_y = (avgUy[fP[0],fP[1]+1] - avgUy[fP[0],fP[1]-1])/(2.0*dX)
        gradUy_x = (avgUy[fP[0]+1,fP[1]] - avgUy[fP[0]-1,fP[1]])/(2.0*dX)
        value = gradUx_y*gradUy_x - gradUx_x*gradUy_y
        print('value = ',value)
        print('point = ',[mx[fP[1],fP[0]],my[fP[1],fP[0]]])
        if value < 0.0:
            #We have a vortex!
            print('vortex: fP = ',fP)
            savedVortexCores.append([mx[fP[0],fP[1]],my[fP[0],fP[1]]])
        else:
            #Check if saddle (stagnation) point
            #Will be a stagnation point if velocity switches direction
            #We will check both x and y directions
            if (avgUx[fP[0]+1,fP[1]]*avgUx[fP[0]-1,fP[1]] <= 0.0 and
                avgUy[fP[0],fP[1]+1]*avgUy[fP[0],fP[1]-1] <= 0.0):
                print('saddle: ',fP)
                #We have a saddle point!
                savedSaddles.append([mx[fP[0],fP[1]],my[fP[0],fP[1]]])

    if cleanedFixedPoints:
        cleanedFixedPoints = np.array(cleanedFixedPoints)
        savedVortexCores = np.array(savedVortexCores)
        savedSaddles = np.array(savedSaddles)
        #Visual Check of vorticity and marked centers
        PlotAvgW(cwd_FIGS,mx,my,avgW,avgUx,avgUy,posData,4,5,savedVortexCores, savedSaddles, 1)
        PlotAvgW(cwd_FIGS,mx,my,avgW,avgUx,avgUy,posData,8,5,savedVortexCores, savedSaddles, 0)
    else:
        print('No fixed Points found')




















