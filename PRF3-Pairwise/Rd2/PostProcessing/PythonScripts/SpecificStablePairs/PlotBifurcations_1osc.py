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
    levels = MaxNLocator(nbins=101).tick_values(-1.0*max(abs(psi2.min()),psi2.max()), max(abs(psi2.min()),psi2.max()))
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
    ax.axis([xmin-0.25,xmax+0.25,ymin-2.0,ymax+2.0])
    #ax.axis([-5.0,5.0,-5.0,5.0])
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
    PlotAvgW(cwd_FIGS,mx,my,avgW,avgUx,avgUy,posData,4,5)

