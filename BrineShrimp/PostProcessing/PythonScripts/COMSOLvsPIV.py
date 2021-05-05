#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 12:32:43 2019

@author: thomas
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import pathlib
from matplotlib.ticker import (MultipleLocator,AutoMinorLocator)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

mpl.rcParams['axes.linewidth'] = 1.5 #set the value globally

#This script will be used to compare the velocity field from PIV and COMSOL data
#COMSOL DATA
#Obtain from ../COMSOL/Results/'mesh dim'/Velocity_'mesh dim'cm_clean.txt
#9 columns of data (x,y,z,x,y,z,u,v,w)
#Only need (u,v,w)
#Reorganize arrays (justU,justV,justW) to be of 3D grid
#Note: Make sure dimension values are in ascending order
#Plot U, V, and W along the XNorm and the ZNorm, each slice.
#Turn slices into a movie where we go from -y to y or -z to z

#CONSTANTS
meshNameList = ['8x15x1','8x15x2','8x15x3','10x20x1','10x20x2','10x20x3','15x30x1','15x30x2','15x30x3','empty']
#meshNameList = ['10x20x2','empty'] 
sliceList = ['btw','on'] 
#x indices where COMSOL data slice corresponds to PIV data
onIdxList = [17,17,17,18,18,18,19,19,19,20] 
btwIdxList = [15,15,15,16,16,16,18,18,18,19]
#onIdxList=[18,20]
#btwIdxList=[16,19]
velStringList = ['V','W','U']

#Obtain python directory    
cwd_PYTHON = os.getcwd()
TINY_NUMBER = 1.0e-15
FIGNUM = 1
PLOTNUM = 0
boolMovie = 0

def getCOMSOLInputData(meshName):
    #Obtain directory of data for timeValue
    cwd_DATA = cwd_PYTHON + '/../COMSOL/Results/'+meshName+'/'
    #Create an array for files U_x, U_y, U_z
    #arrays will be of size (40,320,40) where it is the data of coords (x,y,z)
    allData = np.loadtxt(cwd_DATA+'Velocity_'+meshName+'cm_clean.txt')
    #Split array up into columns
    arrayX = allData[:,0]
    arrayY = allData[:,1]
    arrayZ = allData[:,2]
    arrayU = allData[:,6]
    arrayV = allData[:,7]
    arrayW = allData[:,8]
    #Reshape arrays to be of 3d grid structure
    #Reshaped to where changes are (z,y,x) values
    arrayX = np.reshape(arrayX,(40,320,40))
    arrayY = np.reshape(arrayY,(40,320,40))
    arrayZ = np.reshape(arrayZ,(40,320,40))
    arrayU = np.reshape(arrayU,(40,320,40))
    arrayV = np.reshape(arrayV,(40,320,40))
    arrayW = np.reshape(arrayW,(40,320,40))
    
    return (arrayU, arrayV, arrayW)

def PlotCOMSOLVelData(meshName,velString,data,Vmin,Vmax):
    global FIGNUM, PLOTNUM
    #View Data with imshow
    #View U_? data
    
    #XNorm
    PLOTNUM = 0
    pathlib.Path('../Images/COMSOL/'+meshName+'/XNorm/'+velString).mkdir(parents=True, exist_ok=True)
    cwd_XNORM = cwd_PYTHON + '/../Images/COMSOL/'+meshName+'/XNorm/'+velString+'/'
    for i in range(40):
        #Create figure and axis
        figX = plt.figure(num=FIGNUM,figsize=(16,2),dpi=500)
        axX = figX.add_subplot(111)
        axX.set_xlabel('Y-Axis (mm)')
        axX.set_ylabel('Z-Axis (mm)')
        #Obtain only XNorm slice data
        tempData = data[:,:,i]
        imgX = axX.imshow(tempData,interpolation='nearest',cmap='viridis',
                        extent=[-159,479,-39,39],
                        vmin = Vmin, vmax = Vmax,
                        origin='lower',aspect='equal')
        #fig.colorbar(img,ax=ax)
        figX.tight_layout()
        figX.savefig(cwd_XNORM + 'image%04d.png'%PLOTNUM)
        FIGNUM += 1
        figX.clf()
        plt.close()
        PLOTNUM += 1
    '''#Make the movie!
    os.chdir(cwd_XNORM)
    os.system("ffmpeg -r 10 -i image%04d.png -vcodec mpeg4 -y movie.mp4")
    os.chdir(cwd_PYTHON)
    print('Video Made! XNORM: '+velString)'''
    
    #YNorm
    PLOTNUM = 0
    pathlib.Path('../Images/COMSOL/'+meshName+'/YNorm/'+velString).mkdir(parents=True, exist_ok=True)
    cwd_YNORM = cwd_PYTHON + '/../Images/COMSOL/'+meshName+'/YNorm/'+velString+'/'
    for i in range(320):
        #Create figure
        figY = plt.figure(num=FIGNUM,figsize=(3,3),dpi=500)
        axY = figY.add_subplot(111)
        axY.set_xlabel('X-Axis (mm)')
        axY.set_ylabel('Z-Axis (mm)')
        #Obtain only Ynorm slice data
        tempData = data[:,i,:]
        imgY = axY.imshow(tempData,interpolation='nearest',cmap='viridis',
                        extent=[-39,39,-39,39],
                        vmin = Vmin, vmax = Vmax,
                        origin='lower',aspect='equal')
        #fig.colorbar(img,ax=ax)
        figY.tight_layout()
        figY.savefig(cwd_YNORM + 'image%04d.png'%PLOTNUM)
        FIGNUM += 1
        figY.clf()
        plt.close()
        PLOTNUM += 1
    '''#Make the movie!
    os.chdir(cwd_YNORM)
    os.system("ffmpeg -r 30 -i image%04d.png -vcodec mpeg4 -y movie.mp4")
    os.chdir(cwd_PYTHON)
    print('Video Made! YNORM: '+velString)'''
    
    #ZNorm
    PLOTNUM = 0
    pathlib.Path('../Images/COMSOL/'+meshName+'/ZNorm/'+velString).mkdir(parents=True, exist_ok=True)
    cwd_ZNORM = cwd_PYTHON + '/../Images/COMSOL/'+meshName+'/ZNorm/'+velString+'/'
    for i in range(40):
        #Create figure
        figZ = plt.figure(num=FIGNUM,figsize=(16,2),dpi=500)
        axZ = figZ.add_subplot(111)
        axZ.set_xlabel('Y-Axis (mm)')
        axZ.set_ylabel('X-Axis (mm)')
        #Obtain only ZNorm slice data
        tempData = data[i,:,:]
        tempData = tempData.T
        imgZ = axZ.imshow(tempData,interpolation='nearest',cmap='viridis',
                        extent=[-159,479,-39,39],
                        vmin = Vmin, vmax = Vmax,
                        origin='lower',aspect='equal')
        #fig.colorbar(img,ax=ax)
        figZ.tight_layout()
        figZ.savefig(cwd_ZNORM + 'image%04d.png'%PLOTNUM)
        FIGNUM += 1
        figZ.clf()
        plt.close()
        PLOTNUM += 1
    '''#Make the movie!
    os.chdir(cwd_ZNORM)
    os.system("ffmpeg -r 10 -i image%04d.png -vcodec mpeg4 -y movie.mp4")
    os.chdir(cwd_PYTHON)
    print('Video Made! ZNORM: '+velString)'''
        
    return

def GenerateCOMSOLMovies(meshName,velString):
    #Create all .mp4 movies for the COMSOL data
    cwd_XNORM = cwd_PYTHON + '/../Images/COMSOL/'+meshName+'/XNorm/'+velString+'/'
    #Make the movie!
    os.chdir(cwd_XNORM)
    os.system("ffmpeg -r 10 -i image%04d.png -vcodec libx264 -pix_fmt yuv420p -y movie1.mp4")
    os.chdir(cwd_PYTHON)
    print('Video Made! XNORM: '+velString)
    cwd_YNORM = cwd_PYTHON + '/../Images/COMSOL/'+meshName+'/YNorm/'+velString+'/'
    os.chdir(cwd_YNORM)
    os.system("ffmpeg -r 30 -i image%04d.png -vcodec libx264 -pix_fmt yuv420p -y movie1.mp4")
    os.chdir(cwd_PYTHON)
    print('Video Made! YNORM: '+velString)
    cwd_ZNORM = cwd_PYTHON + '/../Images/COMSOL'+meshName+'/ZNorm/'+velString+'/'
    os.chdir(cwd_ZNORM)
    os.system("ffmpeg -r 10 -i image%04d.png -vcodec libx264 -pix_fmt yuv420p -y movie1.mp4")
    os.chdir(cwd_PYTHON)
    print('Video Made! ZNORM: '+velString)
    return

def getPIVInputData(meshName,sliceLoc):
    #Obtain directory of data for timeValue
    cwd_DATA = cwd_PYTHON + '/../PIV/'+sliceLoc+'/'+meshName+'/'
    #Create an array for files U_x, U_y, U_z
    #arrays will be of size (40,320,40) where it is the data of coords (x,y,z)
    allData = np.loadtxt(cwd_DATA+'PIV_'+meshName+'_'+sliceLoc+'.txt')
    #Split array up into columns
    arrayY = allData[:,0]
    arrayZ = allData[:,1]
    arrayV = allData[:,2]
    arrayW = allData[:,3]
    #Reshape arrays to be of 3d grid structure
    #Reshaped to where changes are (z,y,x) values
    arrayY = np.reshape(arrayY,(47,64))
    arrayZ = np.reshape(arrayZ,(47,64))
    arrayV = np.reshape(arrayV,(47,64))
    arrayW = np.reshape(arrayW,(47,64))
    
    return (arrayV, arrayW)

def PlotPIVVelData(meshName,sliceLoc,velString,data,Vmin,Vmax):
    global FIGNUM, PLOTNUM
    #View Data with imshow
    #View U_? data
    
    #XNorm
    PLOTNUM = 0
    pathlib.Path('../Images/PIV/'+meshName+'/XNorm/'+velString).mkdir(parents=True, exist_ok=True)
    cwd_XNORM = cwd_PYTHON + '/../Images/PIV/'+meshName+'/XNorm/'+velString+'/'
    #Create figure and axis
    figX = plt.figure(num=FIGNUM,figsize=(6,4),dpi=500)
    axX = figX.add_subplot(111)
    axX.set_xlabel('Y-Axis (mm)')
    axX.set_ylabel('Z-Axis (mm)')
    #Obtain only XNorm slice data
    imgX = axX.imshow(data,interpolation='nearest',cmap='viridis',
                    extent=[-55.3845,55.3845,3.551119316,74.43389895],
                    vmin = Vmin, vmax = Vmax,
                    origin='upper',aspect='equal')
    #fig.colorbar(img,ax=ax)
    figX.tight_layout()
    figX.savefig(cwd_XNORM + 'PIV'+meshName+'_'+velString+'_'+sliceLoc+'.png')
    FIGNUM += 1
    figX.clf()
    plt.close()
    PLOTNUM += 1
        
    return

def PlotSliceVelData(meshName,sliceLoc,velString,PIVdata,COMSOLdata,Vmin,Vmax):
    global FIGNUM, PLOTNUM
    #View Data with imshow
    #View U_? data
    
    '''#XNorm
    PLOTNUM = 0
    pathlib.Path('../Images/Slice/'+meshName+'/XNorm/'+velString).mkdir(parents=True, exist_ok=True)
    cwd_SLICE = cwd_PYTHON + '/../Images/Slice/'+meshName+'/XNorm/'+velString+'/'
    #Create figure and axis
    figX = plt.figure(num=FIGNUM,figsize=(7,4),dpi=500)
    axX = figX.add_subplot(111)
    axX.set_title(COMSOLorPIV+': '+meshName)
    axX.set_xlabel('Y-Axis (mm)')
    axX.set_ylabel('Z-Axis (mm)')
    #Obtain only XNorm slice data
    if(COMSOLorPIV == 'COMSOL'):
        imgX = axX.imshow(data,interpolation='nearest',cmap='viridis',
                        extent=[-55,55,0.0,78],
                        vmin = Vmin, vmax = Vmax,
                        origin='lower',aspect='equal')
    else:
        imgX = axX.imshow(data,interpolation='nearest',cmap='viridis',
                          extent=[-55.38455,55.38455,3.551119316,74.43389895],
                          vmin = Vmin, vmax = Vmax,
                          origin='upper',aspect='equal')
    #fig.colorbar(img,ax=ax)
    figX.tight_layout()
    figX.savefig(cwd_SLICE + COMSOLorPIV+meshName+'_'+velString+'_'+sliceLoc+'.png')
    FIGNUM += 1
    figX.clf()
    plt.close()
    PLOTNUM += 1'''
    
    #XNorm
    PLOTNUM = 0
    #Not Normalized
    pathlib.Path('../Images/Slice/Raw/'+velString).mkdir(parents=True, exist_ok=True)
    cwd_SLICE = cwd_PYTHON + '/../Images/Slice/Raw/'+velString+'/'
    #Normalized
    #pathlib.Path('../Images/Slice/Normalized/'+velString).mkdir(parents=True, exist_ok=True)
    #cwd_SLICE = cwd_PYTHON + '/../Images/Slice/Normalized/'+velString+'/'
    #Create figure and axis
    fig = plt.figure(num=FIGNUM,figsize=(14,4),dpi=500)
    csfont = {'fontname':'Times New Roman'}
    axP = fig.add_subplot(121)
    axP.set_title('PIV: '+meshName+': '+sliceLoc+' pegs',fontsize=16,**csfont)
    axP.set_xlabel('Y-Axis (mm)',fontsize=14,**csfont)
    axP.set_ylabel('Z-Axis (mm)',fontsize=14,**csfont)
    axC = fig.add_subplot(122)
    axC.set_title('COMSOL: '+meshName+': '+sliceLoc+' pegs',fontsize=16,**csfont)
    axC.set_ylabel('Z-Axis (mm)',fontsize=14,**csfont)
    axC.set_xlabel('Y-Axis (mm)',fontsize=14,**csfont)
    #axB = fig.add_subplot(1,4,2)
    #Obtain only XNorm slice data
    imgC = axC.imshow(COMSOLdata,interpolation='nearest',cmap='viridis',
                    extent=[-55,55,0.0,78],
                    vmin = Vmin, vmax = Vmax,
                    origin='lower',aspect='equal')
    imgP = axP.imshow(PIVdata,interpolation='nearest',cmap='viridis',
                      extent=[-55.38455,55.38455,3.551119316,74.43389895],
                      vmin = Vmin, vmax = Vmax,
                      origin='upper',aspect='equal')
    axinsC = inset_axes(axC,
                   width="5%",  # width = 5% of parent_bbox width
                   height="50%",  # height : 50%
                   loc='center left',
                   bbox_to_anchor=(1.05, 0., 1, 1),
                   bbox_transform=axC.transAxes,
                   borderpad=0,
                   )
    axinsP = inset_axes(axP,
                   width="5%",  # width = 5% of parent_bbox width
                   height="50%",  # height : 50%
                   loc='center left',
                   bbox_to_anchor=(1.05, 0., 1, 1),
                   bbox_transform=axP.transAxes,
                   borderpad=0,
                   )
    cbarP = fig.colorbar(imgP,cax=axinsP)
    cbarP.set_label('V (m/s)')
    cbar = fig.colorbar(imgC,cax=axinsC)
    cbar.set_label('V (m/s)')
    fig.tight_layout()
    fig.savefig(cwd_SLICE + meshName+'_'+velString+'_'+sliceLoc+'.png')
    FIGNUM += 1
    fig.clf()
    plt.close()
    PLOTNUM += 1
        
    return

def PlotLineoutData(meshName,sliceLoc,PIVdata,COMSOLdata):
    global FIGNUM
    #Plot velocity data in the flow direction along a vertical line in the center
    #of the mesh
    #Not Normalized
    pathlib.Path('../Images/Lineout/Raw').mkdir(parents=True, exist_ok=True)
    cwd_LINE = cwd_PYTHON + '/../Images/Lineout/Raw/'
    #Normalized
    #pathlib.Path('../Images/Lineout/Normalized').mkdir(parents=True, exist_ok=True)
    #cwd_LINE = cwd_PYTHON + '/../Images/Lineout/Normalized/'
    #Create figure and axis
    fig = plt.figure(num=FIGNUM,figsize=(7,4),dpi=500)
    ax = fig.add_subplot(111)
    csfont = {'fontname':'Times New Roman'}
    ax.set_title(meshName+': '+sliceLoc+' pegs',**csfont,fontsize=16)
    ax.set_xlabel('Height (mm)',fontsize=14,**csfont)
    ax.set_ylabel('V (m/s)',fontsize=14,**csfont)
    COMSOLz = np.arange(0,80,2)
    PIVz = np.linspace(74.43389895,3.551119316,47) #In reverse order
    #Create databases
    COMSOLdict = {'z':COMSOLz,'V':COMSOLdata}
    COMSOLdb = pd.DataFrame(data=COMSOLdict)
    PIVdict = {'z':PIVz,'V':PIVdata}
    PIVdb = pd.DataFrame(data=PIVdict)
    #Remove NaN from the data
    COMSOLdb = COMSOLdb.dropna()
    #print(COMSOLdb)
    PIVdb = PIVdb.dropna()
    #print('\n')
    #print(PIVdb)
    #Plot PIV and COMSOL vel data
    ax.plot([0.0,80],[0.0,0.0],color='k',linewidth=2)
    ax.plot(COMSOLdb.z,COMSOLdb.V,color='r',linewidth=2,label='COMSOL')
    ax.scatter(COMSOLdb.z,COMSOLdb.V,color='r',s=9,label='')
    ax.plot(PIVdb.z,PIVdb.V,color='b',linewidth=2,label='PIV')
    ax.scatter(PIVdb.z,PIVdb.V,color='b',s=9,label='')
    #Vmin = np.amin(COMSOLdata)
    Vmax = np.nanmax(COMSOLdata) + 0.01
    #Not Normalized
    ax.axis([0,80,-0.01,0.075])
    #Normalized
    #ax.axis([0,80,0.0,1.0])
    SetAxesParameters(ax,0.0,80.0,10.0,5.0,0.1)
    ax.legend(loc='best',fontsize='small')
    fig.tight_layout()
    fig.savefig(cwd_LINE + meshName+'_'+sliceLoc+'.png')
    print('Lineout Plot Created: '+meshName+': '+sliceLoc)
    FIGNUM += 1
    fig.clf()
    plt.close()
    
    return

def PlotVelSpectrum(meshName,sliceLoc,PIVdata,COMSOLdata):
    global FIGNUM
    #Plot a density histogram of the velocity in the flow direction for COMSOL and PIV
    #Not Normnalized
    pathlib.Path('../Images/VelDist/Raw').mkdir(parents=True, exist_ok=True)
    cwd_HIST = cwd_PYTHON + '/../Images/VelDist/Raw/'
    #Normalized
    #pathlib.Path('../Images/VelDist/Normalized').mkdir(parents=True, exist_ok=True)
    #cwd_HIST = cwd_PYTHON + '/../Images/VelDist/Normalized/'
    #Create figure and axis
    fig = plt.figure(num=FIGNUM,figsize=(8,4),dpi=500)
    axP = fig.add_subplot(121)
    axC = fig.add_subplot(122)
    csfont = {'fontname':'Times New Roman'}
    axP.set_title('PIV: '+meshName+': '+sliceLoc+' pegs',**csfont,fontsize=16)
    axP.set_xlabel('V (m/s)',fontsize=14,**csfont)
    axP.set_ylabel('Frequency %',fontsize=14,**csfont)
    axC.set_title('COMSOL: '+meshName+': '+sliceLoc+' pegs',**csfont,fontsize=16)
    axC.set_xlabel('V (m/s)',fontsize=14,**csfont)
    axC.set_ylabel('Frequency %',fontsize=14,**csfont)
    #Create databases
    COMSOLdict = {'V':COMSOLdata}
    COMSOLdb = pd.DataFrame(data=COMSOLdict)
    PIVdict = {'V':PIVdata}
    PIVdb = pd.DataFrame(data=PIVdict)
    #Remove/Fill NaN from the data
    COMSOLdb = COMSOLdb.fillna(0.0)
    PIVdb = PIVdb.fillna(0.0)
    #Plot PIV and COMSOL vel data
    weights = np.ones_like(PIVdb.V)/float(len(PIVdb.V))       
    #axP.hist(PIVdb.V,bins=50,range=(-0.01,0.15),weights=weights,color='b')
    axP.hist(PIVdb.V,bins=50,weights=weights,color='b')
    weights = np.ones_like(COMSOLdb.V)/float(len(COMSOLdb.V)) 
    #axC.hist(COMSOLdb.V,bins=50,range=(-0.01,0.15),weights=weights,color='r')
    axC.hist(COMSOLdb.V,bins=50,weights=weights,color='r')
    #Axes Parameters
    axP.tick_params(which='major',axis='both',direction='in',length=6,width=1)
    axP.tick_params(which='minor',axis='both',direction='in',length=4,width=0.75)
    axC.tick_params(which='major',axis='both',direction='in',length=6,width=1)
    axC.tick_params(which='minor',axis='both',direction='in',length=4,width=0.75)
    #Not Normalized
    axP.axis([-0.01,0.075,0.0,0.35])
    axC.axis([-0.01,0.075,0.0,0.35])
    #Normalized
    #axP.axis([0.0,1.0,0.0,0.25])
    #axC.axis([0.0,1.0,0.0,0.25])
    fig.tight_layout()
    fig.savefig(cwd_HIST + meshName+'_'+sliceLoc+'.png')
    print('Histogram Plot Created: '+meshName+': '+sliceLoc)
    FIGNUM += 1
    fig.clf()
    plt.close()
    
    return  

def SetAxesParameters(ax,xMin,xMax,xStep,xMinor,yMinor):
    #Axes Parameters
    ax.tick_params(which='major',axis='both',direction='in',length=6,width=1)
    ax.tick_params(which='minor',axis='both',direction='in',length=4,width=0.75)
    ax.set_axisbelow(False)
    ax.set_xticks(np.arange(xMin,xMax+0.1,step=xStep))
    ax.xaxis.set_minor_locator(MultipleLocator(xMinor))
    ax.yaxis.set_minor_locator(MultipleLocator(yMinor))

def CompareCOMSOLandPIV(meshName,sliceLoc,PIVdata, COMSOLdata):
    #Break up the PIV and COMSOL data into its components
    PIVdataV, PIVdataW = PIVdata[0], PIVdata[1]
    COMSOLdataV, COMSOLdataW, COMSOLdataU = COMSOLdata[0], COMSOLdata[1], COMSOLdata[2]
    #Calculate velocity magnitudes (squared)
    PIVdataMag2 = np.multiply(PIVdataV,PIVdataV) + np.multiply(PIVdataW,PIVdataW)
    COMSOLdataMag2 = (np.multiply(COMSOLdataU,COMSOLdataU) + np.multiply(COMSOLdataV,COMSOLdataV) 
                     + np.multiply(COMSOLdataW,COMSOLdataW))
    #Subset data to Lineout data around y = -1 mm
    #Only do it for the flow direction
    loutPIV = PIVdataV[:,31]
    loutCOMSOL = COMSOLdataV[:,27]
    print(loutCOMSOL.shape)
    PlotLineoutData(meshName,sliceLoc,loutPIV, loutCOMSOL)
    #Using all 2D X-Norm slice data, obtain velocity spectrums for PIV and COMSOL
    histPIV = PIVdataV.reshape(47*64)
    histCOMSOL = COMSOLdataV.reshape(40*56)
    PlotVelSpectrum(meshName,sliceLoc,histPIV,histCOMSOL)
    
    return

if __name__ == '__main__':
    count = 0
    for meshName in meshNameList:
        #COMSOL
        arrayU, arrayV, arrayW = getCOMSOLInputData(meshName)
        #Normalize COMSOL data
        #arrayU = arrayU/(np.nanmax(arrayU))# - np.nanmin(arrayU))
        #rrayV = arrayV/(np.nanmax(arrayV))# - np.nanmin(arrayV))
        #arrayW = arrayW/(np.nanmax(arrayW))# - np.nanmin(arrayW))
        COMSOLvelData = [arrayV, arrayW, arrayU]
        #Get Slice Data (btw, on) for COMSOL
        #Also subset to just be 73.846% of mesh (52:108) in y-dir
        btwCOMSOLvelData = [arrayV[:,52:108,btwIdxList[count]],arrayW[:,52:108,btwIdxList[count]],arrayU[:,52:108,btwIdxList[count]]]
        onCOMSOLvelData = [arrayV[:,52:108,onIdxList[count]],arrayW[:,52:108,onIdxList[count]],arrayU[:,52:108,onIdxList[count]]]
        count += 1
        #PIV
        for sliceLoc in sliceList:
            arrayV, arrayW = getPIVInputData(meshName,sliceLoc)
            #Normalize PIV data
            #arrayV = arrayV/(np.nanmax(arrayV))# - np.nanmin(arrayV))
            #arrayW = arrayW/(np.nanmax(arrayW))# - np.nanmin(arrayW))
            PIVvelData = [arrayV, arrayW]
            #Plot Both COMSOL and PIV data slices
            for i in range(2):
                velString = velStringList[i]
                Vmin = np.nanmin(PIVvelData[i])
                Vmax = np.nanmax(PIVvelData[i])
                #PlotSliceVelData(meshName,sliceLoc,velString,PIVvelData[i],Vmin,Vmax,'PIV')
                if(sliceLoc == 'btw'):
                    #Vmin = np.nanmin(btwCOMSOLvelData[i])
                    #Vmax = np.nanmax(btwCOMSOLvelData[i])
                    print(meshName+', '+sliceLoc+ ': Vmin = %.3e: Vmax =  Vmin, Vmax = %.3e'%(Vmin,Vmax))
                    PlotSliceVelData(meshName,sliceLoc,velString,PIVvelData[i],btwCOMSOLvelData[i],Vmin,Vmax)
                elif(sliceLoc == 'on'):
                    Vmin = np.nanmin(onCOMSOLvelData[i])
                    Vmax = np.nanmax(onCOMSOLvelData[i])
                    PlotSliceVelData(meshName,sliceLoc,velString,PIVvelData[i],onCOMSOLvelData[i],Vmin,Vmax)
            #Identify Slice in COMSOL data that corresponds to PIV
            if(sliceLoc == 'btw'):
                #Comparing Velocity data in between pegs
                #8x15: pegs 3-4
                #10x20: pegs 4-5
                #15x30: pegs 7-8
                CompareCOMSOLandPIV(meshName,sliceLoc,PIVvelData,btwCOMSOLvelData)
            elif(sliceLoc == 'on'):
                #Comparing Velocity data on a central peg
                #8x15: peg 4
                #10x20: peg 5
                #15x30: peg 8
                CompareCOMSOLandPIV(meshName,sliceLoc,PIVvelData,onCOMSOLvelData)
        
        
    '''COMSOL Images/Movies
    for meshName in meshNameList:
        if(boolMovie == 1):
            break
        arrayU, arrayV, arrayW = getCOMSOLInputData(meshName)
        velData = [arrayU, arrayV, arrayW]
        velString = ['U','V','W']
        for i in range(len(velData)):
            Vmin = np.nanmin(velData[i])
            Vmax = np.nanmax(velData[i])
            print(Vmin, Vmax)
            PlotCOMSOLVelData(meshName,velString[i],velData[i],Vmin,Vmax)
            
    velString = ['U','V','W']
    for meshName in meshNameList:
        for i in range(3):
            GenerateCOMSOLMovies(meshName,velString[i])'''
    
    '''PIV Images'''
    '''#PIV data and stuff goes here
    for meshName in meshNameList:
        for sliceLoc in sliceList: 
            arrayV, arrayW = getPIVInputData(meshName,sliceLoc)
            velData = [arrayV, arrayW]
            velString = ['V','W']
            for i in range(len(velData)):
                Vmin = np.nanmin(velData[i])
                Vmax = np.nanmax(velData[i])
                print(Vmin,Vmax)
                PlotPIVVelData(meshName,sliceLoc,velString[i],velData[i],Vmin,Vmax)'''
    
    '''COMPARE'''
    #Compare COMSOL data to PIV data
