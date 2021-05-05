#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 17:01:32 2019

@author: thomas
"""

#In this script, we will be using the expansion and compression displacement
# data created by PRQuantitative.py
#We will attempt to find a fit for the data SSL and LSL
#The fit will be based on 3 parameters: Re = A_r*r/delta^2, St_r = r/A, and d_0?

#MODULES
import os,sys
import re
import numpy as np
import pandas as pd
from scipy import stats
from sklearn import linear_model
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import (MultipleLocator,AutoMinorLocator)
from scipy.signal import savgol_filter

#CONSTANTS
cwd_PYTHON = os.getcwd()
PERIOD = 0.1
FREQ = 10.0
DENS = 2.0
MAXAMP = 0.3
RSMALL = 0.15
RLARGE = 0.3
RSL_DEFAULT = 0.75

rslList = [1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.8,2.0]
WholeExpFit = np.zeros((2,5))
WholeComFit = np.zeros((2,5))

def StoreData(cwd_DATA):
    global RSL_DEFAULT
    #Load position data
    #Columns
    #Re rsl amp y_e y_c y_net
    data = pd.read_csv(cwd_DATA+'/Shifted/ExpCom.txt',delimiter = ' ')
    data = data[data.Re <= 150.0].copy()
    data = data.sort_values(by=['Re','rsl','amp'])
    data = data.reset_index(drop=True)
    #Renormalize
    data['y_e'] /= RSMALL
    data['y_c'] /= RSMALL
    #Let's create 5 new variables (log values)
    data['logRe'] = np.log10(data.Re)
    data['logy_e'] = np.log10(abs(data['y_e']))
    data['logy_c'] = np.log10(abs(data['y_c']))
    data['epsilon'] = data.amp/RSMALL
    data['logeps'] = np.log10(data['epsilon'])
    #data['d_0'] = RSL_DEFAULT*data.rsl
    data['d_0'] = RSL_DEFAULT*data.rsl/RSMALL
    data['logd_0'] = np.log10(data['d_0'])
    data['Diff'] = abs(data['y_e'] + data['y_c'])
    #Attempt2
    data['M2'] = data['Re']/data['epsilon']
    data['logM2'] = np.log10(data['M2'])
    PlotAllStrokes(data)
    data.to_csv(cwd_DATA + '/v2_ExpComAllVar.txt',index=None,header=True)
    #sys.exit(0)
    return data

def PlotAllStrokes(data):
    ampList = [0.12,0.15,0.18,0.21,0.24]
    rslList = [1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.8,2.0]
    #Color Palette for plotting
    R1=[255/255,255/255,153/255,153/255,204/255]
    G1=[153/255,204/255,255/255,204/255,153/255]
    B1=[204/255,153/255,153/255,255/255,255/255]
    
    ShadeValue = 1 #used to scale line color so there is a gradient as rsl changes

    #Create Figure for each amp value
    fig = plt.figure(figsize=(9,4),dpi=200,num=100)
    axPos = fig.add_subplot(121)
    axNeg = fig.add_subplot(122)
    #Labels and Such
    axPos.set_xlabel('Re',fontsize=14)
    axNeg.set_xlabel('Re',fontsize=14)
    axPos.set_ylabel(r'$\Delta \hat{y}_{exp}$ (m)',fontsize=14)
    axNeg.set_ylabel(r'$\Delta \hat{y}_{com} - \Delta \hat{y}_{min}$ (m)',fontsize=14)

    for idxAmp in range(5):
        ampValue = ampList[idxAmp]
        ampData = data[data.amp == ampValue].copy()
        '''#Create Figure for each amp value
        fig = plt.figure(figsize=(9,4),dpi=200,num=100)
        axPos = fig.add_subplot(121)
        axNeg = fig.add_subplot(122)
        #Labels and Such
        axPos.set_xlabel('Re',fontsize=14)
        axNeg.set_xlabel('Re',fontsize=14)
        axPos.set_ylabel(r'$\Delta y_{exp}$ (m)',fontsize=14)
        axNeg.set_ylabel(r'$\Delta y_{com} - \Delta y_{min}$ (m)',fontsize=14)'''
        #axPos.set_title(r'Expansion and Compression: A%.2fm'%ampValue)
        ShadeValue = 1.0
        for idxRSL in range(len(rslList)):
            rslValue = rslList[idxRSL]
            rslData = ampData[ampData.rsl == rslValue].copy()
            rslData = rslData.sort_values(by=['Re'])
            rslData = rslData.reset_index(drop=True)
            #ShadeValue = 1.0 - 0.05*idxRSL/(len(rslList))
            #Select RGB Color
            R=R1[idxAmp]*ShadeValue
            G=G1[idxAmp]*ShadeValue
            B=B1[idxAmp]*ShadeValue
            
            seriesCMin = rslData[['y_c']].min()
            miny_c = seriesCMin[0]
            rslData['y_c'] -= miny_c
            rslData['logy_c'] = np.log10(rslData['y_c'])
            rslData = rslData.loc[1:]
            
            #Plot Re_c
            axNeg.plot([20.0,20.0],[-10,10],color='gray',ls=':')
            axPos.plot([20.0,20.0],[-10,10],color='gray',ls=':')
            #Expansion
            axPos.plot(rslData['Re'],rslData['y_e']/RSMALL,color=(R,G,B),zorder=5)
            axPos.scatter(rslData['Re'],rslData['y_e']/RSMALL,color=(R,G,B),s=9,zorder=5)
            #Compression
            axNeg.plot([0.05,150.0],[0.0,0.0],color='k')
            axNeg.plot(rslData['Re'],rslData['y_c']/RSMALL,color=(R,G,B),zorder=5)
            axNeg.scatter(rslData['Re'],rslData['y_c']/RSMALL,color=(R,G,B),s=9,zorder=5)
            ShadeValue -=0.05
        
        axPos.set_xlim(0.5,200.0)
        axNeg.set_xlim(0.5,200.0)
        axNeg.set_ylim(1.0e-2,5.0e0)
        axPos.set_ylim(1.0e-2,5.0e0)
        axPos.set_xscale('log')
        axPos.set_yscale('log')
        axNeg.set_xscale('log')
        axNeg.set_yscale('log')
        
        #Axes Parameters
        axPos.tick_params(which='major',axis='both',direction='in',length=6,width=1)
        axPos.tick_params(which='minor',axis='both',direction='in',length=4,width=0.75)
        axNeg.tick_params(which='major',axis='both',direction='in',length=6,width=1)
        axNeg.tick_params(which='minor',axis='both',direction='in',length=4,width=0.75)
    
    fig.tight_layout()
    #figName = '../PR/Shifted/Fits/TestCompression/ExpComStroke_A'+str(ampValue)+'.png'
    figName = cwd_PYTHON+'/../Version2/PaperFigures/Images/ExpComStroke_RawData_All_v2.svg'
    fig.savefig(figName)
    fig.clf()
    
    return

def FindRoot(data):
    #Identify 2 locations
    #1) where v goes from + to -
    #2) where v goes from - to +
    #+ to - v will always happen first
    nRe = len(data)
    idxRe = 1
    root = None
    
    #Find where v goes from - to +
    #Save as root
    while(idxRe < nRe-1):
        if(data[idxRe]*data[idxRe-1] <= 0.0 and 
           data[idxRe]*data[idxRe+1] >= 0.0 and
           data[idxRe] >= 0.0):
            #- to + root has been found
            print('- to +: b4y_c = %.3f\ty_c = %.3f\ta4y_c = %.5f'%(data[idxRe-1],
                                                            data[idxRe],
                                                            data[idxRe+1]))
            #Save + to - index value
            root = idxRe
        idxRe += 1
        
    return root

def FindWholeFit(data,y_name): 
    #In this function, we will use sklearn.linear_model to find
    # a mutivariable linear regression of the expansion or compression data
    
    '''#Remove Outliers
    filteredData = RejectOutliers(data,y_name)'''
    #Variables to be fit
    xData = data[['logRe','logeps','logd_0']]
    #Data to be power fit
    yData = data[y_name]    
        
    #Construct linear model for Expansion/Compression
    lm = linear_model.LinearRegression()
    model = lm.fit(xData,yData)
    #Calculate predicted values from model
    predictions = lm.predict(xData)
    #Calculate R^2 value
    score = lm.score(xData,yData)
    #Return Coefficients of linear model
    coef = lm.coef_
    #Return intercept of linear model
    intercept = lm.intercept_
    
    #Print Stats
    print("="*40)
    print('Mult Vairable Whole Model:')
    print(y_name)
    print('M2 exponent = %.3f'%coef[0])
    print('eps exponent = %.3f'%coef[1])
    print('d0 exponent = %.3f'%coef[2])
    print('intercept = %.3f'%intercept)
    print('R^2 value = %.3f'%score)

    return (coef[0],coef[1],coef[2],intercept, score)

def PlotWholeFit(ax,data,fitList,var,fitColor):
    #Color Palette for plotting
    R1=[255/255,255/255,153/255,153/255,204/255]
    G1=[153/255,204/255,255/255,204/255,153/255]
    B1=[204/255,153/255,153/255,255/255,255/255]
    ShadeValue = 1.0 #used to scale line color so there is a gradient as rsl changes
    
    #Calculate Fit Expression
    powRe, poweps, powd_0, intercept = fitList[0], fitList[1], fitList[2], fitList[3]
    xFit = np.linspace(0.5,250.0,1000)
    epsilon = 0.18/RSMALL
    #d_0 = rsl*RSL_DEFAULT
    d_0 = 1.0*RSL_DEFAULT/RSMALL
    yFit = 10.0**(intercept)*xFit**(powRe)*(epsilon)**(poweps)*(d_0)**(powd_0)
    
    #Plot Re_c line
    ax.plot([20.0,20.0],[-1,1],color='gray',ls=':')
    #Plot Fit
    ax.plot(xFit,yFit/((epsilon)**(poweps)*(d_0)**(powd_0)),color=fitColor,lw=2,ls='--')
    
    for rsl in rslList:
        rslData = data[data.rsl == rsl].copy()
        rslData = rslData.reset_index(drop=True)
    
        for idxAmp in range(0,5):
            ampValue = 0.12 + 0.03*idxAmp
            ampData = rslData[rslData.amp == ampValue].copy()
        
            #Select RGB Color
            R=R1[idxAmp]*ShadeValue
            G=G1[idxAmp]*ShadeValue
            B=B1[idxAmp]*ShadeValue

            #ax.plot(xFit,yFit/((epsilon)**(poweps)),color=(R,G,B),lw=2,ls='--')
            #Plot Raw Data
            ax.scatter(ampData['Re'],abs(ampData[var])/(ampData['epsilon']**(poweps)*ampData['d_0']**(powd_0)),color=(R,G,B),s=20,zorder=5)#,edgecolor='k',zorder=5,linewidth=0.5)
            #ax.scatter(ampData['Re'],abs(ampData[var])/(ampData['epsilon']**(poweps)),color=(R,G,B),s=20,zorder=5)#,edgecolor='k',zorder=5,linewidth=0.5)
        ShadeValue -= 0.05
    '''if(name=='exp'):
        ax.legend(title=r'%.3flog(Re)+%.3flog(St)+%.3flog($d_0$)+%.3f'%(powRe,powSt,powd_0,intercept),loc='best',fontsize='x-small')
    else:
        ax.legend(title=r'%.3flog(Re)+%.3flog(St)+%.3flog($d_0$)+%.3f'%(powCRe,powCSt,powCd_0,interceptC),loc='best',fontsize='x-small')'''
    
    #Axes Parameters
    ax.tick_params(which='major',axis='both',direction='in',length=6,width=1)
    ax.tick_params(which='minor',axis='both',direction='in',length=4,width=0.75)
    
    return ax

def PlotAllCompression(ax,data):
    #Color Palette for plotting
    R1=[255/255,255/255,153/255,153/255,204/255]
    G1=[153/255,204/255,255/255,204/255,153/255]
    B1=[204/255,153/255,153/255,255/255,255/255]
    ShadeValue = 1.0 #used to scale line color so there is a gradient as rsl changes
    for rsl in rslList:
        rslData = data[data.rsl == rsl].copy()
        rslData = rslData.reset_index(drop=True)
    
        for idxAmp in range(0,5):
            ampValue = 0.12 + 0.03*idxAmp
            ampData = rslData[rslData.amp == ampValue].copy()
        
            #Select RGB Color
            R=R1[idxAmp]*ShadeValue
            G=G1[idxAmp]*ShadeValue
            B=B1[idxAmp]*ShadeValue
        
            #Plot Re_c line
            ax.plot([20.0,20.0],[-1,1],color='gray',ls=':')
            #Plot Raw Data
            ax.plot(ampData['Re'],abs(ampData['y_c']),color=(R,G,B),lw=1)
            ax.scatter(ampData['Re'],abs(ampData['y_c']),color=(R,G,B),s=20,zorder=5)
        ShadeValue -= 0.05
    '''if(name=='exp'):
        ax.legend(title=r'%.3flog(Re)+%.3flog(St)+%.3flog($d_0$)+%.3f'%(powRe,powSt,powd_0,intercept),loc='best',fontsize='x-small')
    else:
        ax.legend(title=r'%.3flog(Re)+%.3flog(St)+%.3flog($d_0$)+%.3f'%(powCRe,powCSt,powCd_0,interceptC),loc='best',fontsize='x-small')'''
    
    #Axes Parameters
    ax.tick_params(which='major',axis='both',direction='in',length=6,width=1)
    ax.tick_params(which='minor',axis='both',direction='in',length=4,width=0.75)
    
    return ax

if __name__ == '__main__':
    #This is where the main part of the code will be conducted
    #Obtain the directory where data is stored
    cwd_DATA = cwd_PYTHON + '/../PR'
    allData = StoreData(cwd_DATA)
    csfont = {'fontname':'Times New Roman'}
    #Here we need to split up the data based on the following criteria
    #1) Expansion: Split based on minimum value
    # Do not include minimum and value to the left or right
    #2) Compression: Split based on when sign switches (root)
    # Do not include root and value to the left or right
    #But first we need to find the root for each set of (A,d_0) data
    #Create DataFrames which will store the 2 regions of info
    vEData = pd.read_csv(cwd_DATA+'/v2_ExpComAllVar.txt',delimiter = ' ',nrows=0)
    iEData = vEData.copy()
    vCData = iEData.copy()
    iCData = vCData.copy()
    #Find root and split data up based on rsl and A.
    #Loop over RSL
    for idxRSL in range(len(rslList)):
        rslValue = rslList[idxRSL]
        rslData = allData[allData.rsl == rslValue].copy()
        #Loop over A
        for idxAmp in range(0,5):
            ampValue = 0.12 + 0.03*idxAmp
            print('A = %.2f: RSL = %.2f'%(ampValue,rslValue))
            ampData = rslData[rslData.amp == ampValue].copy()
            #Sort based on Re
            ampData = ampData.sort_values(by=['Re'])
            ampData = ampData.reset_index(drop=True)
            print(ampData['y_e'])
            #Now we should have data for 20 Re
            print('# of Re = ',len(ampData.Re))
            #Expansion Data
            #Split data up based on minimum value
            seriesMin = ampData[['y_e']].idxmin()
            idxMin = seriesMin[0]
            print('idxMin = ',idxMin)
            vEData = pd.concat([vEData,ampData.iloc[:idxMin+1]], ignore_index=True, sort=True)
            iEData = pd.concat([iEData,ampData.iloc[idxMin+1:]], ignore_index=True, sort=True)
            #Compression Data
            seriesCMin = ampData[['y_c']].min()
            miny_c = seriesCMin[0]
            ampData['y_c'] -= miny_c
            ampData['logy_c'] = np.log10(ampData['y_c'])
            ampData = ampData.loc[1:]
            vCData = pd.concat([vCData,ampData.iloc[:idxMin+1]], ignore_index=True,sort=True)
            iCData = pd.concat([iCData,ampData.iloc[idxMin+1:]], ignore_index=True,sort=True)

    
    vEData = vEData.sort_values(by=['Re','amp','rsl'])
    vEData = vEData.reset_index(drop=True)
    vCData = vCData.sort_values(by=['Re','amp','rsl'])
    vCData = vCData.reset_index(drop=True)
    iEData = iEData.sort_values(by=['Re','amp','rsl'])
    iEData = iEData.reset_index(drop=True)
    iCData = iCData.sort_values(by=['Re','amp','rsl'])
    iCData = iCData.reset_index(drop=True)
    #Data has been split up appropriately
    #For each dataframe, we will perform 4 linear regressions
    #1) logRe only
    #2) logRe and logSt
    #3) logRe and logd_0
    #4) logRe, logSt, and logd_0
    
    #There will be 5 sets of figures (1 for viscous region, 1 for inertial region)
    #1) logRe only (sorted by 9 d_0) 5 amps on each
    #2) logRe only (sorted by 5 A) 9 RSL each
    #3) logRe and logSt (sorted by 9 d_0) 
    #4) logRe and log_d_0 (sorted by 5 A)
    #5) logRe, logSt, and logd_0 (1 Figure of all)
    
    #Make a figure for each Whole Fit
    #First Find fits
    #Variables of interest: Re and St and d_0
    #Plot Fits for Whole Data set, (Re, St, and d_0)
    WholeExpFit[0] = FindWholeFit(vEData,'logy_e')
    WholeComFit[0] = FindWholeFit(vCData[vCData.Re <= 2.0],'logy_c')
    WholeExpFit[1] = FindWholeFit(iEData,'logy_e')
    WholeComFit[1] = FindWholeFit(vCData[vCData.Re > 2.0],'logy_c')
    #Important Figure Data
    cwd_PLOT = cwd_PYTHON + '/../Version2/PaperFigures/Images'
    #Create a Figure for ReStd_0 Fit Data
    fig, ax = plt.subplots(nrows=1, ncols=2, num=1, 
                             figsize=(9,4),dpi=200)
    #Axes Labels for ReSt Fit
    ax[0].set_xlabel(r'Re',fontsize=14,**csfont)
    ax[0].set_ylabel(r'$\Delta \hat{y}_{exp}$ / $\epsilon^a \hat{d}_0^b$',fontsize=14,**csfont)
    ax[1].set_xlabel(r'Re',fontsize=14,**csfont)
    ax[1].set_ylabel(r'($\Delta \hat{y}_{com} - \Delta \hat{y}_{min}$) / $\epsilon^a \hat{d}_0^b$',fontsize=14,**csfont)
    #vEData
    imgName = 'Final'
    ax[0] = PlotWholeFit(ax[0],vEData,WholeExpFit[0],'y_e','k')
    #iEData
    ax[0] = PlotWholeFit(ax[0],iEData,WholeExpFit[1],'y_e','tab:red')
    ax[0].set(xlim=(0.5,200.0),ylim=(5.0e-5/RSMALL,5.0e-2/RSMALL))
    #ax[0].set(xlim=(0.5,200.0),ylim=(1.0e-4,1.0e-1))
    ax[0].set_xscale('log')
    ax[0].set_yscale('log')
    #ax[0].set_aspect(1.0)
    #cData
    #Reg1Test
    ax[1] = PlotWholeFit(ax[1],vCData[vCData.Re <= 2.0],WholeComFit[0],'y_c','k')
    #Reg2Test
    ax[1] = PlotWholeFit(ax[1],vCData[vCData.Re > 2.0],WholeComFit[1],'y_c','tab:blue')
    #Reg3 Test
    WholeComFit[1] = FindWholeFit(iCData,'logy_c')
    ax[1] = PlotWholeFit(ax[1],iCData,WholeComFit[1],'y_c','tab:red')
    #Plot All of Compression Unfit
    #ax[1] = PlotAllCompression(ax[1],vCData)
    ax[1].set(xlim=(0.5,200.0),ylim=(5.0e-5/RSMALL,5.0e-2/RSMALL))
    #ax[1].set(xlim=(0.5,150.0),ylim=(1.0e-4,8.0e-2))
    ax[1].set_xscale('log')
    ax[1].set_yscale('log')
    #ax[1].set_aspect(0.5)
    #ax[1].invert_yaxis()
    fig.tight_layout()
    figName = cwd_PLOT+'/'+imgName+'_WholeFitv2_l.svg'
    fig.savefig(figName)
    fig.clf()
    plt.close()
    
    

            
    
    