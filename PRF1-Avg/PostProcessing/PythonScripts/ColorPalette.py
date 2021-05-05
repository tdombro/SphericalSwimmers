#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 10:43:31 2018

@author: thomas
"""

import numpy as np
import matplotlib.pyplot as plt

def main():
    #Color Palette for plotting
    #Colors are in the following order: Red/Pink, Orange, Green, Blue, Purple
    R1=[255/255,255/255,153/255,153/255,204/255]
    G1=[153/255,204/255,255/255,204/255,153/255]
    B1=[204/255,153/255,153/255,255/255,255/255]
    
    #Sample Data (Just using numpy for ease)
    var1List = [0.2,0.4,0.6]
    var2List = [0.5,1.0,1.5]
    time = np.linspace(0.0,10.0,10)
    velocity = np.zeros((3,3,10))
    k = 1.0
    for i in range(len(var1List)):
        velocity[i,0,:] = np.linspace(0.0,10.0*k,10)
        velocity[i,1,:] = np.linspace(0.0,15.0*k,10)
        velocity[i,2,:] = np.linspace(0.0,17.0*k,10)
        k += 0.2
        print('='*40+'\n')
        print('i = ',i)
        print(velocity[i,:,:])
        print('-'*40+'\n')
        
    fig = plt.figure(num=0, figsize=(4,4),dpi=120)
    ax = fig.add_subplot(111)
    ax.set_title('V vs Time: Color Preview')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Vel (m/s)')

    
    j = 0 #initial color (Red/Pink)
    for k in range(len(var1List)):
        i = 1 #used to scale line color so there is a gradient as rsl changes (Darkness)
        for m in range(len(var2List)):
            R = R1[j]*i
            G = G1[j]*i
            B = B1[j]*i
            #Plot each (var1, var2) line
            ax.plot(time,velocity[k,m,:],
                       label='$var1$='+str(var1List[k])+' $var2$='+str(var2List[m]),color=(R,G,B),linewidth=2)
            ax.scatter(time,velocity[k,m,:],
                          label=None,color=(R,G,B),s=18)
            i -= 0.66/(len(var2List)+1) #Increase Darkness (-0.66) value can be changed
        j += 1 #Increase color index (change color based on var1)
    
    lgd = ax.legend(loc=2, bbox_to_anchor=(1.05,1),borderaxespad=0,ncol=1,fontsize='x-small')
    fig.savefig('ColorPreview.png',bbox_extra_artists=(lgd,),bbox_inches='tight')
    plt.show()
    return
    
#------------------__END MAIN__-----------------------------------
main()