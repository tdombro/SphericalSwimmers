#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 15:39:29 2019

@author: thomas
"""

import os, shutil
import pathlib
from shutil import copyfile
import zipfile

#Copy IBAMR files into every simulation directory
#1)Loop over R,Theta,Configuration
#2)Create SSL and LSL directories
#3)Copy all files in Structures directory over to SSL/LSL directory
#4)Copy all files in Para/Anti directory over to SSL/LSL directory
#5)Zip Structures folder

#Lists
ThetaList=['0.0','22.5','45.0','67.5','90.0','112.5','135.0','157.5','180.0',
           '202.5','225.0','247.5','270.0','292.5','315.0','337.5']
HxList=['-0.5','0.5','1.5','2.5','3.5','4.5','5.5','6.5','7.5','8.5','9.5','10.5','11.5','12.5']
HyList=['-13','-11','-9','-7','-5','-3','-1','1','3','5','7','9']

cwd_PYTHON = os.getcwd()

def zipFiles(src,dst):
    zf = zipfile.ZipFile('%s.zip' % (dst), 'w', zipfile.ZIP_DEFLATED)
    abs_src = os.path.abspath(src)
    for dirname, subdirs, files in os.walk(src):
        for filename in files:
            absname = os.path.abspath(os.path.join(dirname, filename))
            arcname = absname[len(abs_src) + 1:]
            print('zipping %s as %s' % (os.path.join(dirname, filename),
                                        arcname))
            zf.write(absname, arcname)
    zf.close()
    
if __name__ == '__main__':
    
    for Theta in ThetaList:
        for Hx in HxList:
            dirNames = [ name for name in os.listdir('./Theta'+Theta+'/Hx'+Hx+'/') if os.path.isdir(os.path.join('./Theta'+Theta+'/Hx'+Hx+'/', name)) ]
            print('dirHys = ',dirNames)
            for dirHy in dirNames:
                pathlib.Path('Theta'+Theta+'/Hx'+Hx+'/'+dirHy+'/').mkdir(parents=True, exist_ok=True)
                cwd_COPY = cwd_PYTHON + '/Theta'+Theta+'/Hx'+Hx+'/'+dirHy+'/'
                #Sim FILES
                copyfile('main.C',cwd_COPY+'main.C')
                copyfile('Makefile',cwd_COPY+'Makefile')
                copyfile('update_springs.C',cwd_COPY+'update_springs.C')
                copyfile('update_springs.h',cwd_COPY+'update_springs.h')
                copyfile('input2dSSL',cwd_COPY+'input2D')
                os.chdir(cwd_PYTHON)

