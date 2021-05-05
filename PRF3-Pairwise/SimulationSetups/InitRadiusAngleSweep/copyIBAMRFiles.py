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
#RList = ['5','7','10']
RList = ['3','5','6','7']
configList = ['Parallel','Anti','PerpL','PerpS']
moveDirList = ['SSL','LSL','Stat']

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
    for R in RList:
        dirNames = [ name for name in os.listdir('./'+R+'/') if os.path.isdir(os.path.join('./'+R+'/', name)) ]
        #dirNames = os.listdir('./'+R+'/')
        print(dirNames)
        for dirName in dirNames:
            configNames = [ name for name in os.listdir('./'+R+'/'+dirName+'/') if os.path.isdir(os.path.join('./'+R+'/'+dirName+'/', name)) ]
            print(configNames)
            for config in configNames:
                for move in moveDirList:
                    pathlib.Path(R+'/'+dirName+'/'+config+'/'+move+'/').mkdir(parents=True, exist_ok=True)
                    cwd_COPY = cwd_PYTHON + '/'+R+'/'+dirName+'/'+config+'/'+move+'/'
                    #Sim FILES
                    copyfile('main.C',cwd_COPY+'main.C')
                    copyfile('Makefile',cwd_COPY+'Makefile')
                    copyfile('update_springs.C',cwd_COPY+'update_springs.C')
                    copyfile('update_springs.h',cwd_COPY+'update_springs.h')
                    if(move == 'SSL'):
                        copyfile('input2dSSL',cwd_COPY+'input2dSSL')
                        #os.system("sed 's/ANGLE/'$visc'/g;s/AMP_VALUE/'$amp'/g;s/RSL_VALUE/'$rsl_value'/g' spherobot2D.bsub > ../Structures/v$visc/amp$amp/1botCIB.bsub")
                        #copyfile('scriptSSL.sh',cwd_COPY+'scriptSSL.sh')
                    elif(move == 'LSL'):
                        copyfile('input2dLSL',cwd_COPY+'input2dLSL')
                        #copyfile('scriptLSL.sh',cwd_COPY+'scriptLSL.sh')
                    else:
                        copyfile('input2dStat',cwd_COPY+'input2dStat')
                        #copyfile('scriptStat.sh',cwd_COPY+'scriptStat.sh')
                    #Mesh FILES
                    cwd_MESH = cwd_PYTHON + '/'+R+'/'+dirName+'/'+config+'/'
                    os.chdir(cwd_MESH)
                    if(move == 'SSL'):
                        os.system("cp bot* skel* SSL/")
                    elif(move == 'LSL'):
                        os.system("cp bot* skel* LSL/")
                    else:
                        os.system("cp bot* skel* Stat/")
                    os.chdir(cwd_PYTHON)
    #zipFiles('../Structures','PairwiseSetup')
    #zipFiles('5/','PairwiseSetup5')
