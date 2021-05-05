#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 14:16:33 2018

@author: thomas
"""

def main():
    
    #nvert_bot1up = 4921
    #nvert_bot1low = 1330
    nvert_boundary = 34878
    Ks = 1.0e6
    
    f = open('boundary.target','w')
    f.write("%i\n" %(nvert_boundary))
    for i in range(nvert_boundary):
        f.write("%i %.3e\n" % (i,Ks))
    f.close()
    
    '''f = open('botup1.target','w')
    f.write("%i\n" %(nvert_bot1up))
    for i in range(nvert_bot1up):
        f.write("%i %.3e\n" % (i,Ks))
    f.close()
    
    f = open('botlow1.target','w')
    f.write("%i\n"%(nvert_bot1low))
    for i in range(nvert_bot1low):
        f.write("%i %.3e\n" %(i,Ks))
    f.close'''
    
#--------------------END OF MAIN---------------------
main()
