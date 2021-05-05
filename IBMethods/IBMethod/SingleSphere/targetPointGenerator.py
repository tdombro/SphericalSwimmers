#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 14:16:33 2018

@author: thomas
"""

def main():
    
    nvert_cylinder = 1746
    Ks = 6.0e5
    
    f = open('cylinder.target','w')
    f.write("%i\n" %(nvert_cylinder))
    for i in range(nvert_cylinder):
        f.write("%i %.3e\n" % (i,Ks))
    f.close()
    
#--------------------END OF MAIN---------------------
main()
