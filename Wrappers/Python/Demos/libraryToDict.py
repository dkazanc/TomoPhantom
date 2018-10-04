# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 11:12:51 2018

@author: ofn77899
"""
import re
from tomophantom.TomoP2D import Objects2D

def modelfile2Dtolist(filepath, model):
    #read all the file content
    f = open(filepath).readlines()
    
    # model 4 is lines 56
    i = 0
    model = model
    components = 0
#    timesteps = 1
    while (True):
        if f[i] == 'Model : {0:02};\n'.format(model):
            components = int (re.search('^Components : (\d+);$',f[i+1]).groups()[0] )
#            timesteps = int( re.search('^TimeSteps : (\d+);$',f[i+2]).groups()[0] )
#            print (i)
            i += 3
#            print (i)
#            print (f[i])
            break
        i += 1
    #print (i)
    
    objectlist = []
    
    for j in range(components):
        descr = f[i+j].split()
        objstr, C0, x0, y0, a, b, phi = descr[2:]
        for oo in Objects2D:
            if oo.value==objstr:
                break
#        print (phi)
        objectlist.append( {'Obj' : oo , 
                            'C0' : float(C0), 
                            'x0' : float(y0),
                            'y0' : float(x0),
                            'a'  : float(b),
                            'b'  : float(a),
                            'phi': float(phi[:-1])} )
    return objectlist

if __name__ == "__main__":
    
    pathTP = '../../functions/models/Phantom2DLibrary.dat'

    print (modelfile2Dtolist(pathTP , 6))