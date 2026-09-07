#!usr/bin/python 

import subprocess as sp
import os 

def run_mbdyn(model, output, path2mbdyn, coupled):
    oldpwd = os.getcwd()
    try:
        os.chdir(path2mbdyn)  
        if coupled:  
            sp.run('mbdyn ' + model + ' -o ' + output + ' 2>&1 &' , shell=True)
        else: 
            sp.run('mbdyn ' + model + ' -o ' + output , shell=True)
        os.chdir(oldpwd)
        cwd = os.getcwd()
    except:
        pass
