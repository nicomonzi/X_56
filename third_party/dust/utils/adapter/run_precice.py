#!usr/bin/python

import tempfile 
import os, sys
# set to path of MBDyn support for python communication
sys.path.append('/usr/local/mbdyn/libexec/mbpy') 

from mbc_py_interface import mbcNodal

import precice
from precice import *
from mbdynInterface import MBDynInterface
from mbdynAdapter import MBDynAdapter   
import subprocess as sp
import numpy as np


def cleanup(path2mbdyn, path2dust, path2precice): 

    if os.path.isfile(os.path.join(path2mbdyn,'precice-MBDyn-iterations.json')):
        events = os.path.join(path2mbdyn,'precice-MBDyn-events.json')
        iteration = os.path.join(path2mbdyn,'precice-MBDyn-iterations.log') 
        sp.run('rm -fv ' + events,      shell=True)
        sp.run('rm -fv ' + iteration,   shell=True)

    if os.path.isdir(os.path.join(path2mbdyn,'__pycache__')):      
        pycache = os.path.join(path2mbdyn,'__pycache__')
        sp.run('rm -rfv ' + pycache, shell=True)

    #remove temporary folder 
    sp.run('rm -rfv ' + os.path.join(path2mbdyn,'mbdyn','.mbdyn_*'),        shell=True)

    if os.path.isfile(os.path.join(path2dust,'precice-dust-convergence.log')): 
        convergence = os.path.join(path2dust,'precice-dust-convergence.log')
        events = os.path.join(path2dust,'precice-dust-events.json')
        iteration =  os.path.join(path2dust,'precice-dust-iterations.log')

        sp.run('rm -fv ' + convergence, shell=True)
        sp.run('rm -fv ' + events,      shell=True)
        sp.run('rm -fv ' + iteration,   shell=True)
        
    # remove dummy precice-run folder
    if os.path.isdir(os.path.join(path2precice,'precice-run')):
        sp.run('rm -rfv '+ os.path.join(path2precice,'precice-run'), shell=True)
        sp.run('rm -fv ' + os.path.join(path2precice,'nnodes.dat'),  shell=True)


def setup_socket(): 
    tmpdir = tempfile.mkdtemp('', '.mbdyn_')
    path = tmpdir + '/mbdyn.sock'
    print(' path: ', path)
    os.environ['MBSOCK'] = path
    return path

def run_precice(path, dt, path2mbdyn): 

    oldpwd = os.getcwd()
    os.chdir(path2mbdyn)
    #> Does the other solver needs MBDyn nodes to build its own mesh?
    writeNodes = False
    #> Construct MBDyn/mbc_py interface
    #> Initialize MBDyn/mbc_py interface: negotiate and recv()
    mbd = MBDynInterface()
    # get node number 
    rr = mbd.refConfigNodes() 
    nnodes = np.size(rr,0)
    # initialize socket for data exchange
    mbd.initialize( path=path, verbose=1, nnodes=nnodes, accels=1, \
                    dumpAuxFile=True )
    
    #> Send MBDyn exposed nodes to the other solver, if needed
    n = mbd.socket.nnodes
    print(' n: ', n)
    print(" Initialize MBDyn adapter ")
    #try: 
    adapter = MBDynAdapter(mbd)
    #except:
    #    pass
    if ( adapter.debug ):
        print(' participant: ', adapter.p["name"] )
        print(' solver     : ', adapter.p["mesh"]["name"] )
        print(' fields     : ', adapter.p["fields"] )
    #> Start coupled simulation with PreCICE
    adapter.runPreCICE(dt) 

    os.environ.pop("MBSOCK")    
    os.chdir(oldpwd)
    cwd = os.getcwd() 
