import os
import sys
import subprocess
import time
import pdb
import argparse
# set to path of MBDyn support for python communication
sys.path.append('/usr/local/mbdyn/libexec/mbpy')

import tempfile
from numpy import *
import numpy as np
# ===========================================================================

#exe_path = '../../../../../build/bin'
exe_path = '../../../../build/bin'
exe_path = os.path.abspath(exe_path)
dust_pre    = exe_path + '/dust_pre'
dust = exe_path + '/dust'
dust_post   = exe_path + '/dust_post'

# ===========================================================================
# Timing parameters (NOTE: need to set also in in dust.in)
dt = 0.1

#> MBDyn model parameters 
nnodes = 1 + 2 + 2

# file name 
file_str = 'test'
dust_log = 'dust_' + file_str + '.log'
mbdyn_file = 'rolling_wing.mbd'

output_mbdyn = 'Output'

if os.path.isdir(output_mbdyn):
    pass
else:
    subprocess.run('mkdir ' + output_mbdyn, shell= True)

# check existens 
tmpdir = tempfile.mkdtemp('', '.mbdyn_')
# ===========================================================================

# clean up
print('\033[0;33m ------------------------------------------ \033[0m')
print('\033[0;33m ◀ Removing __pycache__ and old log files ▶ \033[0m')
print('\033[0;33m ------------------------------------------ \033[0m')

if os.path.isfile('precice-MBDyn-events.json'):
    subprocess.run('rm -fv precice-MBDyn-events.json' , shell=True)
    subprocess.run('rm -fv precice-MBDyn-iterations.log' , shell=True)

if os.path.isdir('__pycache__'):            
    subprocess.run('rm -rfv __pycache__' , shell=True)

if os.path.isdir('/tmp/'  + tmpdir  ):
    subprocess.run('rm -rfv /tmp' + tmpdir , shell=True)

if os.path.isfile('../dust/precice-dust-convergence.log'):
    subprocess.run('rm -fv ../dust/precice-dust-convergence.log' , shell=True)
    subprocess.run('rm -fv ../dust/precice-dust-events.json' , shell=True)
    subprocess.run('rm -fv ../dust/precice-dust-iterations.log' , shell=True)
    subprocess.run('rm -fv log_dust' , shell=True)

# remove dummy precice-run folder
if os.path.isdir('../precice-run'):
    subprocess.run('rm -rfv ../precice-run' , shell=True)
    subprocess.run('rm -fv ../nnodes.dat' , shell=True)
        
tmpdir = tempfile.mkdtemp('', '.mbdyn_')
path = tmpdir + '/mbdyn.sock'
print(' path: ', path)
os.environ['MBSOCK'] = path

# ===========================================================================
# launch DUST
print('\033[0;33m ------------------------------------ \033[0m')
print('\033[0;33m ◀ DUST prepocessor and DUST append ▶ \033[0m')
print('\033[0;33m ------------------------------------ \033[0m')
dust = subprocess.run("cd ../dust &&" + dust_pre + "> log_pre.log &&" + dust + ">" + dust_log + " &", shell=True)

# set launch dust && MBDyn command
str_mbdyn = 'mbdyn -f ' + mbdyn_file + ' -o ' + output_mbdyn + '/' + file_str + ' > ' + output_mbdyn + '/' + file_str + '.txt 2>&1 &'
str_run   =  str_mbdyn  + '> log.mbd &'    
mbdyn = subprocess.run(str_run , shell=True)


#> Import stuff for precice 
from mbc_py_interface import mbcNodal
import precice
from precice import *
from mbdynInterface import MBDynInterface
from mbdynAdapter import MBDynAdapter    

#> Does the other solver needs MBDyn nodes to build its own mesh?
writeNodes = False
#> Construct MBDyn/mbc_py interface
#> Initialize MBDyn/mbc_py interface: negotiate and recv()

mbd = MBDynInterface()

mbd.initialize( path=path, verbose=1, nnodes=nnodes, accels=1, \
                dumpAuxFile=True )

#> ==============================================================
#> Send MBDyn exposed nodes to the other solver, if needed
#> Build reference "Lagrangian" grid for preCICE-MBDyn solver,
#  from the position negotiated through MBDyn-mbc_py interface
print('\033[0;33m ---------------------------- \033[0m')
print('\033[0;33m ◀ Initialize MBDyn adapter ▶ \033[0m')
print('\033[0;33m ---------------------------- \033[0m')
adapter = MBDynAdapter( mbd )

#> ==============================================================
print('\033[0;32m ------------------------------- \033[0m')
print('\033[0;32m ◀ Starting coupled simulation ▶ \033[0m')
print('\033[0;32m ------------------------------- \033[0m')
#> Start coupled simulation with PreCICE
adapter.runPreCICE(dt)
        
os.environ.pop("MBSOCK")

