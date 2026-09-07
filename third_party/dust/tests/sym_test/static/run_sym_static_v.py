import subprocess as sbprc
import numpy as np
import sys
import fileinput
import argparse
import os
import math as mth
import h5py

pars = argparse.ArgumentParser()
pars.add_argument('-t','--tolerance', nargs='?', help='tolerance for the test',
                  required=False, default=1e-10, type=float)
pars.add_argument('-r','--release', nargs='?', help='release to check with',
                  required=False, default='0.6.0')
args = pars.parse_args()

args = pars.parse_args()
exe_path = 'build/bin'
exe_path = os.path.abspath(exe_path)
pre    = exe_path + '/dust_pre'
solver = exe_path + '/dust'
post   = exe_path + '/dust_post'
tol = args.tolerance

pwd = os.getcwd()
sets_descr = ['Basic wing with symmetry']

#run the first set of simulations, set a
os.chdir("./tests/sym_test/static/v")
pre_ret = sbprc.call([pre])
runs = ['dust.in']
descr = ['symmetry static vortex lattice']
runs_names = ['symmetry static vortex lattice']
sol_dsets = ['/ParticleWake/WakePoints', '/ParticleWake/WakeVort', 
      '/Components/Comp001/Solution/Vort', '/Components/Comp001/Solution/Pres']

sol_descr = ['Particles position   ', 'Particles Intensity  ', 
             'Component 1 Intensity', 'Component 1 Pressure ']

errors = np.zeros([len(runs),len(sol_dsets)])
for i, run in enumerate(runs):
  solv_ret = sbprc.call([solver,run])

  ref_name = 'ref'
  check_name = 'test'
  suffix = '.h5'
  folder = './Output/'
  res = '_res_0011'

  ref_filename = folder+ref_name+res+suffix
  check_filename = folder+check_name+res+suffix
  ref_file = h5py.File(ref_filename, 'r')
  check_file = h5py.File(check_filename, 'r')
  for j, dset_name in enumerate(sol_dsets):
    ref_data   = np.array(ref_file[dset_name])
    check_data = np.array(check_file[dset_name])
    err = np.linalg.norm(check_data-ref_data)/np.linalg.norm(ref_data)
    errors[i,j] = err

#delete all the produced files
files = os.listdir('Output/')
for file in files:
  if file.startswith('test'):
    os.remove(os.path.join('Output/',file))
os.remove('geo_input.h5')

#print the errors:
print('Difference w.r.t. reference:')
for i, run in enumerate(descr):
  print('In run ',run,':')
  for j, dset in enumerate(sol_descr):
    print('Difference on ',dset,': ',errors[i,j])

for i, run in enumerate(runs):
  for j, dset_name in enumerate(sol_dsets):
    errors[i,j] = err
    if err > tol:
      raise TypeError('Test Failed!')

#whatever happened, get back to the previous folder
os.chdir(pwd)
