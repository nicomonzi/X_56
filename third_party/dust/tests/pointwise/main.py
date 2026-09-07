import os

from run_dust import *
from run_mbdyn import *
from run_precice import *

# ===========================================================================

# > MBDyn model parameters
dt = 0.1
path2mbdyn = os.path.join('structure')
path2dust = os.path.join('dust')
path2precice = os.path.join('')
mbdyn_model = 'rolling_wing.mbd'


class default():
    def __init__(self):
        # Define common files
        self.dust_pre_input = 'dust_pre.in'
        self.dust_input = 'dust.in'
        self.dust_post = 'dust_post.in'


d = default()
# cleanup previous .log and cached files
cleanup(path2mbdyn, path2dust, path2precice)
path = setup_socket()
run_dust_pre(path2dust, d)
run_dust(path2dust, d)
run_mbdyn(mbdyn_model, 'Output', path2mbdyn, 1)
run_precice(path, dt, path2mbdyn)
run_dust_post()
