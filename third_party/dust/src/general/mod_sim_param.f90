!./\\\\\\\\\\\...../\\\......./\\\..../\\\\\\\\\..../\\\\\\\\\\\\\.
!.\/\\\///////\\\..\/\\\......\/\\\../\\\///////\\\.\//////\\\////..
!..\/\\\.....\//\\\.\/\\\......\/\\\.\//\\\....\///.......\/\\\......
!...\/\\\......\/\\\.\/\\\......\/\\\..\////\\.............\/\\\......
!....\/\\\......\/\\\.\/\\\......\/\\\.....\///\\...........\/\\\......
!.....\/\\\......\/\\\.\/\\\......\/\\\.......\///\\\........\/\\\......
!......\/\\\....../\\\..\//\\\...../\\\../\\\....\//\\\.......\/\\\......
!.......\/\\\\\\\\\\\/....\///\\\\\\\\/..\///\\\\\\\\\/........\/\\\......
!........\///////////........\////////......\/////////..........\///.......
!!=========================================================================
!!
!! Copyright (C) 2018-2024 Politecnico di Milano,
!!                           with support from A^3 from Airbus
!!                    and  Davide   Montagnani,
!!                         Matteo   Tugnoli,
!!                         Federico Fonte
!!
!! This file is part of DUST, an aerodynamic solver for complex
!! configurations.
!!
!! Permission is hereby granted, free of charge, to any person
!! obtaining a copy of this software and associated documentation
!! files (the "Software"), to deal in the Software without
!! restriction, including without limitation the rights to use,
!! copy, modify, merge, publish, distribute, sublicense, and/or sell
!! copies of the Software, and to permit persons to whom the
!! Software is furnished to do so, subject to the following
!! conditions:
!!
!! The above copyright notice and this permission notice shall be
!! included in all copies or substantial portions of the Software.
!!
!! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
!! OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
!! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
!! HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
!! WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
!! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
!! OTHER DEALINGS IN THE SOFTWARE.
!!
!! Authors:
!!          Federico Fonte
!!          Davide Montagnani
!!          Matteo Tugnoli
!!=========================================================================


!> Module to define the structrues containing the simulation parameters
!!
module mod_sim_param

use mod_param, only: &
  wp, max_char_len

use mod_handling, only: &
  error, warning, printout

use mod_basic_io, only: &
  read_real_array_from_file  
  
use mod_hdf5_io, only: &
  h5loc, write_hdf5_attr, open_hdf5_file, close_hdf5_file, read_hdf5

use mod_parse, only: &
  t_parse, &
  countoption , &
  getstr, getlogical, getreal, getint, getrealarray, getintarray, &
  ignoredParameters, finalizeParameters

implicit none

public :: t_sim_param, sim_param, create_param_main, create_param_pre, &
          create_param_post, create_param_test_particle, &  
          init_sim_param, init_sim_param_particle

private

type t_sim_param

  !For debugging purpose
  character(max_char_len) :: basename_debug

  !Time:
  !> Start time
  real(wp) :: t0
  !> Time step
  real(wp) :: dt
  !> Final time
  real(wp) :: tend
  !> Number of timesteps
  integer  :: n_timesteps
  !> Vector of time instants
  real(wp) , allocatable :: time_vec(:)
  character(len=max_char_len) :: integrator 
  !> Actual time
  real(wp) :: time
  !> Previous time
  real(wp) :: time_old
  !> ndt between 2 wake updates
  integer :: ndt_update_wake
  !> file particles
  character(len=max_char_len) :: particles_file 

  !> Output detailed geometry each timestep
  logical :: output_detailed_geo

  !Physical parameters:
  !> Altitude
  real(wp) :: altitude 
  !> Units of the input data
  character(len=max_char_len) :: units 
  !> Free stream pressure
  real(wp) :: P_inf
  !> Free stream density
  real(wp) :: rho_inf
  !> Free stream velocity
  real(wp) :: u_inf(3)
  !> Reference velocity (magnitude of u_inf unless specified)
  real(wp) :: u_ref
  !> Free stream speed of sound
  real(wp) :: a_inf
  !> Free stream dynamic viscosity
  real(wp) :: mu_inf
  !> Free stream kinematic viscosity
  real(wp) :: nu_inf

  !Wake
  !> Scaling of the first implicit panel
  real(wp) :: first_panel_scaling
  !> Minimum velocity at the trailing edge
  real(wp) :: min_vel_at_te
  !> Is the wake rigid?
  logical :: rigid_wake
    !> Velocity of the rigid wake
    real(wp)  :: rigid_wake_vel(3)
  !> Number of wake panels
  integer :: n_wake_panels
  !> Number of wake particles
  integer :: n_wake_particles
  !> Minimum and maximum of the particles box
  real(wp) :: particles_box_min(3)
  real(wp) :: particles_box_max(3)
  !> Minimun particle magnitude allowed (suppress if lower)
  real(wp) :: mag_threshold
  !> Minimum particle core radius allowed (suppress if lower)
  real(wp) :: sigma_threshold
  !> Maximum particle magnitude allowed (clamp if higher, <=0 disables the check)
  real(wp) :: mag_max
  !> Maximum particle core radius allowed in rVPM (<=0 disables the check)
  real(wp) :: sigma_max
  !> Maximum relative particle intensity increment per wake update (<=0 disables the check)
  real(wp) :: alpha_rel_increment_max

    !> Wake initial condition
  integer  :: part_n0 
  real(wp), allocatable :: part_pos0(:,:)
  real(wp), allocatable :: part_vort0_dir(:,:)
  real(wp), allocatable :: part_vort0_mag(:)
  real(wp), allocatable :: part_vol(:) 
  
  !> Join close trailing edges
  logical :: join_te
  !> All trailing edges closer than join_te_factor will be joined
  real(wp) :: join_te_factor
  !> Wake refinement with subparticles
  logical :: refine_wake  
  !> k_refine
  integer :: k_refine
  !> refine tolerance
  real(wp) :: tol_refine
  !> Wake interpolation
  logical :: interpolate_wake
  !> Trailing edge autoscaling
  logical :: autoscale_te

  !Method parameters:  
  !> Multiplier for far field threshold computation on doublet
  real(wp) :: FarFieldRatioDoublet
  !> Multiplier for far field threshold computation on sources
  real(wp) :: FarFieldRatioSource
  !> Thresold for considering the point in plane in doublets
  real(wp) :: DoubletThreshold
  !> Rankine Radius for vortices
  real(wp) :: RankineRad
  !> Octree vortex particle size core
  real(wp) :: octree_vortex_rad
  !> Vortex Radius for vortex particles
  real(wp) :: VortexRad
  !> Vortex Radius coefficient for vortex particles
  !> if too low or negative reverts to original behaviour (VortexRad)
  real(wp) :: KVortexRad  
  !> Particle volume computation coefficient (applies to rVol)
  real(wp) :: KVol 
  !> Complete cutoff radius
  real(wp) :: CutoffRad
  !> use the vortex stretching or not
  logical :: use_vs
  !> use the vortex stretching from elements
  logical :: vs_elems
  !> use the divergence filtering
  logical :: use_divfilt
  !> Pedrizzetti relaxation coefficient for divergence filtering
  real(wp):: alpha_divfilt
  !> use the vorticity diffusion or not
  logical :: use_vd
  !> use turbulent viscosity or not
  logical :: use_tv
  !> use the penetration avoidance
  logical :: use_pa
  !> check radius for penetration avoidance
  real(wp) :: pa_rad_mult
  !> element impact radius for penetration avoidance
  real(wp) :: pa_elrad_mult
  !> simulate viscosity effects or not
  logical :: use_ve
  !> use reformulated formulation rVPM (Alvarez 2023)
  logical   :: use_reformulated
    !> rVPM coefficients
    real(wp)  :: f                   
    real(wp)  :: g
  !> suppress the initial wake
  logical :: suppress_wake
  !> suppress the initial wake after n timesteps
  integer :: suppress_wake_nsteps 
  !> suppress the wake below y_fountain before n timesteps
  integer :: suppress_fountain_nsteps 

  !Lifting Lines
  character(len=max_char_len) :: llSolver
  !> Reynolds corrections of .c81 tables
  logical  :: llReynoldsCorrections
  !> n factor for Reynolds corrections of .c81 tables: (Re/Re_T)^n
  real(wp) :: llReynoldsCorrectionsNfact
  !> Maximum number of iteration in LL algorithm
  integer  :: llMaxIter
  !> Tolerance for the relative error in fixed point iteration for LL
  real(wp) :: llTol 
  !> Damping param in fixed point iteration for LL used to avoid oscillations
  real(wp) :: llDamp
  !> Avoid "unphysical" separations in inner sections of LL? :: llTol
  logical  :: llStallRegularisation
  !> Number of "unphysical" separations that can be removed
  integer  :: llStallRegularisationNelems
  !> Number of iterations between two regularisation processes
  integer  :: llStallRegularisationNiters
  !> Reference stall AOA for regularisation
  real(wp) :: llStallRegularisationAlphaStall
  !> Constant Artificial Viscosity for regularisation
  real(wp) :: llArtificialViscosity
  !> Adaptive Artificial Viscosity algorithm
  logical  :: llArtificialViscosityAdaptive
  !> Adaptive Artificial Viscosity algorithm, reference AOA
  real(wp) :: llArtificialViscosityAdaptive_Alpha
  !> Adaptive Artificial Viscosity algorithm, blending interval
  ! between full regularisation and no regularisation ( AOA )
  real(wp) :: llArtificialViscosityAdaptive_dAlpha
  !> Use AVL expression for inviscid load computation ( ~ VL )
  logical  :: llLoadsAVL

  
  !FMM parameters
  !> Employing the FMM method
  logical :: use_fmm
    !> Employing the FMM method also for panels
    logical :: use_fmm_pan
    !> Size of the Octree box
    real(wp) :: BoxLength
    !> Number of boxes in each direction
    integer :: NBox(3)
    !> Origin of the Octree system of boxes
    real(wp) :: OctreeOrigin(3)
    !> Number of Octree levels
    integer :: NOctreeLevels
    !> Minimum number of particles for each box
    integer :: MinOctreePart
    !> Multipole expansion degree
    integer :: MultipoleDegree
    !> Use dynamic levels
    logical :: use_dyn_layers
      !> Maximum number of octree levels
      integer :: NMaxOctreeLevels
      !> Time ratio that triggers the increase of levels
      real(wp) :: LeavesTimeRatio
    !> use particles redistribution
    logical :: use_pr
      !> Level at which is checked the presence of panels
      integer :: lvl_solid
      real(wp) :: part_redist_ratio

  !HCAS parameters
  !> Use hcas
  logical :: hcas
    !> Time of deployment of the hcas
    real(wp) :: hcas_time
    !> Velocity of the hcas
    real(wp) :: hcas_vel(3)


  !Handling parameters:
  !> Debug level
  integer :: debug_level
  !> Output interval
  real(wp) :: dt_out
  !> Basename
  character(len=max_char_len) :: basename
  !> Geometry file
  character(len=max_char_len) :: GeometryFile
  !> References file
  character(len=max_char_len) :: ReferenceFile
  !> Restart from file
  logical :: restart_from_file
  !> Restart file
  character(len=max_char_len) :: restart_file
  !> Reset the time after restart
  logical :: reset_time

  !Variable wind
  !> Use time-varying free stream velocity
  logical  :: time_varying_u_inf
  !> Use non-uniform free stream velocity (along one direction)
  logical  :: non_uniform_u_inf
    !> Filename for time-varying or non-uniform free stream velocity
    character(len=max_char_len) :: u_inf_filen
    real(wp), allocatable       :: u_inf_mat(:,:)
    real(wp), allocatable       :: u_inf_comps(:,:)
    real(wp), allocatable       :: u_inf_time(:)
    real(wp), allocatable       :: u_inf_coord(:)
    integer                     :: non_uniform_u_inf_dir
    integer                     :: nt_u_inf
    integer                     :: nc_u_inf
  !> Gust
  logical :: use_gust
  !> Gust type
  character(len=max_char_len) :: GustType
  !> Gust parameters
  real(wp) :: gust_origin(3)
  real(wp) :: gust_front_direction(3)
  real(wp) :: gust_front_speed
  real(wp) :: gust_u_des
  real(wp) :: gust_perturb_direction(3)
  real(wp) :: gust_gradient
  real(wp) :: gust_time

  !> Vl correction 
  logical   :: vl_correction = .false. 
  character(len=max_char_len) :: vlSolver
  real(wp)  :: vl_tol
  real(wp)  :: vl_relax
  real(wp)  :: vl_alpha_step_max
  real(wp)  :: vl_dclda_min
  integer   :: vl_maxiter
  integer   :: vl_startstep
  integer   :: vl_iter_ave 
  logical   :: vl_dynstall = .false.
  logical   :: rel_aitken  = .false. 
  logical   :: vl_ave      = .false.

  !> kutta parameters 
  logical   :: kutta_correction = .false.
  real(wp)  :: kutta_beta 
  real(wp)  :: kutta_tol
  integer   :: kutta_maxiter
  integer   :: kutta_startstep
  integer   :: kutta_update_jacobian

  !> PreCICE
#if USE_PRECICE
  character(len=max_char_len) :: precice_config
#endif
  
contains
  procedure, pass(this) :: save_param => save_sim_param
end type t_sim_param

type(t_sim_param) :: sim_param
character(len=*), parameter     :: this_mod_name = 'mod_sim_param' 

real(wp), allocatable :: particlesMat(:,:)

!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------

!> Subroutines for main dust, dust_pre and dust_post prms creation 
!  to avoid clutter in source files
!> Create parameter object for dust main parameters
subroutine create_param_main(prms)
  type(t_parse), intent(inout) :: prms
  
  !> Define the parameters to be read
  !> Time
  call prms%CreateRealOption('tstart', 'Starting time')
  call prms%CreateRealOption('tend',   'Ending time')
  call prms%CreateRealOption('dt',     'time step')
  call prms%CreateIntOption ('timesteps', 'number of timesteps')
  call prms%CreateRealOption('dt_out', 'output time interval')
  call prms%CreateRealOption('dt_debug_out', 'debug output time interval')
  call prms%CreateIntOption ('ndt_update_wake', 'n. dt between two wake updates', '1')
  
  !> Reformulated formulation                                         
  call prms%CreateLogicalOption('reformulated','Employ rVPM by Alvarez','T')
  call prms%CreateRealOption('f','rVPM coefficient f','0.0')
  call prms%CreateRealOption('g','rVPM coefficient g','0.2')
  !> treatment of the initial wake 
  call prms%CreateLogicalOption('suppress_wake','Suppress the initial wake','F') 
  call prms%CreateIntOption('suppress_wake_nsteps','Suppress the initial wake after n timesteps', '100') 
  !> remove fountain effect 
  call prms%CreateIntOption('suppress_fountain_nsteps','Suppress the wake below y_fountain before n timesteps', '1')  
  !> Integration 
  call prms%CreateStringOption('integrator', &
                            &'integrator solver: Euler or low storage RK','euler')  
  !> Input
  call prms%CreateStringOption('geometry_file','Main geometry definition file')
  call prms%CreateStringOption('reference_file','Reference frames file','no_set')
  
  !> Output
  call prms%CreateStringOption('basename','output basename','./')
  call prms%CreateStringOption('basename_debug','output basename for debug','./')
  call prms%CreateLogicalOption('output_start', "output values at starting &
                                                            & iteration", 'F')
  call prms%CreateLogicalOption('output_detailed_geo', "output at each &
                      &timestep the detailed geometry in the results file", 'F')
  call prms%CreateIntOption('debug_level', 'Level of debug verbosity/output', '1')
  
  !> Restart
  call prms%CreateLogicalOption('restart_from_file','restarting from file?','F')
  call prms%CreateStringOption('restart_file','restart file name')
  call prms%CreateLogicalOption('reset_time','reset the time from previous execution?','F')
  
  !> Parameters: reference conditions 
  call prms%CreateRealOption('altitude', "Altitude in meters", '0.0') 
  call prms%CreateStringOption('units', "units of the input data", 'SI')  
  call prms%CreateRealArrayOption('u_inf', "free stream velocity", '(/1.0, 0.0, 0.0/)')
  call prms%CreateRealOption('u_ref', "reference velocity")             
  call prms%CreateRealOption('P_inf', "free stream pressure", '101325')    
  call prms%CreateRealOption('rho_inf', "free stream density", '1.225')   
  call prms%CreateRealOption('a_inf', "Speed of sound", '340.0')        ! m/s   for dimensional sim
  call prms%CreateRealOption('mu_inf', "Dynamic viscosity", '0.000018') ! kg/ms

  !> Wake 
  call prms%CreateIntOption('n_wake_panels', 'number of wake panels','1')
  call prms%CreateIntOption('n_wake_particles', 'number of wake particles', '10000')
  call prms%CreateRealArrayOption('particles_box_min', 'min coordinates of the &
                                  &particles bounding box', '(/-10.0, -10.0, -10.0/)')
  call prms%CreateRealArrayOption('particles_box_max', 'max coordinates of the &
                                  &particles bounding box', '(/10.0, 10.0, 10.0/)')
  call prms%CreateRealOption('mag_threshold', "Minimum particle magnitude allowed", '1.0e-9')
  call prms%CreateRealOption('sigma_threshold', "Minimum particle core radius allowed", '1.0e-9')
  call prms%CreateRealOption('mag_max', "Maximum particle magnitude allowed (clamp)", '1.0e+10')
  call prms%CreateRealOption('sigma_max', 'Maximum particle core radius allowed in rVPM (-1 disables the check)', '-1.0')
  call prms%CreateRealOption('alpha_rel_increment_max', 'Maximum relative &
            &particle intensity increment per wake update (-1 disables the check)', '-1.0')
  
  call prms%CreateRealOption('implicit_panel_scale', "Scaling of the first implicit wake panel", '0.3')
  call prms%CreateRealOption('implicit_panel_min_vel', "Minimum velocity at the trailing edge", '1.0e-8')
  call prms%CreateLogicalOption('rigid_wake','rigid wake?','F')
  call prms%CreateRealArrayOption('rigid_wake_vel', "rigid wake velocity" )
  call prms%CreateLogicalOption('join_te','join trailing edge','F')
  call prms%CreateRealOption('join_te_factor', "join the trailing edges when closer than factor*te element size",'1.0' )
  call prms%CreateLogicalOption('refine_wake','refined wake with subparticles','F')
  call prms%CreateIntOption('k_refine','refine factor for wake subdivision with subparticles','1')
  call prms%CreateRealOption('tol_refine','tolerance for wake refinement','0.2')
  call prms%CreateLogicalOption('interpolate_wake','interpolate wake subparticles','F')
  call prms%CreateLogicalOption('autoscale_te','Autoscale the first te panel as dx / (Vlocal*dt) & 
                                where dx is the panel width (spanwise direction) and Vlocal is   &
                                norm2(Vinf-vte)','F') 

  !> Regularisation 
  call prms%CreateStringOption('particles_file', 'file with particles initial condition', 'particles.dat') 
  call prms%CreateRealOption('far_field_ratio_doublet', &
        "Multiplier for far field threshold computation on doublet", '10.0')
  call prms%CreateRealOption('far_field_ratio_source', &
        "Multiplier for far field threshold computation on sources", '10.0')
  call prms%CreateRealOption('doublet_threshold', &
        "Thresold for considering the point in plane in doublets", '1.0e-6')
  call prms%CreateRealOption('rankine_rad', &
        "Radius of Rankine correction for vortex induction near core", '0.1')
  call prms%CreateRealOption('octree_vortex_rad', &
        "Vortex particle core size for octree initialization", '0.1')
  call prms%CreateRealOption('vortex_rad', &
        "Radius of vortex core, for particles", '0.1')
  call prms%CreateRealOption('k_vortex_rad', &
        "Radius coefficient of vortex core, for particles", '1.0') ! default is ON
  call prms%CreateRealOption('k_vol', &
        "Radius coefficient of volume, for particles", '1.0') 
  call prms%CreateRealOption('cutoff_rad', &
        "Radius of complete cutoff  for vortex induction near core", '0.001')
  
  !> Lifting line elements
  call prms%CreateStringOption('ll_solver','Solver for the LL elements', &
                            &'gamma_method')
  call prms%CreateLogicalOption('ll_reynolds_corrections', &
                            &'Use Reynolds corrections for the .c81 tables?', 'F')
  call prms%CreateRealOption('ll_reynolds_corrections_nfact', &
                            &'Exponent in (Re/Re_T)^n correction', '0.2')
  call prms%CreateIntOption('ll_max_iter', &
                            &'Maximum number of iteration in LL algorithm', '100')
  call prms%CreateRealOption('ll_tol', 'Tolerance for the relative error in &
                            &fixed point iteration for LL','1.0e-6' )
  call prms%CreateRealOption('ll_damp', 'Damping param in fixed point iteration &
                            &for LL used to avoid oscillations','25.0')
  call prms%CreateLogicalOption('ll_stall_regularisation', &
                            &'Avoid "unphysical" separations in inner sections of LL?','T')
  call prms%CreateIntOption('ll_stall_regularisation_nelems', &
                            &'Number of "unphysical" separations to be removed', '1' )
  call prms%CreateIntOption('ll_stall_regularisation_niters', &
                            &'Number of timesteps between two regularisations', '1' )
  call prms%CreateRealOption('ll_stall_regularisation_alpha_stall', &
                            &'Stall angle used as threshold for regularisation [deg]', '15.0' )
  call prms%CreateRealOption('ll_artificial_viscosity', &
                            &'Constant artificial viscosity for regularisation', '0.0' )
  call prms%CreateLogicalOption('ll_artificial_viscosity_adaptive', &
                            &'Adaptive artificial viscosity algorithm', 'F' )
  call prms%CreateRealOption('ll_artificial_viscosity_adaptive_alpha', &
                            &'Adaptive Artificial Viscosity algorithm, reference AOA [deg]')
  call prms%CreateRealOption('ll_artificial_viscosity_adaptive_dalpha', &
                            &'Adaptive Artificial Viscosity algorithm, reference AOA [deg]')
  call prms%CreateLogicalOption('ll_loads_avl', &
                            &'Use AVL expression for inviscid load computation','F')
  
  !> VL correction parameter 
  call prms%CreateStringOption('vl_solver', &
                            &'Solver for VL correction: gamma_method or alpha_method', 'gamma_method')
  call prms%CreateRealOption('vl_relax', 'Relaxation factor for rhs update','0.3')
  call prms%CreateIntOption('vl_maxiter', &
                            &'Maximum number of iteration in VL algorithm', '100')
  call prms%CreateRealOption('vl_tol', 'Tolerance for the absolute error on lift coefficient in &
                            &fixed point iteration for VL','1.0e-2' )
  call prms%CreateRealOption('vl_alpha_step_max', 'Maximum absolute alpha increment per VL iteration [deg] (alpha_method)', '2.0')
  call prms%CreateRealOption('vl_dclda_min', 'Minimum absolute dCl/dalpha used in alpha_method [1/rad]', '3.0')
  call prms%CreateIntOption('vl_start_step', &
                            &'Step in which the VL correction start', '0')
  call prms%CreateLogicalOption('vl_dynstall', 'Dynamic stall on corrected VL', 'F')
  call prms%CreateLogicalOption('aitken_relaxation', 'Employ aitken acceleration method during &   
                                &the fixed point iteration', 'T')  
  call prms%CreateLogicalOption('vl_average', 'Average panel intensity between the last iterations', 'F')  
  call prms%CreateIntOption('vl_average_iter', 'Number of iterations to average', '10')  
  
  !> kutta correction parameters
  call prms%CreateLogicalOption('kutta_correction', 'Employ kutta condition', 'F') 
  call prms%CreateRealOption('kutta_beta', 'perturbation factor for kutta condition', '0.01') 
  call prms%CreateRealOption('kutta_tol', 'tolerance for kutta condition', '1.0e-4') 
  call prms%CreateIntOption('kutta_maxiter', 'maximum number of iterations for kutta condition', '100')
  call prms%CreateIntOption('kutta_start_step', 'step in which the kutta condition starts', '1') 
  call prms%CreateIntOption('kutta_update_jacobian', 'step frequency where the Jacobian is updated', '1')   
  
  !> Octree and multipole data 
  call prms%CreateLogicalOption('fmm','Employ fast multipole method?','T')
  call prms%CreateLogicalOption('fmm_panels','Employ fast multipole method &
                                &also for panels?','F')
  call prms%CreateRealOption('box_length','length of the octree box')
  call prms%CreateIntArrayOption('n_box','number of boxes in each direction')
  call prms%CreateRealArrayOption( 'octree_origin', "rigid wake velocity" )
  call prms%CreateIntOption('n_octree_levels','number of octree levels')
  call prms%CreateIntOption('min_octree_part','minimum number of octree particles')
  call prms%CreateIntOption('multipole_degree','multipole expansion degree')
  call prms%CreateLogicalOption('dyn_layers','Use dynamic layers','F')
  call prms%CreateIntOption('nmax_octree_levels','maximum number of octree levels')
  call prms%CreateRealOption('leaves_time_ratio','Ratio that triggers the &
                                            &increase of the number of levels')
  
  !> Models options
  call prms%CreateLogicalOption('vortstretch','Employ vortex stretching','T')
  call prms%CreateLogicalOption('vortstretch_from_elems','Employ vortex stretching&
                                & from geometry elements','F')
  call prms%CreateLogicalOption('divergence_filtering','Employ divergence filtering', 'T')
  call prms%CreateRealOption('alpha_divfilt','Pedrizzetti relaxation coefficient', '0.3')
  call prms%CreateLogicalOption('diffusion','Employ vorticity diffusion','T')
  call prms%CreateLogicalOption('turbulent_viscosity','Employ turbulent &
                                &viscosity','F')
  call prms%CreateLogicalOption('penetration_avoidance','Employ penetration avoidance','F')
  call prms%CreateRealOption('penetration_avoidance_check_radius', &
        'Check radius for penetration avoidance','5.0')
  call prms%CreateRealOption('penetration_avoidance_element_radius', &
        'Element impact radius for penetration avoidance','1.5')
  call prms%CreateLogicalOption('viscosity_effects','Simulate viscosity &
                                                                & effects','F')
  call prms%CreateLogicalOption('particles_redistribution','Employ particles &
                                                          &redistribution','F')
  call prms%CreateIntOption('octree_level_solid','Level at which the panels &
                            & are considered for particles redistribution')
  call prms%CreateRealOption('particles_redistribution_ratio','How many times &
            &a particle need to be smaller than the average of the cell to be&
            & eliminated','3.0')
  
  !> HCAS
  call prms%CreateLogicalOption('HCAS','Hover Convergence Augmentation System', 'F')
  call prms%CreateRealOption('HCAS_time','HCAS deployment time')
  call prms%CreateRealArrayOption('HCAS_velocity','HCAS velocity')
  
  !> Variable wind
  call prms%CreateLogicalOption('time_varying_u_inf', "Use time-varying free stream velocity", 'F')
  call prms%CreateLogicalOption('non_uniform_u_inf', "Use non-uniform free stream velocity", 'F')
  call prms%CreateIntOption('non_uniform_u_inf_dir','Direction of variation of non-uniform u_inf')
  call prms%CreateStringOption('u_inf_file', "file .dat containing time/coordinates + variable Ux, Uy, Uz")
  call prms%CreateLogicalOption('gust','Gust perturbation','F')
  call prms%CreateStringOption('gust_type','Gust model','AMC')
  call prms%CreateRealArrayOption('gust_origin','Gust origin point')
  call prms%CreateRealArrayOption('gust_front_direction','Gust front direction vector')
  call prms%CreateRealOption('gust_front_speed','Gust front speed')
  call prms%CreateRealOption('gust_u_des','Design gust velocity')
  call prms%CreateRealArrayOption('gust_perturbation_direction','Gust perturbation &
                                &direction vector','(/0.0, 0.0, 1.0/)')
  call prms%CreateRealOption('gust_gradient','Gust gradient')
  call prms%CreateRealOption('gust_start_time','Gust starting time','0.0')
  
  !> preCICE
#if USE_PRECICE
  call prms%CreateStringOption('precice_config','PreCICE configuration file','./../precice-config.xml')
#endif

end subroutine create_param_main

!> Create parameter object for dust preprocessor parameters
subroutine create_param_pre(prms)
  type(t_parse), intent(inout) :: prms
    
  call prms%CreateStringOption('comp_name','Component Name', multiple=.true.)
  call prms%CreateStringOption('geo_file','Geometry definition files', multiple=.true.)
  call prms%CreateStringOption('ref_tag','Reference Tag of the component', multiple=.true.)
  call prms%CreateRealOption('tol_se_wing','Global parameter for closing gaps','0.001')
  call prms%CreateRealOption('inner_product_te','Global parameter for edge identification','-0.5')
  call prms%CreateStringOption('file_name','Preprocessor output file')

end subroutine create_param_pre

!> Create parameter object for dust postprocessor parameters
subroutine create_param_post(prms, sbprms, bxprms)
  type(t_parse), intent(inout) :: prms
  type(t_parse), pointer, intent(inout) :: sbprms , bxprms

  call prms%CreateStringOption('basename','Base name of the processed data')
  call prms%CreateStringOption('data_basename','Base name of the data to be &
                                &processed')
  
  call prms%CreateRealOption( 'far_field_ratio_doublet', &
        "Multiplier for far field threshold computation on doublet", '10.0')
  call prms%CreateRealOption( 'far_field_ratio_source', &
        "Multiplier for far field threshold computation on sources", '10.0')
  call prms%CreateRealOption( 'doublet_threshold', &
        "Thresold for considering the point in plane in doublets", '1.0e-6')
  call prms%CreateRealOption( 'rankine_rad', &
        "Radius of Rankine correction for vortex induction near core", '0.1')
  call prms%CreateRealOption('octree_vortex_rad', &
        "Vortex particle core size for octree initialization", '0.1')
  call prms%CreateRealOption( 'vortex_rad', &
        "Radius of vortex core, for particles", '0.1')
  call prms%CreateRealOption( 'cutoff_rad', &
        "Radius of complete cutoff  for vortex induction near core", '0.001')
  
  call prms%CreateSubOption('analysis','Definition of the motion of a frame', &
                            sbprms, multiple=.true.)
  call sbprms%CreateStringOption('type','type of analysis')
  call sbprms%CreateStringOption('name','specification of the analysis')
  call sbprms%CreateIntOption('start_res', 'Starting result of the analysis')
  call sbprms%CreateIntOption('end_res', 'Final result of the analysis')
  call sbprms%CreateIntOption('step_res', 'Result stride of the analysis')
  call sbprms%CreateLogicalOption('average', 'Perform time averaging','F')
  call sbprms%CreateIntOption('avg_res', 'Result for the average visualization')
  call sbprms%CreateLogicalOption('wake', 'Output also the wake for &
                                  &visualization','T')
  call sbprms%CreateLogicalOption('separate_wake', 'Output the wake in a separate &
                                  &way','F')
  call sbprms%CreateStringOption('format','Output format')
  call sbprms%CreateStringOption('component','Component to analyse', &
                                multiple=.true.)
  call sbprms%CreateStringOption('hinge_tag','Hinge to analyse', &
                                multiple=.true.)                              
  call sbprms%CreateStringOption('variable','Variables to be saved: velocity, pressure or&
                                & vorticity', multiple=.true.)
  
  ! probe output -------------
  call sbprms%CreateStringOption('input_type','How to specify probe coordinates',&
                                multiple=.true.)
  call sbprms%CreateRealArrayOption('point','Point coordinates in dust_post.in',&
                                multiple=.true.)
  call sbprms%CreateStringOption('file','File containing the coordinates of the probes',&
                                multiple=.true.)
  ! flow field output --------
  call sbprms%CreateIntArrayOption( 'n_xyz','number of points per coordinate',&
                                multiple=.true.)
  call sbprms%CreateRealArrayOption('min_xyz','lower bounds of the box',&
                                multiple=.true.)
  call sbprms%CreateRealArrayOption('max_xyz','upper bounds of the box',&
                                multiple=.true.)
  ! loads --------------------
  call sbprms%CreateStringOption('comp_name','Components where loads are computed',&
                                multiple=.true.)
  call sbprms%CreateStringOption('reference_tag','Reference frame where loads&
                              & are computed',multiple=.true.)
  ! sectional loads ----------
  call sbprms%CreateRealArrayOption('axis_dir','Direction of the axis defined the reference&
                              & points for sectional loads analisys', multiple=.true.)
  call sbprms%CreateRealArrayOption('axis_nod','Node belonging to the axis used for sectional&
                              & loads analisys', multiple=.true.)
  call sbprms%CreateLogicalOption('lifting_line_data', 'Output lifting line data&
                              & alongside sectional loads','F')
  call sbprms%CreateLogicalOption('vortex_lattice_data', 'Output of corrected vortex lattice data&
                              & alongside sectional loads','F')
  ! chordwise loads
  call sbprms%CreateIntOption('n_station','Number of stations where the loads are extracted' &
                              , multiple=.true.)
  call sbprms%CreateRealArrayOption('span_station','Spanwise coordinates in the component reference frame&
                              & where they are extracted', multiple=.true.)
  
  
  ! sectional loads: box -----
  call sbprms%CreateSubOption('box_sect','Definition of the box for sectional loads', &
                            bxprms)
  call bxprms%CreateRealArrayOption('ref_node','reference node to build the box')
  call bxprms%CreateRealArrayOption('face_vec','vector identifying the direction of the base side &
                            &of the sections')
  call bxprms%CreateRealArrayOption('face_bas','dimension along faceVec of the first and last sections')
  call bxprms%CreateRealArrayOption('face_hei','dimension orthogonal to faceVec of the first and last sections')
  call bxprms%CreateRealArrayOption('span_vec','vector defining the out-of plane direction of the box')
  call bxprms%CreateRealOption('span_len','dimension along the spanVec direction')
  call bxprms%CreateIntOption('num_sect','number of sections')
  call bxprms%CreateLogicalOption('reshape_box','logical input to reshape the box if &
                            &it is "too large"')
  call sbprms%CreateRealArrayOption('axis_mom','axis for the computation of the moment. Perpendicular to sections')
  
  sbprms=>null()

end subroutine create_param_post

!> Subroutines for main dust, dust_pre and dust_post prms creation 
!  to avoid clutter in source files
!> Create parameter object for dust main parameters
subroutine create_param_test_particle(prms)
  type(t_parse), intent(inout) :: prms
  
  !> Define the parameters to be read
  !> Time
  call prms%CreateRealOption('tstart', "Starting time")
  call prms%CreateRealOption('tend',   "Ending time")
  call prms%CreateRealOption('dt',     "time step")
  call prms%CreateIntOption ('timesteps', "number of timesteps")
  call prms%CreateRealOption('dt_out', "output time interval")
  call prms%CreateRealOption('dt_debug_out', "debug output time interval")
  call prms%CreateIntOption ('ndt_update_wake', 'n. dt between two wake updates', '1')

  call prms%CreateStringOption('basename','output basename','./')
  call prms%CreateStringOption('basename_debug','output basename for debug','./')
  call prms%CreateLogicalOption('output_start', "output values at starting &
                                                            & iteration", 'F')
  call prms%CreateIntOption('debug_level', 'Level of debug verbosity/output', '1') 
  call prms%CreateStringOption('particles_file', 'file with particles initial condition', 'particles.dat') 
  
  !> Restart
  call prms%CreateLogicalOption('restart_from_file','restarting from file?','F')
  call prms%CreateStringOption('restart_file','restart file name')
  call prms%CreateLogicalOption('reset_time','reset the time from previous execution?','F')

  !> Parameters: reference conditions 
  call prms%CreateRealArrayOption('u_inf', "free stream velocity", '(/1.0, 0.0, 0.0/)')
  call prms%CreateRealOption('u_ref', "reference velocity")             
  call prms%CreateRealOption('P_inf', "free stream pressure", '101325')    
  call prms%CreateRealOption('rho_inf', "free stream density", '1.225')   
  call prms%CreateRealOption('a_inf', "Speed of sound", '340.0')        ! m/s   for dimensional sim
  call prms%CreateRealOption('mu_inf', "Dynamic viscosity", '0.000018') ! kg/ms

  call prms%CreateIntOption('n_wake_particles', 'number of wake particles', '100000')
  call prms%CreateRealArrayOption('particles_box_min', 'min coordinates of the &
                                  &particles bounding box', '(/-10.0, -10.0, -10.0/)')
  call prms%CreateRealArrayOption('particles_box_max', 'max coordinates of the &
                                  &particles bounding box', '(/10.0, 10.0, 10.0/)')
  call prms%CreateRealOption('mag_threshold', "Minimum particle magnitude allowed", '1.0e-9')
  call prms%CreateRealOption('sigma_threshold', "Minimum particle core radius allowed", '1.0e-9')
  call prms%CreateRealOption('mag_max', "Maximum particle magnitude allowed (clamp)", '1.0e+10')
  call prms%CreateRealOption('sigma_max', 'Maximum particle core radius allowed in rVPM (-1 disables the check)', '-1.0')
  call prms%CreateRealOption('alpha_rel_increment_max', 'Maximum relative &
            &particle intensity increment per wake update (-1 disables the check)', '-1.0')

  !> Regularisation 
  call prms%CreateRealOption('rankine_rad', &
        "Radius of Rankine correction for vortex induction near core", '0.1')
  call prms%CreateRealOption('octree_vortex_rad', &
        "Vortex particle core size for octree initialization", '0.1')
  call prms%CreateRealOption('vortex_rad', &
        "Radius of vortex core, for particles", '0.1')
  call prms%CreateRealOption('k_vortex_rad', &
        "Radius coefficient of vortex core, for particles", '1.0') ! default is ON
  call prms%CreateRealOption('cutoff_rad', &
        "Radius of complete cutoff  for vortex induction near core", '0.001')

  !> Octree and multipole data 
  call prms%CreateLogicalOption('fmm','Employ fast multipole method?','T')
  call prms%CreateLogicalOption('fmm_panels','Employ fast multipole method &
                                &also for panels?','F')
  call prms%CreateRealOption('box_length','length of the octree box')
  call prms%CreateIntArrayOption('n_box','number of boxes in each direction')
  call prms%CreateRealArrayOption( 'octree_origin', "rigid wake velocity" )
  call prms%CreateIntOption('n_octree_levels','number of octree levels')
  call prms%CreateIntOption('min_octree_part','minimum number of octree particles')
  call prms%CreateIntOption('multipole_degree','multipole expansion degree')
  call prms%CreateLogicalOption('dyn_layers','Use dynamic layers','F')
  call prms%CreateIntOption('nmax_octree_levels','maximum number of octree levels')
  call prms%CreateRealOption('leaves_time_ratio','Ratio that triggers the &
                                            &increase of the number of levels')

    !> Models options
  call prms%CreateLogicalOption('vortstretch','Employ vortex stretching','T')
  call prms%CreateLogicalOption('vortstretch_from_elems','Employ vortex stretching&
                                & from geometry elements','F')
  call prms%CreateLogicalOption('divergence_filtering','Employ divergence filtering','T')
  call prms%CreateRealOption('alpha_divfilt','Pedrizzetti relaxation coefficient','0.3')
  call prms%CreateLogicalOption('diffusion','Employ vorticity diffusion','T')
  call prms%CreateLogicalOption('turbulent_viscosity','Employ turbulent &
                                &viscosity','F') 
  call prms%CreateLogicalOption('viscosity_effects','Simulate viscosity &
                                                                & effects','F')
  call prms%CreateLogicalOption('particles_redistribution','Employ particles &
                                                          &redistribution','F')
  call prms%CreateIntOption('octree_level_solid','Level at which the panels &
                            & are considered for particles redistribution')
  call prms%CreateRealOption('particles_redistribution_ratio','How many times &
            &a particle need to be smaller than the average of the cell to be&
            & eliminated','3.0')

  !> Reformulated formulation                                         
  call prms%CreateLogicalOption('reformulated','Employ rVPM by Alvarez','T')
  call prms%CreateRealOption('f','rVPM coefficient f','0.0')
  call prms%CreateRealOption('g','rVPM coefficient g','0.2')

  !> Integrators
  call prms%CreateStringOption('integrator', 'integrator solver: Euler or low storage RK', &
                              'euler') 


end subroutine create_param_test_particle   


!> Initialize all the parameters reading them from the the input file
subroutine init_sim_param(sim_param, prms, nout, output_start)
  class(t_sim_param)          :: sim_param
  type(t_parse)               :: prms
  integer, intent(inout)      :: nout
  logical, intent(inout)      :: output_start

  character(len=*), parameter :: this_sub_name = 'init_sim_param'
  
  !> Timing
  sim_param%t0                  = getreal(prms,     'tstart')
  sim_param%tend                = getreal(prms,     'tend')
  sim_param%dt_out              = getreal(prms,     'dt_out')
  sim_param%debug_level         = getint(prms,      'debug_level')
  sim_param%output_detailed_geo = getlogical(prms,  'output_detailed_geo')
  sim_param%ndt_update_wake     = getint(prms,      'ndt_update_wake')
  
  !> Integration
  sim_param%integrator          = getstr(prms,  'integrator') 
  !> Reference environment values 
  sim_param%altitude            = getreal(prms, 'altitude') 
  sim_param%units               = getstr(prms,  'units')
  if (sim_param%altitude .ne. 0.0_wp) then
    call standard_atmosphere(sim_param) 
  else
    sim_param%P_inf               = getreal(prms, 'P_inf')
    sim_param%rho_inf             = getreal(prms, 'rho_inf')
    sim_param%a_inf               = getreal(prms, 'a_inf')
    sim_param%mu_inf              = getreal(prms, 'mu_inf')
  endif 

  sim_param%nu_inf              = sim_param%mu_inf/sim_param%rho_inf
  sim_param%u_inf               = getrealarray(prms, 'u_inf', 3)
  
  !> Check on reference velocity
  if ( countoption(prms,'u_ref') .gt. 0 ) then
    sim_param%u_ref = getreal(prms, 'u_ref')
  else
    sim_param%u_ref = norm2(sim_param%u_inf)
    if (sim_param%u_ref .le. 1e-10_wp) then
      call error(this_mod_name, this_sub_name,'No reference velocity u_ref provided but &
      &zero free stream velocity. Provide a non-zero reference velocity. &
      &Stopping now before producing invalid results')
    endif
  end if
  
  !> Wake parameters
  sim_param%n_wake_panels         = getint(prms,      'n_wake_panels')
  sim_param%n_wake_particles      = getint(prms,      'n_wake_particles')
  sim_param%particles_box_min     = getrealarray(prms,'particles_box_min',3)
  sim_param%particles_box_max     = getrealarray(prms,'particles_box_max',3)
  sim_param%mag_threshold         = getreal(prms,     'mag_threshold')
  sim_param%sigma_threshold       = getreal(prms,     'sigma_threshold')
  sim_param%mag_max               = getreal(prms,     'mag_max')
  sim_param%sigma_max             = getreal(prms,     'sigma_max')
  sim_param%alpha_rel_increment_max = getreal(prms,   'alpha_rel_increment_max')
  sim_param%rigid_wake            = getlogical(prms,  'rigid_wake')
  sim_param%rigid_wake_vel        = sim_param%u_inf   !> initialisation
  sim_param%refine_wake           = getlogical(prms,  'refine_wake')
  sim_param%k_refine              = getint(prms,      'k_refine')
  sim_param%tol_refine            = getreal(prms,     'tol_refine')
  sim_param%interpolate_wake      = getlogical(prms,  'interpolate_wake')
  sim_param%autoscale_te          = getlogical(prms,  'autoscale_te') 
  !> Check on wake refinement
  if (sim_param%interpolate_wake .and. .not. sim_param%refine_wake) then
        !call warning('dust', 'dust', 'Wake interpolation is selected, but wake refinement &
        !           is not, overriding refinement to True. ')
        sim_param%refine_wake = .true.
  endif

  !> Check on wake panels
  if(sim_param%n_wake_panels .lt. 1) then
    sim_param%n_wake_panels = 1
    call warning(this_mod_name, this_sub_name,'imposed a number of wake panels rows &
                &LOWER THAN 1. At least one row of panels is mandatory, &
                &the simulation will proceed with "n_wake_panels = 1"')
  endif
  if ( sim_param%rigid_wake ) then
    if ( countoption(prms,'rigid_wake_vel') .eq. 1 ) then
      sim_param%rigid_wake_vel    = getrealarray(prms, 'rigid_wake_vel',3)
    else if ( countoption(prms,'rigid_wake_vel') .le. 0 ) then
      call warning(this_mod_name, this_sub_name,'no rigid_wake_vel parameter set, &
            &with rigid_wake = T; rigid_wake_vel = u_inf')
      sim_param%rigid_wake_vel    = sim_param%u_inf
    end if
  end if

  !> Trailing edge 
  sim_param%join_te               = getlogical(prms, 'join_te')
  if (sim_param%join_te) sim_param%join_te_factor=getreal(prms,'join_te_factor')

  !> Names
  sim_param%basename              = getstr(prms, 'basename')
  sim_param%GeometryFile          = getstr(prms, 'geometry_file')
  sim_param%ReferenceFile         = getstr(prms, 'reference_file')

  !> Method parameters
  sim_param%FarFieldRatioDoublet  = getreal(prms,    'far_field_ratio_doublet')
  sim_param%FarFieldRatioSource   = getreal(prms,    'far_field_ratio_source')
  sim_param%DoubletThreshold      = getreal(prms,    'doublet_threshold')
  sim_param%RankineRad            = getreal(prms,    'rankine_rad')
  sim_param%octree_vortex_rad     = getreal(prms,    'octree_vortex_rad')
  sim_param%VortexRad             = getreal(prms,    'vortex_rad')
  sim_param%KVortexRad            = getreal(prms,    'k_vortex_rad')
  sim_param%KVol                  = getreal(prms,    'k_vol')
  sim_param%CutoffRad             = getreal(prms,    'cutoff_rad')
  sim_param%first_panel_scaling   = getreal(prms,    'implicit_panel_scale')
  sim_param%min_vel_at_te         = getreal(prms,    'implicit_panel_min_vel')
  sim_param%alpha_divfilt         = getreal(prms,    'alpha_divfilt')
  sim_param%use_vs                = getlogical(prms, 'vortstretch')
  sim_param%vs_elems              = getlogical(prms, 'vortstretch_from_elems')
  sim_param%use_vd                = getlogical(prms, 'diffusion')
  sim_param%use_divfilt           = getlogical(prms, 'divergence_filtering')
  sim_param%use_tv                = getlogical(prms, 'turbulent_viscosity')
  sim_param%use_ve                = getlogical(prms, 'viscosity_effects')
  sim_param%use_pa                = getlogical(prms, 'penetration_avoidance')
  !> Check on penetration avoidance 
  if(sim_param%use_pa) then
    sim_param%pa_rad_mult = getreal(prms, 'penetration_avoidance_check_radius')
    sim_param%pa_elrad_mult = getreal(prms,'penetration_avoidance_element_radius')
  endif
    
  !> Reformulated formulation (Alvarez rVPM 2023)
  sim_param%use_reformulated      = getlogical(prms, 'reformulated')

  if(sim_param%use_reformulated) then
    sim_param%f                   = getreal(prms, 'f')
    sim_param%g                   = getreal(prms, 'g')
  endif
  !> suppress wake 
  sim_param%suppress_wake         = getlogical(prms, 'suppress_wake') 
  !> n time steps for wake suppression
  sim_param%suppress_wake_nsteps  = getint(prms, 'suppress_wake_nsteps') 
  !> n time steps for fountain suppression (start releasing wake from inboard blade) 
  sim_param%suppress_fountain_nsteps = getint(prms, 'suppress_fountain_nsteps') 

  !> Lifting line elements
  sim_param%llSolver                        = getstr(    prms, 'll_solver')
  sim_param%llReynoldsCorrections           = getlogical(prms, 'll_reynolds_corrections')
  sim_param%llReynoldsCorrectionsNfact      = getreal(   prms, 'll_reynolds_corrections_nfact')
  sim_param%llMaxIter                       = getint(    prms, 'll_max_iter'            )
  sim_param%llTol                           = getreal(   prms, 'll_tol'                )
  sim_param%llDamp                          = getreal(   prms, 'll_damp'               )
  sim_param%llStallRegularisation           = getlogical(prms, 'll_stall_regularisation')
  sim_param%llStallRegularisationNelems     = getint(    prms, 'll_stall_regularisation_nelems')
  sim_param%llStallRegularisationNiters     = getint(    prms, 'll_stall_regularisation_niters')
  sim_param%llStallRegularisationAlphaStall = getreal(   prms, 'll_stall_regularisation_alpha_stall')
  sim_param%llArtificialViscosity           = getreal(   prms, 'll_artificial_viscosity')
  sim_param%llArtificialViscosityAdaptive   = getlogical(prms, 'll_artificial_viscosity_adaptive')
  sim_param%llLoadsAVL                      = getlogical(prms, 'll_loads_avl')

  select case(trim(sim_param%llSolver))
    case('GammaMethod', 'gamma_method', 'gamma')
      sim_param%llSolver = 'GammaMethod'
    case('AlphaMethod', 'alpha_method', 'alpha')
      sim_param%llSolver = 'AlphaMethod'
  end select

  !> check LL inputs
  if (trim(sim_param%llSolver) .ne. 'GammaMethod') then
    write(*,*) ' sim_param%llSolver : ' , trim(sim_param%llSolver)
    call warning(this_mod_name, this_sub_name,' Wrong string for LLsolver. &
                &This parameter is set equal to "GammaMethod" (default) &
                &in init_sim_param() routine.')
    sim_param%llSolver = 'GammaMethod'
  end if

  if ( trim(sim_param%llSolver) .eq. 'GammaMethod' ) then
    if ( sim_param%llArtificialViscosity .gt. 0.0_wp ) then
      call warning(this_mod_name, this_sub_name,'LLartificialViscoisty set as an input, &
          & different from zero, but ll regularisation available only if&
          & LLsolver = AlphaMethod. LLartificialViscosity = 0.0')
      sim_param%llArtificialViscosity = 0.0_wp
      sim_param%llArtificialViscosityAdaptive = .false.
      sim_param%llArtificialViscosityAdaptive_Alpha  = 0.0_wp
      sim_param%llArtificialViscosityAdaptive_dAlpha = 0.0_wp
    end if
    if ( sim_param%llArtificialViscosityAdaptive ) then
      call warning(this_mod_name, this_sub_name,'LLartificialViscosityAdaptive set&
          & as an input, but ll adaptive regularisation available only if&
          & LLsolver = AlphaMethod')
    end if
  end if

  !> The user is required to set _Alpha and _dAlpha for ll adaptive regularisation
  if (trim(sim_param%llSolver) .eq. 'AlphaMethod') then
    if (sim_param%llArtificialViscosityAdaptive) then
      if ((countoption(prms,'ll_artificial_viscosity_adaptive_alpha')  .eq. 0) .or. &
          (countoption(prms,'ll_artificial_viscosity_adaptive_dalpha') .eq. 0)) then
        call error(this_mod_name, this_sub_name,'LLartificialViscosityAdaptive_Alpha or&
          & LLartificialViscosity_dAlpha not set as an input, while LLartificialViscosityAdaptive&
          & is set equal to T. Set these parameters [deg].')
      else
        sim_param%llArtificialViscosityAdaptive_Alpha = getreal(prms,'ll_artificial_viscosity_adaptive_alpha')
        sim_param%llArtificialViscosityAdaptive_dAlpha= getreal(prms,'ll_artificial_viscosity_adaptive_dalpha')
      end if
    end if
  end if
  !> VL correction 
  sim_param%vlSolver                    = getstr(prms, 'vl_solver')
  sim_param%vl_tol                        = getreal(prms, 'vl_tol')
  sim_param%vl_relax                      = getreal(prms, 'vl_relax')
  sim_param%vl_alpha_step_max             = getreal(prms, 'vl_alpha_step_max')
  sim_param%vl_dclda_min                  = getreal(prms, 'vl_dclda_min')
  sim_param%vl_maxiter                    = getint(prms, 'vl_maxiter')
  sim_param%vl_startstep                  = getint(prms, 'vl_start_step')
  sim_param%rel_aitken                    = getlogical(prms, 'aitken_relaxation')
  sim_param%vl_ave                        = getlogical(prms,'vl_average')
  sim_param%vl_iter_ave                   = getint(prms,'vl_average_iter') 

  select case(trim(sim_param%vlSolver))
    case('GammaMethod', 'gamma_method', 'gamma')
      sim_param%vlSolver = 'gamma_method'
    case('AlphaMethod', 'alpha_method', 'alpha')
      sim_param%vlSolver = 'alpha_method'
    case default
      call warning(this_mod_name, this_sub_name,'Wrong string for vl_solver. Setting gamma_method (default).')
      sim_param%vlSolver = 'gamma_method'
  end select

  if (sim_param%vl_alpha_step_max .le. 0.0_wp) then
    call warning(this_mod_name, this_sub_name,'vl_alpha_step_max <= 0.0: setting to 2.0 deg')
    sim_param%vl_alpha_step_max = 2.0_wp
  end if

  if (sim_param%vl_dclda_min .le. epsilon(1.0_wp)) then
    call warning(this_mod_name, this_sub_name,'vl_dclda_min too small: setting to 3.0 1/rad')
    sim_param%vl_dclda_min = 3.0_wp
  end if
  
  !>  VL Dynamic stall
  sim_param%vl_dynstall                   = getlogical(prms, 'vl_dynstall')  
  
  !> Kutta condition 
  sim_param%kutta_correction              = getlogical(prms, 'kutta_correction') 
  sim_param%kutta_beta                    = getreal(prms, 'kutta_beta')
  sim_param%kutta_tol                     = getreal(prms, 'kutta_tol')
  sim_param%kutta_maxiter                 = getint(prms, 'kutta_maxiter')
  sim_param%kutta_startstep               = getint(prms, 'kutta_start_step') 
  sim_param%kutta_update_jacobian         = getint(prms,'kutta_update_jacobian')
  !> Octree and FMM parameters
  sim_param%use_fmm                       = getlogical(prms, 'fmm')

  if(sim_param%use_fmm) then
    sim_param%use_fmm_pan                 = getlogical(prms, 'fmm_panels')
    sim_param%BoxLength                   = getreal(prms, 'box_length')
    sim_param%NBox                        = getintarray(prms, 'n_box',3)
    sim_param%OctreeOrigin                = getrealarray(prms, 'octree_origin',3)
    sim_param%NOctreeLevels               = getint(prms, 'n_octree_levels')
    sim_param%MinOctreePart               = getint(prms, 'min_octree_part')
    sim_param%MultipoleDegree             = getint(prms,'multipole_degree')
    sim_param%use_dyn_layers              = getlogical(prms,'dyn_layers')

    if(sim_param%use_dyn_layers) then
      sim_param%NMaxOctreeLevels          = getint(prms, 'nmax_octree_levels')
      sim_param%LeavesTimeRatio           = getreal(prms, 'leaves_time_ratio')
    else
      sim_param%NMaxOctreeLevels          = sim_param%NOctreeLevels
    endif

    sim_param%use_pr                      = getlogical(prms, 'particles_redistribution')

    if(sim_param%use_pr) then
      sim_param%part_redist_ratio         = getreal(prms,'particles_redistribution_ratio')
      if ( countoption(prms,'octree_level_solid') .gt. 0 ) then
        sim_param%lvl_solid               = getint(prms, 'octree_level_solid')
      else
        sim_param%lvl_solid               = max(sim_param%NOctreeLevels-2,1)
      endif
    endif
  else
    sim_param%use_fmm_pan = .false.
  endif

  !> HCAS
  sim_param%hcas                          = getlogical(prms,'HCAS')
  if(sim_param%hcas) then
    sim_param%hcas_time                   = getreal(prms,'HCAS_time')
    sim_param%hcas_vel                    = getrealarray(prms,'HCAS_velocity',3)
  endif

  !> Variable_wind
  sim_param%time_varying_u_inf      = getlogical(prms,'time_varying_u_inf')
  sim_param%non_uniform_u_inf       = getlogical(prms,'non_uniform_u_inf')
  
  if ( sim_param%time_varying_u_inf .and. sim_param%non_uniform_u_inf ) then
    call error(this_mod_name, this_sub_name,'u_inf can not be both non-uniform and time-varying')
  end if

  if ( sim_param%time_varying_u_inf) then
    if ( countoption(prms,'u_inf_file') .eq. 0 ) then
      call error(this_mod_name, this_sub_name,'Time-varying u_inf file .dat not defined')
    end if
    sim_param%u_inf_filen = trim(getstr(prms,'u_inf_file'))
    !>  Read time and u_inf components
    call read_real_array_from_file ( 4 , trim(sim_param%u_inf_filen) , sim_param%u_inf_mat )
    sim_param%nt_u_inf = size(sim_param%u_inf_mat,1)
    allocate(sim_param%u_inf_comps(3,sim_param%nt_u_inf))
    allocate(sim_param%u_inf_time(sim_param%nt_u_inf)) 
    
    sim_param%u_inf_time  = sim_param%u_inf_mat(:,1)
    sim_param%u_inf_comps = transpose(sim_param%u_inf_mat(:,2:4))

  elseif (sim_param%non_uniform_u_inf) then
    if ( countoption(prms,'u_inf_file') .eq. 0 ) then
      call error(this_mod_name, this_sub_name, 'Non-uniform u_inf file .dat not defined')
    end if 
    sim_param%u_inf_filen = trim(getstr(prms,'u_inf_file'))
    if ( countoption(prms,'non_uniform_u_inf_dir') .eq. 0 ) then
      call error(this_mod_name, this_sub_name, 'Non uniform u_inf direction of variation not provided in dust.in')
    end if
    sim_param%non_uniform_u_inf_dir = getint(prms, 'non_uniform_u_inf_dir')
    if (sim_param%non_uniform_u_inf_dir .ne. 1 .and. &
        sim_param%non_uniform_u_inf_dir .ne. 2 .and. &
        sim_param%non_uniform_u_inf_dir .ne. 3 ) then
      call error(this_mod_name, this_sub_name, 'Non uniform u_inf direction of variation provided is neither 1, 2 or 3')
    end if
    call read_real_array_from_file ( 4 , trim(sim_param%u_inf_filen) , sim_param%u_inf_mat )
    sim_param%nc_u_inf = size(sim_param%u_inf_mat,1)
    allocate(sim_param%u_inf_comps(3,sim_param%nc_u_inf))
    allocate(sim_param%u_inf_coord(sim_param%nc_u_inf)) 
    ! Read coordinates and u_inf components
    sim_param%u_inf_coord  = sim_param%u_inf_mat(:,1)
    sim_param%u_inf_comps = transpose(sim_param%u_inf_mat(:,2:4))
    if ( sim_param%u_inf_coord(1) .gt. sim_param%u_inf_coord(sim_param%nc_u_inf)) then
      call error(this_mod_name, this_sub_name,'Please, provide u_inf file&
                & with coordinate column spanning from lower to greater values')
    end if
    if ( sim_param%u_inf_coord(1) .gt. sim_param%particles_box_min(sim_param%non_uniform_u_inf_dir)) then
      call error(this_mod_name, this_sub_name, 'Velocity profile coordinates list must start& 
                & before or at the relative particle_box_min coordinate')
    else if (sim_param%u_inf_coord(sim_param%nc_u_inf) .lt. &
            & sim_param%particles_box_max(sim_param%non_uniform_u_inf_dir)) then
      call error( this_mod_name, this_sub_name, 'Velocity profile coordinates list must end&
                & after or at the relative particle_box_max coordinate')
    end if
  end if

  sim_param%use_gust                      = getlogical(prms, 'gust')

  if(sim_param%use_gust) then
    sim_param%GustType                    = getstr(prms,'gust_type')
    sim_param%gust_origin                 = getrealarray(prms, 'gust_origin',3)

    if(countoption(prms,'gust_front_direction') .gt. 0) then
      sim_param%gust_front_direction      = getrealarray(prms, 'gust_front_direction',3)
    else
      sim_param%gust_front_direction      = sim_param%u_inf
    end if
    
    if(countoption(prms,'gust_front_speed') .gt. 0) then
      sim_param%gust_front_speed          = getreal(prms, 'gust_front_speed')
    else
      sim_param%gust_front_speed          = norm2(sim_param%u_inf)
    end if
    sim_param%gust_u_des                  = getreal(prms,'gust_u_des')
    sim_param%gust_perturb_direction      = getrealarray(prms,'gust_perturbation_direction',3)
    sim_param%gust_time                   = getreal(prms,'gust_start_time')
    sim_param%gust_gradient               = getreal(prms,'gust_gradient')
  end if
  
  !> PreCICE
#if USE_PRECICE
    sim_param%precice_config              = getstr(prms,'precice_config')
#endif
  
  !> Manage restart
  sim_param%restart_from_file             = getlogical(prms,'restart_from_file')
  if (sim_param%restart_from_file) then

    sim_param%reset_time                  = getlogical(prms,'reset_time')
    sim_param%restart_file                = getstr(prms,'restart_file')
    
    !> Removing leading "./" if present to avoid issues when restarting
    if(sim_param%basename(1:2) .eq. './') sim_param%basename = sim_param%basename(3:)
    
    if(sim_param%restart_file(1:2) .eq. './') sim_param%restart_file = sim_param%restart_file(3:)
    
    call printout('RESTART: restarting from file: '//trim(sim_param%restart_file))
    sim_param%GeometryFile = sim_param%restart_file(1:len(trim(sim_param%restart_file))-11) //'geo.h5'

    !restarting the same simulation, advance the numbers
    if(sim_param%restart_file(1:len(trim(sim_param%restart_file))-12).eq. &
                                                trim(sim_param%basename)) then
    read(sim_param%restart_file(len(trim(sim_param%restart_file))-6:len(trim(sim_param%restart_file))-3),*) nout
      call printout('Identified restart from the same simulation, keeping the&
                    & previous output numbering')
      !> avoid rewriting the same timestep
      output_start = .false.
    endif
    if(.not. sim_param%reset_time) call load_time(sim_param%restart_file, sim_param%t0)
  endif

  !> Check the number of timesteps
  if(CountOption(prms,'dt') .gt. 0) then
    if( CountOption(prms,'timesteps') .gt. 0) then
      call error(this_mod_name, this_sub_name, 'Both number of timesteps and dt are&
                                              & set, but only one of the two can be specified')
    else
      !> get dt and compute number of timesteps
      sim_param%dt     = getreal(prms, 'dt')
      sim_param%n_timesteps = ceiling((sim_param%tend-sim_param%t0)/sim_param%dt) + 1
                              !(+1 for the zero time step)
    endif
  else
    !> get number of steps, compute dt
    sim_param%n_timesteps = getint(prms, 'timesteps')
    sim_param%dt =  (sim_param%tend-sim_param%t0)/&
                      real(sim_param%n_timesteps,wp)
    sim_param%n_timesteps = sim_param%n_timesteps + 1
                            !add one for the first step
  endif
  !> If use variable u_inf is active, check the .dat file time consistency with simulation time
  if ( sim_param%time_varying_u_inf) then
    if ( sim_param%u_inf_time(1) .gt. sim_param%t0 ) then
      write(*,*) ' beginning of the time of variable u_inf : ' , sim_param%u_inf_time(1)
      write(*,*) ' beginning of the simulation time : ' , sim_param%t0
      call error(this_mod_name, this_sub_name,'Error in variable u_inf. &
        & Initial time value of the .dat file greater than initial simulation time.')
    end if
    if ( sim_param%u_inf_time(size(sim_param%u_inf_time)) .lt. sim_param%tend ) then
        write(*,*) ' end of the time of variable u_inf' , sim_param%u_inf_time(size(sim_param%u_inf_time))
        write(*,*) ' end of the time in simulation time : ' , sim_param%tend
        call error(this_mod_name, this_sub_name, 'Error in variable u_inf. &
          & Final time value of the .dat file lower than final simulation time.')
    end if
  end if

end subroutine init_sim_param

!> Initialize all the parameters reading them from the the input file
subroutine init_sim_param_particle(sim_param, prms, nout, output_start)
  class(t_sim_param)          :: sim_param
  type(t_parse)               :: prms
  integer, intent(inout)      :: nout
  logical, intent(inout)      :: output_start
  
  character(len=*), parameter :: this_sub_name = 'init_sim_param'
  !> Timing
  sim_param%t0                  = getreal(prms, 'tstart')
  sim_param%tend                = getreal(prms, 'tend')
  if(CountOption(prms,'dt') .gt. 0) then
    if( CountOption(prms,'timesteps') .gt. 0) then
      call error(this_mod_name, this_sub_name, 'Both number of timesteps and dt are&
                                              & set, but only one of the two can be specified')
    else
      !> get dt and compute number of timesteps
      sim_param%dt     = getreal(prms, 'dt')
      sim_param%n_timesteps = ceiling((sim_param%tend-sim_param%t0)/sim_param%dt) + 1
                              !(+1 for the zero time step)
    endif
  else
    !> get number of steps, compute dt
    sim_param%n_timesteps = getint(prms, 'timesteps')
    sim_param%dt =  (sim_param%tend-sim_param%t0)/&
                      real(sim_param%n_timesteps,wp)
    sim_param%n_timesteps = sim_param%n_timesteps + 1
                            !add one for the first step
  endif 

  sim_param%dt_out              = getreal(prms,'dt_out')
  sim_param%n_timesteps         = ceiling((sim_param%tend-sim_param%t0)/sim_param%dt) + 1
  sim_param%debug_level         = 1
  sim_param%debug_level         = getint(prms, 'debug_level')  
  
  !> file dat 
  sim_param%particles_file      = getstr(prms, 'particles_file') 
  !> Reference environment values
  sim_param%P_inf               = getreal(prms,'P_inf')
  sim_param%rho_inf             = getreal(prms,'rho_inf')
  sim_param%a_inf               = getreal(prms,'a_inf')
  sim_param%mu_inf              = getreal(prms,'mu_inf')
  sim_param%nu_inf              = sim_param%mu_inf/sim_param%rho_inf
  sim_param%u_inf               = getrealarray(prms, 'u_inf', 3)
  
  !> Check on reference velocity
  if ( countoption(prms,'u_ref') .gt. 0 ) then
    sim_param%u_ref = getreal(prms, 'u_ref')
  else
    sim_param%u_ref = norm2(sim_param%u_inf)
    if (sim_param%u_ref .le. 0.0_wp) then
      call error(this_mod_name, this_sub_name,'No reference velocity u_ref provided but &
      &zero free stream velocity. Provide a non-zero reference velocity. &
      &Stopping now before producing invalid results')
    endif
  end if
  
  !> Wake parameters
  sim_param%n_wake_particles      = getint(prms, 'n_wake_particles')
  sim_param%particles_box_min     = getrealarray(prms, 'particles_box_min',3)
  sim_param%particles_box_max     = getrealarray(prms, 'particles_box_max',3)
  sim_param%mag_threshold         = getreal(prms, 'mag_threshold')
  sim_param%sigma_threshold       = getreal(prms, 'sigma_threshold')
  sim_param%mag_max               = getreal(prms, 'mag_max')
  sim_param%sigma_max             = getreal(prms, 'sigma_max')
  sim_param%alpha_rel_increment_max = getreal(prms, 'alpha_rel_increment_max')

  sim_param%basename              = getstr(prms, 'basename')
  !Integrators 
  sim_param%integrator            = getstr(prms, 'integrator')
  !> Manage restart
  sim_param%restart_from_file             = getlogical(prms,'restart_from_file')
  if (sim_param%restart_from_file) then

    sim_param%reset_time                  = getlogical(prms,'reset_time')
    sim_param%restart_file                = getstr(prms,'restart_file')
    
    !> Removing leading "./" if present to avoid issues when restarting
    if(sim_param%basename(1:2) .eq. './') sim_param%basename = sim_param%basename(3:)
    
    if(sim_param%restart_file(1:2) .eq. './') sim_param%restart_file = sim_param%restart_file(3:)
    
    !call printout('RESTART: restarting from file: '//trim(sim_param%restart_file))
    !sim_param%GeometryFile = sim_param%restart_file(1:len(trim(sim_param%restart_file))-11) //'geo.h5'

    !restarting the same simulation, advance the numbers
    if(sim_param%restart_file(1:len(trim(sim_param%restart_file))-12).eq. &
                                                trim(sim_param%basename)) then
    read(sim_param%restart_file(len(trim(sim_param%restart_file))-6:len(trim(sim_param%restart_file))-3),*) nout
      call printout('Identified restart from the same simulation, keeping the&
                    & previous output numbering')
      !> avoid rewriting the same timestep
      output_start = .false.
    endif
    if(.not. sim_param%reset_time) call load_time(sim_param%restart_file, sim_param%t0)
  endif

  !> Method parameters
  sim_param%RankineRad            = getreal(prms,    'rankine_rad')
  sim_param%octree_vortex_rad     = getreal(prms,    'octree_vortex_rad')
  sim_param%VortexRad             = getreal(prms,    'vortex_rad')
  sim_param%CutoffRad             = getreal(prms,    'cutoff_rad')
  sim_param%use_vs                = getlogical(prms, 'vortstretch')
  sim_param%use_vd                = getlogical(prms, 'diffusion')
  sim_param%use_divfilt           = getlogical(prms, 'divergence_filtering')
  sim_param%use_tv                = getlogical(prms, 'turbulent_viscosity')
  sim_param%alpha_divfilt         = getreal(prms,    'alpha_divfilt')

  !> Reformulated formulation (Alvarez rVPM 2023)
  sim_param%use_reformulated      = getlogical(prms, 'reformulated')

  if(sim_param%use_reformulated) then
    sim_param%f                   = getreal(prms, 'f')
    sim_param%g                   = getreal(prms, 'g')
  endif

  !> Octree and FMM parameters
  sim_param%use_fmm                       = getlogical(prms, 'fmm')

  if(sim_param%use_fmm) then
    sim_param%use_fmm_pan                 = getlogical(prms, 'fmm_panels')
    sim_param%BoxLength                   = getreal(prms, 'box_length')
    sim_param%NBox                        = getintarray(prms, 'n_box',3)
    sim_param%OctreeOrigin                = getrealarray(prms, 'octree_origin',3)
    sim_param%NOctreeLevels               = getint(prms, 'n_octree_levels')
    sim_param%MinOctreePart               = getint(prms, 'min_octree_part')
    sim_param%MultipoleDegree             = getint(prms,'multipole_degree')
    sim_param%use_dyn_layers              = getlogical(prms,'dyn_layers')

    if(sim_param%use_dyn_layers) then
      sim_param%NMaxOctreeLevels          = getint(prms, 'nmax_octree_levels')
      sim_param%LeavesTimeRatio           = getreal(prms, 'leaves_time_ratio')
    else
      sim_param%NMaxOctreeLevels          = sim_param%NOctreeLevels
    endif

    sim_param%use_pr                      = getlogical(prms, 'particles_redistribution')

    if(sim_param%use_pr) then
      sim_param%part_redist_ratio         = getreal(prms,'particles_redistribution_ratio')
      if ( countoption(prms,'octree_level_solid') .gt. 0 ) then
        sim_param%lvl_solid               = getint(prms, 'octree_level_solid')
      else
        sim_param%lvl_solid               = max(sim_param%NOctreeLevels-2,1)
      endif
    endif
  else
    sim_param%use_fmm_pan = .false.
  endif

  sim_param%ndt_update_wake       = 1

end subroutine init_sim_param_particle 

subroutine standard_atmosphere(sim_param)
  class(t_sim_param), intent(inout) :: sim_param 

  real(wp) :: P0    ! Sea-level standard atmospheric pressure (Pa) or (lbf/ft^2)
  real(wp) :: L     ! Temperature lapse rate (K/m) or (R/ft)
  real(wp) :: T0    ! Sea-level standard temperature (Kelvin) or (Rankine)
  real(wp) :: g     ! Acceleration due to gravity (m/s^2) or (ft/s^2)
  real(wp) :: M     ! Molar mass of Earth's air (kg/mol) or (lb/mol)
  real(wp) :: R     ! Universal gas constant (J/(mol·K)) or (ft·lbf/(lb·mol·R))
  real(wp) :: gamma ! Adiabatic index for air
  real(wp) :: T, h

  select case(sim_param%units)
    case('SI')
      P0    = 101325.0_wp 
      L     = 6.5e-3_wp 
      T0    = 288.15_wp   
      g     = 9.80665_wp  
      M     = 0.0289644_wp
      R     = 8.31447_wp  
    case('imperial')
      P0    = 2116.22_wp  
      L     = 3.6e-3_wp 
      T0    = 518.67_wp   
      g     = 32.1740_wp  
      M     = 0.0289644_wp
      R     = 1545.35_wp  
  end select

  gamma = 1.4_wp      ! Adiabatic index for air

  !> altitude (m) or (ft)
  h = sim_param%altitude

  ! Calculate temperature at altitude h (K) or (R)
  T = T0 - L * h

  ! Calculate pressure using the standard atmosphere model equation (Pa) or (lbf/ft^2)
  sim_param%P_inf = P0 * (1.0_wp - (L * h) / T0)**(g * M / (R * L))

  ! Calculate density using the ideal gas law (kg/m^3) or (slug/ft^3)
  sim_param%rho_inf = sim_param%P_inf * M / (R * T)

  ! Calculate sound speed (m/s) or (ft/s)
  sim_param%a_inf = sqrt(gamma * R * T)

  ! Calculate viscosity (approximately proportional to temperature)
  select case (sim_param%units)
    case('SI')
      sim_param%mu_inf = T / T0 * 1.7894e-5_wp ! Viscosity at sea level (Pa·s)
    case('imperial')
      sim_param%mu_inf = T / T0 * 3.737e-7_wp ! Viscosity at sea level (lb/(ft·s))
  end select

end subroutine standard_atmosphere

subroutine save_sim_param(this, loc)
  class(t_sim_param) :: this
  integer(h5loc), intent(in) :: loc

  call write_hdf5_attr(this%t0,                  't0'                 , loc)
  call write_hdf5_attr(this%dt,                  'dt'                 , loc)
  call write_hdf5_attr(this%tend,                'tend'               , loc)
  call write_hdf5_attr(this%output_detailed_geo, 'output_detailed_geo', loc)
  call write_hdf5_attr(this%P_inf,               'P_inf'              , loc)
  call write_hdf5_attr(this%rho_inf,             'rho_inf'            , loc)
  call write_hdf5_attr(this%u_inf,               'u_inf'              , loc)
  call write_hdf5_attr(this%u_ref,               'u_ref'              , loc)
  call write_hdf5_attr(this%a_inf,               'a_inf'              , loc)
  call write_hdf5_attr(this%mu_inf,              'mu_inf'             , loc)
  call write_hdf5_attr(this%first_panel_scaling, 'first_panel_scaling', loc)
  call write_hdf5_attr(this%min_vel_at_te,       'min_vel_at_te'      , loc)
  call write_hdf5_attr(this%rigid_wake,          'rigid_wake'         , loc)
  if(this%rigid_wake) &
    call write_hdf5_attr(this%rigid_wake_vel,    'rigid_wake_vel'     , loc)
  call write_hdf5_attr(this%n_wake_panels,       'n_wake_panels'      , loc)
  call write_hdf5_attr(this%n_wake_particles,    'n_wake_particles'   , loc)
  call write_hdf5_attr(this%particles_box_min,   'particles_box_min'  , loc)
  call write_hdf5_attr(this%particles_box_max,   'particles_box_max'  , loc)
  call write_hdf5_attr(this%mag_threshold,       'mag_threshold'      , loc)
  call write_hdf5_attr(this%sigma_threshold,     'sigma_threshold'    , loc)
  call write_hdf5_attr(this%mag_max,             'mag_max'            , loc)
  call write_hdf5_attr(this%sigma_max,           'sigma_max'          , loc)
  call write_hdf5_attr(this%alpha_rel_increment_max, 'alpha_rel_increment_max', loc)

  call write_hdf5_attr(this%join_te,             'join_te'            , loc)
  if(this%join_te) &
    call write_hdf5_attr(this%join_te_factor,    'join_te_factor'     , loc)


  call write_hdf5_attr(this%FarFieldRatioDoublet, 'FarFieldRatioDoublet', loc)
  call write_hdf5_attr(this%FarFieldRatioSource,  'FarFieldRatioSource' , loc)
  call write_hdf5_attr(this%DoubletThreshold,     'DoubletThreshold'    , loc)
  call write_hdf5_attr(this%RankineRad,           'RankineRad'          , loc)
  call write_hdf5_attr(this%VortexRad,            'VortexRad'           , loc)
  call write_hdf5_attr(this%KVortexRad,           'KVortexRad'          , loc)
  call write_hdf5_attr(this%CutoffRad,            'CutoffRad'           , loc)
  call write_hdf5_attr(this%use_vs,               'Vortstretch'         , loc)
  if(this%use_vs) then
    call write_hdf5_attr(this%vs_elems,           'VortstretchFromElems', loc)
    call write_hdf5_attr(this%use_divfilt,        'DivergenceFiltering' , loc)
    call write_hdf5_attr(this%alpha_divfilt,      'AlphaDivFilt'        , loc)
  endif
  call write_hdf5_attr(this%use_vd,               'vortdiff'            , loc)
  call write_hdf5_attr(this%use_tv,               'turbvort'            , loc)
  call write_hdf5_attr(this%use_pa,               'PenetrationAvoidance', loc)
  if(this%use_pa) then
    call write_hdf5_attr(this%pa_rad_mult,        'PenetrationAvoidanceCheckRadius', loc)
    call write_hdf5_attr(this%pa_elrad_mult,      'PenetrationAvoidanceElementRadius', loc)
  endif
  call write_hdf5_attr(this%use_ve,               'ViscosityEffects', loc)
  call write_hdf5_attr(this%use_fmm,              'use_fmm', loc)
  if(this%use_fmm) then
    call write_hdf5_attr(this%use_fmm_pan,         'use_fmm_panels', loc)
    call write_hdf5_attr(this%BoxLength,           'BoxLength', loc)
    call write_hdf5_attr(this%Nbox,                'Nbox', loc)
    call write_hdf5_attr(this%OctreeOrigin,        'OctreeOrigin', loc)
    call write_hdf5_attr(this%NOctreeLevels,       'NOctreeLevels', loc)
    call write_hdf5_attr(this%MinOctreePart,       'MinOctreePart', loc)
    call write_hdf5_attr(this%MultipoleDegree,     'MultipoleDegree', loc)
    call write_hdf5_attr(this%use_dyn_layers,      'use_dyn_layers', loc)
    if(this%use_dyn_layers) then
      call write_hdf5_attr(this%NMaxOctreeLevels,  'NMaxOctreeLevels', loc)
      call write_hdf5_attr(this%LeavesTimeRatio,   'LeavesTimeRatio', loc)
    endif
    call write_hdf5_attr(this%use_pr,              'ParticlesRedistribution', loc)
    if(this%use_pr) then
      call write_hdf5_attr(this%lvl_solid,         'OctreeLevelSolid', loc)
      call write_hdf5_attr(this%part_redist_ratio, &
                                                  'ParticlesRedistributionRatio', loc)
    endif
  endif

  call write_hdf5_attr(this%use_reformulated,     'use_reformulated', loc)
  if(this%use_reformulated) then
    call write_hdf5_attr(this%f, 'f', loc)
    call write_hdf5_attr(this%g, 'g', loc)
  endif

  call write_hdf5_attr(this%debug_level,          'debug_level', loc)
  call write_hdf5_attr(this%dt_out,               'dt_out', loc)
  call write_hdf5_attr(this%basename,             'basename', loc)
  call write_hdf5_attr(this%GeometryFile,         'GeometryFile', loc)
  call write_hdf5_attr(this%ReferenceFile,        'ReferenceFile', loc)
  call write_hdf5_attr(this%restart_from_file,    'restart_from_file', loc)
  if(this%restart_from_file) then
    call write_hdf5_attr(this%restart_file,       'restart_file', loc)
    call write_hdf5_attr(this%reset_time,         'reset_time', loc)
  endif
  call write_hdf5_attr(this%hcas,                 'HCAS', loc)
  if(this%hcas) then
    call write_hdf5_attr(this%hcas_time,          'HCAS_time', loc)
    call write_hdf5_attr(this%hcas_vel,           'HCAS_velocity', loc)
  endif

end subroutine save_sim_param

!----------------------------------------------------------------------
! moved here from mod_dust_io
!> Load the time value from a result file
subroutine load_time(filename, time)
  character(len=*), intent(in) :: filename
  real(wp), intent(out)        :: time

  integer(h5loc)               :: floc

  call open_hdf5_file(filename, floc)
  call read_hdf5(time,'time',floc)
  call close_hdf5_file(floc)

end subroutine load_time

end module mod_sim_param
