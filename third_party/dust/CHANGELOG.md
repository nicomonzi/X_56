#### 0.7.2-b
- Added unsteady contribution to pressure in flowfield and probe postprocessing
- Added ubuntu precompiled packages 
- Fixed some bugs

#### 0.7.1-b

- Fixed some bugs

#### 0.7.0-b

- Added non linear vortex lattice element with .c81 table correction
- Added variable vortex radius to the simulation
- Added multiple regression test and ci/cd gitlab integration
- Added hinge loads postprocessing
- Added adaptive mesh generation for hinged surface  
- Added utils for parsing sectional load in Matlab/Octave
- Added user manual documentation
- Bugfix for ll load postprocessing
- Improved profile generation from .dat file
- Improved coupled simulation (cleanup adapter, improved velocity calculation)
- Improved cmake compilation for intel-mkl library
- Fixed some bugs

#### 0.6.1-b

- Fix trailing edge orientation for hinged wake first panel
- Fixed some bugs when compiling in debug mode
- Improved examples

#### 0.6.0-b

- Added the aeroelastic interface for communication with preCICE
- Added hinged surfaces
- Added variable wind
- Fixed some bugs

#### 0.5.17-b

- Added divergence filtering in the particles wake evolution
- Small improvements

#### 0.5.16-b

- Minor improvements on warnings and parameters

#### 0.5.15-b

- Minor improvements
- Updated examples and quick start guide

#### 0.5.14-b

- Fixed some bugs

#### 0.5.13-b

- Series of updates in lifting lines methods

#### 0.5.12-b

- Added body of revolution type of input
- Updates in the lifting lines models
- Several bugs fixed

#### 0.5.11-b

- Added angle interpolation option in parametric
- Added the possibility to specify timesteps instead of dt
- Added rudimentary regression tests
- Fixed some bugs

#### 0.5.10-b

- Enabled the use of single precision in all the code
- Enabled the interaction of particles on panels with FMM
- Possibility to suppress the trailing edge of single components
- Possibility to set penetration avoidance parameters
- Introduced a turbulent viscosity for particles evolution
- New loads computation on vortex lattice components
- Fixing of minor bugs
- Some minor optimizations

#### 0.5.9-b

- introduced pointwise parametric geometry 

#### 0.5.8-b

- additional output for sectional loads of lifting lines
- added additional output and postprocessing for aeroacoustics

#### 0.5.7-b

- Introduced turbulent viscosity
- Introduced HCAS
- Updated particles redistribution
- Some bug fixed and minor improvements

#### 0.5.6-b

- Minor bugs fixed
- More output in the lifting lines sectional loads

#### 0.5.5-b

- More efficient linear system treatment
- Joining of trailing edges
- Some bugs fixed

#### 0.5.4-b

- Minor bug fixe bug fixess 

#### 0.5.3-b

- Some improvements in the lifting lines
- Minor bugs fixed
- Better output and user interface

#### 0.5.2-a

- Minor improvements
- Some bugs fixed


#### 0.5.1-a

- Minor bugfixes
- Updated and extended examples

#### 0.5.0-a

- Introduced Fast Multipole for particles evolution
- Introduced vortex stretching and diffusion for vortex particles intensity evolution
- Completely new calculation of pressure based on Poisson solution
- Introduced vorticity effects: separation simulation

#### 0.4.1-a

- Added the possibility of defining wings and blades mirrored

#### 0.4.0-a

- Some minor features introduced
- Extensive bug fixing

#### 0.3.0-a

- Added the generation and advection of vortex particles
- Added postprocessing of such particles
- Added postprocessing features, averaging etc.
- Bug fixes

#### 0.2.1-a

- Extended models for multiple geometries (i.e. rotors)
- Added restart capabilities
- Some minor fixes
- Some documentation and examples added

#### 0.2.0-a

- Added the lifting line model for lifting surfaces, to enable table lookup of aerodynamic coefficients
- Added actuator disk model with ring wake for a simple representation of rotors
- Added extensive postprocessing, in binary tecplot, binary vtk and .dat format
  - History of loads of any defined component in any of the defined reference frames
  - History of variables in probing location
  - Volume/slice flowfield visualization
  - Surface visualization on geometry and wake


#### 0.1.0-a

- Added a more comprehensive motion interface
  - New input
  - Possibility of motions of the pole
  - Possibility of sinusoidal motions
  - Possibility of motion history from file
- Added automatic creation of multiple reference frames and components
- First parallel optimizations
