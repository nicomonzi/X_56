# Binary packages
Pre-compiled binary packages are available for Ubuntu 20.04 (focal) and 22.04 (jammy). They include only the standalone DUST version, without the preCICE-MBDyn coupling.

For Ubuntu, you can download and install DUST as follows:
- Download the related package from the [release page](https://public.gitlab.polimi.it/DAER/dust/-/releases)
- Move to the directory where you download the package and run
  ```bash
    sudo apt install ./dust_0.7.2b-1_amd64_focal.deb
  ```
Change `focal` to `jammy` in the line above according to your version.

# Run the installation.sh script (experimental)
  ```bash 
  sed -i 's/\r$//' installation.sh && chmod +x installation.sh
  ```
  ```bash 
  ./installation.sh
  ```
By default only Dust will be installed. Otherwise, if you want Dust coupled with MBDyn through preCICE run the installation script as:
  ```bash 
  ./installation.sh --with-precice YES 
  ```

# Build DUST from source
If you need the coupled preCICE-MBDyn version of DUST, or if you wish to develop the software, you have to build it from source.
## Requirements:

To build DUST are required:

- a Fortran compiler
- CMake 
- a Lapack/BLAS implementation
- an HDF5 library
<br/>supported and tested versions: <ins>1.10</ins> , <ins> 1.12</ins><br/>
- a CGNS library
<br/>supported and tested versions: <ins>3.4.0</ins> , <ins>4.3.0</ins><br/>



A Fortran compiler, CMake and a Lapack implementation can be found pre-packed 
in most Linux distributions.

## Ubuntu 20.04 or above

<details>
  <summary markdown="span"><b>Compilers</b></summary>

  ```bash
  sudo apt install gcc g++ gfortran
  ```
</details>

<details>
  <summary markdown="span"><b>Libraries</b></summary>

  ```bash
  sudo apt install cmake liblapack-dev libblas-dev libopenblas-dev libopenblas0 libcgns-dev libhdf5-dev
  ```

</details>

<details>
  <summary markdown="span"><b>Installation</b></summary>
  
- Create a build folder inside this folder (can be "build" or anything else) and move into it:

  ```bash
  mkdir build && cd build
  ```

- **Configure** cmake with standard options:

  ```bash
  cmake -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE -DWITH_PRECICE=$WITH_PRECICE ../
  ```
  where:
  - **$CMAKE_BUILD_TYPE** can be **Release** or **Debug**
  - **$WITH_PRECICE** can be **YES** or **NO**

  For example: 

  ```bash
  cmake -DCMAKE_BUILD_TYPE=Release -DWITH_PRECICE=NO ../
  ```

- If you want to compile with intel MKL libraries first (optional):
  ```bash
  source /opt/intel/oneapi/setvars.sh
  ```
  then: 
  ```bash
  cmake -DDUST_MKL=YES ../ 
  ``` 
  
- **Build DUST**:

  ```bash
  make
  ```

- Install DUST (with root privileges if needed):

  ```bash
  sudo make install
  ```
  The default install folder should be /usr/local/bin

  Other install folders can be set by setting

  ```bash
  cmake -D CMAKE_INSTALL_PREFIX=/path/to/install/folder ../
  ```

</details>

## Apple Silicon

<details>
  <summary markdown="span"><b>Compilers</b></summary>

  ```bash
  brew install install gcc g++ gfortran
  ```

</details>

<details>
  <summary markdown="span"><b>Libraries</b></summary>

  ```bash
  brew install cmake lapack blas openblas cgns hdf5 llvm
  ```

Then export llvm libraries:

  ```bash
  export LDFLAGS="-L/opt/homebrew/opt/llvm/lib"
  export LDFLAGS="-I/opt/homebrew/opt/llvm/include" 
  echo 'export PATH="/opt/homebrew/opt/llvm/bin:$PATH"' >> ~/.zshrc 
  ```

</details>
<details>
  <summary markdown="span"><b>Installation</b></summary>
  
- Create a build folder inside this folder (can be "build" or anything else) and move into it:

  ```bash
  mkdir build && cd build
  ```
- CGNS folder is located in /opt/homebrew, therefore the cmake command will be:

  ```bash
  cmake -DCMAKE_BUILD_TYPE=Release -DWITH_PRECICE=NO -DCGNS_LIB=/opt/homebrew/lib -DCGNS_INC=/opt/homebrew/include ../
  ```
- Build DUST:

  ```bash
  make
  ```

- Install DUST (with root privileges if needed):

  ```bash
  sudo make install
  ```
</details>

## Coupling with preCICE-MBDyn

Compile DUST with **$WITH_PRECICE**=**ON** and include the MBDyn Adapter and Interface to your Python path.
The adapters are located in the *utils/adapter* folder. To add to the system path, for example, add these line to your ~/.bashrc file:

  ```bash
  export PYTHONPATH="/path/to/dust/utils/adapter":$PYTHONPATH
  ```

<details>
  <summary markdown="span">preCICE</summary>

#### preCICE
For some systems, preCICE is available in form of a pre-build package or a package recipe.
Visit:

<https://precice.org/installation-packages.html>

After installing the library, you need to get the python binding for preCICE with:

  ```bash
  $ pip install pyprecice
  ```

</details>

<details>
  <summary markdown="span">MBDyn</summary>

#### MBDyn
Visit <https://www.mbdyn.org/Installation.html>. 

MBDyn must be compiled on branch **develop**. During the configuration phase enable the following options:
 ```bash
  ./configure --enable-netcdf --with-lapack --enable-python
  ```
</details>

