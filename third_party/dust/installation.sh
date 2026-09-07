#!/bin/bash

# Function to display usage
usage() {
  echo "Usage: $0 [--with-precice YES|NO] [--help]"
  echo
  echo "Options:"
  echo "  --with-precice  Specify whether to build Dust with preCICE support. Default is NO."
  echo "  --help          Display help message" 
  exit 1
}

# Parse arguments for the WITH_PRECICE option
WITH_PRECICE="NO"  # Default value
while [[ "$1" != "" ]]; do
    case $1 in
        --with-precice ) shift
                         WITH_PRECICE="$1"
                         ;;
        --help )         usage
                         ;;
        * )              usage
                         ;;
    esac
    shift
done

# Ensure WITH_PRECICE is either YES or NO
if [[ "$WITH_PRECICE" != "YES" && "$WITH_PRECICE" != "NO" ]]; then
    usage
fi

# Update and upgrade system
sudo apt update && sudo apt upgrade -y

# Detect the Ubuntu version (jammy, focal, etc.)
UBUNTU_VERSION=$(lsb_release -cs)

# Install essential compilers and libraries
sudo apt install -y gcc g++ gfortran cmake cmake-curses-gui liblapack-dev libblas-dev libopenblas-dev libopenblas0 libcgns-dev libhdf5-dev \
libltdl-dev libsuitesparse-dev libnetcdf-dev libnetcdf-c++4-dev python-is-python3 python3-pip python3-numpy python3-venv swig autoconf automake libtool autotools-dev

cd ..
# Install preCICE library if WITH_PRECICE is set to YES
if [ "$WITH_PRECICE" == "YES" ]; then 
  
  wget https://github.com/precice/precice/releases/download/v3.1.2/libprecice3_3.1.2_${UBUNTU_VERSION}.deb
  sudo apt install -y ./libprecice3_3.1.2_${UBUNTU_VERSION}.deb
  # Create and activate a Python virtual environment
  python3 -m venv .venv
  source .venv/bin/activate

  # Install pyprecice in the virtual environment
  pip install pyprecice

  # Clone and build the MBDyn project
  if [ ! -d "mbdyn" ]; then
    git clone https://public.gitlab.polimi.it/DAER/mbdyn.git
  fi

  cd mbdyn || exit
  git checkout develop
  git pull 

  # Run the bootstrap script for MBDyn 
  sh bootstrap.sh

  # Configure the build options for MBDyn
  ./configure --enable-netcdf --with-lapack --enable-python

  # Build and install MBDyn
  make -j8
  sudo make install
  cd ..
  
  # Add MBDyn to the system PATH by updating .bashrc
  echo 'export PATH="/usr/local/mbdyn/bin:$PATH"' >> ~/.bashrc

  # Source .bashrc to update the current session
  source ~/.bashrc

fi

# Clone and build the Dust project
if [ ! -d "dust" ]; then
  git clone https://public.gitlab.polimi.it/DAER/dust.git
fi

DUST_DIR=$(pwd)/dust  # Store the Dust folder path
cd dust || exit
git checkout develop
git pull 
# Check if the build directory exists, if not, create it
if [ ! -d "build" ]; then
  mkdir build
fi

# Set WITH_PRECICE option in cmake based on the input
cd build
if [ "$WITH_PRECICE" == "YES" ]; then
  cmake -DWITH_PRECICE=YES ..
else
  cmake ..
fi

sudo make install -j8
cd ../..

# Add Dust utils/adapter to PYTHONPATH by updating .bashrc
echo "export PYTHONPATH=\"$DUST_DIR/utils/adapter:\$PYTHONPATH\"" >> ~/.bashrc

# Source .bashrc to update the current session
source ~/.bashrc
