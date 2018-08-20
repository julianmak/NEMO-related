#!/bin/bash

set -e

mkdir -pv build install source

export BD=$PWD

#export CC=/usr/bin/gcc-4.9
#export CXX=/usr/bin/g++-4.9
#export FC=/usr/bin/gfortran-4.9
#export F77=/usr/bin/gfortran-4.9
#export CPP=/usr/bin/cpp-4.9

echo CC  compiler = $CC
echo FC  compiler = $FC
echo F77 compiler = $F77
echo CXX compiler = $CXX
echo CPP compiler = $CPP

read -n1 -p "Continue? [y,n]" doit 
case $doit in  
  [yY]) 
    echo "continuing..."
    ;;
  *) 
    echo "set compilers as appropriate in script" 
    exit 1 
    ;;
esac

# forces the programs to look for binaries and libraries in the installation 
# folder first

export C_INCLUDE_PATH=$BD/install/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$BD/install/include:$CPLUS_INCLUDE_PATH
export LIBRARY_PATH=$BD/install/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=$BD/install/lib:$LD_LIBRARY_PATH
export PATH=$BD/install/bin:$PATH

# set up the MPI first
cd $BD/source/
wget http://www.mpich.org/static/downloads/3.0.4/mpich-3.0.4.tar.gz
cd $BD/build/
tar -xvzf $BD/source/mpich-3.0.4.tar.gz
cd mpich-3.0.4
./configure prefix=$BD/install/
make -j 2
make check install

cd $BD/source/
wget http://www.mpich.org/static/downloads/3.0.4/hydra-3.0.4.tar.gz
cd $BD/build/
tar -xvzf $BD/source/hydra-3.0.4.tar.gz
cd hydra-3.0.4
./configure prefix=$BD/install/
make -j 2
make check install

# zlib first
cd $BD/source/
wget http://www.zlib.net/zlib-1.2.11.tar.gz
cd $BD/build/
tar -xvzf $BD/source/zlib-1.2.11.tar.gz
cd zlib-1.2.11
./configure --prefix=$BD/install/
make -j 2
make check install

# HDF5
cd $BD/source
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.19/src/hdf5-1.8.19.tar.gz
cd $BD/build/
tar -xvzf $BD/source/hdf5-1.8.19.tar.gz
cd hdf5-1.8.19
./configure --disable-shared --enable-fortran --enable-cxx --prefix=$BD/install/
make -j 2
make check install

# NetCDF4
cd $BD/source/
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.4.1.1.tar.gz
cd $BD/build/
tar -xvzf $BD/source/netcdf-4.4.1.1.tar.gz
cd netcdf-4.4.1.1
./configure --enable-netcdf4 --disable-shared --prefix=$BD/install/
make -j 2
make check install

# NetCDF4 fortran
cd $BD/source/
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-fortran-4.4.4.tar.gz
cd $BD/build/
tar -xvzf $BD/source/netcdf-fortran-4.4.4.tar.gz
cd netcdf-fortran-4.4.4
./configure --disable-shared --prefix=$BD/install/
make -j 2
make check install


echo "files installed in : $BD/install/" ;;
echo "add it to the PATH variable by adding" ;;
echo " "
echo "export PATH=$BD/install:\$PATH"
echo " "
echo "to the ~/.bashrc file"
echo " "
echo "exiting ..."


