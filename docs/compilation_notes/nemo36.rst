.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

NEMO 3.6 (stable) + XIOS1.0
===========================

System settings: ``gcc-4.9``, Ubuntu 16.04, fresh ``bashrc``

The assumption here is that the compiler is fixed and the packages (HDF5,
NetCDF4 and a MPI binding) are configured to be consistent with the compilers.
See :ref:`here <sec:other-pack>` to check whether the binaries exist, where they
are, and how they might be installed.

XIOS1.0 (svn v703)
------------------

To use NEMO you probably do need `XIOS <http://forge.ipsl.jussieu.fr/ioserver>`_
to do the I/O. The instructions follow the one given in the `XIOS instructions
<http://forge.ipsl.jussieu.fr/ioserver/wiki/documentation>`_ with any errors
that arise. Here XIOS1.0 is used with NEMO3.6 for compatibility reasons.

I recommend creating a folder called ``XIOS``, going into it, and using the
following command creates a folder called ``xios1.0``:

.. code-block:: bash

  mkdir XIOS
  cd XIOS
  svn co http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-1.0@703
  
To get XIOS to compile, the compilers and packages need to be pointed to first,
via modifying files in ``arch``. Since I am using ``gcc``, I did the following
just to make a fresh copy:

.. code-block:: bash

  cd xios1.0/arch
  cp arch-GCC_LINUX.env arch-GCC_local.env
  cp arch-GCC_LINUX.fcm arch-GCC_local.fcm
  cp arch-GCC_LINUX.path arch-GCC_local.path
  
The ``*.env`` file specifies where HDF5 and NetCDF4 binaries live. The ``*.fcm``
file specifies which compilers and options to use. The ``*.path`` file specifies
which paths and options to include. My files look like the following:

.. code-block:: bash

  arch-GCC_local.env (don't copy this line)

  export HDF5_INC_DIR=/usr/local/include
  export HDF5_LIB_DIR=/usr/local/lib

  export NETCDF_INC_DIR=/usr/local/include
  export NETCDF_LIB_DIR=/usr/local/lib
  
Here my HDF5 and NetCDF4 binaries are in ``/usr/local/`` (because I have
``sudo`` access). If they are somewhere else just change it (see :ref:`other
packages <sec:other-pack>` page if you need to install binaries elsewhere for
whatever reason).

.. code-block:: bash

  arch-GCC_local.fcm (don't copy this line)

  ################################################################################
  ###################                Projet XIOS               ###################
  ################################################################################

  %CCOMPILER      /usr/local/bin/mpicc
  %FCOMPILER      /usr/local/bin/mpif90
  %LINKER         /usr/local/bin/mpif90  

  %BASE_CFLAGS    -ansi -w
  %PROD_CFLAGS    -O3 -DBOOST_DISABLE_ASSERTS
  %DEV_CFLAGS     -g -O2 
  %DEBUG_CFLAGS   -g 

  %BASE_FFLAGS    -D__NONE__ 
  %PROD_FFLAGS    -O3
  %DEV_FFLAGS     -g -O2
  %DEBUG_FFLAGS   -g 

  %BASE_INC       -D__NONE__
  %BASE_LD        -lstdc++

  %CPP            cpp-4.9
  %FPP            cpp-4.9 -P
  %MAKE           make
  
Here I have again sudo access to install my MPI things into ``/usr/local``, and
my ``mpicc`` and ``mpif90`` are configured to ``gcc-4.9`` and ``gfortran-4.9``
respectively (check with the command ``mpicc --version`` say). MPI bindings are
used here to avoid a possible error that may pop up in relation to the build
trying to find ``mpi.h``. The ``gmake`` command was swapped out by the ``make``
command (I don't have ``cmake``).

.. code-block:: bash

  arch-GCC_local.path (don't copy this line)

  NETCDF_INCDIR="-I$NETCDF_INC_DIR"
  NETCDF_LIBDIR="-Wl,'--allow-multiple-definition' -L$NETCDF_LIB_DIR"
  NETCDF_LIB="-lnetcdff -lnetcdf"

  MPI_INCDIR=""
  MPI_LIBDIR=""
  MPI_LIB=""

  HDF5_INCDIR="-I$HDF5_INC_DIR"
  HDF5_LIBDIR="-L$HDF5_LIB_DIR"
  HDF5_LIB="-lhdf5_hl -lhdf5 -lhdf5 -lz"

The above has all the OASIS (the atmopshere / ocean coupler) keys removed. (I
added the ``-Wl,'--allow-multiple-definition'`` key for reasons I don't remember
anymore...)

Now it should be ready to compile. Assuming the current directory is
``xios1.0/arch``:

.. code-block:: bash

  cd ../
  ./make_xios --full --prod --arch GCC_local -j2 |& tee compile_log.txt
  
The ``-j2`` option uses two processors to build. The ``tee`` command is to keep
logs of potential errors (the ``|&`` is short for ``2>&1 |``) for debugging the
compiler issues that may arise.

.. note ::

  The following error may show up:
  
  .. code-block:: bash
  
    /home/julian/testing/nemo-6800/xios-703/xios-1.0/inc/netcdf.hpp:20:26: fatal error: netcdf_par.h: No such file or directory
     #  include <netcdf_par.h>
                              ^
    compilation terminated.
    fcm_internal compile failed (256)
    /home/julian/testing/nemo-6800/xios-703/xios-1.0/Makefile:1620: recipe for target 'inetcdf4.o' failed
    
  If you type the command ``find . -type f -iname "netcdf_par.h"`` you will find
  that it lives in ``./extern/src_netcdf4/netcdf_par.h`` so it is not being
  pointed to correctly. The culprit is in ``bld.cfg``:
  
  .. code-block:: bash
  
    bld::tool::cflags    %CFLAGS %CBASE_INC -I${PWD}/extern/src_netcdf -I${PWD}/extern/boost/include -I${PWD}/extern/rapidxml/include -I${PWD}/extern/blitz/include
    
  Where a ``4`` needs to be added to ``src_netcdf``.
  
It should work and takes around 5 mins for me. The main end result is a binary
``xios_server.exe`` in ``xios1.0/bin/`` which NEMO will call.

NEMO3.6 (svn v6800)
-------------------


