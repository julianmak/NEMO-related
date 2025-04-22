.. _sec:nemo42:

NEMO 4.2 + XIOS 2.5
==========================

Tested with

* ``gcc8.3.0`` on a computer cluster (HPC3, with in-built parallel HDF5 and NetCDF4)

The new official page is `here <https://forge.nemo-ocean.eu/nemo/nemo/-/blob/4.2.0/README.rst>`_ and `here <https://sites.nemo-ocean.io/user-guide/install.html#download-and-install-the-nemo-code>`_. Following the instruction there largely works; below details minor things I needed to fix on the particular machine I tested on.

XIOS 2.5 (svn v2462)
--------------------

According to the NEMO `install guide <https://sites.nemo-ocean.io/user-guide/install.html#download-and-install-the-nemo-code>`_ we should use the trunk of XIOS, so

.. code-block :: bash
  
  svn co http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/trunk
  
(the version I happen to get is ``v2462``). Previously I had issues with newer GCC versions that seems to have been circumvented somehow, and XIOS builds with adding ``-D_GLIBCXX_USE_CXX11_ABI=0`` and ``-std=c++11`` to the ``BASE_CFLAGS``.

.. code-block:: none

  # arch-HKUST_HPC3.fcm

  ################################################################################
  ###################                Projet XIOS               ###################
  ################################################################################

  %CCOMPILER      mpicc                # CHANGE ME
  %FCOMPILER      mpif90               # CHANGE ME
  %LINKER         mpif90               # CHANGE ME

  %BASE_CFLAGS    -ansi -w -D_GLIBCXX_USE_CXX11_ABI=0 -std=c++11
  %PROD_CFLAGS    -O3 -DBOOST_DISABLE_ASSERTS
  %DEV_CFLAGS     -g -O2 
  %DEBUG_CFLAGS   -g 

  %BASE_FFLAGS    -D__NONE__ -ffree-line-length-none
  %PROD_FFLAGS    -O3
  %DEV_FFLAGS     -g -O2
  %DEBUG_FFLAGS   -g 

  %BASE_INC       -D__NONE__
  %BASE_LD        -lstdc++

  %CPP            cpp                             # CHANGE ME
  %FPP            cpp -P                          # CHANGE ME
  %MAKE           make

On HPC3 the various NetCDF4 folders are build and dumped in separate folders, so something needed to be done to the ``env`` and ``path`` variable entries, as follows (making sure to load the modules accordingly):

.. code-block:: none

  # arch-HKUST_HPC3.path
  
  NETCDF_INCDIR="-I $NETCDF_INC_DIR -I $NETCDFF_INC_DIR"
  NETCDF_LIBDIR="-L $NETCDF_LIB_DIR -L $NETCDFF_LIB_DIR"
  NETCDF_LIB="-lnetcdff -lnetcdf"

  HDF5_INCDIR="-I $HDF5_INC_DIR"
  HDF5_LIBDIR="-L $HDF5_LIB_DIR"
  HDF5_LIB="-lhdf5_hl -lhdf5 -lhdf5 -lz"
  
.. code-block:: none

  # arch-HKUST_HPC3.env
  
  export HDF5_INC_DIR=/opt/ohpc/pub/libs/gnu8/openmpi3/hdf5/1.10.5/include
  export HDF5_LIB_DIR=/opt/ohpc/pub/libs/gnu8/openmpi3/hdf5/1.10.5/lib

  export NETCDF_INC_DIR=/opt/ohpc/pub/libs/gnu8/openmpi3/netcdf/4.7.1/include
  export NETCDF_LIB_DIR=/opt/ohpc/pub/libs/gnu8/openmpi3/netcdf/4.7.1/lib

  export NETCDFF_INC_DIR=/opt/ohpc/pub/libs/gnu8/openmpi3/netcdf-fortran/4.5.2/include
  export NETCDFF_LIB_DIR=/opt/ohpc/pub/libs/gnu8/openmpi3/netcdf-fortran/4.5.2/lib

I went into ``bld.cfg``, found the line
  
  .. code-block:: none
  
    bld::tool::cflags    %CFLAGS %CBASE_INC -I${PWD}/extern/src_netcdf -I${PWD}/extern/boost/include -I${PWD}/extern/rapidxml/include -I${PWD}/extern/blitz/include
    
and changed ``src_netcdf`` to ``src_netcdf4`` (see :ref:`XIOS1.0 stuff
<sec:nemo36>` for the reason). Then compile as usual:

.. code-block:: none

  cd ../
  ./make_xios --full --prod --arch GCC_local -j2 |& tee compile_log.txt

NEMO 4.2 (Git SHA ``216c746957a674552de5bf02c17d22fa37f2a0d4``)
---------------------------------------------------------------

NEMO is as of writing no longer using SVN, and managing code through Git instead. So I downloaded it by

.. code-block :: bash

  git clone https://forge.nemo-ocean.eu/nemo/nemo.git nemo_4.2.0
  
I downloaded the whole thing and then looked to switch branches. To get only the official release, add the flag ``-b 4.2.0`` (or download the whole thing and then switch using ``git switch --detach 4.2.0``). After some trial and error I basically did

.. code-block :: none

  # arch-HKUST_HPC3.fcm

  %XIOS_HOME           /scratch/PI/jclmak/XIOS_mpi/xios-2.5-r2462

  %CPP                 cpp
  %CPPFLAGS            -P -traditional

  %XIOS_INC            -I%XIOS_HOME/inc
  %XIOS_LIB            -L%XIOS_HOME/lib -lxios

  %NCDF_INC            -I/opt/ohpc/pub/libs/gnu8/openmpi3/netcdf-fortran/4.5.2/include -I/opt/ohpc/pub/libs/gnu8/openmpi3/netcdf/4.7.1/include
  %NCDF_LIB            -L/opt/ohpc/pub/libs/gnu8/openmpi3/netcdf/4.7.1/lib -L/opt/ohpc/pub/libs/gnu8/openmpi3/netcdf-fortran/4.5.2/lib -lnetcdf -lnetcdff -lstdc++
  %FC                  mpif90
  %FCFLAGS             -fdefault-real-8 -O3 -funroll-all-loops -fcray-pointer -cpp -ffree-line-length-none
  %FFLAGS              %FCFLAGS
  %LD                  %FC
  %LDFLAGS
  %FPPFLAGS            -P -C -traditional
  %AR                  ar
  %ARFLAGS             -rs
  %MK                  make
  %USER_INC            %XIOS_INC %NCDF_INC
  %USER_LIB            %XIOS_LIB %NCDF_LIB
  
and everything built fine. The tricky bit was the combination of module loads, and I went one step further and brute force pointed to the relevant ``include`` and ``lib`` folders.

.. note ::

  I had some issues with using older compilers and/or OpenMPI. XIOS will compile fine, but when compiling NEMO experiments will lead to something like
  
  .. code-block :: bash

    There is no specific subroutine for the generic 'mpi_dist_graph_create_adjacent'
    
  Hence the new test compile with newer compilers (because this was the one that already interfaces with the newer OpenMPI3).
  
The usage is as in :ref:`NEMO 4.0 <sec:nemo40>`.

.. note ::

  The `zenodo <https://zenodo.org/record/3767939>`_ repository when I went to check (``ORCA2_ICE_v4.2.tar``) for the inputs when testing ``ORCA2`` was missing stuff (e.g. ``iwd``, internal wave dissipation probably), so I just went into ``namelist_cfg`` and switched it off, and it run as usual.

