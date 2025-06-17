.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   
.. _sec:nemo50:

NEMO 5.0.(1) + XIOS 2.5/3.0
===========================

Tested with

* ``gcc11.4`` on a linux system

Largely similar NEMO 4.2. In this case what happened was I updated my Linux operating system and the compiler also updated. So I am using newer libraries:

.. code-block:: none
  
  mpich-4.3.0
  zlib-1.3.1
  hdf5-1.12.1
  netcdf-c-4.9.3
  netcdf-fortran-4.6.1
  
The thing I did do different here is ``--enable-parallel`` and remove ``--enable-cxx`` for HDF5. XIOS 2.5 probably would still work with serial build, but I found it was bugging out with something in XIOS 3.0 that I couldn't get over (``nc_def_var_filter`` not declared... I don't know the exact reason for the problem so I don't know where to look to fix...)

.. warning::

  When I built the libraries above I skipped some of the checks, because it was bugging out with the MPI tests. Some of it seems to be inconsistency of HDF5 testing protocol with newer MPI behaviour (in the HDF5, preceeded with ``NPROCS=4`` flag). The resulting NEMO outputs seem fine, but proceed at own risk is skipping checks...!
  
XIOS
----

Behaviour between XIOS 2.5 (?; r2628) and XIOS 3.0 (r2763) in the compilation are largely the same except for that error in XIOS 3.0 above.
  
Do the following:

.. code-block:: none

  mkdir XIOS
  cd XIOS
  svn checkout -r 2628 http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-2.5 xios-2.5
  
To get XIOS to compile, the compilers and packages need to be pointed to first, via modifying files in ``arch``. Since I am using ``gcc``, I did the following just to make a fresh copy:

.. code-block:: none

  cd xios2.5/arch
  cp arch-GCC_LINUX.env arch-GCC_local.env
  cp arch-GCC_LINUX.fcm arch-GCC_local.fcm
  cp arch-GCC_LINUX.path arch-GCC_local.path
  
The ``*.env`` file specifies where HDF5 and NetCDF4 libraries live. The ``*.fcm`` file specifies which compilers and options to use. The ``*.path`` file specifies which paths and options to include. My files look like the following:

.. code-block:: none

  # arch-GCC_local.env

  export HDF5_INC_DIR=/usr/local/include       # CHANGE ME
  export HDF5_LIB_DIR=/usr/local/lib           # CHANGE ME

  export NETCDF_INC_DIR=/usr/local/include     # CHANGE ME
  export NETCDF_LIB_DIR=/usr/local/lib         # CHANGE ME
  
You could get an idea where the HDF5 and NetCDF4 directories are by doing ``which h5copy`` and ``which nc-config`` (assuming these are on ``$PATH``), which should give you a ``directory/bin``, and it is the ``directory`` part you want. If you did install the libraries somewhere else as in :ref:`other packages <sec:other-pack>`, say, then make sure the ``which`` commands are pointing to the right place.

.. code-block:: none

  # arch-GCC_local.fcm

  ################################################################################
  ###################                Projet XIOS               ###################
  ################################################################################

  %CCOMPILER      /usr/local/bin/mpicc                # CHANGE ME (if building own MPI bindings)
  %FCOMPILER      /usr/local/bin/mpif90               # CHANGE ME
  %LINKER         /usr/local/bin/mpif90               # CHANGE ME

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

  %CPP            cpp                            # CHANGE ME
  %FPP            cpp -P                         # CHANGE ME
  %MAKE           make
  
Check the MPI locations and versions by doing ``which mpicc`` and ``mpicc --version`` say. If they are the right ones you could just have ``mpicc`` instead of the full path as given above. MPI bindings are used here to avoid a possible error that may pop up in relation to the build trying to find ``mpi.h``. The ``gmake`` command was swapped out by the ``make`` command (I don't have ``cmake`` on the laptop).

.. note ::
  
  I added ``-D_GLIBCXX_USE_CXX11_ABI=0`` and ``-std=c++11`` to ``%BASE_CFLAGS`` for reasons documented in another page.

.. code-block:: none

  # arch-GCC_local.path

  NETCDF_INCDIR="-I$NETCDF_INC_DIR"
  NETCDF_LIBDIR="-Wl,'--allow-multiple-definition' -L$NETCDF_LIB_DIR"
  NETCDF_LIB="-lnetcdff -lnetcdf"

  MPI_INCDIR=""
  MPI_LIBDIR=""
  MPI_LIB=""

  HDF5_INCDIR="-I$HDF5_INC_DIR"
  HDF5_LIBDIR="-L$HDF5_LIB_DIR"
  HDF5_LIB="-lhdf5_hl -lhdf5 -lhdf5 -lz"

The above has all the OASIS (the atmosphere / ocean coupler) keys removed. I added the ``-Wl,'--allow-multiple-definition'`` key for reasons I don't remember anymore...

For a HDF5 build with parallel capabilities, I did not need to modify anything in ``bld.cfg`` (for serial builds I changed ``src_netcdf`` to ``src_netcdf4``).

Now it should be ready to compile. Assuming the current directory is ``xios2.5/arch``:

.. code-block:: none

  cd ../
  ./make_xios --full --prod --arch GCC_local -j2 |& tee compile_log.txt
  
The ``-j2`` option uses two processors to build. The ``tee`` command is to keep logs of potential errors (the ``|&`` is short for ``2>&1 |``) for debugging errors that may arise.


NEMO 5.0 (also for 5.0.1 tag)
-----------------------------

I checked out NEMO with (change to ``5.0.1`` tag as appropriate)

.. code-block :: bash

  git clone --branch 5.0 https://forge.nemo-ocean.eu/nemo/nemo.git nemo5.0
  
A similar procedure to specify compilers and where XIOS lives needs to be done for NEMO. Again, because of the compilers I am using:

.. code-block :: bash
  
  cd nemo5.0/arch
  cp arch-linux_gfortran.fcm ./gfortran_local.fcm
  
None of the fcm files associated with gfortran actually worked for me out of the box so here is my build of it (click :ref:`HERE <sec:nemo-fcm-log>` for a detailed log of how I got to the following):

.. code-block :: none

  # gfortran_local.fcm
  
  # generic gfortran compiler options for linux
  # NCDF_INC    netcdf include file
  # NCDF_LIB    netcdf library
  # FC          Fortran compiler command
  # FCFLAGS     Fortran compiler flags
  # FFLAGS      Fortran 77 compiler flags
  # LD          linker
  # LDFLAGS     linker flags, e.g. -L<lib dir> if you have libraries in a
  # FPPFLAGS    pre-processing flags
  # AR          assembler
  # ARFLAGS     assembler flags
  # MK          make
  # USER_INC    additional include files for the compiler,  e.g. -I<include dir>
  # USER_LIB    additional libraries to pass to the linker, e.g. -l<library>

  %NCDF_HOME           /usr/local                                        # CHANGE ME

  %XIOS_HOME           /home/julian/testing/gcc4.9-builds/XIOS/xios-2.5  # CHANGE ME

  %CPP	               cpp                                               # CHANGE ME
  %CPPFLAGS            -P -traditional

  %XIOS_INC            -I%XIOS_HOME/inc
  %XIOS_LIB            -L%XIOS_HOME/lib -lxios

  %NCDF_INC            -I%NCDF_HOME/include
  %NCDF_LIB            -L%NCDF_HOME/lib -lnetcdf -lnetcdff -lstdc++
  %FC                  mpif90                                            # CHANGE ME
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

Go into the configuration folder by

.. code-block :: bash
  
  cd ../cfgs

I tend to add ``key_nosignedzero`` into the configuration `cpp-*` file. If XIOS3 wants to be used, add in ``key_xios3``. Compile with

.. code-block :: bash
  
    ../makenemo -r GYRE_PISCES -m gfortran_local -j 2 -v 1

where the ``-v 1`` flag turns on verbosity and outputs things to screen.

XIOS 2.5 behaves as I expect, but XIOS 3.0 has some quirks I haven't figured out; it's probably some extra options in XIOS3 that I don't know about yet.
