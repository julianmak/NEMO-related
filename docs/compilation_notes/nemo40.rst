.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   
.. _sec:nemo40:

NEMO 4.0 (beta) + XIOS 2.5
==========================

The code structure in NEMO 4.0 and the use of some commands are slightly
different (at least in v9925) and will be documented below (please see the
`official NEMO annoucement
<http://forge.ipsl.jussieu.fr/nemo/wiki/Users/Agenda/2018-07-11>`_ for details).
If you get errors that are not documented here, see if :ref:`the XIOS1.0 NEMO3.6
<sec:nemo36>` page contains the relevant errors.

System settings: ``gcc-4.9``, Ubuntu 16.04, ``bashrc`` as configured in the page
one level up.

The assumption here is that the compiler is fixed and the packages (e.g.,
NetCDF4 and a MPI bindings) are configured to be consistent with the compilers.
See :ref:`here <sec:other-pack>` to check whether the binaries exist, where they
are, and how they might be installed.

The instructions below assumes ``gcc4.9`` compilers but works for ``gcc5.4``
compilers too (with one extra flag required in XIOS compilation). I additionally
defined some extra variables

.. code-block:: bash

  export $BD=/home/julian/testing/gcc4.9-builds # CHANGE ME

  export C_INCLUDE_PATH=$BD/install/include:$C_INCLUDE_PATH
  export CPLUS_INCLUDE_PATH=$BD/install/include:$CPLUS_INCLUDE_PATH
  export LIBRARY_PATH=$BD/install/lib:$LIBRARY_PATH
  export LD_LIBRARY_PATH=$BD/install/lib:$LD_LIBRARY_PATH
  
otherwise I have found the resulting libraries and binaries are not necessarily
linked to the right ones (I have a few versions of libraries at different places
as a result of the testing recorded here).

XIOS 2.5 (svn v1566)
--------------------

.. note ::

  Looks like you could use XIOS 2.0 with NEMO 4.0, so if the following doesn't
  work for you, try compiling :ref:`XIOS 2.0 <sec:nemo37>` instead.
  
Do the following:

.. code-block:: none

  mkdir XIOS
  cd XIOS
  svn co http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-2.5@1566 xios-2.5
  
To get XIOS to compile, the compilers and packages need to be pointed to first,
via modifying files in ``arch``. Since I am using ``gcc``, I did the following
just to make a fresh copy:

.. code-block:: none

  cd xios2.5/arch
  cp arch-GCC_LINUX.env arch-GCC_local.env
  cp arch-GCC_LINUX.fcm arch-GCC_local.fcm
  cp arch-GCC_LINUX.path arch-GCC_local.path
  
The ``*.env`` file specifies where HDF5 and NetCDF4 binaries live. The ``*.fcm``
file specifies which compilers and options to use. The ``*.path`` file specifies
which paths and options to include. My files look like the following:

.. code-block:: none

  # arch-GCC_local.env

  export HDF5_INC_DIR=/usr/local/include       # CHANGE ME
  export HDF5_LIB_DIR=/usr/local/lib           # CHANGE ME

  export NETCDF_INC_DIR=/usr/local/include     # CHANGE ME
  export NETCDF_LIB_DIR=/usr/local/lib         # CHANGE ME
  
You could check where the HDF5 and NetCDF4 directories are by doing ``which
h5copy`` and ``which nc-config``, which should give you a ``directory/bin``, and
it is the ``directory`` part you want. If you did install the libraries
somewhere else as in :ref:`other packages <sec:other-pack>`, say, then make sure
the ``which`` commands are pointing to the right place.

.. code-block:: none

  # arch-GCC_local.fcm

  ################################################################################
  ###################                Projet XIOS               ###################
  ################################################################################

  %CCOMPILER      /usr/local/bin/mpicc                # CHANGE ME
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

  %CPP            cpp-4.9                             # CHANGE ME
  %FPP            cpp-4.9 -P                          # CHANGE ME
  %MAKE           make
  
Check the MPI locations by doing ``which mpicc`` and ``mpicc --version`` say. If
they are the right ones you could just have ``mpicc`` instead of the full path
as given above. MPI bindings are used here to avoid a possible error that may
pop up in relation to the build trying to find ``mpi.h``. The ``gmake`` command
was swapped out by the ``make`` command (I don't have ``cmake``).

.. note ::

  For ``gcc5.4`` and maybe newer versions, doing just the above when compiling
  leads to a whole load of errors about clashing in C++:
  
  .. code-block:: bash
    
    .../include/boost/functional/hash/extensions.hpp:69:33: error: ‘template<class T, class A> std::size_t boost::hash_value’ conflicts with a previous declaration
     std::size_t hash_value(std::list<T, A> const& v)
                                 ^
  
  Adding ``-D_GLIBCXX_USE_CXX11_ABI=0`` to ``%BASE_CFLAGS`` fixes these.

  A difference I've found between XIOS 2.5 and other XIOS versions is that doing
  just the above might lead to an error like the following:
  
  .. code-block:: bash
  
    This file requires compiler and library support for the ISO C++ 2011 standard. This support is currently experimental, and must be enabled with the -std=c++11 or -std=gnu++11 compiler options.

  Adding ``-std=c++11`` to ``%BASE_CFLAGS`` seems to fix this.
  
  You might also get the following:
  
  .. code-block:: bash
  
    SUBROUTINE cxios_set_interpolate_domain_read_write_convention(interpolate_domain_hdl, read_write_convention, read_write_conventi
                                                                                                                                    1
    Error: Unexpected junk in formal argument list at (1)
    
  The Fortran lines are too long, so fix this by adding
  ``-ffree-line-length-none`` to ``%BASE_FFLAGS``.

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

The above has all the OASIS (the atmosphere / ocean coupler) keys removed. I
added the ``-Wl,'--allow-multiple-definition'`` key for reasons I don't remember
anymore...

I went into ``bld.cfg``, found the line
  
  .. code-block:: none
  
    bld::tool::cflags    %CFLAGS %CBASE_INC -I${PWD}/extern/src_netcdf -I${PWD}/extern/boost/include -I${PWD}/extern/rapidxml/include -I${PWD}/extern/blitz/include
    
and changed ``src_netcdf`` to ``src_netcdf4`` (see :ref:`XIOS1.0 stuff
<sec:nemo36>` for the reason).

Now it should be ready to compile. Assuming the current directory is
``xios2.5/arch``:

.. code-block:: none

  cd ../
  ./make_xios --full --prod --arch GCC_local -j2 |& tee compile_log.txt
  
The ``-j2`` option uses two processors to build. The ``tee`` command is to keep
logs of potential errors (the ``|&`` is short for ``2>&1 |``) for debugging the
compiler issues that may arise. It should work and takes around 5 mins to
compile for me. The main end result is are binaries in ``xios2.5/bin/`` which
NEMO will call.


NEMO 4.0 (svn v9925)
--------------------

There is a restructuring of folders (see the `official annoucement
<http://forge.ipsl.jussieu.fr/nemo/wiki/Users/Agenda/2018-07-11>`_ for details)
so the commands below will reflect this.

Check out a version of NEMO. I have another folder separate to the XIOS folders
to contain the NEMO codes and binaries:

.. code-block :: bash

  mkdir NEMO
  cd NEMO
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/trunk@9925 nemo4.0-9925
  
This checks out version 9925 (NEMO 4.0 beta) and dumps it into a folder called
``nemo4.0-9925`` (change the target path to whatever you like). A similar
procedure to specify compilers and where XIOS lives needs to be done for NEMO.
Again, because I of the compilers I am using:

.. code-block :: bash
  
  cd nemo4.0-9925/arch
  cp arch-linux_gfortran.fcm ./gfortran_local.fcm
  
None of the fcm files associated with gfortran actually worked for me out of the
box so here is my build of it (click :ref:`HERE <sec:nemo-fcm-log>` for a
detailed log of how I got to the following):

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

  %CPP	               cpp-4.9                                           # CHANGE ME
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

The main changes are (again, see :ref:`here <sec:nemo-fcm-log>` for an attempt
at the reasoning and a log of errors that motivates the changes):

* added ``%NCDF_HOME`` to point to where NetCDF lives
* added ``%XIOS_*`` keys to point to where XIOS lives
* added ``%CPP`` and flags, consistent with using ``gcc4.9``
* added the ``-lnetcdff`` and ``-lstdc++`` flags to NetCDF flags
* using ``mpif90`` which is a MPI binding of ``gfortran-4.9``
* added ``-cpp`` and ``-ffree-line-length-none`` to Fortran flags
* swapped out ``gmake`` with ``make``

Go into the configuration folder by

.. code-block :: bash
  
  cd ../cfgs

One of the things I noticed is that ``makenemo`` now seems to work slightly
differently (at least with this version). Normally you can do ``makenemo -r GYRE
-n GYRE_testing -j0 -m gcc_fortran_local``, which copies a configuration but
does not compile it, so you can edit the ``cpp`` flags before compiling (and
note that it adds an entry into ``works_cfgs.txt``). However now it seems you
have to specify a ``-r`` flag or a ``-d`` flag (which specifies what NEMO
modules the configuration should have), whereas before just a ``-n`` flag would
work by itself. 

You could just compile as usual with ``makenemo`` (see :ref:`NEMO 3.6
<sec:nemo36>` for syntax). The slightly untidy way to circumvent errors that I
know will come up was to do the following:

1. Open ``refs_cfg.txt``, copy the ``GYRE_PISCES OCE TOP`` line and paste it at the bottom, but then change the configuration name (``GYRE_PISCES`` to ``GYRE_testing`` in my case), save and close it;

2. Then do

  .. code-block :: bash
  
    mkdir GYRE_testing
    rsync -arv GYRE_PISCES/* GYRE_testing/
    
3. I opened ``/GYRE_testing/cpp_GYRE_testing.fcm`` and replaced ``key_top`` with ``key_nosignedzero`` (does not compile TOP for speed speeds, and make sure zeros are not signed), save it;

4. Compile with

  .. code-block :: bash
  
    ./makenemo -j2 -n GYRE_testing -m gfortran_local |& tee compile_log.txt
  
  which should take a few minutes. 

.. warning ::

  See if this untidy feature of ``makenemo`` has been fixed in the trunk?

Note the executable ``opa`` is now called ``nemo`` (so make sure you change
those submission scripts on the relevant clusters if you use NEMO on them).
Check that it does run with the following:

.. code-block :: bash

  cd GYRE_testing/EXP00
  mpiexec -n 1 ./nemo
  

Note that what used to be ``solver.stat`` is now called ``run.stat``, and there
is an extra ``run.stat.nc`` for whatever reason. The ``ocean.output`` file is
still the same.
