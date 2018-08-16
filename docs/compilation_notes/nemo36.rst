.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

NEMO 3.6 (stable) + XIOS 1.0
============================

System settings: ``gcc-4.9``, Ubuntu 16.04, ``bashrc`` as configured in the page
one level up.

The assumption here is that the compiler is fixed and the packages (e.g.,
NetCDF4 and a MPI bindings) are configured to be consistent with the compilers.
See :ref:`here <sec:other-pack>` to check whether the binaries exist, where they
are, and how they might be installed.

XIOS 1.0 (svn v703)
-------------------

To use NEMO you probably do need `XIOS <http://forge.ipsl.jussieu.fr/ioserver>`_
to do the I/O. The instructions follow the one given in the `XIOS instructions
<http://forge.ipsl.jussieu.fr/ioserver/wiki/documentation>`_ with any errors
that arise. Here XIOS1.0 is used with NEMO3.6 for compatibility reasons.

For the purposes here I created a folder called ``XIOS`` and used ``svn`` to get
XIOS1.0 (which is going to be ``XIOS/xios1.0``):

.. code-block:: none

  mkdir XIOS
  cd XIOS
  svn co http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-1.0@703
  
To get XIOS to compile, the compilers and packages need to be pointed to first,
via modifying files in ``arch``. Since I am using ``gcc``, I did the following
just to make a fresh copy:

.. code-block:: none

  cd xios1.0/arch
  cp arch-GCC_LINUX.env arch-GCC_local.env
  cp arch-GCC_LINUX.fcm arch-GCC_local.fcm
  cp arch-GCC_LINUX.path arch-GCC_local.path
  
The ``*.env`` file specifies where HDF5 and NetCDF4 binaries live. The ``*.fcm``
file specifies which compilers and options to use. The ``*.path`` file specifies
which paths and options to include. My files look like the following:

.. code-block:: none

  arch-GCC_local.env (don't copy this line)

  export HDF5_INC_DIR=/usr/local/include
  export HDF5_LIB_DIR=/usr/local/lib

  export NETCDF_INC_DIR=/usr/local/include
  export NETCDF_LIB_DIR=/usr/local/lib
  
Here my HDF5 and NetCDF4 binaries are in ``/usr/local/`` (because I have
``sudo`` access). If they are somewhere else then specify another path (see
:ref:`other packages <sec:other-pack>` page if you need to install binaries
elsewhere for whatever reason).

.. code-block:: none

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
  
Here I have again my MPI things are in ``/usr/local``, with my ``mpicc`` and
``mpif90`` bound to ``gcc-4.9`` and ``gfortran-4.9`` respectively (check with
the command ``mpicc --version`` say). MPI bindings are used here to avoid a
possible error that may pop up in relation to the build trying to find
``mpi.h``. The ``gmake`` command was swapped out by the ``make`` command (I
don't have ``cmake``).

.. code-block:: none

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

.. code-block:: none

  cd ../
  ./make_xios --full --prod --arch GCC_local -j2 |& tee compile_log.txt
  
The ``-j2`` option uses two processors to build. The ``tee`` command is to keep
logs of potential errors (the ``|&`` is short for ``2>&1 |``) for debugging the
compiler issues that may arise.

.. note ::

  The following error may show up:
  
  .. code-block:: none
  
    /home/julian/testing/nemo-6800/xios-703/xios-1.0/inc/netcdf.hpp:20:26: fatal error: netcdf_par.h: No such file or directory
     #  include <netcdf_par.h>
                              ^
    compilation terminated.
    fcm_internal compile failed (256)
    /home/julian/testing/nemo-6800/xios-703/xios-1.0/Makefile:1620: recipe for target 'inetcdf4.o' failed
    
  The command ``find . -type f -iname "netcdf_par.h"`` shows that there is a
  copy of the file in ``./extern/src_netcdf4/netcdf_par.h`` and it is not being
  pointed to correctly. The culprit is in ``bld.cfg``:
  
  .. code-block:: none
  
    bld::tool::cflags    %CFLAGS %CBASE_INC -I${PWD}/extern/src_netcdf -I${PWD}/extern/boost/include -I${PWD}/extern/rapidxml/include -I${PWD}/extern/blitz/include
    
  Where ``src_netcdf`` needs to be changed to ``src_netcdf4``.
  
It should work and takes around 5 mins to compile for me. The main end result is
are binaries in ``xios1.0/bin/`` which NEMO will call.

.. note ::
  
  ``xios_server.exe`` is one of the other binaries built from compiling but is
  not required for small runs on a laptop. For its use on a cluster see for
  example the instructions on the `NOCL ARCHER guide
  <https://nemo-nocl.readthedocs.io/en/latest/work_env/archer.html>`_.

NEMO 3.6 (svn v6800)
--------------------

Check out a version of NEMO. I have another folder separate to the XIOS folders
to contain the NEMO codes and binaries:

.. code-block :: bash

  mkdir NEMO
  cd NEMO
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/trunk@6800 nemo3.6-6800
  
This checks out version 6800 (NEMO 3.6) and dumps it into a folder called
``nemo3.6-6800`` (change the target path to whatever you like). A similar
procedure to specify compilers and where XIOS lives needs to be done for NEMO.
Again, because I am using the ``gcc4.9`` compilers:

.. code-block :: bash
  
  cd nemo3.6-6800/NEMOGCM/ARCH
  cp OLD/gfortran_linux.fcm ./gfortran_local.fcm
  
None of the fcm files associated with gfortran actually worked for me out of the
box so here is my build of it (click :ref:`HERE <sec:nemo-fcm-log>` for a
detailed log of how I got to the following):

.. code-block :: none

  gfortran_local.fcm (don't copy this line)
  
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

  %NCDF_HOME           /usr/local

  %XIOS_HOME           /home/julian/testing/nemo-6800/xios-703/xios-1.0

  %CPP	               cpp-4.9
  %CPPFLAGS            -P -traditional

  %XIOS_INC            -I%XIOS_HOME/inc 
  %XIOS_LIB            -L%XIOS_HOME/lib -lxios

  %NCDF_INC            -I%NCDF_HOME/include
  %NCDF_LIB            -L%NCDF_HOME/lib -lnetcdf -lnetcdff -lstdc++
  %FC	                 mpif90
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

.. note::

  Before doing the following, it might be worthwhile doing
  
  .. code-block :: bash
  
    cd ../CONFIG/
    ./makenemo -j0 -r GYRE -n GYRE_testing -m gfortran_local
    
  and then editing ``/GYRE_testing/cpp_GYRE_testing.fcm`` to add the
  ``key_nosignedzero`` key to the end. ``-j0`` doesn't do the compile but does
  the folder creation and initial file copying. See the note at the bottom of
  the page.

To compile a configuration (using the GYRE config):
  
.. code-block :: bash
  
  cd ../CONFIG/
  ./makenemo -j2 -r GYRE -n GYRE_testing -m gfortran_local |& tee compile_log.txt
  
This uses two processors, with ``GYRE`` as a reference, builds a new folder
called ``GYRE_testing``, with the specified architecture file, and outputs a
log.

.. note ::

  The ``-r GYRE`` flag here only needs to be done once to create an extra folder
  and add GYRE_testing to ``cfg.txt``. The subsequent compilations should then
  read, e.g., ``./makenemo -n GYRE_testing -m gfortran_local``.
  
Check that it does run with the following:

.. code-block :: bash

  cd GYRE_testing/EXP00
  mpiexec -n 1 ./opa
  
This may be ``mpirun`` instead of ``mpiexec``, and ``-n 1`` just runs it as a
single core process. Change ``nn_itend = 4320`` in ``nn_itend = 120`` to only
run it for 10 days (``rdt = 7200`` which is 2 hours). With all the defaults as
is, there should be some ``GYRE_5d_*.nc`` data in the folder. You can read this
with ``ncview`` (see the ncview `page
<http://cirrus.ucsd.edu/~pierce/software/ncview/index.html>`_ or, if you have
``sudo`` access, you can install it through ``sudo apt-get install ncview``),
bearing in mind that this is actually a rotated gyre configuration (see the
following `NEMO forge page
<http://forge.ipsl.jussieu.fr/nemo/doxygen/node109.html?doc=NEMO>`_ or search
for ``gyre`` in the `NEMO book
<https://www.nemo-ocean.eu/wp-content/uploads/NEMO_book.pdf>`_).

.. note ::

  My run actually crashed immediately. By looking into ``ocean.output`` and
  searching for ``E R R O R`` shows that ``key_nosignedzero`` needs to be added
  to ``/GYRE_testing/cpp_GYRE_testing.fcm``. Rebuilding with the key then works
  fine.

