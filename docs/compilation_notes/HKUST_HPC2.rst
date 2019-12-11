.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   
.. _sec:hkust:

HKUST HPC2 compilation
======================

.. warning::

  To add: SLURM submission script, post processing script

The build uses NEMO 3.7/4.0 + XIOS 2.0 as the example. For installing other
versions, extrapolate from the other notes.

HKUST HPC2 is a cluster with SLURM. Modules are loaded on an per basis via
sourcing some shell scripts. The following is going to use ``gcc`` and
``openmpi``, but in theory the corresponding intel compilers should work too
(not tested):

.. code-block:: bash

  source /usr/local/setup/gcc-g++-4.9.2.sh
  source /usr/local/setup/openmpi-2.0.0.sh

The notes are **(psuedo)-chronological** (complete with errors) rather than
the final product to highlight some pitfalls and workarounds to do with HDF5 and
NetCDF4 compatibility (the system itself does not have parallel HDF5 or NetCDF4
and it was a mystery which compiler the library was built with).

XIOS (1st try that doesn't quite work)
--------------------------------------

.. warning::

  Doing whatever is detailed here in this subsection will get XIOS compiled, but
  then it turns out when compiling NEMO that the system NetCDF4 is incompatible
  with the chosen compiler (I still have no idea which compiler was used for the
  system NetCDF4). The final working solution is to compile (a much more
  up-to-date) HDF5 and NetCDF4 separately; this means the final
  ``arch-HKUST_HPC2.env`` will be different.

I did the usual things of downloading XIOS and copying the arch files in

.. code-block:: bash

  cd $PI_HOME # <--- this is the "work" directory (which is generically not ~/)
  mkdir XIOS
  cd XIOS
  svn co http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/trunk@1322 xios-2.0
  cd xios2.0/arch
  cp arch-GCC_LINUX.env arch-HKUST_HPC2.env
  cp arch-GCC_LINUX.fcm arch-HKUST_HPC2.fcm
  cp arch-GCC_LINUX.path arch-HKUST_HPC2.path
  
HDF5 and NetCDF4 is not included in a load so it is there by default. Doing for
example ``locate libnetcdf`` tells me that I should have the following:

.. code-block:: none

  # arch-HKUST_HPC2.env

  export HDF5_INC_DIR=/usr/include
  export HDF5_LIB_DIR=/usr/lib64

  export NETCDF_INC_DIR=/usr/include
  export NETCDF_LIB_DIR=/usr/lib64
  
Following this I did

.. code-block:: none

  export LD_LIBRARY_PATH="/usr/lib64:$LD_LIBRARY_PATH"
  
as it seems to help the programs find the libraries.
  
I must have done something a bit weird to the ``arch-HKUST_HPC2.path``:

.. code-block:: none

  NETCDF_INCDIR="-I$NETCDF_INC_DIR"
  #NETCDF_LIBDIR='-Wl,"--allow-multiple-definition" -Wl,"-Bstatic" -L$NETCDF_LIB_DIR'
  NETCDF_LIBDIR='-Wl,"--allow-multiple-definition" -L$NETCDF_LIB_DIR'
  NETCDF_LIB="-lnetcdff -lnetcdf"

  HDF5_LIBDIR="-L$HDF5_LIB_DIR"
  HDF5_LIB="-lhdf5_hl -lhdf5 -lz"

Not sure where I got the ``-Bstatic`` flag from initially (maybe from the ARCHER
compilation). If that flag is there when doing the compiling then I get the
error

.. code-block:: bash

  ### ERROR ###
  linker error: ld cannot locate lnetcdf etc.
  
but doing something like ``ld [-L/usr/lib64] -lnetcdf --verbose`` or using
whatever the ``ld`` is actually called because of the modified ``$PATH`` clearly
shows success. The same happens when the intel compilers are used. Anyway, using
the following (the system had ``gmake`` so I left it; ``make`` should work too)

.. code-block:: none

  # arch-HKUST_HPC2.fcm

  ################################################################################
  ###################                Projet XIOS               ###################
  ################################################################################

  %CCOMPILER      mpicc
  %FCOMPILER      mpif90
  %LINKER         mpif90  

  %BASE_CFLAGS    -ansi -w
  %PROD_CFLAGS    -O3 -DBOOST_DISABLE_ASSERTS
  %DEV_CFLAGS     -g -O2 
  %DEBUG_CFLAGS   -g 

  %BASE_FFLAGS    -D__NONE__ -ffree-line-length-none 
  %PROD_FFLAGS    -O3
  %DEV_FFLAGS     -g -O2
  %DEBUG_FFLAGS   -g 

  %BASE_INC       -D__NONE__
  %BASE_LD        -lstdc++

  %CPP            cpp
  %FPP            cpp -P
  %MAKE           gmake

followed by

.. code-block:: bash

  cd ../
  [CPPFLAGS=-I/usr/include LDFLAGS=-L/usr/lib64] ./make_xios --full --prod --arch HKUST_HPC2 -j4 |& tee compile_log.txt

seems to do the job. I think I did go into ``bld.cfg`` and changed
``src_netcdf`` to ``src_netcdf4`` for safety; don't remember needing this in
ARCHER (did need it when doing a local compilation).

NEMO (1st try that doesn't quite work)
--------------------------------------

.. warning::

  Again this doesn't quite work because of NetCDF4 Fortran compiler
  compatibility. The final working ``arch-HKUST_HPC2.fcm`` has a modified
  ``%NCDF_INC`` and ``%NCDF_LIB``.

As advertised, when doing the following

.. code-block:: bash

  cd $PI_HOME
  mkdir NEMO
  cd NEMO
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/NEMO/trunk@8666 nemo3.7-8666
  cd nemo3.7-8666/NEMOGCM/ARCH
  cp OLD/arch-gfortran_linux.fcm ./arch-HKUST_HPC2.fcm
  
using
  
.. code-block :: none

  # arch-HKUST_HPC2.fcm
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

  %XIOS_HOME           /home/share/jclmak/XIOS/xios-2.0

  %CPP                 cpp
  %CPPFLAGS            -P -traditional

  %XIOS_INC            -I%XIOS_HOME/inc
  %XIOS_LIB            -L%XIOS_HOME/lib -lxios

  %NCDF_INC            -I/usr/include
  %NCDF_LIB            -L/usr/lib64 -lnetcdf -lnetcdff -lstdc++
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
  
When building with

.. code-block:: bash

  cd ../CONFIG
  ./makenemo -r GYRE_PISCES -n GYRE_testing -m HKUST_HPC2 -j0
  nano GYRE_testing/cpp_GYRE_testing.fcm # (have key_top -> key_nosignedzero)
  ./makenemo -n GYRE_tesitng -m HKUST_HPC2 -j4
  
throws up the error that NetCDF4 being called was built with a different
gfortran compiler. So the workaround here is build the dependencies
separately...

zlib, HDF5 and NetCDF4
----------------------

I have not figured out how to get the parallel builds of HDF5 and NetCDF4 done
successfully. Without it NEMO still works fine it just means each processor
spits out the data associated with the tile it is assigned to: the ``one_file``
option in ``file_def_nemo.xml`` doesn't work without parallel NetCDF4 and only
``multiple_file`` is allowed (it will crash the first time step it tries to
write). The workaround here is to at the post-processing stage rely on the NEMO
``TOOLS/REBUILD_NEMO`` to recombine the files if required.

I built everything as follows (see :ref:`here <sec:other-pack>` for more details
on the commands maybe):

.. warning::
  ``LD_LIBRARY_FLAG`` definitely does not point to ``/usr/lib64`` now, though I
  don't remember if I strictly needed to set it to ``$PI_HOME/custom_libs/lib``

.. code-block:: bash

  ### initialise
  cd $PI_HOME
  mkdir custom_libs
  cd custom_libs
  mkdir sources
  cd sources
  
  # zlib
  wget http://www.zlib.net/zlib-1.2.11.tar.gz
  tar -xvzf $BD/source/zlib-1.2.11.tar.gz
  cd zlib-1.2.11
  CFLAGS=-fPIC ./configure --prefix=$PI_HOME/custom_libs # -fPIC for shared libraries
  make -j 4
  make check install
  
  # HDF5
  cd $PI_HOME/custom_libs/sources
  wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.19/src/hdf5-1.8.19.tar.gz
  tar -xvzf $BD/source/hdf5-1.8.19.tar.gz
  cd hdf5-1.8.19
  CPPFLAGS=-I$PI_HOME/custom_libs/include LDFLAGS=-L$PI_HOME/custom_libs/lib \
    CFLAGS=-fPIC ./configure --enable-shared --enable-fortran --prefix=$PI_HOME/custom_libs
  make -j 4
  make check install # <---- this step takes a while
  
  # NetCDF (C)
  cd $PI_HOME/custom_libs/sources
  wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.4.1.1.tar.gz
  tar -xvzf $BD/source/netcdf-4.4.1.1.tar.gz
  cd netcdf-4.4.1.1
  CPPFLAGS=-I$PI_HOME/custom_libs/include LDFLAGS=-L$PI_HOME/custom_libs/lib \
    ./configure --enable-netcdf4 --enable-shared --prefix=$PI_HOME/custom_libs
  make -j 4
  make check install # <---- this step takes a while
  
  # NetCDF (Fortran)
  cd $PI_HOME/custom_libs/sources
  wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-fortran-4.4.4.tar.gz
  tar -xvzf $BD/source/netcdf-fortran-4.4.4.tar.gz
  cd netcdf-fortran-4.4.4
  CPPFLAGS=-I$PI_HOME/custom_libs/include LDFLAGS=-L$PI_HOME/custom_libs/lib \
    ./configure --enable-shared --prefix=$PI_HOME/custom_libs
  make -j 4
  make check install
  
My written notes says I made sure ``LD_LIBRARY_PATH`` pointed to
``$PI_HOME/custom_libs/libs`` for the NetCDF4-fortran ``./configure`` part.

Building XIOS and NEMO again
----------------------------

I rebuilt XIOS after changing ``arch-HKUST_HPC2.env`` to (probably added to
``LD_LIBRARY_PATH``):

.. code-block:: none

  # arch-HKUST_HPC2.env

  export HDF5_INC_DIR=$PI_HOME/custom_libs/include
  export HDF5_LIB_DIR=$PI_HOME/custom_libs/lib

  export NETCDF_INC_DIR=$PI_HOME/custom_libs/include
  export NETCDF_LIB_DIR=$PI_HOME/custom_libs/lib
  
For the NEMO part, ``arch-HKUST_HPC2.fcm`` now has the following:

.. code-block:: none

  %NCDF_INC            -I/$PI_HOME/custom_libs/include
  %NCDF_LIB            -L$PI_HOME/custom_libs/lib -lnetcdf -lnetcdff -lstdc++
  
Then finally everything works. I'm going to make use of the NEMO
``TOOLS/REBUILD_NEMO`` to have a single NetCDF file so I additionally do the
following (starting from the ``CONFIG`` folder):

.. code-block:: bash

  cd ../TOOLS
  ./maketools -n REBUILD_NEMO -m HKUST_HPC2
  
which results in a ``TOOLS/REBUILD_NEMO/rebuild_nemo.exe`` that I am going to
use in my post-processing script later.

Running NEMO on the HPC2
------------------------

The system uses SLURM and the key commands are

* ``sbatch [submit_nemo]``: submits the job detailed in ``submit_nemo`` (see below) 
* ``scancel [job ID]``: cancel the job
* ``sinfo``: check status of queues available
* ``squeue -u $USER``: check job info for ``$USER``

``sbatch`` could be used with arguments but I am going to have everything within
``submit_nemo`` itself. The generic one I use is as follows (I have some ASCII
art in there because I got bored at some point):

.. code-block:: bash

  #!/bin/bash

  # NOTE: Lines starting with "#SBATCH" are valid SLURM commands or statements,
  #       while those starting with "#" and "##SBATCH" are comments.  Uncomment
  #       "##SBATCH" line means to remove one # and start with #SBATCH to be a
  #       SLURM command or statement.

  #===============================================================
  # DEFINE SOME JUNK FOR THE SUBMISSION (??? make this more flexible with e.g. queues?)
  #===============================================================

  #SBATCH -J gyre04       # job name 
  #SBATCH -o stdouterr    # output and error file name
  #SBATCH -n 24           # total number of mpi tasks requested
  #SBATCH -N 1            # total number of nodes requested
  #SBATCH -p ssci         # queue (partition) -- standard, development, etc.
  #SBATCH -t 12:00:00     # maximum runtime

  # Enable email notificaitons when job begins and ends, uncomment if you need it
  ##SBATCH --mail-user=user_name@ust.hk #Update your email address
  ##SBATCH --mail-type=begin
  ##SBATCH --mail-type=end

  # Setup runtime environment if necessary
  # For example, setup MPI environment
  source /home/jclmak/nemo_env.sh
  # or you can source ~/.bashrc or ~/.bash_profile

  #===============================================================
  # LAUNCH JOB
  #===============================================================

  echo " _ __   ___ _ __ ___   ___         "
  echo "| '_ \ / _ \ '_ ' _ \ / _ \        "
  echo "| | | |  __/ | | | | | (_) |       "
  echo "|_| |_|\___|_| |_| |_|\___/  v3.7  "

  # Go to the job submission directory and run your application
  cd $PI_HOME/NEMO/nemo3.7-8666/NEMOGCM/CONFIG/GYRE_testing/EXP00/
  mpirun -n 24 ./opa

  #===============================================================
  # POSTPROCESSING
  #===============================================================

  # kills the daisy chain if there are errors

  if grep -q 'E R R O R' ocean.output ; then

    echo "E R R O R found, exiting..."
    echo "  ___ _ __ _ __ ___  _ __  "
    echo " / _ \ '__| '__/ _ \| '__| "
    echo "|  __/ |  | | | (_) | |    "
    echo " \___|_|  |_|  \___/|_|    "
    echo "check out ocean.output or stdouterr to see what the deal is "

    exit
  else 
    echo "going into postprocessing stage..."
    # cleans up files, makes restarts, moves files, resubmits this pbs

    bash ./postprocess.sh
    exit
  fi

Here because I am not using ``xios_server.exe`` I don't strictly need the ``-n
24`` after ``mpirun`` (it will then just use however many cores that's given in
``#SBATCH -n``). Maybe see the :ref:`Oxford ARC <sec:oxford>` one to see how it
might work when ``xios_server.exe`` is run alongside NEMO to do the I/O (see why
you might want to do this on the `NEMO page <https://www.nemo-ocean.eu/framework/components/interfaces/>`_.
