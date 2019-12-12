.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   
.. _sec:oxford:

Oxford ARC compilation
======================

.. note::

  To add: post processing script

The build uses NEMO 3.7/4.0 + XIOS 2.0 as the example. For installing other
versions, extrapolate from the other notes.

Annoyingly (!) everything basically works out of the box because all the
dependency modules have been built already! This has not been the usual
experience I have with XIOS and NEMO...

Building NEMO and XIOS
----------------------

Log on first using:

.. code-block:: bash

  ssh [-X] phys????@arcus-b.arc.ox.ac.uk
  
On doing ``module list`` we should for the first time see no modules are loaded
by default. If there are do ``module purge`` just for safety. Doing ``module
avail`` shows the list of modules available. I'm going to use the ``gcc`` one,
and by doing

.. code-block:: bash

  module load /netcdf-parallel/4.4__mvapich2__gcc
  
this loads NetCDF4 as well as its dependencies (which should be HDF5, gcc4.9.2
and the relevant mvapich. Once I've done this I did

.. code-block:: bash

  echo $LD_LIBRARY_PATH
  > /system/software/arcus-b/lib/netcdf/4.4/mvapich2-2.1.0__gcc-4.9.2/lib:/system/software/arcus-b/lib/hdf5/1.8.12/mvapich2-2.1.0__gcc-4.9.2/lib ...
  
which tells me where the NetCDF and HDF5 libraries live. Then for XIOS I do

.. code-block:: bash

  cd $DATA # <--- this is the "work" directory (which is generically not ~/)
  mkdir XIOS
  cd XIOS
  svn co http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/trunk@1322 xios-2.0
  cd xios2.0/arch
  cp arch-GCC_LINUX.env arch-GCC_ARC.env
  cp arch-GCC_LINUX.fcm arch-GCC_ARC.fcm
  cp arch-GCC_LINUX.path arch-GCC_ARC.path
  
with

.. code-block:: none

  # arch-GCC_ARC.env

  export HDF5_INC_DIR=/system/software/arcus-b/lib/hdf5/1.8.12/mvapich2-2.1.0__gcc-4.9.2/include
  export HDF5_LIB_DIR=/system/software/arcus-b/lib/hdf5/1.8.12/mvapich2-2.1.0__gcc-4.9.2/lib

  export NETCDF_INC_DIR=/system/software/arcus-b/lib/netcdf/4.4/mvapich2-2.1.0__gcc-4.9.2/include
  export NETCDF_LIB_DIR=/system/software/arcus-b/lib/netcdf/4.4/mvapich2-2.1.0__gcc-4.9.2/lib
  
and the other two as default. Running

.. code-block:: bash

  cd ../
  ./make_xios --full --prod --arch HKUST_HPC2 -j4 |& tee compile_log.txt

seems to do the job. I think I did go into ``bld.cfg`` and changed
``src_netcdf`` to ``src_netcdf4`` for safety; don't remember needing this in
ARCHER (did need it when doing a local compilation). If that doesn't work
consider adding ``CPPFLAGS`` and ``LDFLAGS`` before the ``./make_xios`` command
to force the program to look in the specified place.

NEMO is then built as follows:

.. code-block:: bash

  cd $DATA
  mkdir NEMO
  cd NEMO
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/NEMO/trunk@8666 nemo3.7-8666
  cd nemo3.7-8666/NEMOGCM/ARCH
  cp OLD/arch-gfortran_linux.fcm ./arch-GCC_ARC.fcm
  
using
  
.. code-block :: none

  # arch-GCC_ARC.fcm
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

  %XIOS_HOME           $DATA/XIOS/xios-2.0

  %CPP                 cpp
  %CPPFLAGS            -P -traditional

  %XIOS_INC            -I%XIOS_HOME/inc
  %XIOS_LIB            -L%XIOS_HOME/lib -lxios

  %NCDF_INC            -I/system/software/arcus-b/lib/netcdf/4.4/mvapich2-2.1.0__gcc-4.9.2/include
  %NCDF_LIB            -L/system/software/arcus-b/lib/netcdf/4.4/mvapich2-2.1.0__gcc-4.9.2/lib -lnetcdf -lnetcdff -lstdc++
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
  
followed by

.. code-block:: bash

  cd ../CONFIG
  ./makenemo -r GYRE_PISCES -n GYRE_testing -m GCC_ARC -j0
  nano GYRE_testing/cpp_GYRE_testing.fcm # (have key_top -> key_nosignedzero)
  ./makenemo -n GYRE_tesitng -m HKUST_HPC2 -j4
  
and it should work. One more thing we will do is to make ``TOOLS/REBUILD_NEMO``:

.. code-block:: bash

  cd ../TOOLS
  ./maketools -n REBUILD_NEMO -m GCC_ARC

Running NEMO on the ARC
-----------------------

The system uses SLURM and the key commands are

* ``sbatch [submit_nemo]``: submits the job detailed in ``submit_nemo`` (see below) 
* ``scancel [job ID]``: cancel the job
* ``sinfo``: check status of queues available
* ``squeue -u $USER``: check job info for ``$USER``

``sbatch`` could be used with arguments but I am going to have everything within
``submit_nemo`` itself. Check balance and budget account names with the
``mybalance`` command. Running ``sinfo`` shows the queue available is called
``compute``. One thing to note is that ARC has 16 cores per node and this is
reflected in the core/node request numbers.

Oxford ARC does have parallel NetCDF so I can use XIOS in detached mode. To do
this I link ``xios_server.exe`` to the folder:

.. code-block:: bash

  cd GYRE_testing/EXP00
  ln -s $DATA/XIOS/xios2.0/bin/xios_server.exe .
  
Modify ``iodef.xml`` so that the user server boolean is ``true``. Additionally I
go into ``file_def_nemo.xml`` and swap out ``multiple_file`` at the top header
to ``one_file``, which then spits out a single NetCDF file. This however only
works for the diagnostic files but not the restart files, so recombining the
restart files we are going to call ``TOOLS/REBUILD_NEMO`` in the post-processing
script.

The generic submission script I use (based on the one given on the `NOCL
page <https://nemo-nocl.readthedocs.io/en/latest/work_env/mobius.html>`_) is as
follows (I have some ASCII art in there because I got bored at some point):

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
  #SBATCH -n 32           # total number of mpi tasks requested
  #SBATCH -N 2            # total number of nodes requested
  #SBATCH -p compute      # queue (partition) -- standard, development, etc.
  #SBATCH -t 12:00:00     # maximum runtime

  # Enable email notificaitons when job begins and ends, uncomment if you need it
  ##SBATCH --mail-user=user_name@ust.hk #Update your email address
  ##SBATCH --mail-type=begin
  ##SBATCH --mail-type=end

  # Setup runtime environment if necessary
  module purge
  module load netcdf-parallel/4.4__mvapich2__gcc

  #===============================================================
  # LAUNCH JOB
  #===============================================================

  echo " _ __   ___ _ __ ___   ___         "
  echo "| '_ \ / _ \ '_ ' _ \ / _ \        "
  echo "| | | |  __/ | | | | | (_) |       "
  echo "|_| |_|\___|_| |_| |_|\___/  v3.7  "

  # Go to the job submission directory and run your application
  cd /data/phys-geometric/phys1342/NEMO/nemo3.7-8666/NEMOGCM/CONFIG/GYRE_testing/EXP00
  mpirun -n 2 ./xios_server.exe : -n 30 ./opa

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

    #bash ./postprocess.sh
    exit
  fi

The ratio of ``XIOScore`` to ``NEMOcore`` I never found to lead to major
differences for the size of runs I do (not larger than 300 cores); vaguely
remember reading somewhere that ``XIOScore`` hovering between 5 to 10 per cent
of ``NEMOcore`` is ok.
