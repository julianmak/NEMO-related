.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   
.. _sec:oxford:

Oxford ARC compilation
======================

The build uses NEMO 3.7/4.0 + XIOS 2.0 as the example. For installing other versions, extrapolate from the other notes.

Annoyingly (!) everything basically works out of the box because all the dependency modules have been built already! This has not been the usual experience I have with XIOS and NEMO...

Building NEMO and XIOS
----------------------

Log on first using:

.. code-block:: bash

  ssh [-X] phys????@arcus-b.arc.ox.ac.uk
  
On first login doing ``module list`` should show no modules; if there are then might want to do ``module purge`` just for safety. Doing ``module avail`` shows the list of modules available. I'm going to use the ``gcc`` one, and by doing

.. code-block:: bash

  module load /netcdf-parallel/4.4__mvapich2__gcc
  
this loads NetCDF4 as well as its dependencies (which should be HDF5, gcc4.9.2 and the relevant mvapich); I added that line to ``~/.bashrc`` so it loads from now on when logging in. Then I did

.. code-block:: bash

  echo $LD_LIBRARY_PATH
  > /system/software/arcus-b/lib/netcdf/4.4/mvapich2-2.1.0__gcc-4.9.2/lib:/system/software/arcus-b/lib/hdf5/1.8.12/mvapich2-2.1.0__gcc-4.9.2/lib ...
  
which tells me where the NetCDF and HDF5 libraries live. So for XIOS I do

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

seems to do the job. I think I did go into ``bld.cfg`` and changed ``src_netcdf`` to ``src_netcdf4`` for safety; don't remember needing this in ARCHER (did need it when doing a local compilation). If that doesn't work consider adding ``CPPFLAGS`` and ``LDFLAGS`` before the ``./make_xios`` command to force the program to look in the specified place.

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
  ./makenemo -n GYRE_tesitng -m GCC_ARC -j4
  
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

``sbatch`` could be used with arguments but I am going to have everything within ``submit_nemo`` itself. Check balance and budget account names with the ``mybalance`` command. Running ``sinfo`` shows the queue available is called ``compute``. One thing to note is that ARC has 16 cores per node and this is reflected in the core/node request numbers.

Oxford ARC does have parallel NetCDF so I can use XIOS in detached mode. To do this I link ``xios_server.exe`` to the folder:

.. code-block:: bash

  cd GYRE_testing/EXP00
  ln -s $DATA/XIOS/xios2.0/bin/xios_server.exe .
  
Modify ``iodef.xml`` so that the user server boolean is ``true``. Additionally I go into ``file_def_nemo.xml`` and swap out ``multiple_file`` at the top header to ``one_file``, which then spits out a single NetCDF file. This however only works for the diagnostic files but not the restart files, so recombining the restart files we are going to call ``TOOLS/REBUILD_NEMO`` in the post-processing script.

The generic submission script I use (based on the one given on the `NOCL page <https://nemo-nocl.readthedocs.io/en/latest/work_env/mobius.html>`_) is as follows (I have some ASCII art in there because I got bored at some point):

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

    bash ./postprocess.sh >& cleanup.log
    exit
  fi

The ratio of ``XIOScore`` to ``NEMOcore`` I never found to lead to major differences for the size of runs I do (not larger than 300 cores); vaguely remember reading somewhere that ``XIOScore`` hovering between 5 to 10 per cent of ``NEMOcore`` is ok.

The following post-processing script requires a few prepping (I make no apologies for the bad code and the script being fickle; feel free to modify as you see fit):

* copying the ``nn_date0`` line into ``namelist_cfg`` from say ``namelist_ref`` if it doesn't exist already, because the time-stamps are modified by modifying ``nn_date0``
* do a search in ``namelist_cfg`` and make sure there is only ever one mention of ``nn_date0`` (otherwise it grabs the wrong lines)
* ``nn_date0`` should not begin with zeros (e.g. ``10101`` rather than ``010101`` in ``yymmdd``)
*  in the experiment folder, do ``mkdir RESTARTS OUTPUTS`` (otherwise there is no folder to copy into)

The ``postprocess.sh`` I cooked up is here:

.. code-block:: bash

  #!/bin/bash
  #! postprocess.sh
  #! Script to clean up the NEMO outputs

  export BASE_DIR=$DATA/NEMO/nemo3.7-8666/NEMOGCM/
  export MODEL=GYRE
  export NUM_CPU=30

  # time-stamp increment, yymmdd
  export DATE_INC=100000

  # when to stop the daisy chaining, yymmdd
  export THRESH=10

  # error catching (only when restart files etc cannot be copied or made)
  export ERR_CATCH=0

  ########################################################
  # 0) recombine files to one netcdf (restarts and/or outputs)
  # restarts: extract the restart file time-step stamp
  #              based on the *0000.nc restart which should (!) always exist
  #           rebuild the restart file in the submission directory
  # outputs:  put them in manually and just do a grab
  #           this assumes only files at the current time-stamp is there,
  #              otherwise it will bug out as it grabs wrong files
  ########################################################

  # restart files
  export RES_TIMESTAMP=$(echo $(ls -d ${MODEL}_*_restart_0000.nc) | awk -F _ '{print $2 }')

  $BASE_DIR/TOOLS/REBUILD_NEMO/rebuild_nemo ${MODEL}_${RES_TIMESTAMP}_restart $NUM_CPU
  if (($? > 0)); then 
    ERR_CATCH=$((ERR_CATCH + 1))
    echo "  ERR: making the restart file in the folder"
  fi
  ##$BASE_DIR/TOOLS/REBUILD_NEMO/rebuild_nemo ${MODEL}_${RES_TIMESTAMP}_restart_ice $NUM_CPU

  # output files (assumes a grid_T always exists)
  #export OUT_FREQ=$(echo $(ls -d ${MODEL}_*_grid_T_0000.nc) | awk -F _ '{print $2 }')
  #export OUT_START=$(echo $(ls -d ${MODEL}_*_grid_T_0000.nc) | awk -F _ '{print $3 }')
  #export OUT_END=$(echo $(ls -d ${MODEL}_*_grid_T_0000.nc) | awk -F _ '{print $4 }')

  #$BASE_DIR/TOOLS/REBUILD_NEMO/rebuild_nemo ${MODEL}_${OUT_FREQ}_${OUT_START}_${OUT_END}_grid_T $NUM_CPU
  #$BASE_DIR/TOOLS/REBUILD_NEMO/rebuild_nemo ${MODEL}_${OUT_FREQ}_${OUT_START}_${OUT_END}_grid_U $NUM_CPU
  #$BASE_DIR/TOOLS/REBUILD_NEMO/rebuild_nemo ${MODEL}_${OUT_FREQ}_${OUT_START}_${OUT_END}_grid_V $NUM_CPU
  #$BASE_DIR/TOOLS/REBUILD_NEMO/rebuild_nemo ${MODEL}_${OUT_FREQ}_${OUT_START}_${OUT_END}_grid_W $NUM_CPU

  # add more things in here if output freqs are different etc

  ########################################################
  # 1) pull out some variables to modify namelist file
  ########################################################

  # pull the number out
  # add the increment to it for new date
  # subtract appropriately to get the date stamp 
  #   (e.g. 110101 - 8871 = 101230) and bulk out zeros

  export OLD_DATE_STR=$(grep -ri "nn_date0" namelist_cfg)
  export OLD_DATE_NUM=$(echo ${OLD_DATE_STR} | sed -e 's/[^0-9 ]//g' | awk '{print $NF}')
  export NEW_DATE_NUM=$((OLD_DATE_NUM + DATE_INC))

  # 8871 for 30 days a month (so the RES_STAMP=yyyy1230)
  # otherwise do 8870        (so the RES_STAMP=yyyy1231)
  # do something else for other time units
  export RES_STAMP=$(printf %08d $((NEW_DATE_NUM - 8871)))

  ########################################################
  # 2) move files around and tidy up
  ########################################################

  cp -pv ${MODEL}_${RES_TIMESTAMP}_restart.nc ./RESTARTS/${MODEL}_${RES_STAMP}_restart.nc
  cp -pv ./output.namelist.dyn ./OUTPUTS/output.namelist.dyn.${RES_STAMP}
  #cp -pv ${MODEL}_${RES_TIMESTAMP}_restart_ice.nc ./RESTARTS/${MODEL}_${RES_STAMP}_restart_ice.nc
  #cp -pv ./output.namelist.ice ./OUTPUTS/output.namelist.ice.${RES_STAMP}
  cp -pv ./ocean.output ./OUTPUTS/ocean.output.${RES_STAMP}
  cp -pv ./solver.stat ./OUTPUTS/solver.stat.${RES_STAMP}
  cp -pv ./stdouterr ./OUTPUTS/stdouterr.${RES_STAMP}
  cp -pv ./namelist_cfg ./OUTPUTS/namelist_cfg.${RES_STAMP}

  #cp -pv ./volume_transport ./OUTPUTS/volume_transport.${RES_STAMP}
  #cp -pv ./salt_transport ./OUTPUTS/salt_transport.${RES_STAMP}
  #cp -pv ./heat_transport ./OUTPUTS/heat_transport.${RES_STAMP}

  rm -v ${MODEL}_${RES_TIMESTAMP}_restart*
  rm -v restart.nc 
  #rm -v restart_ice.nc
  rm -v ${MODEL}_*_????.nc
  mv ${MODEL}*.nc ./OUTPUTS

  cp -pv RESTARTS/${MODEL}_${RES_STAMP}_restart.nc ./restart.nc
  if (($? > 0)); then
    ERR_CATCH=$((ERR_CATCH + 1))
    echo "  ERR: copying restart file into folder"
  fi

  #cp -pv RESTARTS/${MODEL}_${RES_STAMP}_restart_ice.nc ./restart_ice.nc
  #if (($? > 0)); then 
  #  ERR_CATCH=$((ERR_CATCH + 1))  
  #  echo "  ERR: copying restart_ice file into folder"
  #fi

  ########################################################
  # 3) if all good, then modify namelist_cfg and resbumit
  ########################################################

  if (($ERR_CATCH > 0)) || ((${NEW_DATE_NUM} > $THRESH)); then
    if (($ERR_CATCH > 0)); then
      echo " "
      echo " "
      echo " "
      echo "ERR: caught a non-zero exit status, check cleanup.log for what the deal was"
      echo "ERR: caught a non-zero exit status, check cleanup.log for what the deal was"
    else
      echo "OK: grabbed time stamp ${NEW_DATE_NUM} larger than threshold ${THRESH}, breaking..."
      echo "OK: grabbed time stamp ${NEW_DATE_NUM} larger than threshold ${THRESH}, breaking..."
      # WARNING: this assumes that OLD_DATE_NUM is the only number within the file, which should
      #          really be true
      sed -i "s/${OLD_DATE_NUM}/${NEW_DATE_NUM}/g" namelist_cfg
    fi
    echo " "
    echo " "
    echo " "
    echo " "
    echo " ... a wild Totoro appeared and blocked your resubmission!"
    echo "         ,--'''',--.__,---[],-------._                               "
    echo "       ,'   __,'            \         \--''''''==;-                  "
    echo "     ,' _,-'  '/---.___     \       ___\   ,-'','                    "
    echo "    /,-'      / ;. ,.--'-.__\  _,-'' ,| ','   /                      "
    echo "   /''''''-._/,-|:\       []\,' '''-/:;-. '. /                       "
    echo "             '  ;:::      ||       /:,;  '-.\                        "
    echo "                =.,'__,---||-.____',.=                               "
    echo "                =(:\_     ||__    ):)=                               "
    echo "               ,'::::'----||::'--':::'._                             "
    echo "             ,':::::::::::||::::::::::::'.                           "
    echo "    .__     ;:::.-.:::::__||___:::::.-.:::\     __,                  "
    echo "       '''-;:::( O )::::>_|| _<::::( O )::::-'''                     "
    echo "   =======;:::::'-':::::::||':::::::'-':::::\=======                 "
    echo "    ,--'';:::_____________||______________::::''----.          , ,   "
    echo "         ; ::'._(    |    |||     |   )_,'::::\_,,,,,,,,,,____/,'_,  "
    echo "       ,;    :::'--._|____[]|_____|_.-'::::::::::::::::::::::::);_   "
    echo "      ;/ /      :::::::::,||,:::::::::::::::::::::::::::::::::::/    "
    echo "     /; ''''''----------/,'/,__,,,,,____:::::::::::::::::::::,'      "
    echo "     ;/                :);/|_;| ,--.. . '''-.:::::::::::::_,'        "
    echo "    /;                :::):__,'//''\\. ,--.. \:::,:::::_,'           "
    echo "   ;/              :::::/ . . . . . . //''\\. \::':__,'              "
    echo "   ;/          :::::::,' . . . . . . . . . . .:'::\                  "
    echo "   ';      :::::::__,'. ,--.. . .,--. . . . . .:'::'                 "
    echo "   ';   __,..--'''-. . //''\\. .//''\\ . ,--.. :':::'                "
    echo "   ;    /  \\ .//''\\ . . . . . . . . . //''\\. :'::'                "
    echo "   ;   /       . . . . . . . . . . . . . . . . .:'::'                "
    echo "   ;   (          . . . . . . . . . . . . . . . ;:::'                "
    echo "   ,:  ;,            . . . . . . . . . . . . . ;':::'                "
    echo "   ,:  ;,             . . . . . . . . . . . . .;':::'                "
    echo "   ,:   ;,             . . . . . . . . . . . . ;'::;'                "
    echo "     :   ;             . . . . . . . . . . . ,':::;                  "
    echo "      :   '.          . . . . . . . .. . . .,':::;'                  "
    echo "       :    '.       . . . . . . . . . . . ;::::;'                   "
    echo "        '.    '-.   . . . . . . . . . . ,-'::::;                     "
    echo "          ':_    ''--..___________..--'':::::;''                     "
    echo "             '._::,.:,.:,:_ctr_:,:,.::,.:_;''                        "
    echo "________________''\/'\/\/''''''\/'\/''\/'____________________________"

  else
  # WARNING: this assumes that OLD_DATE_NUM is the only number within the file, which should
    #          really be true
    sed -i "s/${OLD_DATE_NUM}/${NEW_DATE_NUM}/g" namelist_cfg
    
    echo "grabbed time stamp ${NEW_DATE_NUM} smaller than threshold ${THRESH}, resubmitting..."
    echo "grabbed time stamp ${NEW_DATE_NUM} smaller than threshold ${THRESH}, resubmitting..."
    echo "grabbed time stamp ${NEW_DATE_NUM} smaller than threshold ${THRESH}, resubmitting..."
    echo "grabbed time stamp ${NEW_DATE_NUM} smaller than threshold ${THRESH}, resubmitting..."
    echo " "
    echo "OK: ...and here is Christopher resbumitting the job for you......"
    echo "                  ,-.____,-.          "
    echo "                  /   ..   \          "
    echo "                 /_        _\         "
    echo "                |'o'      'o'|        "
    echo "               / ____________ \       "
    echo "             , ,'    '--'    '. .     "
    echo "            _| |              | |_    "
    echo "          /  ' '              ' '  \  "
    echo "         (    ',',__________.','    ) "
    echo "          \_    ' ._______, '     _/  "
    echo "             |                  |     "
    echo "             |    ,-.    ,-.    |     "
    echo "              \      ).,(      /      "
    echo "         gpyy   \___/    \___/        "
    sbatch submit_nemo

  fi

  exit
  
The output recombination steps have bene commented out because ARC does have parallel NetCDF4 and so the ``one_file`` option in ``field_def_nemo.xml`` already takes care of the outputs.
