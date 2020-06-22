.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   
.. _sec:nemo37:

NEMO 3.7/4.0 + XIOS 2.0
=======================

Tested with

* ``gcc4.9``, ``gcc5.4`` on a linux system
* ``gcc4.8`` on a Mac (El Capitan OSX 10.11)

This is the version I first implemented GEOMETRIC in, which is a development
version I guess (?) that eventually led to NEMO 4.0. The code structure largely
follows NEMO 3.6 but the commands are slightly different. 

If you get errors that are not documented here, see if :ref:`the XIOS1.0 NEMO3.6
<sec:nemo36>` page contain the relevant errors.

The assumption here is that the compiler is fixed and the packages (e.g.,
NetCDF4 and a MPI bindings) are configured to be consistent with the compilers.
See :ref:`here <sec:other-pack>` to check whether the binaries exist, where they
are, and how they might be installed separately if need be. All the ``#CHANGE
ME`` highlighted below needs be modified to point to the appropriate paths or
binaries (soft links with ``ln -s`` are ok). 

The instructions below uses ``gcc4.9`` for demonstration (modifications with
``gcc5.4`` as appropriate). I defined some extra variables on a Linux machine:

.. code-block:: bash

  export $BD=/home/julian/testing/gcc4.9-builds # CHANGE ME

  export C_INCLUDE_PATH=$BD/install/include:$C_INCLUDE_PATH
  export CPLUS_INCLUDE_PATH=$BD/install/include:$CPLUS_INCLUDE_PATH
  export LIBRARY_PATH=$BD/install/lib:$LIBRARY_PATH
  export LD_LIBRARY_PATH=$BD/install/lib:$LD_LIBRARY_PATH
  
You shouldn't need to do the above if the packages are forced to look at the
right place (e.g. via ``-L`` and/or ``-I`` flags with path to libraries and
include files respectively). Not all of these are necessary depending on whether
you choose to build/have static or dynamic libraries, and the
``LD_LIBRARY_PATH`` seems to sort out a lot of problems with linking libraries.

On a Mac done through anaconda the above was not necessary. My understanding is
that setting these variables might not actually do anything unless an option is
specifically enabled in Xcode.

XIOS 2.0 (svn v1322)
--------------------

Do the following:

.. code-block:: none

  mkdir XIOS
  cd XIOS
  svn checkout -r 1322 http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/trunk xios-2.0
  
.. note ::

  ``svn checkout
  http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-2.0 xios-2.0``
  also works with instructions below.
  
To get XIOS to compile, the compilers and packages need to be pointed to first,
via modifying files in ``arch``. Since I am using ``gcc``, I did the following
just to make a fresh copy:

.. code-block:: none

  cd xios2.0/arch
  cp arch-GCC_LINUX.env arch-GCC_local.env
  cp arch-GCC_LINUX.fcm arch-GCC_local.fcm
  cp arch-GCC_LINUX.path arch-GCC_local.path
  
The ``*.env`` file specifies where HDF5 and NetCDF4 libraries live. The
``*.fcm`` file specifies which compilers and options to use. The ``*.path`` file
specifies which paths and options to include. My files look like the following:

.. code-block:: none

  # arch-GCC_local.env

  export HDF5_INC_DIR=/usr/local/include       # CHANGE ME
  export HDF5_LIB_DIR=/usr/local/lib           # CHANGE ME

  export NETCDF_INC_DIR=/usr/local/include     # CHANGE ME
  export NETCDF_LIB_DIR=/usr/local/lib         # CHANGE ME
  
You could get an idea where the HDF5 and NetCDF4 directories are by doing
``which h5copy`` and ``which nc-config`` (assuming these are on ``$PATH``),
which should give you a ``directory/bin``, and it is the ``directory`` part you
want. If you did install the libraries somewhere else as in :ref:`other packages
<sec:other-pack>`, say, then make sure the ``which`` commands are pointing to
the right place.

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
  
Check the MPI locations and versions by doing ``which mpicc`` and ``mpicc
--version`` say. If they are the right ones you could just have ``mpicc``
instead of the full path as given above. MPI bindings are used here to avoid a
possible error that may pop up in relation to the build trying to find
``mpi.h``. The ``gmake`` command was swapped out by the ``make`` command (I
don't have ``cmake`` on the laptop).

.. note ::

  For ``gcc5.4`` and maybe newer versions, doing just the above when compiling
  leads to a whole load of errors about clashing in C++:
  
  .. code-block:: bash
    
    .../include/boost/functional/hash/extensions.hpp:69:33: error: ‘template<class T, class A> std::size_t boost::hash_value’ conflicts with a previous declaration
     std::size_t hash_value(std::list<T, A> const& v)
                                 ^
  
  Adding ``-D_GLIBCXX_USE_CXX11_ABI=0`` to ``%BASE_CFLAGS`` fixes these.

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
``xios2.0/arch``:

.. code-block:: none

  cd ../
  ./make_xios --full --prod --arch GCC_local -j2 |& tee compile_log.txt
  
The ``-j2`` option uses two processors to build. The ``tee`` command is to keep
logs of potential errors (the ``|&`` is short for ``2>&1 |``) for debugging the
compiler issues that may arise. It should work and takes around 5 mins to
compile for me. The main end result is are binaries in ``xios2.0/bin/`` which
NEMO will call.

NEMO 3.7/4.0 (svn v8666)
------------------------

Check out a version of NEMO. I have another folder separate to the XIOS folders
to contain the NEMO codes and binaries:

.. code-block :: bash

  mkdir NEMO
  cd NEMO
  svn checkout -r 8666 http://forge.ipsl.jussieu.fr/nemo/svn/NEMO/trunk nemo3.7-8666
  
This checks out version 8666 (NEMO 3.7/4.0) and dumps it into a folder called
``nemo3.7-8666`` (change the target path to whatever you like). A similar
procedure to specify compilers and where XIOS lives needs to be done for NEMO.
Again, because of the compilers I am using:

.. code-block :: bash
  
  cd nemo3.7-8666/NEMOGCM/ARCH
  cp OLD/arch-gfortran_linux.fcm ./arch-gfortran_local.fcm
  
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

  %XIOS_HOME           /home/julian/testing/gcc4.9-builds/XIOS/xios-2.0  # CHANGE ME

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

The main changes are (see :ref:`here <sec:nemo-fcm-log>` for an attempt
at the reasoning and a log of errors that motivates the changes):

* added ``%NCDF_HOME`` to point to where NetCDF lives
* added ``%XIOS_*`` keys to point to where XIOS lives
* added ``%CPP`` and flags, consistent with using ``gcc4.9``
* added the ``-lnetcdff`` and ``-lstdc++`` flags to NetCDF flags
* using ``mpif90`` which is a MPI binding of ``gfortran-4.9``
* added ``-cpp`` and ``-ffree-line-length-none`` to Fortran flags
* swapped out ``gmake`` with ``make``

Then, I did (see :ref:`NEMO 3.6 <sec:nemo36>` for the reason):
  
.. code-block :: bash

  cd ../CONFIG/
  ./makenemo -j0 -r GYRE_PISCES -n GYRE_testing -m gfortran_local
    
Edit ``/GYRE_testing/cpp_GYRE_testing.fcm`` and replaced ``key_top`` with
``key_nosignedzero`` (does not compile TOP for speed reasons, and make sure
zeros are not signed). Then
  
.. code-block :: bash
  
  ./makenemo -j2 -n GYRE_testing -m gfortran_local |& tee compile_log.txt
  
This uses two processors, with ``GYRE`` as a reference, builds a new folder
called ``GYRE_testing``, with the specified architecture file, and outputs a
log.

.. note ::

  The ``-r GYRE`` flag here only needs to be done once to create an extra folder
  and add ``GYRE_testing`` to ``cfg.txt``. The subsequent compilations should
  then read, e.g., ``./makenemo -n GYRE_testing -m gfortran_local``.
  
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

  If your installation compiles but does not run with the following error
  
  .. code-block :: bash

    dyld: Library not loaded: @rpath/libnetcdff.6.dylib
    Referenced from: /paths/./nemo
    Reason: no suitable image found.  Did find:
    /usr/local/lib/libnetcdff.6.dylib: stat() failed with errno=13

  then it is not finding the right libraries. These could be fixed by adding the
  ``-Wl,-rpath,/fill me in/lib`` flag to the relevant flags bit in the ``*.fcm``
  files (or possibly in XIOS the ``path`` and/or ``env`` ) to specify exactly
  where the libraries live. This can happen for example on a Mac or if the
  libraries are installed not at the usual place.
  
.. note ::

  One infuriating problem I had specifically with a Mac (though it might be a
  ``gcc4.8`` issue) is that the run does not get beyond the initialisation
  stage. Going into ``ocean.output`` and searching for ``E R R O R`` shows that
  it complained about a misspelled namelist item (in my case it was in the
  ``namberg`` namelist). If you go into ``output.namelist.dyn`` and look for the
  offending namelist is that it might be reading in nonsense. This may happen if
  the comment character ``!`` is right next to a variable, e.g.

  ::
  
    ln_icebergs = .true.!this is a comment
    
  Fix this by adding a white space, i.e.
  
  ::
  
    ln_icebergs = .true. !this is a comment
