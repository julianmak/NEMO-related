.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   
.. _sec:nemo-fcm-log:

detailed log of generating arch-gfortran_local.fcm
==================================================

Log is longer than it needs to be to highlight some possible errors that could
be thrown up. The following is sort of adapting a version of the fcm file that
works on ARCHER (e.g. `NOCL ARCHER guide
<https://nemo-nocl.readthedocs.io/en/latest/work_env/archer.html>`_). 

Starting with ``arch-gfortran_local.fcm`` as copied from
``OLD/arch-gfortran_linux.fcm``:

.. code-block :: none

  %NCDF_INC            -I/usr/local/netcdf/include
  %NCDF_LIB            -L/usr/local/netcdf/lib -lnetcdf
  %FC	             gfortran
  %FCFLAGS             -fdefault-real-8 -O3 -funroll-all-loops -fcray-pointer 
  %FFLAGS              %FCFLAGS
  %LD                  gfortran
  %LDFLAGS
  %FPPFLAGS            -P -C -traditional
  %AR                  ar
  %ARFLAGS             -rs
  %MK                  make
  %USER_INC            %NCDF_INC
  %USER_LIB            %NCDF_LIB

To be consistent with where I installed NetCDF4 and what compilers I am using,
the above is modified to the following:

.. code-block :: none

  %NCDF_HOME           /usr/local

  %NCDF_INC            -I%NCDF_HOME/include
  %NCDF_LIB            -L%NCDF_HOME/lib -lnetcdf
  %FC	                 gfortran-4.9
  %FCFLAGS             -fdefault-real-8 -O3 -funroll-all-loops -fcray-pointer 
  %FFLAGS              %FCFLAGS
  %LD                  %FC
  %LDFLAGS
  %FPPFLAGS            -P -C -traditional
  %AR                  ar
  %ARFLAGS             -rs
  %MK                  make
  %USER_INC            %NCDF_INC
  %USER_LIB            %NCDF_LIB

CPP flags
---------

Compile with ``./makenemo -j2 -r GYRE -n GYRE_testing -m gfortran_local |& tee compile_log.txt`` gives

.. code-block :: none
  
  %CPP: variable not expanded
  
C pre-processing is not defined, so add in the following keys:

.. code-block :: none

  %CPP	               cpp-4.9
  %CPPFLAGS            -P -traditional
  
``cpp-4.9`` is my C-preprocessing associated with ``gcc4.9``; do something like
``whereis cpp-4.9`` to find if it is the right binary to call. (Note C++ appears
as ``g++`` or ``mpicxx`` for example). 

Fortran talking with CPP
------------------------

Compile with ``./makenemo -n GYRE_testing clean_config && ./makenemo -r GYRE -n GYRE_testing -m gfortran_local -j2`` (note the clean compile) gives something like:

.. code-block :: none

  /home/julian/testing/nemo-6800/nemo6800/NEMOGCM/CONFIG/GYRE_testing/BLD/ppsrc/nemo/bdylib.f90:1.1:

  /* Copyright (C) 1991-2016 Free Software Foundation, Inc.
   1
  Error: Invalid character in name at (1)

This is saying that the Fortran interpreter is not recognising the formatting.
This is fixed by adding a ``-cpp`` flag:

.. code-block :: none

  %FCFLAGS             -fdefault-real-8 -O3 -funroll-all-loops -fcray-pointer -cpp
  
Fortran line length
-------------------

Same again, gets rid of that error but then something like the below appears:

.. code-block :: none

  /home/julian/testing/nemo-6800/nemo6800/NEMOGCM/CONFIG/GYRE_testing/BLD/ppsrc/nemo/lib_mpp.f90:3011.132:

  4,num_fields), zfoldwk(jpi,4,num_fields), znorthgloio(jpi,4,num_fields,jpni

This indicates that the line is too short (there is a default limit on how many
characters a line should have in Fortran). Fix this by adding
``-ffree-line-length-none``:

.. code-block :: none
  
  %FCFLAGS             -fdefault-real-8 -O3 -funroll-all-loops -fcray-pointer -cpp -ffree-line-length-none

XIOS calling
------------

Deals with that, but now it stops with:

.. code-block :: none

  /home/julian/testing/nemo-6800/nemo6800/NEMOGCM/CONFIG/GYRE_testing/BLD/ppsrc/nemo/iom.f90:76.7:

   USE xios
   
XIOS path is not added so add in the following keys:

.. code-block :: none

  %XIOS_HOME           /home/julian/XIOS/xios1.0

  %XIOS_INC            -I%XIOS_HOME/inc 
  %XIOS_LIB            -L%XIOS_HOME/lib -lxios
  
  %USER_INC            %XIOS_INC %NCDF_INC
  %USER_LIB            %XIOS_LIB %NCDF_LIB

Change the ``%XIOS_HOME`` to where XIOS1.0 is installed.

C++ linking
-----------

Passes the XIOS flag and now something like this pops up:

.. code-block :: none

  operator_expr.cpp:(.text._ZN4xios13COperatorExprC2Ev[_ZN4xios13COperatorExprC5Ev]+0xaf3): undefined reference to `std::string::_Rep::_S_empty_rep_storage'
  std::allocator<char> >::basic_string(char const*, std::allocator<char> const&)'
  ...
  /home/julian/testing/nemo-6800/xios-703/xios-1.0/lib/libxios.a(operator_expr.o):(.eh_frame+0x5a7): undefined reference to `__gxx_personality_v0'
  collect2: error: ld returned 1 exit status
  
This is one where having a log is useful. A whole load of error pops up to say
the C++ files are not being interpreted, to do with a linker error (probably
easiest to scroll down from top rather than up from bottom actually). Normally
one might expect that adding the ``-lstdc++`` flag to ``%LDFLAGS`` would work
but it doesn't for whatever reason. The ``-lstdc++`` flag seems to need to go
**right at the end** of the command line, which meant I did the following:

.. code-block :: none

  %NCDF_LIB            -L%NCDF_HOME/lib -lnetcdf -lstdc++
  
If someone could explain to me why this works do send me an e-mail!

MPI errors
----------

Turns out actually there are two errors that get thrown up previously, with the
linker error dominating the output. Sorting the linker one gives something like
the following errors:

.. code-block :: none

  client.cpp:(.text+0x1cd8): undefined reference to `MPI_Intercomm_create'
  client.cpp:(.text+0x1d6d): undefined reference to `MPI_Intercomm_merge'
  client.cpp:(.text+0x1d76): undefined reference to `MPI_Barrier'
  client.cpp:(.text+0x1da8): undefined reference to `MPI_Comm_dup'
  /home/julian/testing/nemo-6800/xios-703/xios-1.0/lib/libxios.a(client.o): In function `xios::CClient::openStream(std::string const&, std::string const&, std::basic_filebuf<char, std::char_traits<char> >*)':
  client.cpp:(.text+0x2076): undefined reference to `MPI_Comm_size'
  collect2: error: ld returned 1 exit status

The compiler used here is the serial rather than the MPI version, so try

.. code-block :: none

  %FC	                 mpif90
  
where ``mpif90`` is my binding of ``gfortran-4.9`` to MPI.

NetCDF errors
-------------

The above procedure gets rid of the MPI errors but throws up a whole load of
NetCDF errors that look like the following:

.. code-block :: none

  obs_fbm.f90:(.text+0x114ef): undefined reference to `__netcdf_MOD_nf90_get_att_text'
  obs_fbm.f90:(.text+0x11536): undefined reference to `__netcdf_MOD_nf90_get_att_text'
  obs_fbm.f90:(.text+0x11591): undefined reference to `__netcdf_MOD_nf90_inquire_dimension'
  obs_fbm.f90:(.text+0x115e1): undefined reference to `__netcdf_MOD_nf90_inquire_dimension'
  collect2: error: ld returned 1 exit status
  
Looks like a NetCDF Fortran error so add the ``-lnetcdff`` flag (notice the
extra `f`):

.. code-block :: none

  %NCDF_LIB            -L%NCDF_HOME/lib -lnetcdf -lnetcdff -lstdc++
  
...success!

