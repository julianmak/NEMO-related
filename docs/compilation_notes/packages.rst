.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _sec:other-pack:

Other packages
==============

Tested with

* ``gcc4.9``, ``gcc5.4`` on a laptop (Ubuntu 16.04)
* ``gcc4.9`` on a modular system (Ubuntu 14.04, Oxford AOPP)

The following packages are needed for NEMO and XIOS and they may need to be
installed or configured accordingly. If you want a script to do all of the
following in one go, then please scroll right to the bottom of this page.

.. note::

  The main issue I have been founding is trying to get the compilations to call
  the correct versions of (static and/or dynamic) libraries that has been
  compiled. If you can get a system manager to install the packages
  (particularly HDF5 and NetCDF4) it would save a lot of time, unless the
  compiled libraries itself are clashing...
  
  As of 21 Aug 2018, the following remains on the agenda:
  
  * intel compilers
  * gcc on Mac OSX
  * reproducing sample compatibility errors

Preliminaries
-------------

The way I went about it was to first choose a set of compilers and use the same
set of compilers to install the dependencies, primarily to avoid errors relating
to compatibility of packages. For example, ``gcc4.9`` was downloaded through
``sudo apt-get install gcc4.9``, or loaded through a network computer through
something like a ``module load`` command. You may have to look it up on the
internet if you don't have either of these.

The order I did them in are:

1. mpich (to bind the set of compilers to a MPI form; I chose ``mpich`` but it should work on ``OpenMP`` too)
2. zlib (1.2.11, for HDF5)
3. hdf5 (1.8.19, for NetCDF)
4. netcdf (4.4.1.1) and netcdf-fortran (4.4.4), for XIOS

Within a folder called ``gcc4.9-builds``, I added an extra ``extra_variables``
file containing the following:

.. code-block:: bash

  export $BD=/home/julian/testing/gcc4.9-builds # CHANGE ME

  export CC=/usr/bin/gcc-4.9
  export CXX=/usr/bin/g++-4.9
  export FC=/usr/bin/gfortran-4.9
  export F77=/usr/bin/gfortran-4.9
  export CPP=/usr/bin/cpp-4.9

  export C_INCLUDE_PATH=$BD/install/include:$C_INCLUDE_PATH
  export CPLUS_INCLUDE_PATH=$BD/install/include:$CPLUS_INCLUDE_PATH
  export LIBRARY_PATH=$BD/install/lib:$LIBRARY_PATH
  export LD_LIBRARY_PATH=$BD/install/lib:$LD_LIBRARY_PATH
  export PATH=$BD/install/bin:$PATH
  
Set this by doing ``source extra_variables``, and upon closing the terminal the
variables will be flushed. Some of these may want to be added to ``~/.bashrc``
for convenience. The instructions below attempts to build shared rather than
static libraries, and somewhat depends ``LD_LIBRARY_PATH`` variable being set
(with the added bonus that the ``ldd`` command provides an extra check whether
the correct libraries are being called). Suggestions on how to build the
packages without setting ``LD_LIBRARY_PATH`` or build static packages are given
below (using ``LD_LIBRARY_PATH`` can be dangerous, see e.g., `here
<http://xahlee.info/UnixResource_dir/_/ldpath.html>`_).

.. note::

  Do for example ``$CC --version`` or ``echo $CC`` to see what the variables are
  set to. If you don't want to set the compiler variables then you need to do
  e.g.
  
  .. code:: bash
  
    CC=/usr/bin/gcc-4.9 FC= something ./configure something
    
  where the path points to where the compiler binary lives. This then only sets
  the variable temporarily for the particular command.
  
Some or all of these may be skipped depending on which ones packages you have
already installed and/or configured. The following installs all the libraries
and binaries to the folder specified in ``$BD``; you have ``sudo`` access you
could always just install it to ``/usr/local``. The sub-directories in the
folder are:

* ``source``, where all the compressed files are going to live;
* ``build``, where all the source file folders are going to live
* ``install``, where all the compiled libraries, binaries and header files are going to live.

``source`` and ``build`` can be deleted later.

.. note::

  The binaries built here will not register by default unless it is added to the
  ``$PATH`` variable. If you are going to add to the ``$PATH`` variable, the one
  that gets registered **first** gets priority, i.e.
  
  .. code:: bash
    
    echo $PATH
    > /home/julian/testing/gcc4.9-builds/install/bin:/usr/local/bin
    
  means any binaries in ``/home/julian/testing/gcc4.9-builds/install/bin`` gets
  used first. Do this by adding to ``~/.bashrc`` the following:
  
  .. code:: bash 
  
    export PATH=/usr/local/bin:$PATH
  
  If you don't do this then it just means when you call the binaries you have to
  provide an explicit call, e.g.,
  ``/home/julian/testing/gcc4.9/build/bin/mpif90``. Do for example ``which
  mpif90`` to check what the ``mpif90`` is linked to; if you did add to
  ``$PATH`` then the ``which`` command above should point to the right binary. 

MPICH
-----

Check if there are any MPI capabilities and which compilers they are bound to:

.. code-block:: bash
  
  mpicc --version
  which mpicc
  
If you have these already they may not need to be installed. If they need to be
installed separately for whatever reason, then you could do the following. I
took the source files from the `MPICH website
<http://www.mpich.org/static/downloads/>`_ itself and chose v3.0.4 here. Being
in the ``$BD`` folder, I did:

.. code-block:: bash

  cd $BD/source/
  wget http://www.mpich.org/static/downloads/3.0.4/mpich-3.0.4.tar.gz
  cd $BD/build/
  tar -xvzf $BD/source/mpich-3.0.4.tar.gz
  cd mpich-3.0.4
  ./configure prefix=$BD/install/
  make -j 2
  make check install
  
Within ``install/`` there should now be some folders that can be pointed to for
the binaries, libraries and header files to include for later installations.
  
.. note::

  The ``./configure prefix=`` step requires an absolute (not relative) path;
  change this to change the installation folder.
  

zlib and DF5
------------

Check whether HDF5 exists first (may still need to be installed again for
compatibility reasons). ``h5copy`` is the command that should exist if HDF5 is
installed:

.. code-block:: bash
  
  which h5copy
  h5copy --version
  
If you still want to install both zlib and HDF5, then do the following
(following the instructions on the `Unidata UCAR website
<https://www.unidata.ucar.edu/software/netcdf/netcdf-4/newdocs/netcdf-install/Quick-Instructions.html>`_).
The raw files are taken from the HDF5 website using HDF5 v1.8.19. Again, with
``$BD`` as defined:

.. code-block:: bash
  
  cd $BD/source/
  wget http://www.zlib.net/zlib-1.2.11.tar.gz
  cd $BD/build/
  tar -xvzf $BD/source/zlib-1.2.11.tar.gz
  cd zlib-1.2.11
  CFLAGS=-fPIC ./configure --prefix=$BD/install/
  make -j 2
  make check install
  
  cd $BD/source/
  wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.19/src/hdf5-1.8.19.tar.gz
  cd $BD/build/
  tar -xvzf $BD/source/hdf5-1.8.19.tar.gz
  cd hdf5-1.8.19
  #CPPFLAGS=-I$BD/install/include LDFLAGS=-L$BD/install/lib \
  CFLAGS=-fPIC ./configure --enable-shared --enable-fortran --enable-cxx \
  # --with-zlib=$BD
  --prefix=$BD/install/
  make -j 2
  make check install
  cd $BD
  
.. note::
  
  If ``LD_LIBRARY_PATH`` is set then accordingly then zlib should be detected by
  the HDF5 install. If not, consider including the commented out ``CPPFLAGS``
  and ``LDFLAGS`` or the ``--with-zlib`` line (or both).
  
  HDF5 checking and installation can take a while. If it's more that 30 mins
  however it probably has crashed.
  
  If a shared build option was on, then you can do ``ldd h5copy`` (or wherever
  ``h5copy`` is installed at if the directory has not been added to ``$PATH``)
  to check that ``libhdf5`` does point to where you think it should point to. If
  it isn't, then try the first point in this note.
  
  If an error shows up saying ``recompile with -fPIC``, then trying doing a
  static build. Replace ``--enable-shared`` with ``--disable-shared`` and do the
  first point in this note, possibly adding ``LIBS="-lz -lhdf5`` etc.; see `here
  <https://www.unidata.ucar.edu/software/netcdf/docs/building_netcdf_fortran.html>`_
  for a guide. I would be tempted to keep the ``CFLAGS=-fPIC`` so shared builds
  of NetCDF4 can still be made.

NetCDF4
-------

Check whether NetCDF4 exists first (may still need to be installed again for
compatibility reasons). ``nc-config`` is the command that should exist if
NetCDF4 is installed, and shows where it is installed and what compilers were
used to build it.

.. code-block:: bash
  
  nc-config all
  
If you still want to install it, then do the following (following the
instructions on the `Unidata UCAR website
<https://www.unidata.ucar.edu/software/netcdf/netcdf-4/newdocs/netcdf-install/Quick-Instructions.html>`_).
The raw files are taken from the the NetCDF4 website, using netcdf v4.4.1.1 and
netcdf-fortran v4.4.4:

.. code-block:: bash

  cd $BD/source/
  wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.4.1.1.tar.gz
  cd $BD/build/
  tar -xvzf $BD/source/netcdf-4.4.1.1.tar.gz
  cd netcdf-4.4.1.1
  #CPPFLAGS=-I$BD/install/include LDFLAGS=-L$BD/install/lib \
  ./configure --enable-netcdf4 --enable-shared --prefix=$BD/install/
  make -j 2
  make check install
  
  cd $BD/source/
  wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-fortran-4.4.4.tar.gz
  cd $BD/build/
  tar -xvzf $BD/source/netcdf-fortran-4.4.4.tar.gz
  cd netcdf-fortran-4.4.4
  #CPPFLAGS=-I$BD/install/include LDFLAGS=-L$BD/install/lib \
  ./configure --enable-shared --prefix=$BD/install/
  make -j 2
  make check install
  cd $BD
  
.. note::
  
  NetCDF4 checking and installation can take a while. If it's more that 30 mins
  however it probably has crashed.
  
  If a shared build option was on, then you can do ``ldd ncdump`` (or wherever
  ``ncdump`` was installed if the directory has not been added to ``$PATH``) and
  check that ``libnetcdf``, ``libhdf5`` and ``libz`` really does point to where
  you think it should point to. If not, consider doing something similar to the
  HDF5 note above.
  
  If an error shows up saying ``recompile with -fPIC``, then trying doing a
  static build (I had this problem on one of the computers where the Fortran
  part is static). See HDF5 note above.

  I had a problem with not having the m4 package, which I just installed as the
  installation commands above, with the binaries found from ``wget
  ftp://ftp.gnu.org/gnu/m4/m4-1.4.10.tar.gz``. This is not in the script below.

This should be it! Try ``./install/bin/nc-config --all`` and/or
``./install/bin/nf-config --all`` to see where everything is configured. The
things in ``build/`` and ``source/`` may now be deleted.

Combined shell script
---------------------

A script that does **all** of the above in one go may be found in the following
commands (use at your own risk):

.. code-block :: bash

  mkdir gcc4.9-builds/               # CHANGE ME
  cd gcc4.9-builds/                  # CHANGE ME
  wget https://raw.githubusercontent.com/julianmak/NEMO-related/master/docs/compilation_notes/compile_dependencies.sh
  chmod +x compile_dependencies.sh
  
Before you execute the shell script with ``./compile_dependencies.sh``, make
sure the compilers are pointed to appropriately. You can do this in
``~/.bashrc`` (see first code block on this page) or within the shell script
itself (it is commented out at the moment). If some packages already exist and
you don't want them installed, comment the appropriate lines.
