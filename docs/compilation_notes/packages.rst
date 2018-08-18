.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _sec:other-pack:

Other packages
==============

The following packages are needed for NEMO and XIOS and they may need to be
installed or configured accordingly. If you want a script to do all of the
following in one go, then please scroll right to the bottom of this page.

Preliminaries
-------------

.. note::

  As of 18 Aug 2018, the following are remains on the agenda:
  
  * testing with ``gcc5.4``
  * testing with the intel compilers (using the Oxford system)
  * reproducing sample compatibility errors (using the Oxford system)
  * compilation of a script file that will do **everything** below

The way I went about it was to first choose a set of compilers and use the same
set of compilers to install some of the dependencies, primarily to avoid errors
relating to compatibility of packages. For example, ``gcc4.9`` was downloaded
through ``sudo apt-get install gcc4.9``, or loaded through a network computer
through something like a ``module load`` command. You may have to look it up on
the internet if you don't have either of these.

The order I did them in are:

1. mpich and hydra (to bind the set of compilers to a MPI form; I chose ``mpich`` but it should work on ``OpenMP`` too)
2. zlib (1.2.11, for HDF5)
3. hdf5 (1.8.19, for NetCDF)
4. netcdf (4.4.1.1) and netcdf-fortran (4.4.4), for XIOS

I put the following in my ``~/.bashrc`` file and called it with ``source
~/.bashrc``; the following codes have it explicitly defined just in case:

.. code-block:: bash

  export CC=/usr/bin/gcc-4.9
  export CXX=/usr/bin/g++-4.9
  export FC=/usr/bin/gfortran-4.9
  export F77=/usr/bin/gfortran-4.9
  export CPP=/usr/bin/cpp-4.9
  
.. note::

  Do for example ``$CC --version`` or ``echo $CC`` to see what the variables are
  set to.
  
  If you don't want to set the variables to anything, then you need to do e.g.
  
  .. code:: bash
  
    CC=/usr/bin/gcc-4.9 FC= something ./configure something
    
  where the path points to where the compiler binary lives. This then only sets
  the variable temporarily for the particular command, but you do need to it for
  all the ``./configure`` commands below.
  
Some or all of these may be skipped depending on which ones packages you have
already installed and/or configured. The following installs all the libraries
and binaries to a folder (``gcc4.9-builds`` for example), but if you have
``sudo`` access you could always just install it to ``/usr/local``. The
sub-directories in the folder are:

* ``source``, where all the compressed files are going to live;
* ``build``, where all the source file folders are going to live
* ``install``, where all the compiled libraries, binaries and header files are going to live.

``source`` and ``build`` can be deleted later.

.. note::

  The binaries built here will not register by default unless you add it to the
  ``$PATH$`` variable in ``~./bashrc`` or equivalent. If you are going to the
  folder to the ``$PATH`` variable, the one that gets registered **first** gets
  priority, i.e.
  
  .. code:: bash
    
    echo $PATH
    > /home/julian/testing/gcc4.9-builds/install/bin:/usr/local/bin
    
  means any binaries in ``/home/julian/testing/gcc4.9-builds/install/bin`` gets
  used first. Do this by adding to ``~/.bashrc`` the following:
  
  .. code:: bash 
  
    export PATH=/usr/local/bin:$PATH
  
  If you don't do this then it just means when you call the binaries you have to
  provide an explicit call, e.g.,
  ``/home/julian/testing/gcc4.9/build/bin/mpif90``. Do for example ``whereis
  mpif90`` to check what the ``mpif90`` is linked to; if you did add to
  ``$PATH`` then the ``whereis`` command above should point to the right binary. 

MPICH
-----

Check if you have any MPI and which compilers they are bound to using, for
example,

.. code-block:: bash
  
  mpicc --version
  whereis mpicc
  
If you have these already they may not need to be installed. If you do want to
install it separately for whatever reason, then you could do the following. I
took the source files from the `MPICH website
<http://www.mpich.org/static/downloads/>`_ itself and chose v3.0.4 here. Being
in the ``gcc.4.9-builds`` folder, I did:

.. code-block:: bash

  cd source/
  wget http://www.mpich.org/static/downloads/3.0.4/mpich-3.0.4.tar.gz
  cd ../build/
  tar -xvzf ../source/mpich-3.0.4.tar.gz
  cd mpich-3.0.4
  ./configure prefix=/home/julian/testing/gcc4.9-builds/install/  # CHANGE ME
  make -j 2
  make check install
  cd ../../
  
  cd source/
  wget http://www.mpich.org/static/downloads/3.0.4/hydra-3.0.4.tar.gz
  cd ../build/
  tar -xvzf ../source/hydra-3.0.4.tar.gz
  cd hydra-3.0.4
  ./configure prefix=/home/julian/testing/gcc4.9-builds/install/  # CHANGE ME
  make -j 2
  make check install
  cd ../../
  
Within ``gcc4.9-builds/install/`` there should now be some folders that can be
pointed to for the binaries, libraries and header files to include for later
installations.
  
.. note::

  The ``./configure prefix=`` step requires an absolute (not relative) path.
  

zlib and HDF5 
-------------

Check whether HDF5 exists first (may still need to be installed again for
compatibility reasons). ``h5copy`` is the command that should exist if HDF5 is
installed.

.. code-block:: bash
  
  whereis h5copy
  h5copy --version
  
If you still want to install it, then do the following (following the
instructions on the `Unidata UCAR website
<https://www.unidata.ucar.edu/software/netcdf/netcdf-4/newdocs/netcdf-install/Quick-Instructions.html>`_).
The raw files are taken from the zlib and HDF5 website, using zlib v1.2.11 and
HDF5 v1.8.19. Again, being in the ``gcc4.9-builds`` directory:

.. code-block:: bash

  cd source/
  wget http://www.zlib.net/zlib-1.2.11.tar.gz
  cd ../build/
  tar -xvzf ../source/zlib-1.2.11.tar.gz
  cd zlib-1.2.11
  CFLAGS=-fPIC ./configure --prefix=/home/julian/testing/gcc4.9-builds/install/  # CHANGE ME
  make -j 2
  make check install
  cd ../../
  
  cd source/
  wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.19/src/hdf5-1.8.19.tar.gz
  cd ../build/
  tar -xvzf ../source/hdf5-1.8.19.tar.gz
  cd hdf5-1.8.19
  CFLAGS=-fPIC ./configure --enable-shared --enable-fortran --enable-cxx --with-zlib=/home/julian/testing/gcc4.9-builds/ --prefix=/home/julian/testing/gcc4.9-builds/install/    # CHANGE ME
  make -j 2
  make check install
  cd ../../
  
.. note::

  HDF5 checking and installation takes a while (~5-10 mins).

  If problems arise, try replacing ``--enable-shared`` with
  ``--disable-shared``. Then it is not necessary to compile zlib or use the
  ``-fPIC`` flags (doesn't matter since the libraries are likely just for
  private consumption).

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
netcdf-fortran v4.4.4. Again, being in the ``gcc4.9-builds`` directory:

.. code-block:: bash

  cd source/
  wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.4.1.1.tar.gz
  cd ../build/
  tar -xvzf ../source/netcdf-4.4.1.1.tar.gz
  cd netcdf-4.4.1.1
  ./configure --enable-netcdf4 --enable-shared --prefix=/home/julian/testing/gcc4.9-builds/install/   # CHANGE ME
  make -j 2
  make check install
  cd ../../
  
  cd source/
  wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-fortran-4.4.4.tar.gz
  cd ../build/
  tar -xvzf ../source/netcdf-fortran-4.4.4.tar.gz
  cd netcdf-fortran-4.4.4
  ./configure --enable-shared --prefix=/home/julian/testing/gcc4.9-builds/install/    # CHANGE ME
  make -j 2
  make check install
  cd ../../
  
This should be it! Try ``./install/bin/nc-config --all`` to see where everything
is configured. The things in ``build/`` and ``source/`` may now be deleted.
  
.. note::

  NetCDF checking and installation takes a while (~5-10 mins). The Fortran
  version however shouldn't take too long.

  If problems arise, try replacing ``--enable-shared`` with
  ``--disable-shared``.


Combined shell script
---------------------

The thing below does **all** of the above in one go (use at your own risk):

.. code-block :: bash

  mkdir gcc4.9-builds/               # CHANGE ME
  cd gcc4.9-builds/                  # CHANGE ME
  wget https://github.com/julianmak/NEMO-related/tree/master/docs/compilation_notes/compile_dependencies.sh
  chmod +x compile_dependencies.sh
  
Before you execute the shell script with ``./compile_dependencies.sh``, make
sure you the compilers are pointed to appropriately. You can do this in
``~/.bashrc`` (see first code block on this page) or within the script itself
(it is commented out at the moment). Some packages may already exist, and if you
don't want them installed you should comment out the appropriate lines in the
script and modify some of the paths accordingly.
