.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _sec:other-pack:

Other packages
==============

Tested with

* ``gcc4.9``, ``gcc5.4`` on a linux system
* ``gcc4.8`` on a Mac (El Capitan OSX 10.11)

The following packages are needed for NEMO and XIOS and they may need to be
installed or configured accordingly. I don't have a windows machine handy (and I
don't really want to try it there either) so for that I would recommend doing
the following through virtualbox or something analogous (which might be another
way to do it on a Mac); I am guessing `cygwin` and the new Windows 10 terminals
might be a possibility.

.. note::

  I would suggest trying the following in reverse order of effort required:

  1. Get someone who knows what they are doing to do it for you! Compiling the following from scratch is not the most interesting activity and is actually quite fiddly (especially the HDF5 and NetCDF4 stuff)...if you don't have access to people who can do that, then try
  2. Doing it through anaconda. There you are somewhat restricted to a certain set of compilers (``gcc 4.8``) but anaconda sorts out the dependencies for you. The only thing then you need to do is to force XIOS and NEMO to use the libraries within the anaconda installation. Failing that...
  3. Do it from scratch. I'm sorry and good luck; see below for some notes to possibly ease your pain.

Anaconda
--------

`Anaconda <https://www.anaconda.com/download/>`_ is a framework mostly for
downloading Python packages, with the added advantage that it resolves the
package dependencies for you (cf. ``apt``, ``yum`` on a Linux machine or
``port`` on a Mac if you have MacPorts). See the `official conda manual
<https://conda.io/docs/index.html>`_ or some of :ref:`my own notes <sec:python>`
on some things to do with installing and managing conda. I used the full
anaconda with Python 3.6 but you could use miniconda or with other pythons
probably.

.. note::

  [20 May 2020] Doing it through anaconda may well only work for Mac, because the gfortran versions does not seem to be available with linux through anaconda...

First I created an environment so all the changes only apply in that
environment:

.. code-block:: bash

  conda create -n nemo python=3.6
  
Accept to install the basic packages for the environment. Then activate the
``nemo`` environment with

.. code-block:: bash

  >> julian@psyduck:~/
  source activate nemo
  >> (nemo) julian@psyduck:~/
  
Now if you have compilers you want to use already then you can skip the compiler
installation. On the Mac I was dealing with there was no ``gcc`` or a Fortran
compiler and I had problems with ``clang``, so I did the following to get a set
of ``gcc`` compilers:

.. code-block:: bash

  conda install gcc
  conda install gfortran_osx-64
  
The second line you should change to ``gfortran_linux-64`` if on a Linux
machine. The command will add some compiler flags that is unset when exiting
from the environment. Check that the compilers are the now default compilers by
doing ``gcc --version`` (which should probably give 4.8) and ``which gcc``
(which should point to the anaconda folder). If not, do something like ``echo
CC`` and ``export CC=/folder/bin`` to force it to point to the right folder
(also do it for ``FC`` and ``CXX``, and maybe put it in the ``$PATH`` variable;
see below).

.. note::

  One thing I found to be an issue is that while ``gfortran`` can compile a
  sample program through ``gfortran hello.f90 -o hi`` with ``hello.f90`` being
  
  .. code-block:: fortran
  
    program hello
      print *, "hi mum"
    end program hello
    
  Executing through ``./hi`` could throw a library complaint:
  
  .. code-block:: bash
  
    dyld: Library not loaded: @rpath/libgfortran.3.dylib
    Referenced from: 
    Reason: no suitable image found.  Did find:
	  /usr/local/lib/libnetcdff.3.dylib: stat() failed with errno=13
	  
  So the problem here is that the computer is looking for the library at the
  wrong place. To force the computer to look at the right place, try
	
  .. code-block:: bash
  
    export FCFLAGS=-Wl,-rpath,${CONDA_PREFIX}/lib
	  
  where ``${CONDA_PREFIX}`` might have been defined by anaconda.

If you already have the MPI capabilities bound to the compilers you will use
then you can skip the following. To make life easier it is advisable to install
either MPICH or OpenMPI. You could try this by

.. code-block:: bash

  conda install -c conda-forge mpich
  conda install -c conda-forge openmpi

and check whether ``which mpicc`` and in particular ``which mpif90``, which
should be pointed to the ``gcc`` compilers. I had a similar problem with
``gfortran`` not being bound properly, which could be fixed with setting
``FCFLAGS``, or to compile it from scratch (see below for the way to do it for
MPICH, which also works for OpenMPI with suitable changes in the hyperlink
address; do a search for this in Google).

To get NetCDF4 and its dependencies I did

.. code-block:: bash

  conda install netcd4
  conda install -c conda-forge netcdf-fortran
  
Do ``which nc-config`` and ``nc-config --all`` to see which paths are being
pointed to. Again, you may need to add the ``FCFLAGS`` detailed above to make
sure it is pointing to the right libraries. Take note of the path where the
libraries and header files live and put those into the XIOS and NEMO files and
that should be it!

Compiling it yourself
---------------------

(Good luck!)

The following has been tried on a Linux machine. I had some problems on a Mac
with ``Clang`` that I don't know how to fix without ``sudo`` access but it is
probably fixable; I have not tried installing things with ``port`` through
MacPorts partly because it requires Xcode to be installed.

A script to do all of the following on a Linux machine in one go can be found at
bottom of this page. The way I went about it was to first choose a set of
compilers and use the same set of compilers to install the dependencies,
primarily to avoid errors relating to compatibility of packages. For example,
``gcc4.9`` was downloaded through ``sudo apt-get install gcc4.9``, or loaded
through a network computer through something like a ``module load`` command. You
may have to look it up on the internet if you don't have either of these.

.. note::

  If you don't have the right compilers you can always try and build your own
  from source, but it takes a while (order of hours) and can be quite fiddly. On
  e.g. HKUST HPC3 I needed some older compilers to play well with XIOS because
  the newer gcc compilers (version after 6) seems to be quite strict with the
  c++ code checking. To do this, I did
  
  .. code-block:: bash
  
    wget http://mirror.koddos.net/gcc/releases/gcc-5.4.0/gcc-5.4.0.tar.gz
    tar -xvzf gcc-5.4.0.tar.gz
    cd gcc-5.4.0/
    ./contrib/download_prerequisites
    cd ..
    mkdir gcc5.4
    cd gcc5.4
    ../gcc-5.4.0/configure --prefix=/scratch/PI/jclmak/custom_libs/gcc5.4/ --enable-languages=c,c++,fortran [--disable-multilib]
    make [-j4]
    make [check] install
  
  The first line grabs a packaged version of gcc, in this case ``5.4.0``; I
  chose the ``x.y.0`` version because I have had problems with the other
  versions with dependency issues with ``flex`` etc. (disclaimer: not checked
  overly rigourously because copmiling take soooo long). After unzipping, the
  4th line downloads the per-requisite libraries into the source folder (gcc
  official website highly recommends you **do not** compile the dependencies
  yourselves manually). 
  
  The 6th and 7th line follows the gcc official recommendation in doing the
  configuring and building **not** in the source directory; change the
  ``--prefix`` to the place where you want to store the libraries, headers and
  binaries. The ``--disable-multilib`` flag forces it to build a 64-bit one only
  (I needed that on the particularly computer). Calling ``make`` will take
  absolutely ages (order of hours, can speed up with giving more CPUs through
  the ``-j`` flag) because it will do a bootstrap build (building needed
  dependencies from existing compiler then using the build tools to build the
  target compiler, then sorting out the dependencies with the newly built
  compilers); can disable but not recommended. 
  
  Once the compilers are built then proceed as usual. Of course if you are on a
  cluster you probably could/should get someone else to do this...
  
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

  # if you want dynamic libraries then have this
  export LD_LIBRARY_PATH=$BD/install/lib:$LD_LIBRARY_PATH
  
  # if you want static libraries then have these
  export C_INCLUDE_PATH=$BD/install/include:$C_INCLUDE_PATH
  export CPLUS_INCLUDE_PATH=$BD/install/include:$CPLUS_INCLUDE_PATH
  export LIBRARY_PATH=$BD/install/lib:$LIBRARY_PATH

  # not strictly required, only for overriding preferences in search for binary
  export PATH=$BD/install/bin:$PATH

For my code testing it doesn't really matter too much whether the libraries are
compiled as static or dynamic because I'm not hugely concerned about performance
and stability, but static is probably safer. Set the above variables by doing
``source extra_variables``; upon closing the terminal the variables will be
flushed. Some of these may want to be added to ``~/.bashrc`` for convenience.
The instructions below attempts to build shared rather than static libraries,
and somewhat depends ``LD_LIBRARY_PATH`` variable being set (with the added
bonus that the ``ldd`` command provides an extra check whether the correct
libraries are being called). Suggestions on how to build the packages without
setting ``LD_LIBRARY_PATH`` or build static packages are given below (using
``LD_LIBRARY_PATH`` can be dangerous, see e.g., `here
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
and binaries to the folder specified in ``$BD``; if you have ``sudo`` access you
install it to ``/usr/local``, although I have found this can be very problematic
if you need to remove the libraries (I've bricked my computer once)... The
sub-directories in the folder are:

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

  The ``./configure prefix=`` step requires an absolute (not relative) path for
  the installation folder.
  

zlib and HDF5
-------------

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
``$BD`` as defined (don't include ``-fPIC`` or ``--enabled-shared`` if you want
the libraries to be static):

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
  CFLAGS=-fPIC ./configure --enable-shared --enable-fortran --enable-cxx
  --prefix=$BD/install/
  make -j 2
  make check install
  cd $BD
  
.. note::
  
  If ``LD_LIBRARY_PATH`` is set then zlib should be detected by the HDF5
  install. If not, consider including the commented out ``CPPFLAGS`` and
  ``LDFLAGS`` line (the ``--with-zlib`` command no longer works in the newer
  HDF5).
  
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
  for a guide.

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
netcdf-fortran v4.4.4 (don't include ``-fPIC`` or ``--enabled-shared`` if you
want the libraries to be static):

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
