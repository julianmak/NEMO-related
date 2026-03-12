.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   
.. _sec:hkusthpc4:

HKUST HPC4 compilation
======================

The build uses NEMO 4.0/4.2 + XIOS 2.5 as the example. For installing other versions, extrapolate from the other notes.

HKUST HPC4 is a cluster with SLURM. The usual ``module load/swap/purge/unload/list/avail`` works here, and the system so far is built by default with ``gcc11.4``. The thing that is different is that it HPC4 uses `spack <https://spack.io/>`_ as a DIY approach for building libraries with specified compilers, intended as a way to avoid clashes and dependency hells (cf. Anaconda).

My TL;DR is that I don't find it dependency hells are completely avoided, but there are certain things it does package up pretty well and can be done in one go, namely all the dependencies under the :ref:`packages <sec:other-pack>` page. The approach I took is to use spack to build the usual OpenMPI, HDF5 and NetCDF4 (with parallel options), then build XIOS manually, then build NEMO as before. Do a minor detour of spack first, then the remaining things are reasonably straightforward.

.. note::

  XIOS is nominally available through spack, but it simply does not work for me for various reasons that I have identified, such as the dependency list having gaps (e.g. subversion, incompatibility with compiler versions etc.), the pull links were wrong (e.g. calling ``svn`` on links that no longer exists), the architecture files do not adapt accordingly etc. I am sure it could work, but I decided it was easier and more worthwhile to build it manually.

Spack
-----

The default HKUST page is `here <https://hkust-hpc-docs.readthedocs.io/latest/software/software-support-overview.html>`_. The page recommends using ``spack-edge`` although the ones below use ``spack``: the former is more updated while the latter is frozen and has some clashes. The usage is similar although behaviour might change.

The default spack is already activated and located at ``${SPACK_ROOT}``. The first run of it sets up a few things with a build by default in ``~/.spack``

.. note::

  To use the "edge" version, do the following

  .. code-block:: bash

    source /opt/shared/.spack-edge/dist/bin/setup-env.sh -y
    
  and maybe put it in the ``.bashrc`` so every time at log in it is picking up the intended version of spack. The default build is then in ``~/.spack-edge``.
  
For me I overrode it by setting

.. code-block:: bash

  export SPACK_USER_CONFIG_PATH=/project/miffy/custom_libs_spack
  export SPACK_USER_CACHE_PATH=/project/miffy/custom_libs_spack
  
in my ``.bashrc``. 

.. note::

  Probably should use ``spack env create -d /project/miffy/custom_libs_spack`` or similar instead. Do that test for later. (Probably call it ``/project/miffy/custom_libs_spack/gcc11.5_xios_libs``, so people can create their own libraries housed in ``custom_libs_spack`` rather than rely on a default.)

Activate environments
^^^^^^^^^^^^^^^^^^^^^
  
In this case I would create an environment called ``xios`` and go into by

.. code-block:: bash

  spack env create xios
  spack env activate xios
  
and you would get out of it by ``spack env deactivate`` (or ``despacktivate``). The building that is done is then done within that environment.

Compilers
^^^^^^^^^

First thing is to see what compilers have been picked up, through

.. code-block:: bash

  spack compiler list
  
I ended up settling on ``gcc@11.4.1`` to have the compiler as close to my laptop as possible, where I know XIOS and NEMO will build. What I do here is to go into where my ``spack.yaml`` lives (in this case ``/project/miffy/custom_libs_spack/environments/xios``) and modify it manually to have something like

.. code-block:: bash

  # This is a Spack Environment file.
  #
  # It describes a set of packages to be installed, along with
  # configuration settings.
  spack:
    packages:
      all:
        compiler: [gcc@11.4.1]
    specs:
    - subversion
    ...
    
to be very specific and force all my packages in this environment to compile with this specific compiler. Could have done this by hand with e.g. ``spack add openmpi%gcc@11.4.1`` (where the ``%`` sign is for specifying the compiler), but that modification to the ``spack.yaml`` file is already supposed to force that.

Adding libraries to environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Note here that adding packages does not mean building them (yet). We can add libraries by

* ``spack add openmpi%gcc@11.4.1`` to specify compiler
* ``spack add netcdf-c ^openmpi@5.0.3%gcc@13.2.0`` to specify dependency on other packages through ``^``
* modify the ``spack.yaml`` file directly (the above commands do that also)

My resulting ``spack.yaml`` that works in the end looks like the following:

.. code-block:: bash

  # This is a Spack Environment file.
  #
  # It describes a set of packages to be installed, along with
  # configuration settings.
  spack:
    packages:
      all:
        compiler: [gcc@11.4.1]
      perl:
        buildable: true
    specs:
    - subversion
    - perl@5.35+cpanm
    - perl-uri
    - openmpi@4.1.2
    - pmix@3.2.1
    - netcdf-c
    - netcdf-fortran
    view: true
    concretizer:
      unify: true
      
The choices I made here are as follows:

* I forced ``perl`` to be buildable, overriding the system default

  - there is a good reason for the default setting as set up by the cluster manager(s), but I needed to have more control over ``perl`` so I overrode that setting

* I need ``subversion`` to copy XIOS in since there is no system ``svn``

* ``perl`` with ``+cpanm`` option to install more packages if needed

* ``perl-uri`` because XIOS does call that

* ``openmpi@4`` to mirror my laptop

* ``pmix@3`` apparently plays better with ``openmpi@4`` then ``openmpi@3``

  - known problems with ``pmix@5`` not playing with ``openmpi@4``

* ``netcdf-c`` and ``netcdf-fortran`` will pull the relevant HDF5 and seems to build with paralell capabilities also
  
.. note::

  ``spack-edge`` may have a more complete version of ``perl`` suitable for the purposes below.
  
  ``spack-edge`` keeps pulling in a clashing compiler for me at various places, which is causing a ton of clashes with HDF5. Probably some default setting that needs to be suppressed, but can't find it; giving up on this for now.

Concretizing environment
^^^^^^^^^^^^^^^^^^^^^^^^

This also doesn't do the building yet, because it needs to sort out the dependency list first. For that we would do something like

.. code-block:: bash

  spack concretize -fU
  
which fixes the exact versions of various libraries to be built and locks it in (cf. containerizing like in docker images). The result is a whole bunch of things as a tree output to screen telling you which hash and versions of things are to be built, as well as generating a ``spack.lock`` where ``spack.yaml`` lives, which can then be used to install the exact same environment somewhere else.

The resulting ``spack.yaml`` and/or ``spack.lock`` can be used to create the exact same versions of the environment somewhere else, either by building directly with it, or intialising an environment, copying those files into that environment, and then doing the next step.

.. note::

  It might refuse to concretize if there are clashes, or it's finding things that are incompatible with settings. I had issues for example when I tried to install ``perl`` without the ``buildable: true`` above, because it then refuses to build the specified version of ``perl``. How these are fixed depends on the reason for the clash though...

Building and installing the libraries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The last step would just be

.. code-block:: bash

  spack install [-v]
  
to install all the packages in some sensible order, or you can manually do it I suppose by doing e.g. ``spack install subversion`` per package for whatever reason.

If it does fail then it will tell you where the sources files and the logs are located, and you would have to go into the relevant folders to try and fix the builds.

.. note::

  Mine broke at ``openssl`` because the installed ``perl`` executable was called ``perl5.35.0`` instead of ``perl``, so it defaulted to ``/usr/bin/perl``. This could be fixed manually by appropriate softlinks and/or adding to the ``$PATH`` variable. I ended up softlinking it and running ``spack install``, which tried skipped the ``perl`` build (because it was already built) and continued with no issues.
  
  If you however concretize it again then it will remove a whole bunch of stuff, which may need to be manually patched again. If there is no concretizing then the ``xios/.spack-env/view/`` gets left alone in that there are no removals at least.
  
I checked the NetCDF was built in the intended place with ``which nc-config``, and whether it was built with parallel capabilities by checking the output of ``nc-config --all``. Can check things accordingly by querying ``$LD_LIBRARY_FLAG``, ``$PATH`` and also doing ``which mpirun`` and ``mpirun --version`` accordingly to make sure things are as intended.

XIOS
----

This next bit basically then proceeds as usual like in the other cases, with some flag modifications rather than anything fundamental (in hindsight; it didn't feel like it was fundamental as I was doing it...)

First make sure to pull a sufficiently updated version of XIOS 2.5 in this case. The updated link I used was

.. code-block:: bash

  svn checkout -r 2628 https://forge.ipsl.jussieu.fr//ioserver/svn/XIOS2/branches/xios-2.5/
  
.. note::

  Older versions the XIOS2.5 (~r15xx) that spack was defaulting to had a whole load of C++ errors with "type mismatch", which seem to all disappear when using a more recent version. Probably a compiler mismatch issue, could not be bothered identifying precisely what the cause was.
  
The changes I to the ``*.fcm`` file were

* add ``/project/miffy/custom_libs_spack/environments/xios/.spack-env/view/[lib,inc,bin]`` to the relevant places

  - the above path is where the environment things are consolidated

* remove the ``-D_GLIBCXX_USE_CXX11_ABI=0`` and ``-std=c++11`` flags

  - didn't seem to need it for ``gcc11``
  - ``-std=c++11`` returns "ambiguous reference" errors
  
* add in ``-ffree-line-length-none`` for ``%BASE_FFLAGS`` for the usual reason

* add in ``-lpmix`` to ``%BASE_LD`` (``-lstdc++`` should already be there)

  - otherwise it crashes at the end at the linker stage because it could find PMIx related things

Could be very specify about the ``mpif90`` and ``mpicc`` versions, but mine were already defined properly it seems. If it builds, then you are probably in business.

NEMO
----

As below. I used NEMO 4.0 r14538 for testing. The arch file is as follows, usual things as far as I can tell.

.. code-block:: bash

  %NCDF_HOME           /project/miffy/custom_libs_spack/environments/xios/.spack-env/view
  %HDF5_HOME           /project/miffy/custom_libs_spack/environments/xios/.spack-env/view
  %XIOS_HOME           /project/miffy/XIOS/xios-2.5-r2628

  %NCDF_INC            -I%NCDF_HOME/include -I%HDF5_HOME/include
  %NCDF_LIB            -L%NCDF_HOME/lib -lnetcdff -lnetcdf
  %XIOS_INC            -I%XIOS_HOME/inc 
  %XIOS_LIB            -L%XIOS_HOME/lib -lxios

  %CPP                 cpp -Dkey_nosignedzero 
  %FC                  mpif90 -c -cpp 
  %FCFLAGS             -fdefault-real-8 -O3 -funroll-all-loops -fcray-pointer -ffree-line-length-none -fallow-argument-mismatch
  %FFLAGS              %FCFLAGS
  %LD                  mpif90
  %LDFLAGS             -lstdc++
  %FPPFLAGS            -P -C -traditional
  %AR                  ar
  %ARFLAGS             -rs
  %MK                  make
  %USER_INC            %XIOS_INC %NCDF_INC
  %USER_LIB            %XIOS_LIB %NCDF_LIB

* ``-Dkeynosignedzero`` in ``%CPP`` may replace the need to put ``key_nosignedzero`` in the ``cpp_GYRE_PISCES.fcm`` or similar? Put it in for safety anyway

* I put in ``-fallow-argument-mismatch`` to suppress a whole load of type errors associated with MPI commands. This is a thing with newer compilers for older code (``gcc`` version > 10 seems have stricter checks).

* I put in ``-lstdc++`` in ``%LDFLAGS`` because it was crashing at the end with linker errors with some stuff relating to XIOS library about undefined references.

Seems to compile and run in the usual way. Did not check performance, although it was clear that the same ``mpirun -np 4 ./nemo`` with and without the XIOS ``one_file`` option has a noticeable (by eye) performance. Could not be bothered to check whether this was because of node allocation , or can be alleviated by allocating some CPUs to ``xios_server``.



