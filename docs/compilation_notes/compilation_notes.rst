.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

NEMO compilation notes
======================

Personally I prefer doing small bits of code testing on a small configuration
(normally GYRE in NEMO) so I have tried to get NEMO working on a local machine,
largely following the instructions from the `NEMO forge
<http://forge.ipsl.jussieu.fr/nemo/wiki/Users/ModelInstall>`_; please consult
that page first as the information on there should be considered the
authoritative version. 

While it is fairly straightforward on a supported cluster/supercomputer (e.g.
first try following `NOCL ARCHER guide
<https://nemo-nocl.readthedocs.io/en/latest/work_env/archer.html>`_) it can be a
bit temperamental on a local machine to do with compiler compatibilities. The
following notes are what I did to get XIOS and NEMO compiling and running, and
will display commands with ``gcc4.9`` compilers (which is my default for other
reasons). Extra things that need to be modified for other compilers I have
tested will be given accordingly (see the top of the individual pages as to
which compilers I have tested the notes with).

I added the following to my ``~/.bashrc``:

.. code-block:: bash

  export CC=/usr/bin/gcc-4.9
  export CXX=/usr/bin/g++-4.9
  export FC=/usr/bin/gfortran-4.9
  export F77=/usr/bin/gfortran-4.9
  export CPP=/usr/bin/cpp-4.9
  
which overrides the default ``gcc5.4`` on my computer.

.. toctree::
   :maxdepth: 2
   :caption: Compilation notes:
   
   nemo36.rst
   nemo37.rst
   nemo40.rst
   packages.rst


