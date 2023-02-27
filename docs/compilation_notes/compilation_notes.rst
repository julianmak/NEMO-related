.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

NEMO compilation notes
======================

These are just my own notes for compiling NEMO on a variety of clusters and
computers in a public place largely so I can look it up as long as I have
internet; if it happens useful for you, great! Please consult the `NEMO forge
page <http://forge.ipsl.jussieu.fr/nemo/wiki/Users/ModelInstall>`_ for the
official details.

While it is fairly straightforward on a supported cluster/supercomputer (e.g.
try `NOCL ARCHER guide
<https://nemo-nocl.readthedocs.io/en/latest/work_env/archer.html>`_) it can be a
bit temperamental on a local machine largely down to library and compiler
compatibility. The following notes are what I did to get XIOS and NEMO compiling
and running, and will display commands with ``gcc4.9`` compilers (which is my
default for other reasons). Extra things that need to be modified for other
compilers I have tested will be given accordingly (see the top of the individual
pages as to which compilers I have tested the notes with).

I added the following to my ``~/.bashrc`` so as to override the default
compilers I had (change these if need be):

.. code-block:: bash

  export CC=/usr/bin/gcc-4.9
  export CXX=/usr/bin/g++-4.9
  export FC=/usr/bin/gfortran-4.9
  export F77=/usr/bin/gfortran-4.9
  export CPP=/usr/bin/cpp-4.9

.. toctree::
   :maxdepth: 2
   :caption: Compilation notes:
   
   nemo36.rst
   nemo37.rst
   nemo40.rst
   nemo42.rst
   Oxford_ARC.rst
   HKUST_HPC2.rst
   HKUST_HPC3.rst
   packages.rst


