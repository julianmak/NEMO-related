.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _sec:other-pack:

Other packages
==============

The following packages are needed for NEMO and XIOS and they may need to be
installed or configured accordingly.

HDF5 
----

Check whether HDF5 exists first (may still need to be installed again for
compatibility reasons). ``h5copy`` is the command that should exist if HDF5 is
installed.

.. code-block:: bash
  
  whereis h5copy
  h5copy --version
  

NetCDF4
-------

Check whether NetCDF4 exists first (may still need to be installed again for
compatibility reasons). ``nc-config`` is the command that should exist if
NetCDF4 is installed, and checks where it is installed and what compilers were
used to build it.

.. code-block:: bash
  
  nc-config all

MPI setup
---------

