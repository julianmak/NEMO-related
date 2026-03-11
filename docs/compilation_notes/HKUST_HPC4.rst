.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   
.. _sec:hkusthpc4:

HKUST HPC4 compilation
======================

The build uses NEMO 4.0/4.2 + XIOS 2.5 as the example. For installing other versions, extrapolate from the other notes.

HKUST HPC4 is a cluster with SLURM. The usual ``module load/swap/purge/unload/list/avail`` works here, and the system so far is built by default with ``gcc11.4``. The thing that is different is that it HPC4 uses `spack <https://spack.io/>`_ as a DIY approach for building libraries with specified compilers, intended as a way to avoid clashes and dependency hells (cf. Anaconda).

My TL;DR is that I don't find it avoids dependency hells completely, but there are certain things it does package up pretty well, namely all the dependencies under the :ref:`packages <sec:other-pack>` page can be done in basically one go. The approach I took is to use spack to build the usual OpenMPI, HDF5 and NetCDF4 (with parallel options), then build XIOS manually, then build NEMO as before. Do a minor detour of spack first, then the remaining things are reasonably straightforward.

.. note::

  XIOS is nominally available through spack, but it simply does not work for me for various reasons that I have identified, such as the dependency list has gaps (e.g. subversion, incompatibility with compiler versions etc.), the pull links were wrong (e.g. calling ``svn`` on a link that no longer exists), the architecture files do not adapt accordingly etc. I am sure it could work, but I decided it was easier to build it manually.

Spack
-----
 
Compilers
---------







