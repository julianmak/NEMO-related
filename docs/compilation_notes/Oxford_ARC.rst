.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   
.. _sec:oxford:

Oxford ARC compilation
======================

The build uses NEMO 3.7/4.0 + XIOS 2.0 as the example. For installing other
versions, extrapolate from the other notes.

.. code-block:: bash

  export $BD=/home/julian/testing/gcc4.9-builds # CHANGE ME

  export C_INCLUDE_PATH=$BD/install/include:$C_INCLUDE_PATH
  export CPLUS_INCLUDE_PATH=$BD/install/include:$CPLUS_INCLUDE_PATH
  export LIBRARY_PATH=$BD/install/lib:$LIBRARY_PATH
  export LD_LIBRARY_PATH=$BD/install/lib:$LD_LIBRARY_PATH

subsection
----------


