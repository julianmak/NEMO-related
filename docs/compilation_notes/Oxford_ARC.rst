.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   
.. _sec:oxford:

Oxford ARC compilation
======================

.. note::

  To add: SLURM submission script, post processing script, ``xios_server.exe`` testing

The build uses NEMO 3.7/4.0 + XIOS 2.0 as the example. For installing other
versions, extrapolate from the other notes.

Annoying (!) everything basically works out of the box! This has not been the
usual experience I have with XIOS and NEMO...

Running NEMO on the ARC
-----------------------

The system uses SLURM and the key commands are

* ``sbatch [submit_nemo]``: submits the job detailed in ``submit_nemo`` (see below) 
* ``scancel [job ID]``: cancel the job
* ``sinfo``: check status of queues available
* ``squeue -u $USER``: check job info for ``$USER``

``sbatch`` could be used with arguments but I am going to have everything within
``submit_nemo`` itself. The generic one I use is as follows (I have some ASCII
art in there because I got bored at some point):

.. code-block:: bash

  #!/bin/bash

  # NOTE: Lines starting with "#SBATCH" are valid SLURM commands or statements,
  #       while those starting with "#" and "##SBATCH" are comments.  Uncomment
  #       "##SBATCH" line means to remove one # and start with #SBATCH to be a
  #       SLURM command or statement.


