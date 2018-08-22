.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _sec:pyCDFTOOLS:

pyCDFTOOLS
==========

For various reasons (mostly personal preference and forcing myself to write in
Python) I made a translation of sorts of `CDFTOOLS
<https://github.com/meom-group/CDFTOOLS>`_ in Python. pyCDFTOOLS I think is:

* slightly more flexible, e.g., no need to recompile if variable name changes between files
* saves on the creation and reading of files
* everything done within Python, rather than Fortran and MATLAB say
* marginally more up-to-date, e.g. dealings with TEOS-10 equation of state

On the other hand, it is

* not as complete, because I only translated ones that I needed (see `here <https://github.com/julianmak/NEMO-related/tree/master/pyCDFTOOLS>`_ for list)...
* not as established and probably slightly error prone
* not as fast (though things that I could not vectorise I used JIT to speed up the looping)
* not NEMO code compliant (CDFTOOLS is designed to conform to NEMO code conventions)

An additional criticism I have is that I wrote pyCDFTOOLS more like
Fortran/MATLAB and not making full use of the Python functionalities (e.g.,
Panda and so forth). I have some idea how I might get it to work but watch this
space...

The routine naming conventions of the programs are basically the same as
CDFTOOLS (see `MEOM page <http://meom-group.github.io/doc/CDFTOOLS/>`_). All
codes with the prefix ``cdf`` are based on CDFTOOLS; all errors are entirely
mine (any things I did change are commented in the code).

Grab it with:

.. code-block:: bash

  git clone https://github.com/julianmak/NEMO-related
  
Some slightly more configuration/model specific Python scripts and notebooks are
in other folders (e.g., ``GYRE`` and ``ORCA``). I tend to just do

.. code-block:: bash

  cd GYRE
  rsync -arv ../pyCDFTOOLS .
  
which then means the scripts and notebooks within the folder have access to the
module, and it separates out a version that I do testing on.

CDFTOOLS itself depends on the following packages (the things I think that come
as standard are omitted):

* numba (for JIT to speed up loops)
* numpy (for tools)
* netCDF4 (for reading)
* scipy (for the occasional times when a MATLAB file is read)

The configuration specific programs depend additionally on Matplotlib and a
whole load of other ones for the ORCA configuration; see the relevant pages. I
installed most of the things through `Anaconda
<https://www.anaconda.com/download/#linux>`_; see the :ref:`Python <sec:python>`
page here for my notes on these.

Use these scripts at your own risk and feel free to modify them (rights etc. as
stated in the license and in line with the CDFTOOLS one). For comparison
purposes you may also want to grab CDFTOOLS to compare results (see the
`CDFTOOLS <https://github.com/meom-group/CDFTOOLS>`_ page):

.. code-block:: bash

  git clone https://github.com/meom-group/CDFTOOLS
  
.. note::
  
  The programs I have uploaded I was satisfied enough with the tests I have
  done, but don't just take my word for it :-)



