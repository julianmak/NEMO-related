.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _sec:analysis tools:

analysis tools
==============

For various reasons (mostly personal preference and forcing myself to write in Python) I previously attemped to make a translation of sorts of `CDFTOOLS <https://github.com/meom-group/CDFTOOLS>`_ in Python, and it ended up being called `pyCDFTOOLS <https://github.com/julianmak/pyCDFTOOLS>`_. This will probably be renamed for reasons detailed later. [22 Apr 2025]

.. note ::

  There is one as part of this repository, but probably don't use that anymore, as the one that will be `here <https://github.com/julianmak/pyCDFTOOLS>`_ and/or one that will replace it is done in a more modern way.
  
There is some talk of a plan for the package. But in the meantime, the idea is that the python version would be:

* slightly more flexible, e.g., no need to recompile if variable name changes between files
* saves on the creation and reading of files
* everything done within Python, rather than Fortran and MATLAB say
* marginally more up-to-date, e.g. dealings with TEOS-10 equation of state

On the other hand, it will

* not NEMO code compliant (CDFTOOLS is designed to conform to NEMO code conventions, hence the eventual renaming)
* not be as complete, because I only translate ones that I needed (and to be honest there are a chunk of CDFTOOL routines that are not arguably a bit niche...)
* not as established and probably slightly error prone

An additional criticism I have is that I write my Python code more like Fortran/MATLAB code. This one is somewhat by design, because I have in mind a set of tools that is more like a calculator, rather than imbuing an object with all the subfunctions etc.

Updates may come on here at some point.



