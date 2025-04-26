.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _sec:python:

Python / Anaconda notes
=======================

At some point I encountered some problem with plotting data in MATLAB (to do
with the tripolar grid meaning the co-ordinate files were not monotonic so
MATLAB hated it), and I went over to Python because the Cartopy and Iris
packages lets me do data projection and plotting in different projects fairly
easily. Here are some notes for Python and Anaconda which may be useful (the
latter might be useful for getting the libraries that NEMO and XIOS need).

Anaconda
--------

Most of these are taken from the `official conda manual
<https://conda.io/docs/index.html>`_ The installation for conda (or the lighter
version miniconda) is somewhat dependent on the OS and the instructions are
`here <https://conda.io/docs/user-guide/install/index.html>`_ You end up
downloading a bash file that you run in the terminal, and from there you can
accept and change some of the settings accordingly. No administrator rights
should be required, though it does mean the installed packages may not be
shareable. The installation will ask if you want to add to your ``$PATH``
variable, which I accepted (it means the some of the anaconda based binaries
take precedence over the system ones).

One conda is installed, I would recommend creating an environment so that if
damage is to occur, it is only within the environment which may be deleted
easily without touching other things. The creation, entering and leaving of the
environment is done by:

.. code-block:: bash

  >> julian@psyduck:~/$ conda create -n nemo python=3.6
  ...
  >> julian@psyduck:~/$
  >> julian@psyduck:~/$ source activate nemo
  >> (nemo) julian@psyduck:~/$ 
  >> (nemo) julian@psyduck:~/$ source deactivate
  >> julian@psyduck:~/$
  
The first command creates and environment called ``nemo`` that uses python 3.6,
and the other commands are self explanatory. An environment may be removed by
issuing the command

.. code-block:: bash

  conda remove --name nemo --all
  
Packages are installed through (make sure you are in an environment first)

.. code-block:: bash

  conda install netcdf
  conda install -c conda-forge netcdf-fortran
  
Some packages need to be searched for in the forge.

Note that while the environment is active some commands take precedence over
others, and a bit of care is needed to make sure the ones you intend to call
really are the ones that are called (e.g. my mercurial command ``hg`` seems to
be overwritten on my machine when I am in my environment). Check with things
like ``which python`` for example which shows which binary the command
``python`` is actually calling. 

Python
------

I mostly develop code in a notebook because I am too heavily influenced by
MATLAB. Notebooks (in particular with Jupyter) lets you write code within cells
that you run and see outputs then and there which is what I am used to. Later on
I do write code in a text editor when I have more specific things I want need to
do.

I normally do the following to get what I need. Within the environment:

.. code-block:: bash

  conda install scipy
  conda install numpy
  conda install matplotlib
  conda install jupyter
  conda install -c conda-forge cartopy
  conda install -c conda-forge iris

I normally install NetCDF as well. ``Numpy`` and ``scipy`` gives the number
crunching stuff I normally need. ``Matplotlib`` gives most of the plotting
capabilities. ``Cartopy`` and ``iris`` are the map and projection packages, and
``jupyter`` is the notebook stuff. To trigger the notebook, I normally do from a
terminal

.. code-block:: bash

  jupyter notebook 2>/dev/null &
  
just to suppress the terminal outputs. The notebook opens in a browser and you
do coding in there (I think there is another software that lets you open and
edit notebooks somewhere else though I've never used it); it's basically
``ipython`` but in a browser. Note that just closing the tabs does not
necessarily close the notebook; you need to do ``files>>close and halt``. Also,
just because the relevant pages are closed in the browser does not mean the
notebook server is shutdown either; you need to click ``logout`` on the top
right corner (assuming you are not using a custom theme which suppresses that).
To kill it in the terminal, either find the job through ``jobs`` and use ``kill
%n`` or do

.. code-block:: bash

  jupyter notebook list
  >> Currently running servers:
  >> http://localhost:8888/?token=7774a1ace4c2a0a1e098a5900f30c67310074a7250bd6c0d :: /home/julian/GitRepo/pydra/wrapper
  >> http://localhost:8889/?token=00b793728b03e2536b5a07a793bbd2a9fc1342469f3cf28d :: /home/julian/Documents/NEMO
  
  jupyter notebook stop 8888
  jupyter notebook list
  >> Currently running servers:
  >> http://localhost:8889/?token=00b793728b03e2536b5a07a793bbd2a9fc1342469f3cf28d :: /home/julian/Documents/NEMO


Some Python banana skins
------------------------

The big banana skin with Python to watch out for is that indexing starts at
``0`` (rather than ``1`` in MATLAB), and index slicing normally omits the last
entry, e.g.

.. code-block:: python

  x_vec = [1, 2, 3, 4, 5, 6]
  x_vec[0:-1]
  >> [1, 2, 3, 4, 5]
  x_vec[1:4]
  >> [2, 3, 4]
  x_vec[0::]
  >> [1, 2, 3, 4, 5, 6]
  x_vec[-1]
  >> 6
  x_vec[-2]
  >> 5

Contrast this to MATLAB which would be

.. code-block:: MATLAB

  x_vec = [1, 2, 3, 4, 5, 6]
  x_vec(0:end-1)
  >> 1, 2, 3, 4, 5
  x_vec(2:4)
  >> 2, 3, 4
  x_vec(:)
  >> 1, 2, 3, 4, 5, 6
  x_vec(end)
  >> 6
  x_vec(end - 1)
  >> 5
  
Another banana skin with python is that data is not necessarily copied when
defining new variables. For example:

.. code-block:: python

  x_vec = [1, 2, 3, 4, 5, 6]
  y_vec = x_vec
  y_vec[0] = 2
  y_vec
  >> [2, 2, 3, 4, 5, 6]
  x_vec
  >> [2, 2, 3, 4, 5, 6]
  
This is especially dangerous if you, like me, do the following in MATLAB:

.. code-block:: MATLAB

  x_vec = zeros(6)
  y_vec = x_vec
  z_vec = x_vec
  
If you really mean to do a copy, do the following:

.. code-block:: python

  from copy import deepcopy
  x_vec = [1, 2, 3, 4, 5, 6]
  y_vec = x_vec
  z_vec = deepcopy(x_vec)
  y_vec[0] = 2
  y_vec
  >> [2, 2, 3, 4, 5, 6]
  x_vec
  >> [2, 2, 3, 4, 5, 6]
  z_vec
  >> [1, 2, 3, 4, 5, 6]
  
Python is really slow with loops, so the more vectorising commands you can use,
the better! If you have routines that you have to use loops in (e.g.
transformation of data from Cartesian co-ordinates to density co-ordinates
through binning into density bins), then consider using ``cypthon`` (write code
in C but call it through Python), ``f2py`` (same but for Fortran), or
``numba``/``JIT`` (compile and run loops, usually on the order of 200 speed up;
restricted to fairly low level commands).
