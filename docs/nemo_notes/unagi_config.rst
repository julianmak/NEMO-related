.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

UNAGI: custom channel model
===========================

Brief overview and sample outputs
---------------------------------

UNAGI (naming based on EEL which was original due to Marina Levy) is a re-entrant :math:`\beta`-plane channel with temperature as the thermodynamic variable that is largely based on Dave Munday's MITgcm channel model reported in :cite:`Munday-et-al15`, as an idealised model to the Antarctic Circumpolar Current. The model at present takes in some to-be specified bathymetry, wind stress profile and initial state, which may be customised accordingly within the ``gen_UNAGI_fields.ipynb``, which may be found at the `host GitHub repository <https://github.com/julianmak/UNAGI_NEMO>`_.

With the choice of SST restoring over the surface layer, to maintain a sensible thermocline the vertical tracer diffusivity is enhanced in a sponge region to the north (see :cite:`Munday-et-al15`). See e.g. :cite:`Abernathey-et-al11` for alternative model formulations. Some model set up choices:

1. relatively long re-entrant zonal channel, no topography except ridge in the middle of the channel extending up to half the depth of the domain
2. fixed sinusoidal wind stress with some peak wind stress value :math:`\tau_0`
3. SST restoring (a relatively hard restoring, the ``rn_dqdt`` value in ``namelist_cfg`` has been amplified by a factor of 2)
4. linearly varying temperature profile at the surface with :math:`e`-folding depth of 1000 metres
5. linear friction
6. linear EOS with only temperature as the thermodynamic variable
7. sponge region to the north where vertical diffusivity is amplified by a factor of 250 from the background value of :math:`10^{-5}\ \mathrm{m}^2\ \mathrm{s}^{-1}`

The diagram below shows the surface relative vorticity (in units of :math:`\mathrm{s}^{-1}`) from the 10km resolution model with biharmonic tracer diffusion and no eddy parameterisation, associated with a rich eddying field. Click `here <https://i.imgur.com/bT37Mo4.gifv>`_ for an animation.

  .. figure:: figs/unagi_R010_xi.png
    :width: 90%
    :align: center
    :alt: unagi_R010_xi
    :name: unagi_R010_xi

How to get the model running
----------------------------

If you just want things to work then try the `zenodo repository <http://dx.doi.org/10.5281/zenodo.8002828>`_, which has all the NEMO modified sources files and model input files required. Some sample analysis code are given in `host GitHub repository <https://github.com/julianmak/UNAGI_NEMO>`_.

It is also fairly quick to recreate the forcing files from scratch, and is likely more informative for making your own models. The relevant notebook is ``gen_UNAGI_fields.ipynb``, given also in `host GitHub repository <https://github.com/julianmak/UNAGI_NEMO>`_. The code can almost be run straight except for one step; see below notes.

.. _sec:build_model:

Building the custom model
-------------------------

The following approach is strictly for NEMO models beyond v3.6, where one can build a customised model through providing a ``domcfg.nc``, which is the main goal here. The details are given below are what I did for the idealised channel model UNAGI; see `here <https://github.com/julianmak/NEMO-related/blob/master/UNAGI/readme_of_sorts.txt>`_ for a step-by-step guide of how I did it.

The biggest obstacle in generating the appropriate ``domcfg.nc`` file for me was in transferring the code that modifies the vertical spacing variables ``e3t/u/v/w`` to have a partial cell description. I first tried to brute force it by writing from scratch a file that provides all the relevant variables needed in the ``domcfg.nc``; see for example the input required in ORCA2. I gave up after a while and fell back to using the NEMO native :cite:`MadecImbard96` grid and the ``TOOLS/DOMAINcfg`` package, as follows:

1. in an external folder (e.g., ``~/Python/NEMO/UNAGI``), create the bathymetry data through a program of your choice (e.g. `Python <https://github.com/julianmak/NEMO-related/blob/master/UNAGI/gen_NEMO_UNAGI_fields.ipynb>`_), and output it as a netCDF file (e.g. ``bathy_meter.nc``)
2. link/copy it as ``bathy_meter.nc`` (the tool requires that specific naming) into the ``TOOLS/DOMAINcfg`` that comes with NEMO 
3. modify the ``namelist_cfg`` file accordingly for the horizontal and vertical grid spacing parameters (see :ref:`here <sec:nemo_packages>` for usage and compiling notes), and the one I used for this model is given as an example in that packages page
4. a ``domcfg.nc`` should result (if not, see ``ocean.output`` for messages), copy it back into the working folder in step 1
5. open ``domcfg.nc`` and use those variables to create the ``state.nc`` and ``forcing.nc`` file again in the program of your choice (this is mostly to keep consistency; I did it in `Python <https://github.com/julianmak/NEMO-related/blob/master/UNAGI/gen_NEMO_UNAGI_fields.ipynb>`_)
6. copy the ``domcfg.nc``, ``state.nc`` and ``forcing.nc`` (I prefixed them with something, e.g. ``UNAGI_domcfg_R010.nc``) and modify the ``namelist_cfg`` accordingly, e.g.

::

  !-----------------------------------------------------------------------
  &namrun        !   parameters of the run
  !-----------------------------------------------------------------------
     cn_exp      =      "UNAGI" !  experience name
     nn_it000    =          1   !  first time step
     nn_itend    =       8640   !  last  time step
     nn_date0    =      10101   !  
     nn_leapy    =         30   !  Leap year calendar (1) or not (0)
     ln_rstart   = .false.   !  start from rest (F) or from a restart file (T)
        nn_euler    =    1            !  = 0 : start with forward time step if ln_rstart=T
        nn_rstctl   =    0            !  restart control ==> activated only if ln_rstart=T
        !                             !    = 0 nndate0 read in namelist
        !                             !    = 1 nndate0 check consistancy between namelist and restart
        !                             !    = 2 nndate0 check consistancy between namelist and restart
     nn_stock    =       8640   !  frequency of creation of a restart file (modulo referenced to 1)
     nn_write    =       8640   !  frequency of write in the output file   (modulo referenced to nn_it000)
  /
  !-----------------------------------------------------------------------
  &namcfg     !   parameters of the configuration   
  !-----------------------------------------------------------------------
     ln_read_cfg = .true.    !  (=T) read the domain configuration file
        !                    !  (=F) user defined configuration  ==>>>  see usrdef(_...) modules
        cn_domcfg = "domcfg_UNAGI" ! domain configuration filename
  /
  ...

That is more or less it. Once you can build the domain variables the model will
at least run and the rest is more to do with experimental design.

Hacking NEMO to get UNAGI
-------------------------

That two main things that needed hacking into NEMO for UNAGI are the vertical tracer diffusion (in the sponge region to the north) and possible combination with the GEOMETRIC parameterisation, the latter could be found `here <https://nemo-related.readthedocs.io/en/latest/GEOMETRIC/geometric.html>`_. For the vertical tracer diffusion given in ``zdfphy.f90``, I hacked an existing variable so that it is dual use to give a specified meridional profile in the vertical diffusivity; search for the variable ``rn_avt_amp`` in `this file <https://github.com/julianmak/UNAGI_NEMO/blob/main/MY_SRC/zdfphy.F90>`_ to see how I did it.

An extra hack I did was to shut off a default warning in ``ldftra.f90`` (with modifications to ``trazdf.f90`` for completeness) that biharmonic tracer diffusion cannot be used with when the GM scheme (``ldfeiv``) is used. Normally if you are using GM you also use isoneutral diffusion rather than biharmonic diffusion, but for my case I do intend on having that specific combination.

.. bibliography:: ../refs.bib
   :filter: docname in docnames



