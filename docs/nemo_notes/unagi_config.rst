.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

UNAGI: custom channel model
===========================

Brief overview and sample outputs
---------------------------------

UNAGI (naming based on EEL which was original due to Marina Levy) is a
re-entrant :math:`\beta`-plane channel with temperature as the thermodynamic
variable that's largely based on Dave Munday's MITgcm channel model reported in
:cite:`Munday-et-al15` as an idealised model to the Antarctic Circumpolar
Current. The model at present takes in some to-be specified bathymetry, wind
stress profile and initial state (which may be customised accordingly within the
``gen_UNAGI_fields.ipynb``).

With the choice of SST restoring over the surface layer, to maintain a sensible
thermocline the vertical tracer diffusivity is enhanced in a sponge region to
the north (see :cite:`Munday-et-al15`). See e.g. :cite:`Abernathey-et-al11` for
alternative model formulations. Some model set up choices:

* relatively long re-entrant zonal channel, no topography except ridge in the middle of the channel extending up to half the depth of the domain
* fixed sinusoidal wind stress with some peak wind stress value :math:`\tau_0`
* SST restoring (a relatively hard restoring, the ``rn_dqdt`` value in ``namelist_cfg`` has been amplified by a factor of 2)
* linearly varying temperature profile at the surface with :math:`e`-folding depth of 1000 metres
* linear friction
* linear EOS with only temperature as the variable
* sponge region to the north where vertical diffusivity is amplified by a factor of 250 from the background value of :math:`10^{-5}\ \mathrm{m}^2\ \mathrm{s}^{-1}`

The diagram below shows an output from the 10km resolution model with biharmonic
tracer diffusion and no eddy parameterisation, showing a rich eddying field.

  .. figure:: figs/unagi_R010_xi.png
    :width: 90%
    :align: center
    :alt: unagi_R010_xi
    :name: unagi_R010_xi
    
    Vertical component of vorticity (in units of :math:`\mathrm{s}^{-1}`) at the surface from UNAGI at 10km resolution. Click `here <https://i.imgur.com/bT37Mo4.gifv>`_ for an animation.

How to get the model running
----------------------------

[TO BE ADDED]

Custom analysis scripts
-----------------------

[TO BE ADDED]

.. _sec:build_model:

Notes 1: building a custom model
--------------------------------



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

Notes 2: hacking NEMO to get UNAGI
----------------------------------

That two main things that needed hacking into NEMO for UNAGI are the vertical
tracer diffusion (in the sponge region to the north) and possible combination
with the GEOMETRIC parameterisation. [TBC, 15 Apr 2019]

.. bibliography:: ../refs.bib
   :filter: docname in docnames



