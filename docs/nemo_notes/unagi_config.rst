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
the north (see :cite:`Munday-et-al15`). The vertical diffusivity is enhanced
through brute force hacking in the ``zdfini.F90`` and ``step.F90`` files.
Analogous approaches using the TKE scheme (e.g. enhancing the minimum energy) or
through providing an extra ``diffkr.nc`` file say could be done but was not
pursued here. See e.g. :cite:`Abernathey-et-al11` for alternative model
formulations.

The default is:

* flat everywhere and a ridge in the middle of the channel extending up to half the depth of the domain
* fixed sinusoidal wind stress with some peak wind stress value :math:`\tau_0`
* SST restoring
* linearly varying temperature profile at the surface with :math:`e`-folding depth of 1000 metres
* linear friction

How to get the model running
----------------------------

Custom analysis scripts
-----------------------

.. _sec:build_model:

Notes: building UNAGI
---------------------

The following approach is strictly for NEMO models beyond v3.6, where one can
build a customised model through providing a ``domcfg.nc``, which is the main
goal of the following text. The details are given below are what I did for the
idealised channel model UNAGI.

The biggest obstacle in generating the appropriate ``domcfg.nc`` file for me was
in transferring the code that modifies the vertical spacing variables
``e3t/u/v/w`` to have a partial cell description (I first tried to brute force
it by writing from scratch a file that provides all the relevant variables
needed in the ``domcfg.nc``; see for example the input required in ORCA2). Which
this in mind, I fell back to using the NEMO native :cite:`MadecImbard96` grid
and the ``TOOLS/DOMAINcfg`` package, as follows:

1. to fill in
2. to fill in
3. to fill in


.. bibliography:: ../refs.bib
   :filter: docname in docnames



