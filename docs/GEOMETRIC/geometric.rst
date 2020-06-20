.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   
.. role:: strike

GEOMETRIC outline
=================

**GEOMETRIC** (*Geometry and Energetics of Ocean Mesoscale Eddies and Their
Rectified Impact on Climate*) is an approach to representing the unresolved
turbulent eddies in ocean climate models, first derived in
:cite:`Marshall-et-al12`. `David Marshall
<https://www.marshallocean.net/geometric>`_'s page has an excellent outline and
summary of GEOMETRIC, so this page will focus on outlining the details
relating to the NEMO implementation.

The implementation of GEOMETRIC was done in NEMO by providing a new module
``ldfeke.f90`` and adding appropriate calls and variables to ``ldftra.f90``,
``step.f90`` ``step_oce.f90`` and ``nemogcm.f90``. This was initially done in
SVN version 8666, which is somewhere between the 3.6 stable and 4.0 beta, by
myself and `Gurvan Madec
<https://scholar.google.com/citations?user=Ewgtd20AAAAJ&hl=fr>`_ back in
November 2017. The current implementation of GEOMETRIC is what may be considered
GM-based :cite:`GentMcWilliams90` and follows the prescription described in
:cite:`Mak-et-al18`. The GEOMETRIC scaling gives :math:`\kappa_{\rm gm} = \alpha
E (N / M^2)` (see below for symbol definitions). While :math:`\alpha` is
prescribed and :math:`M` and :math:`N` are given by the coarse resolution ocean
model, information relating to :math:`E` is provided by a parameterised eddy
energy budget. The recipe for GEOMETRIC then is as follows:

1. time-step the parameterised eddy energy budget to get :math:`E` with info provided by the GCM

2. calculate the new :math:`\kappa_{\rm gm}`

3. use the existing GM routines with new :math:`\kappa_{\rm gm}` and time-step the GCM. Cycle as appropriate.
  
The current NEMO implementation considers an eddy energy field that varies in
longitude, latitude and time (and so :math:`\kappa_{\rm gm}` inherits this
spatio-temporal dependence), given by

.. math::
   \frac{\mathrm{d}}{\mathrm{d} t} \int E\; \mathrm{d}z +
   \nabla \cdot \left( (\tilde{\boldsymbol{u}} - c\boldsymbol{e}_1 ) \int E\; \mathrm{d}z \right) =
   \int \kappa_{\rm gm} \frac{M^4}{N^2}\; \mathrm{d}z -
   \lambda \int E\; \mathrm{d}z +
   \nu_E \nabla^2 \int E\; \mathrm{d}z,
   
(respecively, the time-evolution, advection, source, dissipation and diffusion
of eddy energy), with :math:`\kappa_{\rm gm}` calculated as

.. math::
   \kappa_{\rm gm} = \alpha \frac{\int E\; \mathrm{d}z}{\int \Gamma (M^2 / N)\; \mathrm{d}z} \Gamma(z).
   
The symbols are as follows:

+-------------------------------+-----------------------------------------+----------------------+
| symbol                        | definition                              | units                |
+===============================+=========================================+======================+
| :math:`\alpha`                | eddy efficiency parameter               | :math:`--`           |
|                               | non-dimensional, :math:`|\alpha| \leq 1`|                      |
+-------------------------------+-----------------------------------------+----------------------+
| :math:`E`                     | total eddy energy                       | :math:`m^2\ s^{-2}`  |
+-------------------------------+-----------------------------------------+----------------------+
| :math:`M, N`                  | mean horizontal and vertical buoyancy   | :math:`s^{-1}`       |
|                               | gradient                                |                      |
+-------------------------------+-----------------------------------------+----------------------+
| :math:`\tilde{\boldsymbol{u}}`| depth-mean flow                         | :math:`m^2\ s^{-1}`  |
+-------------------------------+-----------------------------------------+----------------------+
| :math:`c`                     | long Rossby phase speed of 1st          | :math:`m^2\ s^{-1}`  |
|                               | baroclinic mode                         |                      |
+-------------------------------+-----------------------------------------+----------------------+
| :math:`\kappa_{\rm gm}`       | Gent--McWilliams coefficient            | :math:`m^2\ s^{-1}`  |
+-------------------------------+-----------------------------------------+----------------------+
| :math:`\lambda`               | linear damping rate of eddy energy      | :math:`s^{-1}`       |
+-------------------------------+-----------------------------------------+----------------------+
| :math:`\nu_{E}`               | Laplacian diffusion of eddy energy      | :math:`m^2\ s^{-1}`  |
+-------------------------------+-----------------------------------------+----------------------+

Advection
---------

The advection of eddy energy is given in flux form and has a contribution from
the depth-mean flow as well as a contribution associated with the westward
propagation of eddies at the long Rossby phase speed (motivated by e.g.
:cite:`Chelton-et-al11` and :cite:`KlockerMarshall14`). The advection is by the
barotropic mean flow already computed in NEMO, with a first order upwind scheme.
The baroclinic Rossby wave speed is obtained by computing :strike:`the eigenvalue
associated with the first baroclinic mode (see e.g.` eq. 6.11.8 of
:cite:`Gill-GFD`:strike:`) and uses two subroutines (eke_rossby and eke_thomas)`
the WKB expression given in :cite:`Chelton-et-al98` (their equation 2.2):

.. math::
    c_n \approx \frac{1}{n\pi} \int^0_{-H} N(z)\; \mathrm{d}z
    
and the long-phase speed that the total eddy energy is to be advected at is
computed as (e.g. eq. 12.3.13 of :cite:`Gill-GFD`)

.. math::
    c_p \approx -\frac{\beta}{f_0}c_n^2 = -c_n^2 \frac{\cos\phi_0}{2\Omega R \sin^2 \phi_0}
    
In practice the expression diverges at the equator and the actual wave
contribution to eddy energy advection as implemented in GEOMETRIC is taken to be
the minimum of the magnitude of the Rossby long-wave phase speed above and the
tropical planetary wave phase speed (e.g. eq. 12.3.14 of :cite:`Gill-GFD`),
i.e.,

.. math::
    c = \min\left(|c_p|, \left|\frac{c_n}{2n + 1}\right|\right) = \min\left(|c_p|, \left|c_1/3\right|\right)

See :ref:`here <sec:nemo-adv>` for usage and implementation details.

.. note ::
  As of Feb 2019 the removal of the routines to solve the tri-diagonal
  eigenvalue problem means the ``nn_wav_cal`` variable in ``namelist_cfg`` has
  been removed.

Source
------

The source of mesoscale eddy energy here is only from the slumping of neutral
surfaces through the eddy induced velocity as parameterised by the GM scheme
(note that it is positive-definite). These are straight-forwardly computed as is
(rather than using the quasi-Stokes streamfunction) using the already limited
slopes computed in NEMO. See :ref:`here <sec:nemo-sou>` for implementation
details.

Dissipation
-----------

The damping of eddy energy is linearly damped and the coefficient is specified
in ``namelist_cfg`` as a time-scale in *days* (which is subsequently converted
to *per seconds* in ``ldf_eke_init``). There is an option to read in an
externally prepared NetCDF file ``geom_diss_2D.nc`` that varies in longitude and
latitude in anticipation of further investigation. See :ref:`here
<sec:nemo-dis>` for usage details, and **here** for a sample Python Notebook to
generate the file.

Diffusion
---------

The diffusion of eddy energy is through a Laplacian (cf.
:cite:`EdenGreatbatch08`), done through relevant copy and pasting of code that
are in other NEMO modules. The GEOMETRIC scheme is actually stable (most likely
because of the upwinding scheme). The diffusion may be switched off by setting
``rn_eke_lap = 0.`` in ``namelist_cfg`` which will bypass the relevant loop in
``ldf_eke``.

NEMO implementation
-------------------

For details relating to the NEMO implementation (e.g. variable names, numerical
treatment, :code:`namelist_cfg` file), please see the following pages.

.. toctree::
   
   code_details.rst
   assorted_files.rst

.. bibliography:: ../refs.bib
   :filter: docname in docnames
