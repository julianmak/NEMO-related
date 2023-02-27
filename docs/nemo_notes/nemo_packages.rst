.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   
.. _sec:nemo_packages:

Other NEMO packages
===================

Some useful tools that come with NEMO are available in the analog of the
``TOOLS`` folder. These are built using the ``ARCH`` files as you would for
building an experiment with e.g.

.. code-block:: bash
  
  ./maketools -n REBUILD_NEMO -m HKUST_HPC2

Notes on the ones I have used may be found here.

REBUILD_NEMO
------------

XIOS can combine the output-per-CPU cells into one global file, but by default
the restart files and ``mesh_mask.nc`` files are output per CPU, so it is useful
to recombine them. This can be done through the ``REBUILD_NEMO`` package.

Build as usual, and the resulting output should be a ``rebuild_nemo.exe`` in the
folder, to be driven by the script ``rebuild_nemo``. The way to use it is to
call the script as, for example,

.. code-block:: bash

  $BASE_DIR/tools/REBUILD_NEMO/rebuild_nemo ${MODEL}_${RES_TIMESTAMP}_restart $NUM_CPU
  
where ``$BASE_DIR`` is wherever the folder lives, the things to be combined look
like ``${MODEL}_${RES_TMESTAMP}_restart_0000.nc`` (e.g. ``mesh_mask_00??.nc``,
``UNAGI_00051840_restart[_ice].nc``, etc.), and ``$NUM_CPU`` are the number of
files to combine (e.g. if we use 96 cores then we get ``mesh_mask_0000.nc`` to
``mesh_mask_0095.nc``, and we should do ``export NUM_CPU=96``).

SECTIONS_DIADCT
---------------

WEIGHTS
-------

DOMAINcfg
---------

This package generates the ``domain_cfg.nc`` file that encodes the grid
locations and spacings, and is recommended for creating new configurations,
mostly because the vertical grid spacings with partial steps correction are a
bit weird to try and do manually. From the ``readme`` file in there, **you seem
to need to use xios1 to compile this**; see the :ref:`README <sec:nemo36:>`_ for
how various things to watch out for.

When compiled the executable is called ``make_domain_cfg.exe``, and it expects
to read a ``bathy_meter.nc`` (links are ok) and a ``namelist_cfg`` file. The
``namelist_cfg`` file should contain the various settings for horizontal and
vertical grid spacing, which should be consistent with the content in
``bathy_meter.nc``. An example ``namelist_cfg`` is the following:

::

  !-----------------------------------------------------------------------
  &namcfg        !   parameters of the configuration
  !-----------------------------------------------------------------------
     !
     ln_e3_dep   = .true.    ! =T : e3=dk[depth] in discret sens. 
     !                       !      ===>>> will become the only possibility in v4.0
     !                       ! =F : e3 analytical derivative of depth function
     !                       !      only there for backward compatibility test with v3.6
     !                       !      
     cp_cfg      =     "UNAGI"           !  name of the configuration
     jp_cfg      =      100              !  resolution of the configuration
     jpidta      =       90              !  1st lateral dimension ( >= jpi )
     jpjdta      =       26              !  2nd    "         "    ( >= jpj )
     jpkdta      =       31              !  number of levels      ( >= jpk )
     jpiglo      =       90              !  1st dimension of global domain --> i =jpidta
     jpjglo      =       26              !  2nd    -                  -    --> j  =jpjdta
     jpizoom     =       1               !  left bottom (i,j) indices of the zoom
     jpjzoom     =       1               !  in data domain indices
     jperio      =       1               !  lateral cond. type (between 0 and 6) [1 is EW periodicity]
  /
  !-----------------------------------------------------------------------
  &namzgr        !   vertical coordinate
  !-----------------------------------------------------------------------
     ln_zps      = .true.    !  z-coordinate - partial steps
     ln_linssh   = .true.    !  linear free surface
  /
  !-----------------------------------------------------------------------
  &namdom        !  
  !-----------------------------------------------------------------------
     jphgr_msh   =       3               !  type of horizontal mesh
     ppglam0     =       0.0             !  longitude of first raw and column T-point (jphgr_msh = 1)
     ppgphi0     =     -50.0             !  latitude  of first raw and column T-point (jphgr_msh = 1)
     ppe1_deg    =  999999.0             !  zonal      grid-spacing (degrees)
     ppe2_deg    =  999999.0             !  meridional grid-spacing (degrees)
     ppe1_m      =  100000.0             !  zonal      grid-spacing (metres)
     ppe2_m      =  100000.0             !  meridional grid-spacing (metres)
     ppsur       =  999999.0             !  ORCA r4, r2 and r05 coefficients
     ppa0        =  999999.0             ! (default coefficients)
     ppa1        =  999999.0             !
     ppkth       =      18.0             !
     ppacr       =      10.0             !
     ppdzmin     =    10.0               !  Minimum vertical spacing
     pphmax      =    3000.0             !  Maximum depth
     ldbletanh   =  .FALSE.              !  Use/do not use double tanf function for vertical coordinates
     ppa2        =  999999.0             !  Double tanh function parameters
     ppkth2      =  999999.0             !
     ppacr2      =  999999.0             !
  /
  
Here, the configuration is called ``UNAGI``. The ``jp[ijk]data`` is the number
of grid cells in :math:`(x,y,z)`, and I chose ``jp[ij]glo`` to be consistent
with the choice of horizontal sizes. The ``jperio`` denotes the periodicities
(see ``src/domcfg.f90`` for the choices). The present model uses a Cartesian
grid on a :math:`\beta`-plane corresponding to ``jphgr_msh = 3`` (see
``src/domhgr.f90`` for choices), and is centred at longitude 0 and latitude 50 S
(see ``ppglam0`` and ``ppgphi0``). The grid spacing here is 100 km,
correpsonding to ``ppe[12]_m``; the values of ``999999.0`` are options that are
not used.

For the vertical grid, ``ln_zps`` switches on the partial step correction and
takes into account ``bathy_meter.nc``. The vertical spacing is governed through
the parameters ``ppkth``, ``ppacr``, ``ppdzmin`` and ``pphmax``
(:cite:`MadecImbard96`; unless you use the double tanh option).

.. note ::

  Note NEMO 4.2 seems to be using different namings and convention (`see here
  <https://sites.nemo-ocean.io/user-guide/migration.html`_), so be careful with
  above.

NESTING (AGRIF)
---------------
