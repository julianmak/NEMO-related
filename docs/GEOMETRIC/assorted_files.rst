.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Jul  4 10:59:03 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

files
=====

Below are sample files or file generators you may need for using GEOMETRIC in
NEMO.

``namelist_cfg``
----------------

::

  !----------------------------------------------------------------------------------
  &namtra_ldfeiv !   eddy induced velocity param.
  !----------------------------------------------------------------------------------
     ln_ldfeiv     =.true.   ! use eddy induced velocity parameterization
     ln_ldfeiv_dia =.false.  ! diagnose eiv stream function and velocities
     rn_aeiv_0     = 1000.   ! eddy induced velocity coefficient   [m2/s]
     nn_aei_ijk_t  =   32    ! space/time variation of the eiv coeficient
     !                                !   =-20 (=-30)    read in eddy_induced_velocity_2D.nc (..._3D.nc) file
     !                                !   =  0           constant 
     !                                !   = 10 F(k)      =ldf_c1d 
     !                                !   = 20 F(i,j)    =ldf_c2d 
     !                                !   = 21 F(i,j,t)  =Treguier et al. JPO 1997 formulation
     !                                !   = 30 F(i,j,k)  =ldf_c2d + ldf_c1d
     !                                !   = 32 F(i,j,t)  = GEOMETRIC parameterization        (=> fill namldf_eke)
     ln_eke_equ    =.true.   ! switch on update of GEOMETRIC eddy energy equation            (=> fill namldf_eke)
                             ! forced to be true if nn_aei_ijk_t = 32
  /
  !----------------------------------------------------------------------------------
  &namldf_eke !   GEOMETRIC param. (total EKE equation)                           (nn_aei_ijk_t = 32)
  !----------------------------------------------------------------------------------
     rn_ekedis      =  180.      ! dissipation time scale of EKE [days]
        nn_eke_dis  =  -20       ! dissipation option
        !                             !   =  0  constant in space
        !                             !   =-20  read in geom_diss_2D.nc file
     rn_geom        =  0.1       ! geometric parameterization master coefficient (>0 & <1)
     rn_eke_init    =  1.e-2     ! initial total EKE value
     rn_eke_min     =  1.e+0     ! background value of total EKE
     rn_ross_min    =  4.e+3     ! tapering of aeiv based on min Rossby radius [m]
     !                           !   set to zero to not taper it
     rn_eke_lap     =  2000.     ! Laplacian diffusion coefficient of EKE
     !                           !   this is in all options below, so set it to zero and nothing is done
     rn_aeiv_min    =  1.e+1     ! minimum bound of eiv coefficient
     rn_aeiv_max    =  1.5e+4    ! maximum bound of eiv coefficient
     rn_SFmin       =  1.0       ! minimum bound of Structure Function
     rn_SFmax       =  1.0       ! maximum bound of Structure Function
     nn_eke_opt     =  1         ! options for terms to include in EKE budget
     !                                !   =  0  PE->EKE conversion, dissipation only 
     !                                !   =  1  as 0 but with advection
     !                                !   =  2  as 1 but with additional KE->EKE conversion
     !                                !   for testing purposes:
     !                                !   = 88  only advection by depth-averaged flow
     !                                !   = 99  only Laplacian diffusion
     ln_adv_wav     =  .true.   ! include advection at long Rossby speed
        nn_wav_cal  =  60             ! number of time steps between eigenvalue calculation
  /


``field_def_nemo-opa.xml``
--------------------------

These following needs to be added into ``field_def_nemo-opa.xml`` if any of the
GEOMETRIC routines in ``ldfeke`` is used, so XIOS does not crash the runs. Call
these in ``file_def_nemo.xml`` as appropriate (see a sample below).

(Note: you may or may not find the ``bn2`` variable (the vertical buoyancy
frequency diagnostic) in the T grid for instead of the W grid. The file below
has ``bn2`` moved to the W grid group.)

.. code-block:: xml

  <!-- T grid -->
      
  <field_group id="grid_T" grid_ref="grid_T_2D" >
  
     <!-- GEOMETRIC fields (requires nn_aei_ijk_t = 32)  -->
     <field id="eke"               long_name="total EKE (EKE+EPE)"                     unit="m3/s2" />
     <field id="trd_eke_adv_ubt"   long_name="ubt advective trend of EKE (LHS)"        unit="m3/s3" />
     <field id="trd_eke_adv_wav"   long_name="wav advective trend of EKE (LHS)"        unit="m3/s3" />
     <field id="trd_eke_lap"       long_name="diffusive trend of EKE (RHS)"            unit="m3/s3" />
     <field id="trd_eke_peS"       long_name="PE to EKE source trend (RHS)"            unit="m3/s3" />
     <field id="trd_eke_keS"       long_name="KE to EKE source trend (RHS)"            unit="m3/s3" />
     <field id="trd_eke_dis"       long_name="dissipation trend of EKE (RHS)"          unit="m3/s3" />

  </field_group>
  
  <!-- W grid -->
      
  <field_group id="grid_W" grid_ref="grid_W_3D">
     <!-- GEOMETRIC fields (requires nn_aei_ijk_t = 32)  -->
     <field id="aeiv_geom"      long_name="3D w-EIV coefficient from GEOMETRIC param." unit="m2/s" />
     <field id="rossby_rad"     long_name="internal Rossby defromation radius"         unit="m"    grid_ref="grid_W_2D"/>
     <field id="bn2"            long_name="squared Brunt-Vaisala frequency"            unit="s-1"  />
     <field id="c1_vert"        long_name="1st baroclinic mode phase speed"            unit="m/s"  grid_ref="grid_W_2D"/>
     <field id="c_ros"          long_name="long Rossby phase speed"                    unit="m/s"  grid_ref="grid_W_2D"/>

  </field_group>


``file_def_nemo.xml``
---------------------


