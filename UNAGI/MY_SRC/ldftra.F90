MODULE ldftra
   !!======================================================================
   !!                       ***  MODULE  ldftra  ***
   !! Ocean physics:  lateral diffusivity coefficients 
   !!=====================================================================
   !! History :       ! 1997-07  (G. Madec)  from inimix.F split in 2 routines
   !!   NEMO     1.0  ! 2002-09  (G. Madec)  F90: Free form and module
   !!            2.0  ! 2005-11  (G. Madec)  
   !!            3.7  ! 2013-12  (F. Lemarie, G. Madec)  restructuration/simplification of aht/aeiv specification,
   !!                 !                                  add velocity dependent coefficient and optional read in file
   !!            4.0  ! 2017-11  (J. Mak, G. Madec) added GEOMETRIC parameterization
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   ldf_tra_init : initialization, namelist read, and parameters control
   !!   ldf_tra      : update lateral eddy diffusivity coefficients at each time step 
   !!   ldf_eiv_init : initialization of the eiv coeff. from namelist choices 
   !!   ldf_eiv      : time evolution of the eiv coefficients (function of the growth rate of baroclinic instability)
   !!   ldf_eiv_trp  : add to the input ocean transport the contribution of the EIV parametrization
   !!   ldf_eiv_dia  : diagnose the eddy induced velocity from the eiv streamfunction
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE ldfslp          ! lateral diffusion: slope of iso-neutral surfaces
   USE ldfc1d_c2d      ! lateral diffusion: 1D & 2D cases 
   USE diaptr
   !
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O module for ehanced bottom friction file
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE wrk_nemo        ! work arrays
   USE timing          ! timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ldf_tra_init   ! called by nemogcm.F90
   PUBLIC   ldf_tra        ! called by step.F90
   PUBLIC   ldf_eiv_init   ! called by nemogcm.F90
   PUBLIC   ldf_eiv        ! called by step.F90
   PUBLIC   ldf_eiv_trp    ! called by traadv.F90
   PUBLIC   ldf_eiv_dia    ! called by traldf_iso and traldf_iso_triad.F90
   
   !                                   !!* Namelist namtra_ldf : lateral mixing on tracers * 
   !                                    != Operator type =!
   LOGICAL , PUBLIC ::   ln_traldf_lap       !: laplacian operator
   LOGICAL , PUBLIC ::   ln_traldf_blp       !: bilaplacian operator
   !                                    != Direction of action =!
   LOGICAL , PUBLIC ::   ln_traldf_lev       !: iso-level direction
   LOGICAL , PUBLIC ::   ln_traldf_hor       !: horizontal (geopotential) direction
!  LOGICAL , PUBLIC ::   ln_traldf_iso       !: iso-neutral direction                    (see ldfslp)
!  LOGICAL , PUBLIC ::   ln_traldf_triad     !: griffies triad scheme                    (see ldfslp)
   LOGICAL , PUBLIC ::   ln_traldf_msc       !: Method of Stabilizing Correction 
!  LOGICAL , PUBLIC ::   ln_triad_iso        !: pure horizontal mixing in ML             (see ldfslp)
!  LOGICAL , PUBLIC ::   ln_botmix_triad     !: mixing on bottom                         (see ldfslp)
!  REAL(wp), PUBLIC ::   rn_sw_triad         !: =1/0 switching triad / all 4 triads used (see ldfslp)
!  REAL(wp), PUBLIC ::   rn_slpmax           !: slope limit                              (see ldfslp)
   !                                    !=  Coefficients =!
   INTEGER , PUBLIC ::   nn_aht_ijk_t        !: choice of time & space variations of the lateral eddy diffusivity coef.
   REAL(wp), PUBLIC ::   rn_aht_0            !:   laplacian lateral eddy diffusivity [m2/s]
   REAL(wp), PUBLIC ::   rn_bht_0            !: bilaplacian lateral eddy diffusivity [m4/s]

   !                                   !!* Namelist namtra_ldfeiv : eddy induced velocity param. *
   !                                    != Use/diagnose eiv =!
   LOGICAL , PUBLIC ::   ln_ldfeiv           !: eddy induced velocity flag
   LOGICAL , PUBLIC ::   ln_ldfeiv_dia       !: diagnose & output eiv streamfunction and velocity (IOM)
   !                                    != Coefficients =!
   INTEGER , PUBLIC ::   nn_aei_ijk_t        !: choice of time/space variation of the eiv coeff.
   REAL(wp), PUBLIC ::   rn_aeiv_0           !: eddy induced velocity coefficient [m2/s]
   
   LOGICAL , PUBLIC ::   l_ldftra_time = .FALSE.   !: flag for time variation of the lateral eddy diffusivity coef.
   LOGICAL , PUBLIC ::   l_ldfeiv_time = .FALSE.   !: flag for time variation of the eiv coef.
   
   LOGICAL , PUBLIC ::   ln_eke_equ                !: flag for having updates to eddy energy equation
   LOGICAL , PUBLIC ::   l_ldfeke      = .FALSE.   !: GEOMETRIC - total EKE flag
   LOGICAL , PUBLIC ::   l_eke_eiv     = .FALSE.   !: GEOMETRIC - aeiw flag

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ahtu, ahtv   !: eddy diffusivity coef. at U- and V-points   [m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   aeiu, aeiv   !: eddy induced velocity coeff.                [m2/s]

   REAL(wp) ::   r1_4  = 0.25_wp          ! =1/4
   REAL(wp) ::   r1_12 = 1._wp / 12._wp   ! =1/12

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.7 , NEMO Consortium (2015)
   !! $Id: ldftra.F90 7753 2017-03-03 11:46:59Z mocavero $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ldf_tra_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_tra_init  ***
      !! 
      !! ** Purpose :   initializations of the tracer lateral mixing coeff.
      !!
      !! ** Method  : * the eddy diffusivity coef. specification depends on:
      !!
      !!    ln_traldf_lap = T     laplacian operator
      !!    ln_traldf_blp = T   bilaplacian operator
      !!
      !!    nn_aht_ijk_t  =  0 => = constant
      !!                  !
      !!                  = 10 => = F(z) : constant with a reduction of 1/4 with depth 
      !!                  !
      !!                  =-20 => = F(i,j)   = shape read in 'eddy_diffusivity.nc' file
      !!                  = 20    = F(i,j)   = F(e1,e2) or F(e1^3,e2^3) (lap or bilap case)
      !!                  = 21    = F(i,j,t) = F(growth rate of baroclinic instability)
      !!                  !
      !!                  =-30 => = F(i,j,k)   = shape read in 'eddy_diffusivity.nc' file
      !!                  = 30    = F(i,j,k)   = 2D (case 20) + decrease with depth (case 10)
      !!                  = 31    = F(i,j,k,t) = F(local velocity) (  |u|e  /12   laplacian operator
      !!                                                          or |u|e^3/12 bilaplacian operator )
      !!              * initialisation of the eddy induced velocity coefficient by a call to ldf_eiv_init 
      !!            
      !! ** action  : ahtu, ahtv initialized once for all or l_ldftra_time set to true
      !!              aeiu, aeiv initialized once for all or l_ldfeiv_time set to true
      !!----------------------------------------------------------------------
      INTEGER  ::   jk                ! dummy loop indices
      INTEGER  ::   ierr, inum, ios   ! local integer
      REAL(wp) ::   zah0              ! local scalar
      !
      NAMELIST/namtra_ldf/ ln_traldf_lap, ln_traldf_blp  ,                   &   ! type of operator
         &                 ln_traldf_lev, ln_traldf_hor  , ln_traldf_triad,  &   ! acting direction of the operator
         &                 ln_traldf_iso, ln_traldf_msc  , rn_slpmax      ,  &   ! option for iso-neutral operator
         &                 ln_triad_iso , ln_botmix_triad, rn_sw_triad    ,  &   ! option for triad operator
         &                 rn_aht_0     , rn_bht_0       , nn_aht_ijk_t          ! lateral eddy coefficient
      !!----------------------------------------------------------------------
      !
      !  Choice of lateral tracer physics
      ! =================================
      !
      REWIND( numnam_ref )              ! Namelist namtra_ldf in reference namelist : Lateral physics on tracers
      READ  ( numnam_ref, namtra_ldf, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtra_ldf in reference namelist', lwp )
      !
      REWIND( numnam_cfg )              ! Namelist namtra_ldf in configuration namelist : Lateral physics on tracers
      READ  ( numnam_cfg, namtra_ldf, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtra_ldf in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namtra_ldf )
      !
      IF(lwp) THEN                      ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ldf_tra_init : lateral tracer physics'
         WRITE(numout,*) '~~~~~~~~~~~~ '
         WRITE(numout,*) '   Namelist namtra_ldf : lateral mixing parameters (type, direction, coefficients)'
         !
         WRITE(numout,*) '      type :'
         WRITE(numout,*) '         laplacian operator                      ln_traldf_lap   = ', ln_traldf_lap
         WRITE(numout,*) '         bilaplacian operator                    ln_traldf_blp   = ', ln_traldf_blp
         !
         WRITE(numout,*) '      direction of action :'
         WRITE(numout,*) '         iso-level                               ln_traldf_lev   = ', ln_traldf_lev
         WRITE(numout,*) '         horizontal (geopotential)               ln_traldf_hor   = ', ln_traldf_hor
         WRITE(numout,*) '         iso-neutral Madec operator              ln_traldf_iso   = ', ln_traldf_iso
         WRITE(numout,*) '         iso-neutral triad operator              ln_traldf_triad = ', ln_traldf_triad
         WRITE(numout,*) '            iso-neutral (Method of Stab. Corr.)  ln_traldf_msc   = ', ln_traldf_msc
         WRITE(numout,*) '            maximum isoppycnal slope             rn_slpmax       = ', rn_slpmax
         WRITE(numout,*) '            pure lateral mixing in ML            ln_triad_iso    = ', ln_triad_iso
         WRITE(numout,*) '            switching triad or not               rn_sw_triad     = ', rn_sw_triad
         WRITE(numout,*) '            lateral mixing on bottom             ln_botmix_triad = ', ln_botmix_triad
         !
         WRITE(numout,*) '      coefficients :'
         WRITE(numout,*) '         lateral eddy diffusivity   (lap case)   rn_aht_0        = ', rn_aht_0
         WRITE(numout,*) '         lateral eddy diffusivity (bilap case)   rn_bht_0        = ', rn_bht_0
         WRITE(numout,*) '         type of time-space variation            nn_aht_ijk_t    = ', nn_aht_ijk_t
      ENDIF
      !
      !                                ! Parameter control
      !
      IF( .NOT.ln_traldf_lap .AND. .NOT.ln_traldf_blp ) THEN
         IF(lwp) WRITE(numout,*) '   No diffusive operator selected. ahtu and ahtv are not allocated'
         l_ldftra_time = .FALSE.
         RETURN
      ENDIF
      !
      IF( ln_traldf_blp .AND. ( ln_traldf_iso .OR. ln_traldf_triad) ) THEN     ! iso-neutral bilaplacian need MSC
         IF( .NOT.ln_traldf_msc )   CALL ctl_stop( 'tra_ldf_init: iso-neutral bilaplacian requires ln_traldf_msc=.true.' )
      ENDIF
      !
      !  Space/time variation of eddy coefficients 
      ! ===========================================
      !                                               ! allocate the aht arrays
      ALLOCATE( ahtu(jpi,jpj,jpk) , ahtv(jpi,jpj,jpk) , STAT=ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'ldf_tra_init: failed to allocate arrays')
      !
      ahtu(:,:,jpk) = 0._wp                           ! last level always 0  
      ahtv(:,:,jpk) = 0._wp
      !
      !                                               ! value of eddy mixing coef.
      IF    ( ln_traldf_lap ) THEN   ;   zah0 =      rn_aht_0        !   laplacian operator
      ELSEIF( ln_traldf_blp ) THEN   ;   zah0 = ABS( rn_bht_0 )      ! bilaplacian operator
      ENDIF
      !
      l_ldftra_time = .FALSE.                         ! no time variation except in case defined below
      !
      IF( ln_traldf_lap .OR. ln_traldf_blp ) THEN     ! only if a lateral diffusion operator is used
         !
         SELECT CASE(  nn_aht_ijk_t  )                   ! Specification of space time variations of ehtu, ahtv
         !
         CASE(   0  )      !==  constant  ==!
            IF(lwp) WRITE(numout,*) '          tracer mixing coef. = constant = ', zah0 ! jm: moved from rn_aht_0 to zah0
            ahtu(:,:,:) = zah0 * umask(:,:,:)
            ahtv(:,:,:) = zah0 * vmask(:,:,:)
            !
         CASE(  10  )      !==  fixed profile  ==!
            IF(lwp) WRITE(numout,*) '          tracer mixing coef. = F( depth )'
            ahtu(:,:,1) = zah0 * umask(:,:,1)                      ! constant surface value
            ahtv(:,:,1) = zah0 * vmask(:,:,1)
            CALL ldf_c1d( 'TRA', r1_4, ahtu(:,:,1), ahtv(:,:,1), ahtu, ahtv )
            !
         CASE ( -20 )      !== fixed horizontal shape read in file  ==!
            IF(lwp) WRITE(numout,*) '          tracer mixing coef. = F(i,j) read in eddy_diffusivity.nc file'
            CALL iom_open( 'eddy_diffusivity_2D.nc', inum )
            CALL iom_get ( inum, jpdom_data, 'ahtu_2D', ahtu(:,:,1) )
            CALL iom_get ( inum, jpdom_data, 'ahtv_2D', ahtv(:,:,1) )
            CALL iom_close( inum )
            DO jk = 2, jpkm1
               ahtu(:,:,jk) = ahtu(:,:,1) * umask(:,:,jk)
               ahtv(:,:,jk) = ahtv(:,:,1) * vmask(:,:,jk)
            END DO
            !
         CASE(  20  )      !== fixed horizontal shape  ==!
            IF(lwp) WRITE(numout,*) '          tracer mixing coef. = F( e1, e2 ) or F( e1^3, e2^3 ) (lap or blp case)'
            IF( ln_traldf_lap )   CALL ldf_c2d( 'TRA', 'LAP', zah0, ahtu, ahtv )    ! surface value proportional to scale factor
            IF( ln_traldf_blp )   CALL ldf_c2d( 'TRA', 'BLP', zah0, ahtu, ahtv )    ! surface value proportional to scale factor
            !
         CASE(  21  )      !==  time varying 2D field  ==!
            IF(lwp) WRITE(numout,*) '          tracer mixing coef. = F( latitude, longitude, time )'
            IF(lwp) WRITE(numout,*) '                              = F( growth rate of baroclinic instability )'
            IF(lwp) WRITE(numout,*) '                              min value = 0.1 * rn_aht_0'
            IF(lwp) WRITE(numout,*) '                              max value = rn_aht_0 (rn_aeiv_0 if nn_aei_ijk_t=21)'
            IF(lwp) WRITE(numout,*) '                              increased to rn_aht_0 within 20N-20S'
            !
            l_ldftra_time = .TRUE.     ! will be calculated by call to ldf_tra routine in step.F90
            !
            IF( ln_traldf_blp ) THEN
               CALL ctl_stop( 'ldf_tra_init: aht=F(growth rate of baroc. insta.) incompatible with bilaplacian operator' )
            ENDIF
            !
         CASE( -30  )      !== fixed 3D shape read in file  ==!
            IF(lwp) WRITE(numout,*) '          tracer mixing coef. = F(i,j,k) read in eddy_diffusivity.nc file'
            CALL iom_open( 'eddy_diffusivity_3D.nc', inum )
            CALL iom_get ( inum, jpdom_data, 'ahtu_3D', ahtu )
            CALL iom_get ( inum, jpdom_data, 'ahtv_3D', ahtv )
            CALL iom_close( inum )
            DO jk = 1, jpkm1
               ahtu(:,:,jk) = ahtu(:,:,jk) * umask(:,:,jk)
               ahtv(:,:,jk) = ahtv(:,:,jk) * vmask(:,:,jk)
            END DO
            !
         CASE(  30  )      !==  fixed 3D shape  ==!
            IF(lwp) WRITE(numout,*) '          tracer mixing coef. = F( latitude, longitude, depth )'
            IF( ln_traldf_lap )   CALL ldf_c2d( 'TRA', 'LAP', zah0, ahtu, ahtv )    ! surface value proportional to scale factor
            IF( ln_traldf_blp )   CALL ldf_c2d( 'TRA', 'BLP', zah0, ahtu, ahtv )    ! surface value proportional to scale factor
            !                                                    ! reduction with depth
            CALL ldf_c1d( 'TRA', r1_4, ahtu(:,:,1), ahtv(:,:,1), ahtu, ahtv )
            !
         CASE(  31  )      !==  time varying 3D field  ==!
            IF(lwp) WRITE(numout,*) '          tracer mixing coef. = F( latitude, longitude, depth , time )'
            IF(lwp) WRITE(numout,*) '                                proportional to the velocity : |u|e/12 or |u|e^3/12'
            !
            l_ldftra_time = .TRUE.     ! will be calculated by call to ldf_tra routine in step.F90
            !
         CASE DEFAULT
            CALL ctl_stop('ldf_tra_init: wrong choice for nn_aht_ijk_t, the type of space-time variation of aht')
         END SELECT
         !
         IF( ln_traldf_blp .AND. .NOT. l_ldftra_time ) THEN
            ahtu(:,:,:) = SQRT( ahtu(:,:,:) )
            ahtv(:,:,:) = SQRT( ahtv(:,:,:) )
         ENDIF
         !
      ENDIF
      !
   END SUBROUTINE ldf_tra_init


   SUBROUTINE ldf_tra( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_tra  ***
      !! 
      !! ** Purpose :   update at kt the tracer lateral mixing coeff. (aht and aeiv)
      !!
      !! ** Method  :   time varying eddy diffusivity coefficients:
      !!
      !!    nn_aei_ijk_t = 21    aeiu, aeiv = F(i,j,  t) = F(growth rate of baroclinic instability)
      !!                                                   with a reduction to 0 in vicinity of the Equator
      !!    nn_aht_ijk_t = 21    ahtu, ahtv = F(i,j,  t) = F(growth rate of baroclinic instability)
      !!
      !!                 = 31    ahtu, ahtv = F(i,j,k,t) = F(local velocity) (  |u|e  /12   laplacian operator
      !!                                                                     or |u|e^3/12 bilaplacian operator )
      !!
      !! ** action  :   ahtu, ahtv   update at each time step   
      !!                aeiu, aeiv      -       -     -    -   (if ln_ldfeiv=T) 
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! time step
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zaht, zahf, zaht_min, z1_f20       ! local scalar
      !!----------------------------------------------------------------------
      !
      IF( ln_ldfeiv .AND. nn_aei_ijk_t == 21 ) THEN       ! eddy induced velocity coefficients
         !                                ! =F(growth rate of baroclinic instability)
         !                                ! max value rn_aeiv_0 ; decreased to 0 within 20N-20S
         CALL ldf_eiv( kt, rn_aeiv_0, aeiu, aeiv )
      ENDIF
      !
      SELECT CASE(  nn_aht_ijk_t  )       ! Eddy diffusivity coefficients
      !
      CASE(  21  )       !==  time varying 2D field  ==!   = F( growth rate of baroclinic instability )
         !                                             !   min value rn_aht_0 / 10 
         !                                             !   max value rn_aht_0 (rn_aeiv_0 if nn_aei_ijk_t=21)
         !                                             !   increase to rn_aht_0 within 20N-20S
         IF( ln_ldfeiv .AND. nn_aei_ijk_t == 21 ) THEN   ! use the already computed aei.
            ahtu(:,:,1) = aeiu(:,:,1)
            ahtv(:,:,1) = aeiv(:,:,1)
         ELSE                                            ! compute aht. 
            CALL ldf_eiv( kt, rn_aht_0, ahtu, ahtv )
         ENDIF
         !
         z1_f20   = 1._wp / (  2._wp * omega * SIN( rad * 20._wp )  )      ! 1 / ff(20 degrees)   
         zaht_min = 0.2_wp * rn_aht_0                                      ! minimum value for aht
         DO jj = 1, jpj
            DO ji = 1, jpi
               !!gm CAUTION : here we assume lat/lon grid in 20deg N/S band (like all ORCA cfg)
               !!     ==>>>   The Coriolis value is identical for t- & u_points, and for v- and f-points
               zaht = ( 1._wp -  MIN( 1._wp , ABS( ff_t(ji,jj) * z1_f20 ) ) ) * ( rn_aht_0 - zaht_min )
               zahf = ( 1._wp -  MIN( 1._wp , ABS( ff_f(ji,jj) * z1_f20 ) ) ) * ( rn_aht_0 - zaht_min )
               ahtu(ji,jj,1) = (  MAX( zaht_min, ahtu(ji,jj,1) ) + zaht  ) * umask(ji,jj,1)     ! min value zaht_min
               ahtv(ji,jj,1) = (  MAX( zaht_min, ahtv(ji,jj,1) ) + zahf  ) * vmask(ji,jj,1)     ! increase within 20S-20N
            END DO
         END DO
         DO jk = 2, jpkm1                             ! deeper value = surface value
            ahtu(:,:,jk) = ahtu(:,:,1) * umask(:,:,jk)
            ahtv(:,:,jk) = ahtv(:,:,1) * vmask(:,:,jk)
         END DO
         !
      CASE(  31  )       !==  time varying 3D field  ==!   = F( local velocity )
         IF( ln_traldf_lap     ) THEN          !   laplacian operator |u| e /12
            DO jk = 1, jpkm1
               ahtu(:,:,jk) = ABS( ub(:,:,jk) ) * e1u(:,:) * r1_12
               ahtv(:,:,jk) = ABS( vb(:,:,jk) ) * e2v(:,:) * r1_12
            END DO
         ELSEIF( ln_traldf_blp ) THEN      ! bilaplacian operator      sqrt( |u| e^3 /12 ) = sqrt( |u| e /12 ) * e
            DO jk = 1, jpkm1
               ahtu(:,:,jk) = SQRT(  ABS( ub(:,:,jk) ) * e1u(:,:) * r1_12  ) * e1u(:,:)
               ahtv(:,:,jk) = SQRT(  ABS( vb(:,:,jk) ) * e2v(:,:) * r1_12  ) * e2v(:,:)
            END DO
         ENDIF
         !
      END SELECT
      !
      CALL iom_put( "ahtu_2d", ahtu(:,:,1) )   ! surface u-eddy diffusivity coeff.
      CALL iom_put( "ahtv_2d", ahtv(:,:,1) )   ! surface v-eddy diffusivity coeff.
      CALL iom_put( "ahtu_3d", ahtu(:,:,:) )   ! 3D      u-eddy diffusivity coeff.
      CALL iom_put( "ahtv_3d", ahtv(:,:,:) )   ! 3D      v-eddy diffusivity coeff.
      !
!!gm  : THE IF below is to be checked (comes from Seb)
!!jm  : added a .NOT. flag here to avoid duplication of iom_put
      IF( ln_ldfeiv .AND. (.NOT. ln_eke_equ) ) THEN
        CALL iom_put( "aeiu_2d", aeiu(:,:,1) )   ! surface u-EIV coeff.
        CALL iom_put( "aeiv_2d", aeiv(:,:,1) )   ! surface v-EIV coeff.
        CALL iom_put( "aeiu_3d", aeiu(:,:,:) )   ! 3D      u-EIV coeff.
        CALL iom_put( "aeiv_3d", aeiv(:,:,:) )   ! 3D      v-EIV coeff.
      ENDIF
      !
   END SUBROUTINE ldf_tra


   SUBROUTINE ldf_eiv_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_eiv_init  ***
      !!
      !! ** Purpose :   initialization of the eiv coeff. from namelist choices.
      !!
      !! ** Method :
      !!
      !! ** Action :   aeiu , aeiv   : EIV coeff. at u- & v-points
      !!               l_ldfeiv_time : =T if EIV coefficients vary with time
      !!----------------------------------------------------------------------
      INTEGER  ::   jk                ! dummy loop indices
      INTEGER  ::   ierr, inum, ios   ! local integer
      !
      NAMELIST/namtra_ldfeiv/ ln_ldfeiv   , ln_ldfeiv_dia,   &    ! eddy induced velocity (eiv)
         &                    nn_aei_ijk_t, rn_aeiv_0    ,   &    ! eiv  coefficient
         &                    ln_eke_equ                          ! GEOMETRIC eddy energy equation
      !!----------------------------------------------------------------------
      !
      REWIND( numnam_ref )              ! Namelist namtra_ldfeiv in reference namelist : eddy induced velocity param.
      READ  ( numnam_ref, namtra_ldfeiv, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtra_ldfeiv in reference namelist', lwp )
      !
      REWIND( numnam_cfg )              ! Namelist namtra_ldfeiv in configuration namelist : eddy induced velocity param.
      READ  ( numnam_cfg, namtra_ldfeiv, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtra_ldfeiv in configuration namelist', lwp )
      IF(lwm)  WRITE ( numond, namtra_ldfeiv )

      IF(lwp) THEN                      ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ldf_eiv_init : eddy induced velocity parametrization'
         WRITE(numout,*) '~~~~~~~~~~~~ '
         WRITE(numout,*) '   Namelist namtra_ldfeiv : '
         WRITE(numout,*) '      Eddy Induced Velocity (eiv) param.      ln_ldfeiv     = ', ln_ldfeiv
         WRITE(numout,*) '      eiv streamfunction & velocity diag.     ln_ldfeiv_dia = ', ln_ldfeiv_dia
         WRITE(numout,*) '      eddy induced velocity coef.             rn_aeiv_0     = ', rn_aeiv_0
         WRITE(numout,*) '      type of time-space variation            nn_aei_ijk_t  = ', nn_aei_ijk_t
         WRITE(numout,*)
      ENDIF
      !
!!! jm 11 Apr 19: switched off the force stop to have option of GM + no Redi + bilaplacian tracer diffusion
!      IF( ln_ldfeiv .AND. ln_traldf_blp )   CALL ctl_stop( 'ldf_eiv_init: eddy induced velocity ONLY with laplacian diffusivity' )
      IF( ln_ldfeiv .AND. ln_traldf_blp ) THEN
         IF(lwp) WRITE(numout,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
         IF(lwp) WRITE(numout,*) '       jm hack: ln_ldfeiv and ln_traldf_blp both on            '
         IF( ln_traldf_iso )   CALL ctl_stop('STOP', 'ldf_eiv: (jm hack) ln_traldf_iso also on, turn one of ln_traldf_blp or ln_traldf_iso off')
      ENDIF
!!!
      !                                 ! Parameter control
      l_ldfeiv_time = .FALSE.    
      !
      IF( ln_ldfeiv ) THEN                         ! allocate the aei arrays
         ALLOCATE( aeiu(jpi,jpj,jpk), aeiv(jpi,jpj,jpk), STAT=ierr )
         IF( ierr /= 0 )   CALL ctl_stop('STOP', 'ldf_eiv: failed to allocate arrays')
         !
         SELECT CASE( nn_aei_ijk_t )               ! Specification of space time variations of eaiu, aeiv
         !
         CASE(   0  )      !==  constant  ==!
            IF(lwp) WRITE(numout,*) '          eddy induced velocity coef. = constant = ', rn_aeiv_0
            aeiu(:,:,:) = rn_aeiv_0
            aeiv(:,:,:) = rn_aeiv_0
            !
         CASE(  10  )      !==  fixed profile  ==!
            IF(lwp) WRITE(numout,*) '          eddy induced velocity coef. = F( depth )'
            aeiu(:,:,1) = rn_aeiv_0                                ! constant surface value
            aeiv(:,:,1) = rn_aeiv_0
            CALL ldf_c1d( 'TRA', r1_4, aeiu(:,:,1), aeiv(:,:,1), aeiu, aeiv )
            !
         CASE ( -20 )      !== fixed horizontal shape read in file  ==!
            IF(lwp) WRITE(numout,*) '          eddy induced velocity coef. = F(i,j) read in eddy_diffusivity_2D.nc file'
            CALL iom_open ( 'eddy_induced_velocity_2D.nc', inum )
            CALL iom_get  ( inum, jpdom_data, 'aeiu', aeiu(:,:,1) )
            CALL iom_get  ( inum, jpdom_data, 'aeiv', aeiv(:,:,1) )
            CALL iom_close( inum )
            DO jk = 2, jpk
               aeiu(:,:,jk) = aeiu(:,:,1)
               aeiv(:,:,jk) = aeiv(:,:,1)
            END DO
            !
         CASE(  20  )      !== fixed horizontal shape  ==!
            IF(lwp) WRITE(numout,*) '          eddy induced velocity coef. = F( e1, e2 ) or F( e1^3, e2^3 ) (lap or bilap case)'
            CALL ldf_c2d( 'TRA', 'LAP', rn_aeiv_0, aeiu, aeiv )    ! surface value proportional to scale factor
            !
         CASE(  21  )       !==  time varying 2D field  ==!
            IF(lwp) WRITE(numout,*) '          eddy induced velocity coef. = F( latitude, longitude, time )'
            IF(lwp) WRITE(numout,*) '                              = F( growth rate of baroclinic instability )'
            !
            l_ldfeiv_time = .TRUE.     ! will be calculated by call to ldf_tra routine in step.F90
            !
         CASE( -30  )      !== fixed 3D shape read in file  ==!
            IF(lwp) WRITE(numout,*) '          eddy induced velocity coef. = F(i,j,k) read in eddy_diffusivity_3D.nc file'
            CALL iom_open ( 'eddy_induced_velocity_3D.nc', inum )
            CALL iom_get  ( inum, jpdom_data, 'aeiu', aeiu )
            CALL iom_get  ( inum, jpdom_data, 'aeiv', aeiv )
            CALL iom_close( inum )
            !
         CASE(  30  )       !==  fixed 3D shape  ==!
            IF(lwp) WRITE(numout,*) '          eddy induced velocity coef. = F( latitude, longitude, depth )'
            CALL ldf_c2d( 'TRA', 'LAP', rn_aeiv_0, aeiu, aeiv )    ! surface value proportional to scale factor
            !                                                 ! reduction with depth
            CALL ldf_c1d( 'TRA', r1_4, aeiu(:,:,1), aeiv(:,:,1), aeiu, aeiv )
            !
         CASE(  32  )       !==  time varying 3D field  ==!
            IF(lwp) WRITE(numout,*) '          eddy induced velocity coef. = F( latitude, longitude, depth, time )'
            IF(lwp) WRITE(numout,*) '                              = F( total EKE )   GEOMETRIC parameterization'
            !
            ln_eke_equ = .TRUE.       ! force the eddy energy equation to be updated
            l_eke_eiv  = .TRUE.
            !
         CASE DEFAULT
            CALL ctl_stop('ldf_tra_init: wrong choice for nn_aei_ijk_t, the type of space-time variation of aei')
         END SELECT
         !
      ELSE
          IF(lwp) WRITE(numout,*) '   eddy induced velocity param is NOT used neither diagnosed'
          ln_ldfeiv_dia = .FALSE.
      ENDIF
      !
      IF( ln_eke_equ ) THEN
         l_ldfeke   = .TRUE.          ! GEOMETRIC param initialization done in nemogcm_init
         IF(lwp) WRITE(numout,*) '   updating GEOMETRIC eddy energy equation    ln_eke_equ    = ', ln_eke_equ
      ENDIF

!!! jm (14 Apr 19): if GM is involved then switch on the slopes by default
      IF ( ln_ldfeiv ) THEN
         l_ldfslp = .TRUE.
         IF(lwp) WRITE(numout,*) 'HACK ln_ldfeiv is on, switching on l_ldfslp = ', l_ldfslp
      ENDIF
!!!

   END SUBROUTINE ldf_eiv_init


   SUBROUTINE ldf_eiv( kt, paei0, paeiu, paeiv )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_eiv  ***
      !!
      !! ** Purpose :   Compute the eddy induced velocity coefficient from the
      !!              growth rate of baroclinic instability.
      !!
      !! ** Method  :   coefficient function of the growth rate of baroclinic instability
      !!
      !! Reference : Treguier et al. JPO 1997   ; Held and Larichev JAS 1996
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in   ) ::   kt             ! ocean time-step index
      REAL(wp)                        , INTENT(inout) ::   paei0          ! max value            [m2/s]
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   paeiu, paeiv   ! eiv coefficient      [m2/s]
      !
      INTEGER  ::   ji, jj, jk    ! dummy loop indices
      REAL(wp) ::   zfw, ze3w, zn2, z1_f20, zaht, zaht_min, zzaei   ! local scalars
      REAL(wp), DIMENSION(:,:), POINTER ::   zn, zah, zhw, zross, zaeiw   ! 2D workspace
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )   CALL timing_start('ldf_eiv')
      !
      CALL wrk_alloc( jpi,jpj,   zn, zah, zhw, zross, zaeiw )
      !      
      zn   (:,:) = 0._wp      ! Local initialization
      zhw  (:,:) = 5._wp
      zah  (:,:) = 0._wp
      zross(:,:) = 0._wp
      !                       ! Compute lateral diffusive coefficient at T-point
      IF( ln_traldf_triad ) THEN
         DO jk = 1, jpk
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  ! Take the max of N^2 and zero then take the vertical sum 
                  ! of the square root of the resulting N^2 ( required to compute 
                  ! internal Rossby radius Ro = .5 * sum_jpk(N) / f 
                  zn2 = MAX( rn2b(ji,jj,jk), 0._wp )
                  zn(ji,jj) = zn(ji,jj) + SQRT( zn2 ) * e3w_n(ji,jj,jk)
                  ! Compute elements required for the inverse time scale of baroclinic
                  ! eddies using the isopycnal slopes calculated in ldfslp.F : 
                  ! T^-1 = sqrt(m_jpk(N^2*(r1^2+r2^2)*e3w))
                  ze3w = e3w_n(ji,jj,jk) * tmask(ji,jj,jk)
                  zah(ji,jj) = zah(ji,jj) + zn2 * wslp2(ji,jj,jk) * ze3w
                  zhw(ji,jj) = zhw(ji,jj) + ze3w
               END DO
            END DO
         END DO
      ELSE
         DO jk = 1, jpk
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  ! Take the max of N^2 and zero then take the vertical sum 
                  ! of the square root of the resulting N^2 ( required to compute 
                  ! internal Rossby radius Ro = .5 * sum_jpk(N) / f 
                  zn2 = MAX( rn2b(ji,jj,jk), 0._wp )
                  zn(ji,jj) = zn(ji,jj) + SQRT( zn2 ) * e3w_n(ji,jj,jk)
                  ! Compute elements required for the inverse time scale of baroclinic
                  ! eddies using the isopycnal slopes calculated in ldfslp.F : 
                  ! T^-1 = sqrt(m_jpk(N^2*(r1^2+r2^2)*e3w))
                  ze3w = e3w_n(ji,jj,jk) * tmask(ji,jj,jk)
                  zah(ji,jj) = zah(ji,jj) + zn2 * ( wslpi(ji,jj,jk) * wslpi(ji,jj,jk)   &
                     &                            + wslpj(ji,jj,jk) * wslpj(ji,jj,jk) ) * ze3w
                  zhw(ji,jj) = zhw(ji,jj) + ze3w
               END DO
            END DO
         END DO
      END IF

      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zfw = MAX( ABS( 2. * omega * SIN( rad * gphit(ji,jj) ) ) , 1.e-10 )
            ! Rossby radius at w-point taken < 40km and  > 2km
            zross(ji,jj) = MAX( MIN( .4 * zn(ji,jj) / zfw, 40.e3 ), 2.e3 )
            ! Compute aeiw by multiplying Ro^2 and T^-1
            zaeiw(ji,jj) = zross(ji,jj) * zross(ji,jj) * SQRT( zah(ji,jj) / zhw(ji,jj) ) * tmask(ji,jj,1)
         END DO
      END DO

      !                                         !==  Bound on eiv coeff.  ==!
      z1_f20 = 1._wp / (  2._wp * omega * sin( rad * 20._wp )  )
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zzaei = MIN( 1._wp, ABS( ff_t(ji,jj) * z1_f20 ) ) * zaeiw(ji,jj)       ! tropical decrease
            zaeiw(ji,jj) = MIN( zzaei , paei0 )                                  ! Max value = paei0
         END DO
      END DO
      CALL lbc_lnk( zaeiw(:,:), 'W', 1. )       ! lateral boundary condition
      !               
      DO jj = 2, jpjm1                          !== aei at u- and v-points  ==!
         DO ji = fs_2, fs_jpim1   ! vector opt.
            paeiu(ji,jj,1) = 0.5_wp * ( zaeiw(ji,jj) + zaeiw(ji+1,jj  ) ) * umask(ji,jj,1)
            paeiv(ji,jj,1) = 0.5_wp * ( zaeiw(ji,jj) + zaeiw(ji  ,jj+1) ) * vmask(ji,jj,1)
         END DO 
      END DO 
      CALL lbc_lnk( paeiu(:,:,1), 'U', 1. )   ;   CALL lbc_lnk( paeiv(:,:,1), 'V', 1. )      ! lateral boundary condition

      DO jk = 2, jpkm1                          !==  deeper values equal the surface one  ==!
         paeiu(:,:,jk) = paeiu(:,:,1) * umask(:,:,jk)
         paeiv(:,:,jk) = paeiv(:,:,1) * vmask(:,:,jk)
      END DO
      !  
      CALL wrk_dealloc( jpi,jpj,   zn, zah, zhw, zross, zaeiw )
      !
      IF( nn_timing == 1 )   CALL timing_stop('ldf_eiv')
      !
   END SUBROUTINE ldf_eiv


   SUBROUTINE ldf_eiv_trp( kt, kit000, pun, pvn, pwn, cdtype )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_eiv_trp  ***
      !! 
      !! ** Purpose :   add to the input ocean transport the contribution of 
      !!              the eddy induced velocity parametrization.
      !!
      !! ** Method  :   The eddy induced transport is computed from a flux stream-
      !!              function which depends on the slope of iso-neutral surfaces
      !!              (see ldf_slp). For example, in the i-k plan : 
      !!                   psi_uw = mk(aeiu) e2u mi(wslpi)   [in m3/s]
      !!                   Utr_eiv = - dk[psi_uw]
      !!                   Vtr_eiv = + di[psi_uw]
      !!                ln_ldfeiv_dia = T : output the associated streamfunction,
      !!                                    velocity and heat transport (call ldf_eiv_dia)
      !!
      !! ** Action  : pun, pvn increased by the eiv transport
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in   ) ::   kt       ! ocean time-step index
      INTEGER                         , INTENT(in   ) ::   kit000   ! first time step index
      CHARACTER(len=3)                , INTENT(in   ) ::   cdtype   ! =TRA or TRC (tracer indicator)
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pun      ! in : 3 ocean transport components   [m3/s]
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pvn      ! out: 3 ocean transport components   [m3/s]
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pwn      ! increased by the eiv                [m3/s]
      !!
      INTEGER  ::   ji, jj, jk                 ! dummy loop indices
      REAL(wp) ::   zuwk, zuwk1, zuwi, zuwi1   ! local scalars
      REAL(wp) ::   zvwk, zvwk1, zvwj, zvwj1   !   -      -
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zpsi_uw, zpsi_vw
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )   CALL timing_start( 'ldf_eiv_trp')
      !
      CALL wrk_alloc( jpi,jpj,jpk,   zpsi_uw, zpsi_vw )

      IF( kt == kit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'ldf_eiv_trp : eddy induced advection on ', cdtype,' :'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   add to velocity fields the eiv component'
      ENDIF

      
      zpsi_uw(:,:, 1 ) = 0._wp   ;   zpsi_vw(:,:, 1 ) = 0._wp
      zpsi_uw(:,:,jpk) = 0._wp   ;   zpsi_vw(:,:,jpk) = 0._wp
      !
      DO jk = 2, jpkm1
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zpsi_uw(ji,jj,jk) = - 0.25_wp * e2u(ji,jj) * ( wslpi(ji,jj,jk  ) + wslpi(ji+1,jj,jk) )   &
                  &                                       * ( aeiu (ji,jj,jk-1) + aeiu (ji  ,jj,jk) ) * umask(ji,jj,jk)
               zpsi_vw(ji,jj,jk) = - 0.25_wp * e1v(ji,jj) * ( wslpj(ji,jj,jk  ) + wslpj(ji,jj+1,jk) )   &
                  &                                       * ( aeiv (ji,jj,jk-1) + aeiv (ji,jj  ,jk) ) * vmask(ji,jj,jk)
            END DO
         END DO
      END DO
      !
      DO jk = 1, jpkm1
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.               
               pun(ji,jj,jk) = pun(ji,jj,jk) - ( zpsi_uw(ji,jj,jk) - zpsi_uw(ji,jj,jk+1) )
               pvn(ji,jj,jk) = pvn(ji,jj,jk) - ( zpsi_vw(ji,jj,jk) - zpsi_vw(ji,jj,jk+1) )
            END DO
         END DO
      END DO
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               pwn(ji,jj,jk) = pwn(ji,jj,jk) + (  zpsi_uw(ji,jj,jk) - zpsi_uw(ji-1,jj  ,jk)   &
                  &                             + zpsi_vw(ji,jj,jk) - zpsi_vw(ji  ,jj-1,jk) )
            END DO
         END DO
      END DO
      !
      !                              ! diagnose the eddy induced velocity and associated heat transport
      IF( ln_ldfeiv_dia .AND. cdtype == 'TRA' )   CALL ldf_eiv_dia( zpsi_uw, zpsi_vw )
      !
      CALL wrk_dealloc( jpi,jpj,jpk,   zpsi_uw, zpsi_vw )
      !
      IF( nn_timing == 1 )   CALL timing_stop( 'ldf_eiv_trp')
      !
    END SUBROUTINE ldf_eiv_trp


   SUBROUTINE ldf_eiv_dia( psi_uw, psi_vw )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_eiv_dia  ***
      !!
      !! ** Purpose :   diagnose the eddy induced velocity and its associated
      !!              vertically integrated heat transport.
      !!
      !! ** Method :
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   psi_uw, psi_vw   ! streamfunction   [m3/s]
      !
      INTEGER  ::   ji, jj, jk    ! dummy loop indices
      REAL(wp) ::   zztmp   ! local scalar
      REAL(wp), DIMENSION(:,:)  , POINTER ::   zw2d   ! 2D workspace
      REAL(wp), DIMENSION(:,:,:), POINTER ::   zw3d   ! 3D workspace
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start( 'ldf_eiv_dia')
      !
      !                                                  !==  eiv stream function: output  ==!
      CALL lbc_lnk( psi_uw, 'U', -1. )                         ! lateral boundary condition
      CALL lbc_lnk( psi_vw, 'V', -1. )
      !
!!gm      CALL iom_put( "psi_eiv_uw", psi_uw )                 ! output
!!gm      CALL iom_put( "psi_eiv_vw", psi_vw )
      !
      !                                                  !==  eiv velocities: calculate and output  ==!
      CALL wrk_alloc( jpi,jpj,jpk,   zw3d )
      !
      zw3d(:,:,jpk) = 0._wp                                    ! bottom value always 0
      !
      DO jk = 1, jpkm1                                         ! e2u e3u u_eiv = -dk[psi_uw]
         zw3d(:,:,jk) = ( psi_uw(:,:,jk+1) - psi_uw(:,:,jk) ) / ( e2u(:,:) * e3u_n(:,:,jk) )
      END DO
      CALL iom_put( "uoce_eiv", zw3d )
      !
      DO jk = 1, jpkm1                                         ! e1v e3v v_eiv = -dk[psi_vw]
         zw3d(:,:,jk) = ( psi_vw(:,:,jk+1) - psi_vw(:,:,jk) ) / ( e1v(:,:) * e3v_n(:,:,jk) )
      END DO
      CALL iom_put( "voce_eiv", zw3d )
      !
      DO jk = 1, jpkm1                                         ! e1 e2 w_eiv = dk[psix] + dk[psix]
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1  ! vector opt.
               zw3d(ji,jj,jk) = (  psi_vw(ji,jj,jk) - psi_vw(ji  ,jj-1,jk)    &
                  &              + psi_uw(ji,jj,jk) - psi_uw(ji-1,jj  ,jk)  ) / e1e2t(ji,jj)
            END DO
         END DO
      END DO
      CALL lbc_lnk( zw3d, 'T', 1. )      ! lateral boundary condition
      CALL iom_put( "woce_eiv", zw3d )
      !
      !      
      !
      CALL wrk_alloc( jpi,jpj,   zw2d )
      !
      zztmp = 0.5_wp * rau0 * rcp 
      IF( iom_use('ueiv_heattr') .OR. iom_use('ueiv_heattr3d') ) THEN
        zw2d(:,:)   = 0._wp 
        zw3d(:,:,:) = 0._wp 
        DO jk = 1, jpkm1
           DO jj = 2, jpjm1
              DO ji = fs_2, fs_jpim1   ! vector opt.
                 zw3d(ji,jj,jk) = zw3d(ji,jj,jk) + ( psi_uw(ji,jj,jk+1)      - psi_uw(ji,jj,jk)          )   &
                    &                            * ( tsn   (ji,jj,jk,jp_tem) + tsn   (ji+1,jj,jk,jp_tem) ) 
                 zw2d(ji,jj) = zw2d(ji,jj) + zw3d(ji,jj,jk)
              END DO
           END DO
        END DO
        CALL lbc_lnk( zw2d, 'U', -1. )
        CALL lbc_lnk( zw3d, 'U', -1. )
        CALL iom_put( "ueiv_heattr"  , zztmp * zw2d )                  ! heat transport in i-direction
        CALL iom_put( "ueiv_heattr3d", zztmp * zw3d )                  ! heat transport in i-direction
      ENDIF
      zw2d(:,:)   = 0._wp 
      zw3d(:,:,:) = 0._wp 
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zw3d(ji,jj,jk) = zw3d(ji,jj,jk) + ( psi_vw(ji,jj,jk+1)      - psi_vw(ji,jj,jk)          )   &
                  &                            * ( tsn   (ji,jj,jk,jp_tem) + tsn   (ji,jj+1,jk,jp_tem) ) 
               zw2d(ji,jj) = zw2d(ji,jj) + zw3d(ji,jj,jk)
            END DO
         END DO
      END DO
      CALL lbc_lnk( zw2d, 'V', -1. )
      CALL iom_put( "veiv_heattr", zztmp * zw2d )                  !  heat transport in j-direction
      CALL iom_put( "veiv_heattr", zztmp * zw3d )                  !  heat transport in j-direction
      !
      IF( ln_diaptr )  CALL dia_ptr_hst( jp_tem, 'eiv', 0.5 * zw3d )
      !
      zztmp = 0.5_wp * 0.5
      IF( iom_use('ueiv_salttr') .OR. iom_use('ueiv_salttr3d')) THEN
        zw2d(:,:) = 0._wp 
        zw3d(:,:,:) = 0._wp 
        DO jk = 1, jpkm1
           DO jj = 2, jpjm1
              DO ji = fs_2, fs_jpim1   ! vector opt.
                 zw3d(ji,jj,jk) = zw3d(ji,jj,jk) * ( psi_uw(ji,jj,jk+1)      - psi_uw(ji,jj,jk)          )   &
                    &                            * ( tsn   (ji,jj,jk,jp_sal) + tsn   (ji+1,jj,jk,jp_sal) ) 
                 zw2d(ji,jj) = zw2d(ji,jj) + zw3d(ji,jj,jk)
              END DO
           END DO
        END DO
        CALL lbc_lnk( zw2d, 'U', -1. )
        CALL lbc_lnk( zw3d, 'U', -1. )
        CALL iom_put( "ueiv_salttr", zztmp * zw2d )                  ! salt transport in i-direction
        CALL iom_put( "ueiv_salttr3d", zztmp * zw3d )                  ! salt transport in i-direction
      ENDIF
      zw2d(:,:) = 0._wp 
      zw3d(:,:,:) = 0._wp 
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zw3d(ji,jj,jk) = zw3d(ji,jj,jk) + ( psi_vw(ji,jj,jk+1)      - psi_vw(ji,jj,jk)          )   &
                  &                            * ( tsn   (ji,jj,jk,jp_sal) + tsn   (ji,jj+1,jk,jp_sal) ) 
               zw2d(ji,jj) = zw2d(ji,jj) + zw3d(ji,jj,jk)
            END DO
         END DO
      END DO
      CALL lbc_lnk( zw2d, 'V', -1. )
      CALL iom_put( "veiv_salttr", zztmp * zw2d )                  !  salt transport in j-direction
      CALL iom_put( "veiv_salttr", zztmp * zw3d )                  !  salt transport in j-direction
      !
      IF( ln_diaptr ) CALL dia_ptr_hst( jp_sal, 'eiv', 0.5 * zw3d )
      !
      CALL wrk_dealloc( jpi,jpj,   zw2d )
      CALL wrk_dealloc( jpi,jpj,jpk,   zw3d )
      !
      IF( nn_timing == 1 )  CALL timing_stop( 'ldf_eiv_dia')      
      !
   END SUBROUTINE ldf_eiv_dia

   !!======================================================================
END MODULE ldftra
