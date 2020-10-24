MODULE ldfeke
   !!======================================================================
   !!                       ***  MODULE  ldfeke  ***
   !! Ocean physics:  Eddy induced velocity coefficient according to the 
   !!                 GEOMETRIC parameterization scheme
   !!=====================================================================
   !! History :  4.0  !  2017-11  (J. Mak, G. Madec) original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   ldf_eke       : time step depth-integrated EKE and update aeiw (and 
   !!                   from that aeiu and aeiv) according to the GEOMETRIC 
   !!                   parameterization scheme
   !!   ldf_eke_init  : initialization, namelist read, and parameters control
   !!   eke_rst       : read/write eke restart in ocean restart file
   !!----------------------------------------------------------------------
   USE oce            ! ocean: dynamics and active tracers variables
   USE phycst         ! physical constants
   USE dom_oce        ! domain: ocean
   USE ldfslp         ! lateral physics: slope of iso-neutral surfaces
   USE ldftra         ! lateral physics: eddy coefficients
   USE dynldf_lap_blp, ONLY : eke_keS   ! source term of eke due to KE dissipation
   USE zdfmxl         ! vertical physics: mixed layer
   
   USE dynspg_ts , ONLY : un_adv, vn_adv
   
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE prtctl         ! Print control
   USE timing         ! Timing
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ldf_eke        ! routine called in step module
   PUBLIC   ldf_eke_init   ! routine called in opa module


   !                                 !!** Namelist  namldf_eke  **
   REAL(wp) ::   rn_ekedis                  !  dissipation time scale of EKE                  [days]
   REAL(wp) ::   rn_geom                    !  geometric parameterization master coefficient     [-]
   REAL(wp) ::   rn_eke_lap                 !  diffusion of EKE                               [m2/s]
   REAL(wp) ::   rn_eke_init                !  initial value of total EKE                    [m2/s2]
   REAL(wp) ::   rn_eke_min                 !  background value of total EKE                 [m2/s2]
   REAL(wp) ::   rn_ross_min                !  tapering based of minimum Rossby radius           [m]
   REAL(wp) ::   rn_aeiv_min, rn_aeiv_max   !  min and max bounds of aeiv coefficient         [m2/s] 
   REAL(wp) ::   rn_SFmin, rn_SFmax         !  min and max bounds of Structure Function          [-]
   REAL(wp) ::   zf0, zbeta                 !  f0 and beta for computing Rossby speed
   !
   INTEGER  ::   nn_eke_opt                 !  option for what term to include in eddy energy budget
   INTEGER  ::   nn_eke_dis                 !  option for taking constant or spatially varying linear dissipation
   !
   LOGICAL  ::   ln_adv_wav                 !  option for having advection by Rossby wave or not
   LOGICAL  ::   ln_beta_plane              !  option for computing long Rossby phase speed
   !
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   eke_b, eke_n, eke_a   ! vertical sum of total Eddy Kinetic Energy [m3/s2]
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   r1_ekedis             ! linear dissipation rate (= 1/rn_ekedis)   [  /s ]
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   zc1, zc_ros, zadv_wav ! 1st baroclinic mode and long Rossby speed [m /s ]

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2017)
   !! $Id: ldfeke.F90 7813 2017-03-20 16:17:45Z clem $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ldf_eke( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE ldf_eke  ***
      !!
      !! ** Purpose :   implements the GEOMETRIC parameterization 
      !!
      !! ** Notes   :   If nn_aei_ijk_t = 32 then eke and aeiw are BOTH updated
      !!                If ln_eke_equ = .true. in namtra_ldfeiv but nn_aei_ijk_t 
      !!              is something else, then ONLY eddy equation is updated
      !!              (but the eddy energy is passive and doesn't do anything)
      !!
      !! ** Method  :   GEOMETRIC calculates the Gent-McWilliams / eddy induced 
      !!              velocity coefficient according to
      !!
      !!                    aeiw = alpha * (\int E dz) / (\int S M^2 / N dz),
      !!
      !!              where (\int E dz) is the depth-integrated eddy energy 
      !!              (at the previous time level), informed by a parameterized
      !!              depth-integrated eddy energy, where
      !!
      !!    nn_eke_opt    =  0 => default: just PE->EKE growth and linear dissipation
      !!                  !
      !!                  =  1 => default + advection by depth-averaged flow
      !!                  !
      !!                  =  2 => default + advection + contribution to EKE from 
      !!                          momentum dissipation
      !!                  !
      !!                  = 88 => ONLY advection
      !!                  !
      !!                  = 99 => ONLY Laplacian diffusion
      !!
      !!              S is a structure function, and M and N are horizontal and 
      !!              vertical buoyancy frequencies
      !!                  
      !!                  linear dissipation may be specified by
      !!
      !!    nn_eke_dus    =  0 => constant
      !!                  !
      !!                  =-20 => read in a geom_diss_2D.nc field (give it in days)
      !!
      !! ** Action  : * time step depth-integrated eddy energy budget
      !!              * calculate aeiw
      !!
      !! References : Marshall, Maddison, Berloff JPO 2012   ; Mak et al JPO 2018
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !
      INTEGER  ::   ji, jj, jk          ! dummy loop arguments
      INTEGER  ::   ik                  ! local integer
      REAL(wp) ::   zn2_min = 1.e-8     ! minimum value of N2 used to compute structure function
      REAL(wp) ::   z1_f20, zfw, ze3w   ! local scalar
      REAL(wp) ::   zfp_ui, zfp_vj      !   -      -
      REAL(wp) ::   zfm_ui, zfm_vj      !   -      -
      REAL(wp) ::   zn_slp2, zn2        !   -      -
      REAL(wp) ::   zmsku, zaeiu_w      !   -      -
      REAL(wp) ::   zmskv, zaeiv_w      !   -      -
      REAL(wp) ::   zen, zed            !   -      -
      REAL(wp) ::   zeke_rhs            !   -      -
      REAL(wp) ::   zck, zwslpi, zwslpj !   -      -  tapering near coasts
      REAL(wp) ::   zc_rosu             !   -      -
      REAL(wp), DIMENSION(jpi,jpj)     ::   zeke_peS, zn_slp                 ! 2D workspace, PE-KE conversion
      REAL(wp), DIMENSION(jpi,jpj)     ::   zadv_ubt, zwx, zwy               !  -     -      barotropic advection
      REAL(wp), DIMENSION(jpi,jpj)     ::   zlap, zaheeu, zaheev, ztu, ztv   !  -     -      diffusion
      REAL(wp), DIMENSION(jpi,jpj)     ::   zdis                             !  -     -      linear dissipation
      REAL(wp), DIMENSION(jpi,jpj)     ::   zn, zross                        !  -     -      tapering near coasts
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zaeiw, zzSF                      ! 3D workspace
      !!--------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('ldf_eke')
      !
      IF( kt == nit000 .AND. lwp) THEN       !* Control print
         WRITE(numout,*)
         WRITE(numout,*) 'ldf_eke : GEOMETRIC parameterization (total EKE time evolution equation)'
         WRITE(numout,*) '~~~~~~~'
      ENDIF

      !                    !==  EIV mean conversion to EKE (ah N^2 slp2) & N slp  ==!
      !
      !                         work out the 3D structure function here
      !
      ! current: Ferreria et al, zzSF = N^2 / N^2_ref, W-points
      !                          capped between rn_SFmax and rn_SFmin
      zzSF(:,:,:) = 0._wp
      DO jk = 2, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( jk <= nmln(ji,jj) ) THEN ! above and at mixed layer
                  zzSF(ji,jj,jk) = rn_SFmax
               ELSE
                  ik   = MIN( nmln(ji,jj), mbkt(ji,jj) ) + 1     ! one level below the mixed layer (MIN in case ML depth is the ocean depth)
                  zzSF(ji,jj,jk) = MAX( 0._wp , rn2b(ji,jj,jk) ) / MAX( zn2_min , rn2b(ji,jj,ik) )      ! Structure Function : N^2 / N^2_ref
                  zzSF(ji,jj,jk) = MAX(  rn_SFmin , MIN( zzSF(ji,jj,jk) , rn_SFmax )  )
               ENDIF
            END DO
         END DO
      END DO
      !
      !                         !*  parametrized PE_EKE conversion due to eddy induced velocity
      zeke_peS(:,:) = 0._wp
      zn_slp  (:,:) = 0._wp
      zn      (:,:) = 0._wp
      DO jk = 2, jpkm1   ! query: index?
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.      ! NB: ah_slp2 is w-masked
               zmsku = wmask(ji,jj,jk) / MAX(   umask(ji  ,jj,jk-1) + umask(ji-1,jj,jk)          &
                  &                           + umask(ji-1,jj,jk-1) + umask(ji  ,jj,jk) , 1._wp  )
               zmskv = wmask(ji,jj,jk) / MAX(   vmask(ji,jj  ,jk-1) + vmask(ji,jj-1,jk)          &
                  &                           + vmask(ji,jj-1,jk-1) + vmask(ji,jj  ,jk) , 1._wp  )
               !
               zaeiu_w = (   aeiu(ji  ,jj,jk-1) + aeiu(ji-1,jj,jk)    &
                  &        + aeiu(ji-1,jj,jk-1) + aeiu(ji  ,jj,jk)  ) * zmsku
               zaeiv_w = (   aeiv(ji,jj  ,jk-1) + aeiv(ji,jj-1,jk)    &
                  &        + aeiv(ji,jj-1,jk-1) + aeiv(ji,jj  ,jk)  ) * zmskv
               !
               zn_slp2 = (   zaeiu_w * wslpi(ji,jj,jk) * wslpi(ji,jj,jk)   &      ! (slope ** 2) * aeiv
                  &        + zaeiv_w * wslpj(ji,jj,jk) * wslpj(ji,jj,jk)   )      ! JM 28 Jun: undo slope reduction here too?
               zn2     = MAX( 0._wp , rn2b(ji,jj,jk) )
               !
               ze3w      = e3w_b(ji,jj,jk) * tmask(ji,jj,jk)
               zn(ji,jj) = zn(ji,jj) + SQRT( zn2 ) * ze3w                         ! for working out taper at small rossby radius regions later
               !
               zeke_peS(ji,jj) = zeke_peS(ji,jj) + ze3w * zn2 * zn_slp2           ! note this is >=0
               !
               zck     =   ( umask(ji,jj,jk) + umask(ji-1,jj,jk) )   &            ! taken from ldfslp, undo the slope reduction
                  &      * ( vmask(ji,jj,jk) + vmask(ji,jj-1,jk) ) * 0.25         !   near topographic features
               zwslpi  = wslpi(ji,jj,jk) / MAX( zck, 0.1_wp)                      ! (just to avoid dividing by zeros)
               zwslpj  = wslpj(ji,jj,jk) / MAX( zck, 0.1_wp)
               !
               zn_slp(ji,jj) = zn_slp(ji,jj) + ze3w * zzSF(ji,jj,jk)   &          ! note this >=0 and structure function weighted
                  &                          * SQRT( zn2 * ( zwslpi * zwslpi + zwslpj * zwslpj )  )
            END DO
         END DO
      END DO

      !                    !*  upstream advection with initial mass fluxes & intermediate update
      !                          !* upstream tracer flux in the i and j direction
      DO jj = 1, jpjm1
         DO ji = 1, fs_jpim1   ! vector opt.
            ! upstream scheme
            zfp_ui = un_adv(ji,jj) + ABS( un_adv(ji,jj) )
            zfm_ui = un_adv(ji,jj) - ABS( un_adv(ji,jj) )
            zfp_vj = vn_adv(ji,jj) + ABS( vn_adv(ji,jj) )
            zfm_vj = vn_adv(ji,jj) - ABS( vn_adv(ji,jj) )
            zwx(ji,jj) = 0.5 * ( zfp_ui * eke_b(ji,jj) + zfm_ui * eke_b(ji+1,jj  ) )
            zwy(ji,jj) = 0.5 * ( zfp_vj * eke_b(ji,jj) + zfm_vj * eke_b(ji  ,jj+1) )
         END DO
      END DO
      !                           !* divergence of ubt advective fluxes
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zadv_ubt(ji,jj) = - (  zwx(ji,jj) - zwx(ji-1,jj  )   &
               &                 + zwy(ji,jj) - zwy(ji  ,jj-1) ) * r1_e1e2t(ji,jj)
         END DO
      END DO
      IF( ln_linssh ) THEN                !* top value   (linear free surf. only as zwz is multiplied by wmask)
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zadv_ubt(ji,jj) = - wn(ji,jj,1) * eke_b(ji,jj) / ( ht_0(ji,jj) + 1._wp-tmask(ji,jj,1) )   ! jm (03 Mar 18): was eke_n
            END DO
         END DO
      ENDIF
      !
      !                          !* same as above but for advection by Rossby waves
      zadv_wav(:,:) = 0._wp
      IF ( ln_adv_wav ) THEN
         !
         DO jj = 2, jpjm1 ! first work out the Rossby phase speeds
            DO ji = 2, fs_jpim1
               ! compute only for deep enough places 
               IF ( ht_0(ji,jj) > 300._wp ) THEN   ! jm: should use something like ht_b really...
                  ! compute vertical mode phase speed on T point
                  zc1(ji,jj) = MIN( 10._wp, zn(ji,jj) / rpi )
                  ! compute long Rossby phase speed on T point (minus sign later)
                  IF ( ln_beta_plane ) THEN
                     zc_ros(ji,jj) = zc1(ji,jj) * zc1(ji,jj) * zbeta / (zf0 * zf0)
                  ELSE
                     zc_ros(ji,jj) =  zc1(ji,jj) * zc1(ji,jj) * COS( rad * gphit(ji,jj) )   &
                                   / (  ra * ff_t(ji,jj) * SIN( rad * gphit(ji,jj) )        &
                                   + rsmall  )
                  END IF
                  ! cap the Rossby phase speeds by the fastest equatorial Rossby wave speed
                  ! the minus sign for westward propagation goes here
                  zc_ros(ji,jj) = -MIN( zc1(ji,jj) / 3._wp, zc_ros(ji,jj) )
               END IF
            END DO
         END DO
         !
         zwx(:,:) = 0._wp ! wipe the advective contribution from above
         !
         DO jj = 2, jpjm1 ! average Rossby phase speeds onto U points
            DO ji = 1, fs_jpim1
               zc_rosu    = 0.5 * ( zc_ros(ji,jj) + zc_ros(ji+1,jj) ) * umask(ji,jj,1)
               zfp_ui     = zc_rosu + ABS( zc_rosu )
               zfm_ui     = zc_rosu - ABS( zc_rosu )
               zwx(ji,jj) = 0.5 * ( zfp_ui * eke_b(ji,jj) + zfm_ui * eke_b(ji+1,jj) )
            END DO
         END DO
         !                    !* divergence of wav advective fluxes (only e1 here)
         z1_f20 = 1._wp / (  2._wp * omega * sin( rad * 20._wp )  )
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zadv_wav(ji,jj) = - ( zwx(ji,jj) - zwx(ji-1,jj) ) * r1_e1t(ji,jj)
               zadv_wav(ji,jj) = zadv_wav(ji,jj) * MIN(  1._wp, ABS( ff_t(ji,jj) * z1_f20 )  )   ! tropical decrease
            END DO
         END DO
      END IF
      !
                                 !* divergence of diffusive fluxes
                                 !  code adapted from traldf_lap_blp.F90
      IF ( rn_eke_lap >= 0._wp ) THEN
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zaheeu(ji,jj) = rn_eke_lap * e2_e1u(ji,jj) * umask(ji,jj,1)   ! rn_eke_lap is constant (for now) and NOT masked
               zaheev(ji,jj) = rn_eke_lap * e1_e2v(ji,jj) * vmask(ji,jj,1)   !      before it is pahu and pahv which IS maked
            END DO
         END DO
         DO jj = 1, jpjm1         !== First derivative (gradient) ==!
            DO ji = 1, fs_jpim1
               ztu(ji,jj) = zaheeu(ji,jj) * ( eke_b(ji+1,jj  ) - eke_b(ji,jj) )
               ztv(ji,jj) = zaheev(ji,jj) * ( eke_b(ji  ,jj+1) - eke_b(ji,jj) )
            END DO
         END DO
         DO jj = 2, jpjm1         !== Second derivative (divergence), form trend  ==!
            DO ji = fs_2, fs_jpim1
               zlap(ji,jj) = (  ztu(ji,jj) - ztu(ji-1,jj  )     &
                  &           + ztv(ji,jj) - ztv(ji  ,jj-1) )   &
                  &          / e1e2t(ji,jj)
            END DO
         END DO
      ELSE
         zlap(:,:) = 0._wp
      END IF
                                 !* form the trend for linear dissipation
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1
            zdis(ji,jj) = - r1_ekedis(ji,jj) * (eke_b(ji,jj) - rn_eke_min) * tmask(ji,jj,1)
         END DO
      END DO
      !                    !==  time stepping of EKE Eq.  ==!
      !
      ! note: the rn_eke_min term is a forcing in eddy equation, thus damping in mean equation
      !       added to prevent overshoots and oscillations of energy from exponential growth/decay
      !       maintains a background (depth-integrated) eddy energy level
      !
      zeke_rhs = 0._wp
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            !
            SELECT CASE( nn_eke_opt )   ! Specification of what to include in EKE budget
            !
            CASE(   0  )  !  default: just PE->EKE growth and linear dissipation
               zeke_rhs =                                       zeke_peS(ji,jj) + zdis(ji,jj) + zlap(ji,jj)
            CASE(   1  )  !  as default but with full advection
               zeke_rhs = - zadv_ubt(ji,jj) + zadv_wav(ji,jj) + zeke_peS(ji,jj) + zdis(ji,jj) + zlap(ji,jj)
            CASE(   2  )  !  full thing with additional KE->EKE growth 
               zeke_rhs = - zadv_ubt(ji,jj) + zadv_wav(ji,jj) + zeke_peS(ji,jj) + zdis(ji,jj) + zlap(ji,jj) + eke_keS(ji,jj)
            CASE(  88  )  !  ONLY advection by mean flow
               zeke_rhs = - zadv_ubt(ji,jj)
            CASE(  99  )  !  ONLY diffusion
               zeke_rhs =   zlap(ji,jj)
            CASE DEFAULT
               CALL ctl_stop('ldf_eke: wrong choice nn_eke_opt, set at least to 0 (default)')
            END SELECT
            !
            ! leap-frog 
            eke_a(ji,jj) = eke_b(ji,jj) + r2dt * zeke_rhs * ssmask(ji,jj)
            !
         END DO
      END DO
      CALL lbc_lnk( eke_a, 'T', 1. )  ! Lateral boundary conditions on zwi  (unchanged sign)

      IF( neuler == 0 .AND. kt == nit000 ) THEN       ! Euler time-stepping at first time-step (only swap)
         eke_n(:,:) = eke_a(:,:)    
         !
      ELSE                                            ! Leap-Frog + Asselin filter time stepping
         DO jj = 1, jpj
            DO ji = 1, jpi
               zen = eke_n(ji,jj)                                    
               zed = eke_a(ji,jj) - 2._wp * zen + eke_b(ji,jj)  ! time laplacian on tracers
               !
               eke_b(ji,jj) = zen + atfp * zed                  ! eke_b <-- filtered eke_n 
               eke_n(ji,jj) = eke_a(ji,jj)                      ! ele_n <-- eke_a
            END DO
        END DO
      ENDIF
      
      ! initialise it here so XIOS stops complaining...
      zross(:,:)   = 0._wp
      zaeiw(:,:,:) = 0._wp
      IF( l_eke_eiv ) THEN
         !                    !==  resulting EIV coefficient  ==!
         !
         !                          ! surface value
         z1_f20 = 1._wp / (  2._wp * omega * sin( rad * 20._wp )  )
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ! Rossby radius at w-point, Ro = .5 * sum_jpk(N) / f
               zfw = MAX( ABS( ff_t(ji,jj) ) , 1.e-10 )
               zross(ji,jj) = 0.5 * zn(ji,jj) / zfw
               !
               zaeiw(ji,jj,1) = rn_geom * eke_b(ji,jj) / MAX( 1.e-10 , zn_slp(ji,jj) ) * tmask(ji,jj,1)   ! zn_slp has zzSF multiplied to it
               zaeiw(ji,jj,1) = MIN(  rn_aeiv_max, zaeiw(ji,jj,1)  )                                  ! bound aeiv from above
               zaeiw(ji,jj,1) = zaeiw(ji,jj,1)      &
                  ! tanh taper to deal with some some large values near coast
                  &           * 0.5 * (   1._wp - TANH(  ( -ht_0(ji,jj) + 800._wp     ) / 300._wp  )   )   &
                  ! tanh taper of aeiv on internal Rossby radius
                  &           * 0.5 * (   1._wp + TANH(  ( zross(ji,jj) - rn_ross_min ) / 2._wp    )   )
               zaeiw(ji,jj,1) = zaeiw(ji,jj,1) * MIN(  1._wp, ABS( ff_t(ji,jj) * z1_f20 )  )              ! tropical decrease
               zaeiw(ji,jj,1) = MAX(  rn_aeiv_min, zaeiw(ji,jj,1)  )                                  ! bound aeiv from below
            END DO
         END DO
         CALL lbc_lnk( zaeiw(:,:,1), 'W', 1. )
         !
         !                          ! inner value
         ! set bottom to zero and use the un-masked zaeiw(ji,jj,1) first...
         zaeiw(:,:,jpk) = 0._wp   
         DO jk = 2, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zaeiw(ji,jj,jk) = zaeiw(ji,jj,1) * zzSF(ji,jj,jk) * wmask(ji,jj,jk)   ! zzSF has already been capped
               END DO
            END DO
         END DO
         !                          ! aei at u- and v-points
         aeiu(:,:,jpk) = 0._wp   
         aeiv(:,:,jpk) = 0._wp   
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  aeiu(ji,jj,jk) =    (  zaeiw(ji,jj,jk  ) + zaeiw(ji+1,jj,jk  )    &
                     &                 + zaeiw(ji,jj,jk+1) + zaeiw(ji+1,jj,jk+1)  ) &
                     &           / MAX(  wmask(ji,jj,jk  ) + wmask(ji+1,jj,jk  )    &
                     &                 + wmask(ji,jj,jk+1) + wmask(ji+1,jj,jk+1) , 1._wp ) * umask(ji,jj,jk)
                  aeiv(ji,jj,jk) =    (  zaeiw(ji,jj,jk  ) + zaeiw(ji,jj+1,jk  )    &
                     &                 + zaeiw(ji,jj,jk+1) + zaeiw(ji,jj+1,jk+1)  ) &
                     &           / MAX(  wmask(ji,jj,jk  ) + wmask(ji,jj+1,jk  )    &
                     &                 + wmask(ji,jj,jk+1) + wmask(ji,jj+1,jk+1) , 1._wp ) * vmask(ji,jj,jk)
               END DO 
            END DO 
         END DO 
         !                    !--  diagnostics  --!
         CALL lbc_lnk( aeiu , 'U', 1. )
         CALL lbc_lnk( aeiv , 'V', 1. )
         CALL lbc_lnk( zross, 'W', 1. )
      END IF
      !
      IF( lrst_oce )   CALL eke_rst( kt, 'WRITE' )
      
      !                    !==  output EKE related variables  ==!
      CALL lbc_lnk( zadv_ubt , 'T', 1. )
      CALL lbc_lnk( zadv_wav , 'T', 1. )
      CALL lbc_lnk( zdis     , 'T', 1. )
      CALL lbc_lnk( zlap     , 'T', 1. )
      CALL lbc_lnk( zeke_peS , 'T', 1. )
      CALL lbc_lnk( zc1      , 'W', 1. )
      CALL lbc_lnk( zc_ros   , 'W', 1. )
      !
      CALL iom_put( "eke"            , eke_n )         ! parameterized total EKE (EPE+ EKE)
      CALL iom_put( "trd_eke_adv_ubt", zadv_ubt )      ! ubt    advective trend of EKE(LHS)
      CALL iom_put( "trd_eke_adv_wav", zadv_wav )      ! rossby advective trend of EKE(LHS)
      CALL iom_put( "trd_eke_dis"    , zdis )          ! dissipative trend of EKE     (RHS)
      CALL iom_put( "trd_eke_lap"    , zlap )          ! diffusive trend of EKE       (RHS)
      CALL iom_put( "trd_eke_peS"    , zeke_peS )      ! PE to EKE source trend       (RHS)
      CALL iom_put( "trd_eke_keS"    ,  eke_keS )      ! KE to EKE source trend       (RHS)
      CALL iom_put( "aeiv_geom"      , zaeiw )         ! eddy induced coefficient from GEOMETRIC param
      CALL iom_put( "rossby_rad"     , zross )         ! internal Rossby deformation radius
      CALL iom_put( "c1_vert"        , zc1   )         ! 1st baroclinic mode phase speed
      CALL iom_put( "c_ros"          , zc_ros   )      ! long Rossby phase speed
!!
!! jm: modified ldf_tra in ldftra.F90 to include the ln_eke_equ flag 
!!     into consideration, otherwise iom_put is called twice aeiu and aeiv
!!
      CALL iom_put( "aeiu_2d", aeiu(:,:,1) )       ! surface u-EIV coeff.
      CALL iom_put( "aeiv_2d", aeiv(:,:,1) )       ! surface v-EIV coeff.
      CALL iom_put( "aeiu_3d", aeiu(:,:,:) )       ! 3D      u-EIV coeff.
      CALL iom_put( "aeiv_3d", aeiv(:,:,:) )       ! 3D      v-EIV coeff.
      !
      IF( nn_timing == 1 )  CALL timing_stop('ldf_eke')
      !
   END SUBROUTINE ldf_eke

   SUBROUTINE ldf_eke_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_eke_init  ***
      !!                     
      !! ** Purpose :   Initialization of the depth-integrated eke, vertical
      !!              structure function and max/min caps of aeiw when using the 
      !!              GEOMETRIC parameterization
      !!
      !! ** Method  :   Read the namldf_eke namelist and check the parameters
      !!              called at the first timestep (nit000)
      !!
      !! ** input   :   Namlist namldf_eke
      !!----------------------------------------------------------------------
      INTEGER ::   ios, inum
      INTEGER ::   ji, jj
      !!
      NAMELIST/namldf_eke/  rn_ekedis  , nn_eke_dis ,   & ! GEOM master params (lambda and option, alpha, e0)
         &                  rn_geom    , rn_eke_init,   &
         &                  rn_eke_lap ,                & ! coeff of Laplacian diffusion of EKE
         &                  rn_eke_min ,                & ! background value of (depth-integrated) EKE
         &                  rn_ross_min,                & ! taper aeiv based on Rossby internal radius
         &                  rn_aeiv_max, rn_aeiv_min,   & ! caps and options on result aeiv
         &                  rn_SFmin,    rn_SFmax,      & ! caps on vertical structure function
         &                  nn_eke_opt ,                & ! option for EKE budget terms
         &                  ln_adv_wav ,                & ! option for advection by Rossby wave
         &                  ln_beta_plane                 ! beta plane option for computing Rossby wave speed
      !!----------------------------------------------------------------------
      !
      !
! note: GEOMETRIC is *** off *** by default (so no ref list) 
!      REWIND( numnam_ref )              ! Namelist namldf_eke in reference namelist : Turbulent Kinetic Energy
!      READ  ( numnam_ref, namldf_eke, IOSTAT = ios, ERR = 901)
!901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namldf_eke in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namldf_eke in configuration namelist : GEOMETRIC parameterization
      READ  ( numnam_cfg, namldf_eke, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namldf_eke in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namldf_eke )
      !
      IF(lwp) THEN                    !* Control print
         WRITE(numout,*)
         WRITE(numout,*)    'ldf_eke_init : GEOMETRIC parameterization (total EKE time evolution equation)'
         WRITE(numout,*)    '~~~~~~~~~~~~'
         WRITE(numout,*)    '   Namelist namldf_eke : set the GEOMETRIC parameters'
         WRITE(numout,*)    '      aeiw updated according to GEOMETRIC             ln_eke_eiv  = ', l_eke_eiv
         WRITE(numout,*)    '      dissipation time scale of EKE in days           rn_ekedis   = ', rn_ekedis, ' days'
         WRITE(numout,*)    '         option for dissipation field                 nn_eke_dis  = ', nn_eke_dis
         SELECT CASE( nn_eke_dis ) 
         !
         CASE(   0  )
            WRITE(numout,*) '         spatially constant                                         '
         CASE(   2  )
            WRITE(numout,*) '         read in a 2d varying file geom_diss_2D.nc                  '
         END SELECT
         !
         WRITE(numout,*)    '      geometric parameterization master coefficient   rn_geom     = ', rn_geom
         WRITE(numout,*)    '      coeff of Lapican diffusion of EKE               rn_eke_lap  = ', rn_eke_lap
         WRITE(numout,*)    '      initial total EKE value                         rn_eke_init = ', rn_eke_init
         WRITE(numout,*)    '      background values of total EKE                  rn_eke_min  = ', rn_eke_min
         WRITE(numout,*)    '      taper aeiv subject to min Rossby radius         rn_ross_min = ', rn_ross_min
         WRITE(numout,*)    '      maximum bound of aeiv coeff                     rn_aeiv_max = ', rn_aeiv_max
         WRITE(numout,*)    '      minimum bound of aeiv coeff                     rn_aeiv_min = ', rn_aeiv_min
         WRITE(numout,*)    '      minimum bound of Structure Function             rn_SFmin    = ', rn_SFmin
         WRITE(numout,*)    '      maximum bound of Structure Function             rn_SFmax    = ', rn_SFmax
         WRITE(numout,*)    '         [set rnSFmin = rnSFmax = 1 to use aeiv = aeiv(x,y,t)]      '
         WRITE(numout,*)    '      option for terms in eddy energy budget          nn_eke_opt  = ', nn_eke_opt
         WRITE(numout,*)    '         default: just PE->EKE growth and linear dissipation        '
         SELECT CASE( nn_eke_opt ) 
         !
         CASE(   0  )
            WRITE(numout,*) '         default                                                     '
         CASE(   1  )
            WRITE(numout,*) '         default + advection                                         '
         CASE(   2  )
            WRITE(numout,*) '         default + advection + KE->EKE                               '
         CASE(  88  )
            WRITE(numout,*) '         ONLY advection by z-avg mean flow (no growth or dissipation)'          
         CASE(  99  )
            WRITE(numout,*) '         ONLY Laplacian diffusion          (no growth or dissipation)'
         END SELECT
         WRITE(numout,*)    '      advection of energy by Rossby waves             ln_adv_wav  =  ', ln_adv_wav
      ENDIF
      !                                ! allocate eke arrays
      ALLOCATE( eke_b(jpi,jpj)     , eke_n(jpi,jpj)     , eke_a(jpi,jpj)     , eke_keS(jpi,jpj) ,   &
         &      r1_ekedis(jpi,jpj) , zc1(jpi,jpj)       , zc_ros(jpi,jpj)    , zadv_wav(jpi,jpj)     )
      !
      SELECT CASE( nn_eke_dis )               ! Specification of linear dissipation
      !
      CASE(   0  )      !==  constant  ==!
         IF(lwp) WRITE(numout,*) '      linear EKE dissipation coef. = constant = ', rn_ekedis, ' days'
         r1_ekedis(:,:) = 1._wp / (rn_ekedis * rday)
         !
      CASE( -20  )      !== fixed horizontal shape read in file  ==!
         IF(lwp) WRITE(numout,*) '      linear EKE dissipation coef. = F(i,j) read in geom_diss_2D.nc file'
         CALL iom_open ( 'geom_diss_2D.nc', inum )
         CALL iom_get  ( inum, jpdom_data, 'r1_ekedis', r1_ekedis ) ! read in as time-scale in days...
         CALL iom_close( inum )
         DO jj = 1, jpj
            DO ji = 1, jpi
               r1_ekedis(ji,jj) = 1._wp / (r1_ekedis(ji,jj) * rday)   ! ...convert rate in per second
            END DO
         END DO
         !
      CASE DEFAULT
         CALL ctl_stop('ldf_eke_init: wrong choice for nn_eke_dis, the option for linear dissipation in GEOMETRIC')
      END SELECT
      !
      CALL eke_rst( nit000, 'READ' )   !* read or initialize all required files
      !
      ! tidy up the initilisation of zc1
      zc_ros(:,:)   = 0._wp
      zadv_wav(:,:) = 0._wp
      IF ( ln_adv_wav ) THEN
         DO jj = 2, jpjm1
            DO ji = 1, fs_jpim1
               ! set the ones that are not deep enough to be zero
               IF ( ht_0(ji,jj) < 300._wp ) THEN   ! jm: should use something like ht_b really...
                  zc1(ji,jj) = 0._wp
               END IF
            END DO
         END DO
      END IF
      ! if using beta_plane, compute beta and f0 locally
      IF ( ln_beta_plane ) THEN
        zf0 = ff_f(1,1)
        zbeta = (ff_f(1,2) - zf0) / e2t(1,1)
        IF(lwp) WRITE(numout,*)    '         beta plane option is    ln_beta_plane =', ln_beta_plane
        IF(lwp) WRITE(numout,*)    '         f0   = ', zf0,   '    s-1'
        IF(lwp) WRITE(numout,*)    '         beta = ', zbeta, 'm-1 s-1'
      END IF
      !
      !
   END SUBROUTINE ldf_eke_init


   SUBROUTINE eke_rst( kt, cdrw )
     !!---------------------------------------------------------------------
     !!                   ***  ROUTINE eke_rst  ***
     !!                     
     !! ** Purpose :   Read or write EKE file (eke) in restart file
     !!
     !! ** Method  :   use of IOM library
     !!              if the restart does not contain EKE, eke is set
     !!              according to rn_eke_init, and aeiu = aeiv = 10 m s^-2
     !!----------------------------------------------------------------------
     INTEGER         , INTENT(in) ::   kt     ! ocean time-step
     CHARACTER(len=*), INTENT(in) ::   cdrw   ! "READ"/"WRITE" flag
     !
     INTEGER ::   jit, jk                   ! dummy loop indices
     INTEGER ::   id1, id2, id3, id4, id5   ! local integer
     !!----------------------------------------------------------------------
     !
     IF( TRIM(cdrw) == 'READ' ) THEN        ! Read/initialise 
        !                                   ! ---------------
        IF( ln_rstart ) THEN                   !* Read the restart file
           id1 = iom_varid( numror, 'eke_b'   , ldstop = .FALSE. )
           id2 = iom_varid( numror, 'eke_n'   , ldstop = .FALSE. )
           id3 = iom_varid( numror, 'aeiu'    , ldstop = .FALSE. )
           id4 = iom_varid( numror, 'aeiv'    , ldstop = .FALSE. )
           IF(lwp) write(numout,*) ' ldfeke init restart'
           !
           IF( MIN( id1, id2, id3, id4 ) > 0 ) THEN        ! all required arrays exist
              IF(lwp) write(numout,*) ' all files for ldfeke exist, loading ...'
              CALL iom_get( numror, jpdom_autoglo, 'eke_b', eke_b )
              CALL iom_get( numror, jpdom_autoglo, 'eke_n', eke_n )
              CALL iom_get( numror, jpdom_autoglo, 'aeiu' , aeiu )
              CALL iom_get( numror, jpdom_autoglo, 'aeiv' , aeiv )
           ELSE                                                 ! one at least array is missing
              IF(lwp) write(numout,*) ' not all files for ldfeke exist '
              IF(lwp) write(numout,*) '    --- initialize from namelist'
              eke_b(:,:)  = rn_eke_init * ht_n(:,:) * ssmask(:,:)
              eke_n(:,:)  = eke_b(:,:) 
              aeiu(:,:,:) = 10._wp * umask(:,:,:)    ! bottom eddy coeff set to zero at last level
              aeiv(:,:,:) = 10._wp * vmask(:,:,:)
           ENDIF
           ! ln_wav_adv restarts + initialisations
           IF( ln_adv_wav ) THEN
!              id5 = iom_varid( numror, 'zc1' , ldstop = .FALSE. )
!              IF( id5 > 0 ) THEN
!                 IF(lwp) write(numout,*) '    ln_adv_wav on and found zc1 variable, loading ...'
!                 CALL iom_get( numror, jpdom_autoglo, 'zc1'  , zc1 )
!              ELSE
                 IF(lwp) write(numout,*) '    ln_adv_wav on but not found zc1 variable, initialising ...'
                 zc1(:,:)    =  4._wp  ! set it to something everywhere but some are set to zero in ldf_eke_init
!              END IF
           END IF
        ELSE                                          !* Start from rest with a non zero value (required)
            IF(lwp) write(numout,*) ' ldfeke init restart from namelist'
            eke_b(:,:)  = rn_eke_init * ht_n(:,:) * ssmask(:,:)
            eke_n(:,:)  = eke_b(:,:) 
            aeiu(:,:,:) = 10._wp * umask(:,:,:)       ! bottom eddy coeff set to zero at last level
            aeiv(:,:,:) = 10._wp * vmask(:,:,:)
            zc1(:,:)    =  4._wp
            !                              
        ENDIF
        !
     ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN   ! Create restart file
        !                                   ! -------------------
        IF(lwp) WRITE(numout,*) '---- eke-rst ----'
        CALL iom_rstput( kt, nitrst, numrow, 'eke_b', eke_b  )
        CALL iom_rstput( kt, nitrst, numrow, 'eke_n', eke_n  )
        CALL iom_rstput( kt, nitrst, numrow, 'aeiu' , aeiu )
        CALL iom_rstput( kt, nitrst, numrow, 'aeiv' , aeiv )
        CALL iom_rstput( kt, nitrst, numrow, 'zc1'  , zc1 )
        !
     ENDIF
     !
   END SUBROUTINE eke_rst

   !!======================================================================
END MODULE ldfeke
