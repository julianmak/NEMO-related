MODULE traldf
   !!======================================================================
   !!                       ***  MODULE  traldf  ***
   !! Ocean Active tracers : lateral diffusive trends 
   !!=====================================================================
   !! History :  9.0  ! 2005-11  (G. Madec)  Original code
   !!  NEMO      3.0  ! 2008-01  (C. Ethe, G. Madec)  merge TRC-TRA 
   !!            3.7  ! 2013-12  (G. Madec) remove the optional computation from T & S anomaly profiles and traldf_bilapg
   !!             -   ! 2013-12  (F. Lemarie, G. Madec)  triad operator (Griffies) + Method of Stabilizing Correction
   !!             -   ! 2014-01  (G. Madec, S. Masson)  restructuration/simplification of lateral diffusive operators
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_ldf       : update the tracer trend with the lateral diffusion trend
   !!   tra_ldf_init  : initialization, namelist read, and parameters control
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE ldftra         ! lateral diffusion: eddy diffusivity & EIV coeff.
   USE ldfslp         ! lateral diffusion: iso-neutral slope
   USE traldf_lap_blp ! lateral diffusion: laplacian iso-level            operator  (tra_ldf_lap/_blp   routines)
   USE traldf_iso     ! lateral diffusion: laplacian iso-neutral standard operator  (tra_ldf_iso        routine )
   USE traldf_triad   ! lateral diffusion: laplacian iso-neutral triad    operator  (tra_ldf_triad      routine )
   USE trd_oce        ! trends: ocean variables
   USE trdtra         ! ocean active tracers trends
   !
   USE prtctl         ! Print control
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distribued memory computing library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE wrk_nemo       ! Memory allocation
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_ldf        ! called by step.F90 
   PUBLIC   tra_ldf_init   ! called by nemogcm.F90 
   !
   INTEGER ::   nldf = 0   ! type of lateral diffusion used defined from ln_traldf_... (namlist logicals)
   
   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.7 , NEMO Consortium (2015)
   !! $Id: traldf.F90 7765 2017-03-07 12:45:51Z timgraham $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_ldf( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf  ***
      !! 
      !! ** Purpose :   compute the lateral ocean tracer physics.
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !!
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  ztrdt, ztrds
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )   CALL timing_start('tra_ldf')
      !
      IF( l_trdtra )   THEN                    !* Save ta and sa trends
         CALL wrk_alloc( jpi,jpj,jpk,   ztrdt, ztrds ) 
         ztrdt(:,:,:) = tsa(:,:,:,jp_tem) 
         ztrds(:,:,:) = tsa(:,:,:,jp_sal)
      ENDIF
      !
      SELECT CASE ( nldf )                     !* compute lateral mixing trend and add it to the general trend
      CASE ( np_lap   )                                  ! laplacian: iso-level operator
         CALL tra_ldf_lap  ( kt, nit000,'TRA', ahtu, ahtv, gtsu, gtsv, gtui, gtvi, tsb,      tsa, jpts,  1   )
      CASE ( np_lap_i )                                  ! laplacian: standard iso-neutral operator (Madec)
         CALL tra_ldf_iso  ( kt, nit000,'TRA', ahtu, ahtv, gtsu, gtsv, gtui, gtvi, tsb, tsb, tsa, jpts,  1   )
      CASE ( np_lap_it )                                 ! laplacian: triad iso-neutral operator (griffies)
         CALL tra_ldf_triad( kt, nit000,'TRA', ahtu, ahtv, gtsu, gtsv, gtui, gtvi, tsb, tsb, tsa, jpts,  1   )
      CASE ( np_blp , np_blp_i , np_blp_it )             ! bilaplacian: iso-level & iso-neutral operators
         CALL tra_ldf_blp  ( kt, nit000,'TRA', ahtu, ahtv, gtsu, gtsv, gtui, gtvi, tsb      , tsa, jpts, nldf )
      END SELECT
      !
      IF( l_trdtra )   THEN                    !* save the horizontal diffusive trends for further diagnostics
         ztrdt(:,:,:) = tsa(:,:,:,jp_tem) - ztrdt(:,:,:)
         ztrds(:,:,:) = tsa(:,:,:,jp_sal) - ztrds(:,:,:)
         CALL trd_tra( kt, 'TRA', jp_tem, jptra_ldf, ztrdt )
         CALL trd_tra( kt, 'TRA', jp_sal, jptra_ldf, ztrds )
         CALL wrk_dealloc( jpi,jpj,jpk,   ztrdt, ztrds ) 
      ENDIF
      !                                        !* print mean trends (used for debugging)
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' ldf  - Ta: ', mask1=tmask,               &
         &                       tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
      !
      IF( nn_timing == 1 )   CALL timing_stop('tra_ldf')
      !
   END SUBROUTINE tra_ldf


   SUBROUTINE tra_ldf_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf_init  ***
      !! 
      !! ** Purpose :   Choice of the operator for the lateral tracer diffusion
      !!
      !! ** Method  :   set nldf from the namtra_ldf logicals
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio, ierr   ! temporary integers 
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN                     ! Namelist print
         WRITE(numout,*)
         WRITE(numout,*) 'tra_ldf_init : lateral tracer diffusive operator'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namtra_ldf: already read in ldftra module'
         WRITE(numout,*) '      see ldf_tra_init report for lateral mixing parameters'
      ENDIF
      !                                   ! use of lateral operator or not
      nldf   = np_ERROR
      ioptio = 0
      IF( ln_traldf_lap )   ioptio = ioptio + 1
      IF( ln_traldf_blp )   ioptio = ioptio + 1
      IF( ioptio >  1   )   CALL ctl_stop( 'tra_ldf_init: use ONE or NONE of the 2 lap/bilap operator type on tracer' )
      IF( ioptio == 0   )   nldf = np_no_ldf     ! No lateral diffusion
      !
      IF( nldf /= np_no_ldf ) THEN        ! direction ==>> type of operator  
         ioptio = 0
         IF( ln_traldf_lev )   ioptio = ioptio + 1
         IF( ln_traldf_hor )   ioptio = ioptio + 1
         IF( ln_traldf_iso )   ioptio = ioptio + 1
         IF( ioptio >  1 )   CALL ctl_stop( 'tra_ldf_init: use only ONE direction (level/hor/iso)' )
         IF( ioptio == 0 )   CALL ctl_stop( 'tra_ldf_init: use one of (level/hor/iso)' )
         !
         !                                ! defined the type of lateral diffusion from ln_traldf_... logicals
         ierr = 0
         IF( ln_traldf_lap ) THEN         ! laplacian operator
            IF ( ln_zco ) THEN               ! z-coordinate
               IF ( ln_traldf_lev   )   nldf = np_lap     ! iso-level = horizontal (no rotation)
               IF ( ln_traldf_hor   )   nldf = np_lap     ! iso-level = horizontal (no rotation)
               IF ( ln_traldf_iso   )   nldf = np_lap_i   ! iso-neutral: standard  (   rotation)
               IF ( ln_traldf_triad )   nldf = np_lap_it  ! iso-neutral: triad     (   rotation)
            ENDIF
            IF ( ln_zps ) THEN               ! z-coordinate with partial step
               IF ( ln_traldf_lev   )   ierr = 1          ! iso-level not allowed 
               IF ( ln_traldf_hor   )   nldf = np_lap     ! horizontal             (no rotation)
               IF ( ln_traldf_iso   )   nldf = np_lap_i   ! iso-neutral: standard     (rotation)
               IF ( ln_traldf_triad )   nldf = np_lap_it  ! iso-neutral: triad        (rotation)
            ENDIF
            IF ( ln_sco ) THEN               ! s-coordinate
               IF ( ln_traldf_lev   )   nldf = np_lap     ! iso-level              (no rotation)
               IF ( ln_traldf_hor   )   nldf = np_lap_i   ! horizontal             (   rotation)
               IF ( ln_traldf_iso   )   nldf = np_lap_i   ! iso-neutral: standard  (   rotation)
               IF ( ln_traldf_triad )   nldf = np_lap_it  ! iso-neutral: triad     (   rotation)
            ENDIF
         ENDIF
         !
         IF( ln_traldf_blp ) THEN         ! bilaplacian operator
            IF ( ln_zco ) THEN               ! z-coordinate
               IF ( ln_traldf_lev   )   nldf = np_blp     ! iso-level = horizontal (no rotation)
               IF ( ln_traldf_hor   )   nldf = np_blp     ! iso-level = horizontal (no rotation)
               IF ( ln_traldf_iso   )   nldf = np_blp_i   ! iso-neutral: standard  (   rotation)
               IF ( ln_traldf_triad )   nldf = np_blp_it  ! iso-neutral: triad     (   rotation)
            ENDIF
            IF ( ln_zps ) THEN               ! z-coordinate with partial step
               IF ( ln_traldf_lev   )   ierr = 1          ! iso-level not allowed 
               IF ( ln_traldf_hor   )   nldf = np_blp     ! horizontal             (no rotation)
               IF ( ln_traldf_iso   )   nldf = np_blp_i   ! iso-neutral: standard  (   rotation)
               IF ( ln_traldf_triad )   nldf = np_blp_it  ! iso-neutral: triad     (   rotation)
            ENDIF
            IF ( ln_sco ) THEN               ! s-coordinate
               IF ( ln_traldf_lev   )   nldf = np_blp     ! iso-level              (no rotation)
               IF ( ln_traldf_hor   )   nldf = np_blp_it  ! horizontal             (   rotation)
               IF ( ln_traldf_iso   )   nldf = np_blp_i   ! iso-neutral: standard  (   rotation)
               IF ( ln_traldf_triad )   nldf = np_blp_it  ! iso-neutral: triad     (   rotation)
            ENDIF
         ENDIF
         IF( ierr == 1 )   CALL ctl_stop( 'iso-level in z-partial step, not allowed' )
      ENDIF
      !
!!! jm 11 Apr 19: switched off the force stop to have option of GM + no Redi + bilaplacian tracer diffusion
!      IF( ln_ldfeiv .AND. .NOT.( ln_traldf_iso .OR. ln_traldf_triad ) )                                    &
!           &            CALL ctl_stop( 'eddy induced velocity on tracers requires iso-neutral laplacian diffusion' )
      IF( ln_isfcav .AND. ln_traldf_triad ) &
           &            CALL ctl_stop( ' ice shelf cavity and traldf_triad not tested' )
           !
      IF(  nldf == np_lap_i .OR. nldf == np_lap_it .OR. &
         & nldf == np_blp_i .OR. nldf == np_blp_it  )   l_ldfslp = .TRUE.    ! slope of neutral surfaces required 
      !
      IF(lwp) THEN
         WRITE(numout,*)
         SELECT CASE( nldf )
         CASE( np_no_ldf )   ;   WRITE(numout,*) '      ===>>   NO lateral diffusion'
         CASE( np_lap    )   ;   WRITE(numout,*) '      ===>>   laplacian iso-level operator'
         CASE( np_lap_i  )   ;   WRITE(numout,*) '      ===>>   Rotated laplacian operator (standard)'
         CASE( np_lap_it )   ;   WRITE(numout,*) '      ===>>   Rotated laplacian operator (triad)'
         CASE( np_blp    )   ;   WRITE(numout,*) '      ===>>   bilaplacian iso-level operator'
         CASE( np_blp_i  )   ;   WRITE(numout,*) '      ===>>   Rotated bilaplacian operator (standard)'
         CASE( np_blp_it )   ;   WRITE(numout,*) '      ===>>   Rotated bilaplacian operator (triad)'
         END SELECT
      ENDIF
      !
   END SUBROUTINE tra_ldf_init

   !!======================================================================
END MODULE traldf
