MODULE zdfini
   !!======================================================================
   !!                      ***  MODULE  zdfini  ***
   !! Ocean physics :   read vertical mixing namelist and check consistancy
   !!======================================================================
   !! History :  8.0  ! 1997-06  (G. Madec)  Original code from inimix
   !!            1.0  ! 2002-08  (G. Madec)  F90 : free form
   !!             -   ! 2005-06  (C. Ethe) KPP scheme
   !!             -   ! 2009-07  (G. Madec) add avmb, avtb in restart for cen2 advection
   !!            3.7  ! 2014-12  (G. Madec) remove KPP scheme
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   zdf_init    : initialization, namelist read, and parameters control
   !!----------------------------------------------------------------------
   USE par_oce         ! mesh and scale factors
   USE zdf_oce         ! TKE vertical mixing          
   USE sbc_oce         ! surface module (only for nn_isf in the option compatibility test)
   USE zdftke          ! TKE vertical mixing
   USE zdfgls          ! GLS vertical mixing
   USE zdfric          ! Richardson vertical mixing   
   USE zdfddm          ! double diffusion mixing      
   USE zdfevd          ! enhanced vertical diffusion  
   USE tranpc          ! convection: non penetrative adjustment
   USE ldfslp          ! iso-neutral slopes
   !
   USE in_out_manager  ! I/O manager
   USE iom             ! IOM library
   USE lib_mpp         ! distribued memory computing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   zdf_init   ! routine called by opa.F90
   
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: zdfini.F90 7646 2017-02-06 09:25:03Z timgraham $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE zdf_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_init  ***
      !! 
      !! ** Purpose :   initializations of the vertical ocean physics
      !!
      !! ** Method  :   Read namelist namzdf, control logicals 
      !!----------------------------------------------------------------------
      INTEGER          ::   ioptio, ios       ! local integers
      REAL(wp)         ::   rn_avt_amp        ! JM hack: rn_avt0 amplification factor
      !!
      NAMELIST/namzdf/ rn_avm0, rn_avt0, nn_avb, nn_havtb, rn_avt_amp,            &
         &        ln_zdfexp, nn_zdfexp, ln_zdfevd, nn_evdm, rn_avevd, ln_zdfnpc,  &
         &        nn_npc, nn_npcp, ln_zdfqiao
      !!----------------------------------------------------------------------

      REWIND( numnam_ref )              ! Namelist namzdf in reference namelist : Vertical mixing parameters
      READ  ( numnam_ref, namzdf, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namzdf in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namzdf in reference namelist : Vertical mixing parameters
      READ  ( numnam_cfg, namzdf, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namzdf in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namzdf )

      IF(lwp) THEN               !* Parameter print
         WRITE(numout,*)
         WRITE(numout,*) 'zdf_init : vertical physics'
         WRITE(numout,*) '~~~~~~~~'
         WRITE(numout,*) '   Namelist namzdf : set vertical mixing mixing parameters'
         WRITE(numout,*) '      vertical eddy viscosity             rn_avm0   = ', rn_avm0
         WRITE(numout,*) '      vertical eddy diffusivity           rn_avt0   = ', rn_avt0
         WRITE(numout,*) '      constant background or profile      nn_avb    = ', nn_avb
         WRITE(numout,*) '      horizontal variation for avtb       nn_havtb  = ', nn_havtb
         WRITE(numout,*) '      time splitting / backward scheme    ln_zdfexp = ', ln_zdfexp
         WRITE(numout,*) '      number of time step                 nn_zdfexp = ', nn_zdfexp
         WRITE(numout,*) '      enhanced vertical diffusion         ln_zdfevd = ', ln_zdfevd
         WRITE(numout,*) '         applied on momentum (=1/0)       nn_evdm   = ', nn_evdm
         WRITE(numout,*) '      vertical coefficient for evd        rn_avevd  = ', rn_avevd
         WRITE(numout,*) '      non-penetrative convection (npc)    ln_zdfnpc = ', ln_zdfnpc
         WRITE(numout,*) '      npc call  frequency                 nn_npc    = ', nn_npc
         WRITE(numout,*) '      npc print frequency                 nn_npcp   = ', nn_npcp
         WRITE(numout,*) '      Qiao formulation flag               ln_zdfqiao=', ln_zdfqiao
      ENDIF

      !                          !* Parameter & logical controls
      !                          !  ----------------------------
      !
      !                               ! ... check of vertical mixing scheme on tracers
      !                                              ==> will be done in trazdf module
      !
      !                               ! ... check of mixing coefficient
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '   vertical mixing option :'
      ioptio = 0
      IF( lk_zdfcst ) THEN
         IF(lwp) WRITE(numout,*) '      constant eddy diffusion coefficients'
         ioptio = ioptio+1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! JM 04 Feb 19: 
!!! Hack for UNAGI, amplify the tracer diffusivity over a particular region
!!!    the value to be used is 1e-5, to be amplifed to 5e-3 -> factor 500
!!!    step.F90 is modified so that avt is multiplied by avtb_2d
         IF( nn_havtb == 1 ) THEN
            IF(lwp) WRITE(numout,*) '      JM hack: imposing a horizontal structure on constant vertical tracer diffusion'
            IF(lwp) WRITE(numout,*) '               amplfication of rn_avt0 = ', rn_avt0, ' over sponge region by factor of rn_avt_amp = ', rn_avt_amp
         END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ENDIF
      IF( lk_zdfric ) THEN
         IF(lwp) WRITE(numout,*) '      Richardson dependent eddy coefficients'
         ioptio = ioptio+1
      ENDIF
      IF( lk_zdftke ) THEN
         IF(lwp) WRITE(numout,*) '      TKE dependent eddy coefficients'
         ioptio = ioptio+1
      ENDIF
      IF( lk_zdfgls ) THEN
         IF(lwp) WRITE(numout,*) '      GLS dependent eddy coefficients'
         ioptio = ioptio+1
      ENDIF
      IF( ioptio == 0 .OR. ioptio > 1 )   &
         &   CALL ctl_stop( ' one and only one vertical diffusion option has to be defined ' )
      IF( ( lk_zdfric .OR. lk_zdfgls ) .AND. ln_isfcav )   &
         &   CALL ctl_stop( ' only zdfcst and zdftke were tested with ice shelves cavities ' )
      !
      !                               ! ... Convection
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '   convection :'
      !
#if defined key_top
      IF( ln_zdfnpc )   CALL ctl_stop( ' zdf_init: npc scheme is not working with key_top' )
#endif
      !
      ioptio = 0
      IF( ln_zdfnpc ) THEN
         IF(lwp) WRITE(numout,*) '      use non penetrative convective scheme'
         ioptio = ioptio+1
      ENDIF
      IF( ln_zdfevd ) THEN
         IF(lwp) WRITE(numout,*) '      use enhanced vertical dif. scheme'
         ioptio = ioptio+1
      ENDIF
      IF( lk_zdftke ) THEN
         IF(lwp) WRITE(numout,*) '      use the 1.5 turbulent closure'
      ENDIF
      IF( lk_zdfgls ) THEN
         IF(lwp) WRITE(numout,*) '      use the GLS closure scheme'
      ENDIF
      IF ( ioptio > 1 )   CALL ctl_stop( ' chose between ln_zdfnpc and ln_zdfevd' )
      IF( ioptio == 0 .AND. .NOT.( lk_zdftke .OR. lk_zdfgls ) )           &
         CALL ctl_stop( ' except for TKE or GLS physics, a convection scheme is',   &
         &              ' required: ln_zdfevd or ln_zdfnpc logicals' )

      !                               !* Background eddy viscosity and diffusivity profil
      IF( nn_avb == 0 ) THEN                ! Define avmb, avtb from namelist parameter
         avmb(:) = rn_avm0
         avtb(:) = rn_avt0                     
      ELSE                                  ! Background profile of avt (fit a theoretical/observational profile (Krauss 1990)
         avmb(:) = rn_avm0
         avtb(:) = rn_avt0 + ( 3.e-4_wp - 2._wp * rn_avt0 ) * 1.e-4_wp * gdepw_1d(:)   ! m2/s
         IF(ln_sco .AND. lwp)   CALL ctl_warn( 'avtb profile not valid in sco' )
      ENDIF
      !
      IF( ln_rstart ) THEN                  !  Read avmb, avtb in restart (if exist)
         ! if ln_traadv_cen, avmb, avtb have been modified in traadv_cen2 module. 
         ! To ensure the restartability, avmb & avtb are written in the restart 
         ! file in traadv_cen2 end read here. 
         IF( iom_varid( numror, 'avmb', ldstop = .FALSE. ) > 0 ) THEN
            CALL iom_get( numror, jpdom_unknown, 'avmb', avmb )
            CALL iom_get( numror, jpdom_unknown, 'avtb', avtb )
         ENDIF
      ENDIF
      !                                     ! 2D shape of the avtb
      avtb_2d(:,:) = 1.e0                        ! uniform 
      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! JM 01 Mar 19: 
!!! Hack for UNAGI, amplify the tracer diffusivity over a particular region of width 300km
!!!    rn_avt0 = 1e-5, to be amplifed by some factor taken from input
!!!    maximum at the 2400km which is beyond the last wet point 
!!!    end of sponge region to be 2400 - 300 = 2100, so choose a rest gphit value just less than that
!!!    step.F90 is modified so that avt is multiplied by avtb_2d
!!!    
      IF( nn_havtb == 1 ) THEN
         avtb_2d(:,:) = 1. + 0.5 * rn_avt_amp * (    1. + COS(  rpi * ( gphit(:,:) - 2400.) / 300.  )   )   &
            &              - 0.5              * (    1. + COS(  rpi * ( gphit(:,:) - 2400.) / 300.  )   )
         WHERE(  gphit <= 2099. )   avtb_2d = 1.e0
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      IF( nn_havtb == 1 ) THEN                   ! decrease avtb in the equatorial band
!           !  -15S -5S : linear decrease from avt0 to avt0/10.
!           !  -5S  +5N : cst value avt0/10.
!           !   5N  15N : linear increase from avt0/10, to avt0
!           WHERE(-15. <= gphit .AND. gphit < -5 )   avtb_2d = (1.  - 0.09 * (gphit + 15.))
!           WHERE( -5. <= gphit .AND. gphit <  5 )   avtb_2d =  0.1
!           WHERE(  5. <= gphit .AND. gphit < 15 )   avtb_2d = (0.1 + 0.09 * (gphit -  5.))
!      ENDIF
      !
   END SUBROUTINE zdf_init

   !!======================================================================
END MODULE zdfini
