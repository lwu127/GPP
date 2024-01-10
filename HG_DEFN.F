
!------------------------------------------------------------------------!
!  The Community Multiscale Air Quality (CMAQ) system software is in     !
!  continuous development by various groups and is based on information  !
!  from these groups: Federal Government employees, contractors working  !
!  within a United States Government contract, and non-Federal sources   !
!  including research institutions.  These groups give the Government    !
!  permission to use, prepare derivative works of, and distribute copies !
!  of their work in the CMAQ system to the public and to permit others   !
!  to do so.  The United States Environmental Protection Agency          !
!  therefore grants similar permission to use the CMAQ system software,  !
!  but users are requested to provide copies of derivative works or      !
!  products designed to operate in the CMAQ system to the United States  !
!  Government without restrictions as to use by others.  Software        !
!  that is used with the CMAQ system but distributed under the GNU       !
!  General Public License or the GNU Lesser General Public License is    !
!  subject to their copyright restrictions.                              !
!------------------------------------------------------------------------!

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE GPP( CGRID, JDATE, JTIME, TSTEP )
C  hgii gas species gas-particle partitioning.

C  Contains:
C     Subroutine extract_soa
C     Subroutine update_orgvapor
C     Function findVapor
C     Function findOrgprod
C     Subroutine orgaer
C-----------------------------------------------------------------------
      USE GRID_CONF, ONLY: NCOLS, NROWS, NLAYS, IO_PE_INCLUSIVE
      USE RXNS_DATA, ONLY: MECHNAME
      USE AERO_DATA, ONLY: COAG_BUDGET, COND_BUDGET, GROWTH_BUDGET, NPF_BUDGET,
     &                     FIXED_sg, AEROMODE_DIAM, AEROMODE_LNSG, AEROMODE_DENS,
     &                     MOMENT0_CONC, MOMENT2_CONC, MOMENT3_CONC,
     &                     AEROSPC, AEROSPC_CONC, N_MODE, AH2O_IDX,
     &                     EXTRACT_AERO, UPDATE_AERO, CALCMOMENTS
      USE PRECURSOR_DATA, ONLY: PRECURSOR_CONC, PHG_RATE,
     &                         EXTRACT_PRECURSOR, UPDATE_PRECURSOR
      USE AEROMET_DATA, ONLY: AIRTEMP, AIRPRES, AIRQV, AIRDENS, AIRRH,
     &                        H2OSATVP, H2OVP, MWWAT, MWAIR 
      USE UTILIO_DEFN, ONLY: OPEN3, INTERPX, WRITE3, FSREAD3, ALLVAR3,
     &                       XSTAT1, TIME2SEC, SEC2TIME
#ifdef twoway
     &                       , INDEX1, XSTAT3 
#endif
      USE CGRID_SPCS, ONLY: NSPCSD
#ifdef twoway
     &                      , N_GC_CONC, GC_CONC, GC_STRT, GC_CONC_MAP
#endif
      USE PMDIAG_DATA_MODULE, ONLY: NUM_PMDIAG_SPC, A_PMDIAG_SPC_MAP, PMDIAG_SPC_RECORD
      USE RUNTIME_VARS, ONLY: LOGDEV, PMDIAG, APMDIAG, APMDIAG_ZHI, END_TIME
#ifdef twoway
      USE twoway_data_module
#endif
      USE CENTRALIZED_IO_MODULE, ONLY: interpolate_var
#ifdef sens
      USE DDM3D_DEFN, ONLY : NP, NPMAX, SENGRID
      Use aero_ddm3d, ONLY : s_so4rate, s_precursor_conc, s_aerospc_conc, s_phg_rate
#endif 
      USE SHG


      IMPLICIT NONE 

C *** Includes:
  
      INCLUDE SUBST_FILES_ID  ! file name parameters

C *** Arguments:

C *** CGRID is conc field (including gas and aerosol variables)
      REAL, POINTER :: CGRID( :,:,:,: )      ! concentrations
      INTEGER      JDATE        ! Current model date, coded YYYYDD
      INTEGER      JTIME        ! Current model time, coded HHMMSS
      INTEGER      TSTEP( 3 )   ! Time step vecto (HHMMSS)
                                ! TSTEP(1) = local output step
                                ! TSTEP(2) = sciproc sync. step (chem)
                                ! TSTEP(3) = twoway model time step w.r.t wrf time 
                                !            step and wrf/cmaq call frequency

C *** Local Variables:

      CHARACTER( 16 ), SAVE :: PNAME = 'HG_DEFN'
      CHARACTER( 16 ) :: VNAME         ! variable name
      CHARACTER( 96 ) :: XMSG = ' ' 

      INTEGER  MDATE, MTIME, MSTEP     ! julian date, time and 
                                       ! timestep in sec 
      INTEGER  C, R, L, V, N           ! loop counters 
      INTEGER  SPC                     ! species loop counter
      INTEGER  STRT, FINI              ! loop induction variables
      INTEGER  ALLOCSTAT               ! memory allocation status

      LOGICAL  LERROR                  ! Error flag

C *** Variables to set time step for writing visibility file
      INTEGER, SAVE :: WSTEP = 0 ! local write counter
      LOGICAL, SAVE :: WRITETIME = .FALSE. ! logical write flag

C *** Meteorological variables 
      REAL PRES  (NCOLS,NROWS,NLAYS)   ! Atmospheric pressure [ Pa ]
      REAL TA    (NCOLS,NROWS,NLAYS)   ! Air temperature [ K ]
      REAL DENS  (NCOLS,NROWS,NLAYS)   ! Air density [ kg/m**-3 ]
      REAL QV    (NCOLS,NROWS,NLAYS)   ! Water vapor mixing ratio [ kg/kg ]

c *** synchronization time step [ s ]
      REAL DT

C *** Statement Function ***************
      REAL ESATL ! arithmetic statement function for vapor pressure [Pa]
      REAL TT

C *** Coefficients for the equation, ESATL defining saturation vapor pressure
      REAL, PARAMETER :: AL = 610.94
      REAL, PARAMETER :: BL = 17.625
      REAL, PARAMETER :: CL = 243.04

      REAL, ALLOCATABLE, SAVE :: TEMP_DATA(:)     
 
C *** values of AL, BL, CL are from :
C     Alduchov and Eskridge, "Improved mMagnus Form Approximations of 
C     Saturation Vapor Pressure,"
C     Jour. of Applied Meteorology, col. 35,
C     pp 601-609, April, 1996

      ESATL( TT ) = AL * EXP( BL * ( TT - 273.15 ) / ( TT - 273.15 + CL ) )

C *** End Statement Function ************

      logical, save :: now = .true.

C *** first pass flag
      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      LOGICAL, SAVE :: FIRST_CTM_PMDIAG_1 = .TRUE.

C *** ratio of molecular weights of water vapor to dry air = 0.622015
      REAL, PARAMETER :: EPSWATER = MWWAT / MWAIR

C-------------------------Begin body of HG_DEFN-----------------------------
      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
        
         ALLOCATE( TEMP_DATA( NUM_PMDIAG_SPC ),
     &             STAT=ALLOCSTAT)

      END IF

      MDATE = JDATE
      MTIME = JTIME
      MSTEP = TIME2SEC( TSTEP( 2 ) )
      CALL NEXTIME ( MDATE, MTIME, SEC2TIME( MSTEP / 2 ) )

      WSTEP = WSTEP + TIME2SEC( TSTEP( 2 ) )
      IF ( WSTEP .GE. TIME2SEC( TSTEP( 1 ) ) ) WRITETIME = .TRUE.

C *** Set floating point synchronization time step:
      DT = FLOAT( MSTEP )  ! set time step in seconds

C *** Get Meteorological Variables

C *** pressure [Pa]
      CALL interpolate_var ('PRES', mdate, mtime, PRES)

C *** temperature [K]
      CALL interpolate_var('TA', mdate, mtime, TA)

C *** specific humidity [g H2O/g air]
      CALL interpolate_var('QV', mdate, mtime, QV)

C *** air density [kg/m3]
      CALL interpolate_var('DENS', mdate, mtime, DENS) 

C--------------------- Begin loops over grid cells ---------------------

C *** initialize conc arrays

      DO L = 1, NLAYS
         DO R = 1, NROWS
            DO C = 1, NCOLS

C *** Grid cell meteorological data.
               AIRTEMP  = TA  ( C,R,L )
               AIRPRES  = PRES( C,R,L )
               AIRQV    = QV  ( C,R,L )
               AIRDENS  = DENS( C,R,L )
               H2OSATVP = ESATL( AIRTEMP )
               H2OVP    = AIRPRES * AIRQV / ( EPSWATER + AIRQV ) 
               AIRRH    = MAX( 0.005, MIN( 0.99, H2OVP / H2OSATVP ) )
C *** Extract grid cell concentrations of aero species from CGRID
C     into aerospc_conc in aero_data module (set minimum)
#ifdef sens
               CALL EXTRACT_AERO( CGRID( C,R,L,: ), .TRUE., SENGRID( C,R,L,:,: ), .TRUE. )
#else
               CALL EXTRACT_AERO( CGRID( C,R,L,: ), .TRUE. )
#endif

C *** Extract hgii concentrations from CGRID 
               CALL EXTRACT_HG( CGRID( C,R,L,: ) )
C        tempcorr = 10 - 2500.0/airtemp    ! log(1/k)     
C        kcorr = 10**tempcorr           ! 1/k, (GOM*PM2.5)/PBM
C        f = kcorr/pm25       ! GOM/PBM
C        fconv = 1.0/f       ! PBM/GOM
C        fgas = 1.0/(fconv+1.0)    ! GOM/(GOM+PBM) 
C        faer = 1.0/(f+1.0)  ! PBM/(PBM+GOM)

C *** Call phgaer to calculate gas particle partitioning of hgii species
               CALL HGAER( DT, L ) 

C *** Update gas hgii concentrations back to CGRID
C               CALL UPDATE_PHG( CGRID( C,R,L,: ) )
C *** Update gas hgii concentrations back to CGRID
               CALL UPDATE_GHG( CGRID( C,R,L,: ) )
C *** Update aerosol variables conc back into CGRID (set minimum) 
C     and convert M2 to dry and save as surface area
#ifdef sens
               CALL UPDATE_AERO( CGRID( C,R,L,: ), .TRUE., SENGRID(C,R,L,:,: ) )
#else
               CALL UPDATE_AERO( CGRID( C,R,L,: ), .TRUE. )
#endif
            END DO
         END DO
      END DO

      RETURN
      END

