
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
      Module shg

C  Defines aerosol species arrays and parameters required in SOA processing.

C  Contains:
C     Subroutine extract_soa
C     Subroutine update_orgvapor
C     Function findVapor
C     Function findOrgprod
C     Subroutine orgaer
C-----------------------------------------------------------------------
      Implicit None

      ! Define Logical values as T and F for the OA table
      Logical, Parameter, Private :: T = .true.
      Logical, Parameter, Private :: F = .false.

      Integer, Parameter :: n_phg_list = 22   ! # of potential partitioning SVOCs
      Integer, Save      :: n_phg             ! # of actual simulated SOA vapors

      Type phg_type
         Character( 16 ) :: pm_name     ! Organic Aerosol Species Name
C         Character(len=:),allocatable:: pm_name
         Character( 16 ) :: gas_name    ! Condensable Vapor Species Name
         Real            :: mw          ! Molar mass of species
         Real            :: alpha       ! Mass-based stoichiometric coefficients [ug/m^3]/[ug/m^3]
         Real            :: cstar       ! Effective saturation concentrations [ug/m^3] at 298 K
         Real            :: enth        ! Enthalphy of Vaporization [J/mol]
         Logical         :: nonvol      ! Flag identifying nonvolatile species
      End Type phg_type

      Type( phg_type ), Allocatable, Save :: phgspc( : )

      Type( phg_type ), Save :: phg_list( n_phg_list ) = (/

      !         PM          Vapor       Molar       Alpha   CStar^    Enth-   Nonvol-     
      !         Name        Name        Mass                          alpy    atile        
      !        ----------  ----------  ----------   ------  --------  -------  -----    
!     & phg_type('APHG'    , 'HG'     ,  200.6,    1.000,   40.807,   59.1E3,   F ), 
     & phg_type('AHGCL2'  , 'HGCL2'  ,  271.6,    1.000,   40.807,   81.0E3,   F ), 
     & phg_type('AHGOOH'  , 'HGOOH'  ,  233.6,    1.000,   40.807,   81.0E3,   F ),
     & phg_type('AHGOH2'  , 'HGOH2'  ,  234.6,    1.000,   40.807,   81.0E3,   F ),
     & phg_type('AOHHGONO', 'OHHGONO',  263.6,    1.000,   40.807,   81.0E3,   F ),
     & phg_type('AOHHGOOH', 'OHHGOOH',  250.6,    1.000,   40.807,   81.0E3,   F ),
     & phg_type('AHGOHBRO', 'HGOHBRO',  313.5,    1.000,   40.807,   81.0E3,   F ),
     & phg_type('AHGOHCLO', 'HGOHCLO',  269.1,    1.000,   40.807,   81.0E3,   F ),
     & phg_type('ACLHGCLO', 'CLHGCLO',  287.6,    1.000,   40.807,   81.0E3,   F ),
     & phg_type('ABRHGCLO', 'BRHGCLO',  332.0,    1.000,   40.807,   81.0E3,   F ),
     & phg_type('ACLHGBRO', 'CLHGBRO',  332.0,    1.000,   40.807,   81.0E3,   F ),
     & phg_type('ACLHGBR' , 'CLHGBR' ,  316.0,    1.000,   40.807,   81.0E3,   F ),
     & phg_type('ABRHGBRO', 'BRHGBRO',  376.4,    1.000,   40.807,   81.0E3,   F ),
     & phg_type('AHGBR2'  , 'HGBR2'  ,  360.4,    1.000,   40.807,   81.0E3,   F ),
     & phg_type('AHGBRO'  , 'HGBRO'  ,  296.5,    1.000,   40.807,   81.0E3,   F ),
     & phg_type('AHGCLO'  , 'HGCLO'  ,  252.1,    1.000,   40.807,   81.0E3,   F ),
     & phg_type('ACLHGONO', 'CLHGONO',  282.1,    1.000,   40.807,   81.0E3,   F ),
     & phg_type('ABRHGONO', 'BRHGONO',  326.5,    1.000,   40.807,   81.0E3,   F ),
     & phg_type('ACLHGOOH', 'CLHGOOH',  269.1,    1.000,   40.807,   81.0E3,   F ),
     & phg_type('ABRHGOOH', 'BRHGOOH',  313.5,    1.000,   40.807,   81.0E3,   F ),
     & phg_type('ACLHGOH' , 'CLHGOH' ,  253.1,    1.000,   40.807,   81.0E3,   F ),
     & phg_type('ABRHGOH' , 'BRHGOH' ,  297.5,    1.000,   40.807,   81.0E3,   F ),
     & phg_type('ABRHGOI' , 'BRHGOI' ,  423.4,    1.000,   40.807,   81.0E3,   F ) /)
      ! Table Comments:      ! *Reaction Counters are only needed if the vapor species is not formed
      !  directly in the gas-phase chemistry mechanism.
      ! ^Nonvolatile species are given a default C* of 1e-10 ug/m3 and
      !  enthalpy of vaporization equal to 1 kJ/mol.

      ! Variables for saving properties and concentrations of organic
      ! compounds
      !   Molecular Weights
       Real, Allocatable    :: phg_mw   ( : )   ! PM species
       Real, Allocatable    :: ghg_mw   ( : )   ! Vapor Species
                                                   ! that are linked to yields
      !   Concentrations
       Real, Allocatable    :: ghg_conc( : )   ! Vapors
       Real, Allocatable    :: phg_conc( : )   ! Aerosols 
       Real, Allocatable    :: ghg_chg( : )   ! Aerosols 
       Real, Allocatable    :: phgi_conc( : )   ! Aerosols 
       Real, Allocatable    :: phgj_conc( : )   ! Aerosols 
       Character(16), Allocatable    :: pmvar( : )   ! Aerosols 
       Real            :: srfa

      ! Variables for mapping OA to aerosol and CGRID arrays
      Integer, Allocatable :: map_PHGtoAERO( : )    ! mercury aerosol pointers to aerospc
      Integer, Allocatable :: map_GHGtoCGRID( : )   ! pointers of vapor species to CGRID
      Integer, Allocatable :: map_PHGItoCGRID( : )  ! pointers of PBM species to CGRID
      Integer, Allocatable :: map_PHGJtoCGRID( : )  ! pointers of PBM species to CGRID
      Integer, Allocatable :: map_PMtoCGRID( : )  ! pointers of PBM species to CGRID
      
      ! Variables for controlling behavior of OA module
      Logical, Private, Save :: mapped              = .False.
      Logical, Private, Save :: RXNS_eflag          = .False. ! error flag for RXNS modules used
      Logical, Private, Save :: phg_eflag            = .False. ! error flag for soa_defn
      Character( 16 ), Private, Save :: pname = 'SHG        '

      Contains

C-----------------------------------------------------------------------
      Subroutine extract_hg( conc )

C  Extracts the required hg data from CGRID into the conc array,gas is in the
C  unit of ppmv, particulate matter is in the unit of ug/m3 and is converted to 
C  ppmv.
C  Gas speices index is mapGHGtoCGRID(spc), particulate spcecies index is 
C  mapPHGtoCGRID(spc). Output arrays are hgii_conc(spc)/phg_mw(spc) and 
C  phg_conc(spc)/phg_mw(spc).

C------------------------------------------------------------------------
      Use rxns_data, only: mechname
      Use aero_data, only: AE_eflag, findAero, aerospc, aerolist, n_aerolist, aerospc_mw, n_aerospc, aerospc_conc
      Use aeromet_data, only: airdens, inv_mwair, min_gasconc
      Use cgrid_spcs, only: n_gc_spc, gc_spc, n_gc_g2ae, gc_g2ae, gc_strt, gc_molwt, gc_g2ae_map, 
     &                      n_ae_spc, ae_spc, ae_strt, ae_molwt
      Use runtime_vars, only: logdev
      Use utilio_defn, only: index1, xstat3

      Implicit None

      ! Arguments:
      Real, Intent( In ) :: conc( : )

      ! Local Variables:
      Character( 300 ):: xmsg
      Real            :: gasconv
      Real            :: ppm2ug 
c      Real            :: conc_ec
      Integer         :: n, a, g, b, n1, n2,n3
      Integer         :: spc, idx_eci, idx_ecj
      Integer         :: map_PHGtoPHGlist( 200 )
      Integer         :: i,idx_i,idx_j, index_n1, index_n2, index_n3, n_eci, n_ecj
      Character(16) :: varnm1, varnm2


      ! Map All OA Species and Load Properties
      If ( .Not. mapped ) Then
         mapped = .TRUE.

         ! Allocate the OA property and mapping variables now that the
         ! number of active OA species has been determined.
         Allocate ( phgspc ( n_phg_list ))
         Allocate ( phg_mw ( n_phg_list ))
         Allocate ( ghg_mw ( n_phg_list ))
         Allocate ( ghg_conc( n_phg_list ))
         Allocate ( map_PHGtoAERO( n_phg_list ))
         Allocate ( map_GHGtoCGRID( n_phg_list ))
         Allocate ( map_PMtoCGRID( n_aerolist ))
         Allocate ( pmvar( n_aerolist ))

         gasconv = airdens * inv_mwair ! Compute gas conversion constant       
         phg_mw = 200.0
         Do spc = 1, n_phg_list
            phgspc( spc ) = phg_list(spc)
            phg_mw( spc ) = phgspc( spc )%mw 
C            Write(logdev, '(2X, I3,E14.7)') spc,phg_mw(spc)
         End Do
         map_PHGtoAERO = 0
 
         Do spc = 1, n_phg_list
            map_PHGtoAERO( spc ) = findAero( phgspc( spc )%pm_name, .False. )
         End Do

         ! Determine the location on CGRID of condensable gases in 
         ! equilibrium with OA particle species. When a match is found, 
         ! populate the properties of those gas species.
         map_GHGtoCGRID = 0
         ghg_mw = 200.0 ! Default Molecular Weight for all Organic Vapor Species
         Do spc = 1,n_phg_list
           If( phgspc( spc )%gas_name == '' ) Cycle
       
           g = index1( phgspc( spc )%gas_name, n_gc_spc, gc_spc )
           If ( g .ge. 1 ) Then
               map_GHGtoCGRID( spc ) = gc_strt - 1 + g
               ghg_mw( spc ) = gc_molwt( g )
           End If
c           Write(logdev, '(2X, I3,I3)') spc,map_GHGtoCGRID(spc)
         End Do

      End if

      ! Copy grid cell concentrations of vapor species
      ghg_conc   = 0.0    ! product hgii gas concentrations

      Do spc = 1, n_phg_list                            
         n = map_GHGtoCGRID( spc )
         If ( n .Ne. 0 ) Then
            ghg_conc( spc ) = Max( conc( n ), min_gasconc )  !Unit:ppmv
         End If
C         Write(logdev, '(2X, I3,E14.7)') spc,ghg_conc(spc)
         
      End Do
      ! Extract accumulation and atk and coarse mode aerosol surface area (m2/m3)
      n1 = index1( "SRFATKN", n_ae_spc, ae_spc ) 
      n2 = index1( "SRFACC", n_ae_spc, ae_spc )
      n3 = index1( "SRFCOR", n_ae_spc, ae_spc )
      index_n1 = ae_strt - 1 + n1
      index_n2 = ae_strt - 1 + n2
      index_n3 = ae_strt - 1 + n3
      srfa = conc( index_n1 ) + conc( index_n2 ) + conc( index_n3 )
C      Write(logdev, '(3X, I3,I3,E14.7)') n1,n2,srfa

C      n_eci = index1( "AECI", n_ae_spc, ae_spc )
C      n_ecj = index1( "AECJ", n_ae_spc, ae_spc )
C      idx_eci = ae_strt - 1 + n_eci
C      idx_ecj = ae_strt - 1 + n_ecj
C      conc_ec = conc( idx_eci ) + conc ( idx_ecj )  !  EC PM2.5 in the unit of ug/m3
C      conc_ec1 = conc_ec*1E-6/12.0  ! EC PM2.5 concentrations in the unit of mol/m3               

C      Write(logdev, '(2X, E14.7, E14.7)') conc_ec, conc_ec1

C 99903 Format(1X,'species name is: ', a /)
99904 Format(1X,'aerosol product from namelist species: ', a /)


      Return
      End Subroutine extract_hg

C-----------------------------------------------------------------------
      Subroutine update_ghg( conc )

C  Populates CGRID from the conc array with udated hgii gas concentrations.

C-----------------------------------------------------------------------

      Use aeromet_data
      Use utilio_defn
      Use runtime_vars, only: logdev
      Use aero_data
      Use cgrid_spcs, only: n_gc_spc, gc_spc, n_gc_g2ae, gc_g2ae, gc_strt, gc_molwt, gc_g2ae_map, 
     &                      n_ae_spc, ae_spc, ae_strt, ae_molwt

      Implicit None

      ! Arguments:
      Real, Intent( Out ) :: conc( : )

      ! Local Variables:
      Character( 80 ) :: xmsg
      Integer         :: m, n, hgcl2_idx
      Integer         :: spc
      Real            :: ppm2ug 
      Real            :: gasconv
      Real            :: hgcl2_chg

      If ( .Not. mapped ) Then
         xmsg = 'CGRID Species has not been mapped'
         Call m3exit( pname, 0, 0, xmsg, xstat3 )
      End If

c      m = index1("HGCL2", n_gc_spc, gc_spc)
c      hgcl2_idx = gc_strt - 1 + m 
      gasconv = rgasuniv*airtemp/airpres ! Compute gas conversion constant       

      ! Copy hg gas conc back to grid cell concentrations
      Do spc = 1, n_phg_list
         n = map_GHGtoCGRID( spc )
         If ( n .Ne. 0 ) Then
            ppm2ug = gasconv/phg_mw( spc )
            conc( n ) = Max ( ghg_conc( spc )*ppm2ug, min_gasconc )
c            Write(logdev, '(2X, E14.7,E14.7)') ghg_conc(spc),conc(n)
         End If
      End Do
c      hgcl2_chg = 0
c      Do spc = 1, n_phg_list
c        ppm2ug = gasconv/phg_mw(spc)
c        hgcl2_chg = hgcl2_chg + ghg_chg( spc )*ppm2ug
C        Write( logdev, '(2X, I3, E14.7)') spc, hgcl2_chg
c      End Do 
c      conc( hgcl2_idx ) = ghg_conc( 1 )*ppm2ug + hgcl2_chg
C      Write( logdev, '(3X, E14.7, E14.7, E14.7)') ghg_conc( 1 ), hgcl2_chg, conc( hgcl2_idx)

      Return
      End Subroutine update_ghg

C-----------------------------------------------------------------------
      Subroutine update_phg( conc )

C  Populates CGRID from the conc array with updated hgii aerosol concentrations.NOT USED

C-----------------------------------------------------------------------

      Use aeromet_data
      Use utilio_defn, only: xstat3
      Use runtime_vars, only: logdev
      Use aero_data

      Implicit None

      ! Arguments:
      Real, Intent( Out ) :: conc( : )

      ! Local Variables:
      Character( 80 ) :: xmsg
      Integer         :: n, m
      Integer         :: spc
      Real            :: ppm2ug
      Real            :: gasconv

      If ( .Not. mapped ) Then
         xmsg = 'CGRID Species has not been mapped'
         Call m3exit( pname, 0, 0, xmsg, xstat3 )
      End If

      gasconv = airpres/(rgasuniv*airtemp) ! Compute gas conversion constant       
      ! Copy hg pariculate matter conc back to grid cell concentrations
      Do spc = 1, n_phg_list
         n = map_PHGJtoCGRID( spc )
         m = map_PHGItoCGRID( spc )
         ppm2ug = gasconv*phg_mw( spc )           
         If ( n .Ne. 0 ) Then
           conc( n ) = Max ( phgj_conc( spc )*ppm2ug, min_gasconc )
         End If
         If ( m .Ne. 0 ) Then
           conc( m ) = Max ( phgi_conc( spc )*ppm2ug, min_gasconc )
         End If
C           Write(logdev, '(2X, I3,E14.7)') spc,phg_conc(spc)
      End Do

      Return
      End Subroutine update_phg

C-----------------------------------------------------------------------
      Subroutine hgaer( dt, layer )

C Updates CGRID via several pathways for secondary organic aerosol (SOA)
C formation, as recommended by Edney et al. (2007).  These include SOA
C formation from isoprene, monoterpenes, sesquiterpenes, long alkanes, and
C aromatics (incl. benzene).

C-----------------------------------------------------------------------

C Key Subroutines/Functions called:  newt, soabisection

      Use aero_data
      Use aeromet_data
      Use utilio_defn
      Use rxns_data, only : MECHNAME

      Implicit None

      ! Arguments:
      Real    :: dt            ! Synchronization time step [ s ]
      Integer :: layer         ! model layer number

      ! Local variables:
      Logical, Save :: first_time = .True.
      Character( 300 ) :: xmsg
      Integer       :: i, im, indx, ispc, n, nsvol

      Real, Allocatable, Save ::  mw_inv( : )     ! Inverse MW of SVOCs [ mol/g ]

      ! Parameters & variables for adjusting cstar to ambient conditions
      Real, Parameter :: tref   = 298.0          ! reference temperature [ K ]
      Real, Parameter :: rgasconst = 8.206E-5
      Real, Parameter :: trefm1 = 1.0 / tref     ! inverse of reference temperature
      Real, Parameter :: rgas1  = 1.0 / rgasuniv ! reciprocal of universal gas constant
      Real, Parameter :: p298 = 3.06E-5  ! 1/RT at 1 atm and 298 K [ mole/m**3 ]
      Real convfac, tt1, tt2, pm25,pm 
      Real tempcorr      ! temperature correction factor for cstar
      Real, Parameter :: diffhg = 9.76E-6  ! Diffusivity of hgii [m2 s-1]
      Real, Parameter :: alphhg = 1.0      ! accomodation coefficient
!      Real, Parameter :: ratio = 3.44E-8   ! test2: 1.72E-6     ! test 1: 3.44E-7
!      Real, Parameter :: ratio = 1.0   ! test2: 1.72E-6     ! test 1: 3.44E-7
!      Real, Parameter :: ratio = 0.5   ! test2: 1.72E-6     ! test 1: 3.44E-7
!      Real, Parameter :: ratio = 3.0   ! April
!      Real, Parameter :: ratio = 4.0   ! May
!      Real, Parameter :: ratio = 2.0   ! June
      Real, Parameter :: ratio = 4.0   ! July

      ! Variables used in equilibrium partitioning calculations
      Real(8) GRtmp( n_mode )   ! Dummy variable for accurate treatment of growth to specific moment
      Real GR3( n_phg_list,n_mode ) ! 3rd Moment Growth for each mode and compound
      Real GR3FRAC( n_phg_list,n_mode ) ! Fraction of each mode growing/shrinking
      Real cbar_hg(n_phg_list) ! On-line molecular speed of each organic
      Real dv_hg           ! On-line gas-phase diffusivity of each organic
      Real totdry           ! dry aerosols 
      Real kp           ! dry aerosols 
      Real fconv           ! dry aerosols 
      Real f, ns, pt, p0, lnp, atsp, conc_oa, c_oa, fom, ptc, kpd, kpb      ! dry aerosols 
      Real faer              ! fraction of total in aerosol, intermediate value
      Real fgas              ! fraction of total in aerosol, intermediate value
      Real            :: ppm2ug 
      Real            :: gasconv
      Real, Allocatable, Save          :: rhg ( : ) ! reactive mercury conc before current time step [ ug/m3 ]
      Real, Allocatable, Save          :: phgi_c0( : ) ! i mode particulate mercury before current time step [ ug/m3 ]
      Real, Allocatable, Save          :: phgj_c0( : ) ! j mode particulate mercury before current time step [ ug/m3 ]
      Real, Allocatable, Save          :: ghg_c0( : ) ! j mode particulate mercury before current time step [ ug/m3 ]
      Real, Allocatable, Save          :: phg_c0( : ) ! j mode particulate mercury before current time step [ ug/m3 ]
      Real, Allocatable, Save          :: phg_conc( : ) ! j mode particulate mercury before current time step [ ug/m3 ]
      Real, Allocatable, Save          :: daer( : )  ! Total (g+p) conc before current time step [ ug/m3 ]
      Real, Allocatable, Save          :: dgas( : )  ! Total (g+p) conc before current time step [ ug/m3 ]
      Real, Allocatable, Save          :: totaer ( : ) ! Total (g+p) conc before current time step [ ug/m3 ]
      Real, Allocatable, Save          :: Phi( : )
      Real, Allocatable, Save          :: conc_om( : )

      ! Variables for computing the budget
      REAL :: CBUDGET0_NUM ( N_MODE )
      REAL :: CBUDGET0_SRF ( N_MODE )
      REAL :: CBUDGET0_MASS( N_AEROSPC,N_MODE )


      ! Variables for updating 2nd and 3rd moments
      Real(8) m0_init( 2 )  ! initial 0 moment, wet [ mom/m3 ]
      Real(8) m1wet_init( 2 )  ! initial 1st moment, wet [ mom/m3 ]
      Real(8) m2wet_initD( 2 )  ! initial 2nd moment, wet [ mom/m3 ]
      Real m3wet_init( 2 )  ! initial 3rd moment, wet [ mom/m3 ]
      Real m2wet_init( 2 )  ! initial 2nd moment, wet [ mom/m3 ]
      Real m3wet_final( 2 )  ! final 3rd moment with updated PHG [ mom/m3 ]
      Real m2wet_final( 2 )  ! final 2nd moment with updated PHG [ mom/m3 ]


C-----------------------------------------------------------------------

      If ( first_time )  Then
         first_time = .False.

         ! Set unit conversion and inverse mw constants. Allocate space
         ! for semivolatile species partitioning variables.
         nsvol = count( .not.phgspc( : )%nonvol )

         Allocate( rhg( nsvol ), daer( nsvol ), phg_c0( nsvol ),
     &     totaer(nsvol), phg_conc( nsvol), Phi( nsvol ),
     &     ghg_c0(nsvol), ghg_chg(nsvol), dgas(nsvol),conc_om(3) )

      End If ! first_time

      ! Equilibrium Partitioning Calculations

      ! Set temperature factors
      tt2 = 1.0 / airtemp - trefm1

      ! moles in the system
      gasconv = rgasuniv*airtemp/airpres ! Compute gas conversion constant       
      pm25 = sum( aerospc_conc( :,1 )+aerospc_conc(:,2)+
     &            aerospc_conc(:,3), mask=aerospc(:)%Dry )   ! pm25 concentration in the unit of ug/m3
!      pm = sum( aerospc_conc( :,1 )+aerospc_conc(:,2),
!     &           mask=aerospc(:)%Dry )   
      conc_oa = sum( aerospc_conc( :,1 )+aerospc_conc(:,2)+
     &               aerospc_conc(:,3), mask=aerospc(:)%om )   ! oa concentration in the unit of ug/m3
      c_oa = 1.0E-6*conc_oa/200.0 ! oa concentration in the unit of mol/m3  
      conc_om(1) = sum( aerospc_conc( :,1 ), mask=aerospc(:)%om)
      conc_om(2) = sum( aerospc_conc( :,2 ), mask=aerospc(:)%om)
      conc_om(3) = sum( aerospc_conc( :,3 ), mask=aerospc(:)%om)
      atsp = 10000*srfa/pm25   ! specific surface area for the particulate matter (cm2/ug)
!      Write(logdev, '(2X, E14.7)') pm25

      fom = conc_oa/pm25     ! om mass fraction ug/ug     
      lnp = 81041*tt2*rgas1    ! ln(P(298)/p)
      p0 = exp(lnp)   ! p(298)/p(t)
      pt = p298/p0    ! vapor pressure at temperature t, unit: torr 
      ptc = 16.0*pt
      ns = c_oa/(10000*srfa)  ! surface concentration of sorption sites ( mol/cm2),srfa in the unit of m2/m3      
      kpd = ns*atsp*airtemp/ptc  
c      kp = kpb + kpd   ! PBM/(GOM*PM2.5)
c      f = kp*pm25    ! PBM/GOM
c      fconv = 1.0/f  ! GOM/PBM
c      fgas = 1.0/(fconv+1.0)
c      faer = 1.0/(f+1.0)
C      Write(logdev, '(4X, E14.7, E14.7, E14.7, E14.7)') pm25, atsp, kp, airtemp 
      ! Compound
      Do i = 1, n_phg_list
        n = map_PHGtoAERO( i )
        kpb = 760*rgasconst*fom*airtemp*1E-6/(phg_mw( i )*pt)   ! rgas1=8.206E-5
c        kpb = 760*rgasuniv*fom*airtemp*1E-6/phg_mw( i )
C        Write(logdev, '(1X, E14.7)') kpb 
        kp = 10*(kpb + kpd)   ! adjust value: 30
        f = kp*pm25    ! PBM/GOM
c        fconv = 1.0/f  ! GOM/PBM
        fgas = 1.0/(f+1.0)
        faer = f/(f+1.0)
        ! Sum particle-phase and total organic mass for this
        ! semivolatile species
        ppm2ug = gasconv/phg_mw( i )
        phg_c0( i ) = (aerospc_conc( n,2 ) + aerospc_conc( n,1 ) +
     &                 aerospc_conc( n,3 )) ! ug/m3 
        ghg_c0( i ) = ghg_conc( i )/ppm2ug   ! unit ug/m3 
        rhg( i ) = ghg_c0( i ) + phg_c0( i )   ! ug/m3
C        Write(logdev, '(2X, I3, E14.7)') i, ctoti(i)
        phg_conc( i ) = rhg( i )*faer  ! unit:ug/m3
        ghg_conc( i ) = rhg( i )*fgas
c        daer( i ) = phg_c0( i ) - phg_conc( i )    ! in the unit of ppmv
C        ghg_conc( i ) = ghg_conc( i ) + daer( i )    ! in the unit of ppmv
c        Write(logdev, '(2X, I3, E14.7)') i,phg_conc(i)
        ! Constrain caer to values between conmin and ctotf

        ! Check that the mass of aerosol left has not gone below a
        ! miminum value
        If ( phg_conc( i ) .Lt. tiny(0.0) ) Then
           Write( logdev,* ) 'ghg less than zero for ghg species: ', 
     &                        phgspc( i )%pm_name, ' ', phg_conc( i ),
     &                        ' reset to evapmin.'
           rhg( i ) = rhg( i ) - phg_conc( i ) + evapmin
           phg_conc( i ) = evapmin
        End If

        ! Check that the mass of aerosol condensed does not exceed
        ! the total
        If ( phg_conc( i ) - rhg( i ) .Gt. tiny(0.0) ) Then
            Write( logdev,* ) 'phg exceeds rhg at i = ', i
            Write( logdev,* ) 'phg: ', phg_conc( i ), 'rhg: ', rhg( i ) 
            Write( logdev,* ) 'phg reset to rhg'
            phg_conc( i ) = rhg( i )
        End If

        ! Calculate total change in hg aerosol concentrations
        daer( i ) = phg_conc( i ) - phg_c0( i )    ! in the unit of ppm
        dgas( i ) = ghg_conc( i ) - ghg_c0( i )    ! in the unit of ppm

        ! Transfer the semivolatile mass
c        totaer( i ) = phg_c0( i )
        ! Sum the total aerosol mass of this compound across all
        ! aerosol modes
        indx = map_PHGtoAERO( i )
c        totaer( i ) = max( sum( aerospc_conc( indx,: ) , mask = 
c     &              .NOT.aero_missing( indx,: ) ), conmin )
     
C        Write(logdev, '(2X, I3,E14.7)') indx, totaer( i ) 
        ! Transfer the semivolatile mass
        Do im = 1,3
           If ( aero_missing( indx,im ) ) Cycle
c           If ( daer( i ) .lt. 0.0 ) Then
              ! Evaporate using mode-dependent mass fraction
              Phi( i ) = conc_om( im ) / conc_oa
              aerospc_conc( indx,im ) = 
     &             aerospc_conc( indx,im ) + daer( i ) * Phi( i )
c           Else
              ! Condense using mode-dependent condensaiton flux 
              ! fraction
c              aerospc_conc( indx,im ) = aerospc_conc( indx,im ) 
c     &             + daer( i ) * GR3FRAC( i,im )
C              Write(logdev, '(2X, E14.7, E14.7)') GR3FRAC(i,1), GR3FRAC(i,2) 
c           End If
        End Do
C        If ( dgas( i ) .lt. 0.0 ) Then
        ghg_conc( i ) = rhg( i ) - phg_conc( i )
C           ghg_chg( i ) = min_gasconc
C        Else
C           ghg_conc( i ) = ghg_c0( i )
C           ghg_chg( i ) = dgas( i ) 
C        End If
C        Write(logdev, '(4X, I3,E14.7,E14.7,E14.7)') i, dgas( i ), ghg_conc( i ), ghg_chg( i )  
C        Write(logdev, '(4X, E14.7,E14.7,E14.7,E14.7)') airtemp,ghg_conc(i),phg_conc(i),ug2ppm
C        Write(logdev, '(4X, E14.7,E14.7,E14.7,E14.7)') airtemp,tempcorr,kcorr,f  
      End Do

C      Write(logdev, '1X,E14.7)') aerospc_conc(73,2) 
      Return
      End Subroutine hgaer
C-----------------------------------------------------------------------
      End Module shg
