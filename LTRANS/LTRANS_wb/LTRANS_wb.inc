! ******************************* LTRANS Include File *******************************

!*** NUMBER OF PARTICLES ***
INTEGER, PARAMETER :: numpar = 10 !13738          ! Number of particles

!*** TIME PARAMETERS ***
!REAL,    PARAMETER :: days   = 5.625        ! Number of days to run the model
REAL,    PARAMETER :: days   = 1.0        ! Number of days to run the model

!INTEGER, PARAMETER :: iprint = 3600         ! Print interval for LTRANS output (s); 3600 = every hour
INTEGER, PARAMETER :: iprint = 21600         ! Print interval for LTRANS output (s); 3600 = every hour

INTEGER, PARAMETER :: dt  = 3600            ! External time step (duration between hydro model predictions) (s) 
INTEGER, PARAMETER :: idt = 120             ! Internal (particle tracking) time step (s)
DOUBLE PRECISION, PARAMETER :: Delay = 0.0  ! Time (s) to delay particle release 

!*** ROMS HYDRODYNAMIC MODULE PARAMETERS ***
INTEGER, PARAMETER :: us = 12               ! Number of Rho grid s-levels in ROMS hydro model
INTEGER, PARAMETER :: ws = 13               ! Number of W grid s-levels in ROMS hydro model
INTEGER, PARAMETER :: tdim = 1              ! Number of time steps per ROMS hydro predictions file
REAL,    PARAMETER :: hc  = 2.5             ! Min Depth - used in ROMS S-level transformations
DOUBLE PRECISION, PARAMETER:: z0 = 0.0005   ! ROMS roughness parameter
DOUBLE PRECISION, PARAMETER :: ConstantHTurb = 1.0  ! Constant value of horizontal turbulence (m2/s)

!*** TURBULENCE MODULE PARAMETERS ***
LOGICAL, PARAMETER :: HTurbOn = .TRUE.      ! Horizontal Turbulence on (.TRUE.) or off (.FALSE.)
LOGICAL, PARAMETER :: VTurbOn = .TRUE.      ! Vertical Turbulence on (.TRUE.) or off (.FALSE.)
INTEGER, PARAMETER :: p2 = ws * 4	        ! Number of proliferated points in vertical turbulence module 
     !Note: Only change p2 with a great deal of caution.

!*** BEHAVIOR MODULE PARAMETERS ***
INTEGER, PARAMETER :: Behavior = 0          ! Behavior type (specify a number)
! Note: The behavior types numbers are: 0 Passive, 1 near-surface, 2 near-bottom, 
!       3 DVM, 4 C.virginica oyster larvae, 5 C.ariakensis oyster larvae, 6 constant)

DOUBLE PRECISION, PARAMETER :: deadage = 1.8*24.*3600.  ! Age at which a particle stops moving (i.e., dies) (s)
!  Note: deadage can be used to stop particle motion for all behavior types (0-6)

DOUBLE PRECISION, PARAMETER :: pediage = 1.2*24.*3600.   ! Age when particle reaches max swim speed and can settle (s)
!  Note: for oyster larvae behavior types (4 & 5), pediage = age at which a particle becomes a pediveliger
!  Note: pediage does not cause particles to settle if the Settlement module is not on
DOUBLE PRECISION, PARAMETER :: swimstart = 0.5*24.*3600. ! Age that swimming or sinking begins (s) 1 day = 1.*24.*3600.
DOUBLE PRECISION, PARAMETER :: swimslow = 0.00025        ! Swimming speed when particle begins to swim (m/s)
DOUBLE PRECISION, PARAMETER :: swimfast = 0.003          ! Maximum swimming speed (m/s)  0.05 m/s for 5 mm/s
!  Note: for constant swimming speed for behavior types 1,2 & 3, set swimslow = swimfast = constant speed

DOUBLE PRECISION, PARAMETER :: Sgradient = 1.0		! Salinity gradient threshold that cues larval behavior (psu/m)
!  Note: This parameter is only used if Behavior = 4 or 5. 

DOUBLE PRECISION, PARAMETER :: constant = -0.0003        ! Sinking velocity for behavior type 6
!  Note: This parameter is only used if Behavior = 6.

!* DVM. The following are parameters for the Diurnal Vertical Migration (DVM) behavior type:
DOUBLE PRECISION, PARAMETER :: twistart = 4.801821   ! Time of twilight start (hr) **
DOUBLE PRECISION, PARAMETER :: twiend = 19.19956     ! Time of twilight end (hr) **
DOUBLE PRECISION, PARAMETER :: daylength = 14.39774  ! Length of day (hr) **
DOUBLE PRECISION, PARAMETER :: Em = 1814.328         ! Irradiance at solar noon (microE m^-2 s^-1) **
DOUBLE PRECISION, PARAMETER :: Kd = 1.07             ! Vertical attenuation coefficient
DOUBLE PRECISION, PARAMETER :: thresh = 0.0166       ! Light threshold that cues behavior (microE m^-2 s^-1)
!  Note: These values were calculated for September 1 at the latitude of 37.0 (Chesapeake Bay mouth)
!  Note: Variables marked with ** were calculated with light_v2BlueCrab.f (not included in LTRANS yet)
!  Note: These parameters are only used if Behavior = 3 

!*** SETTLEMENT MODULE PARAMETERS ***
LOGICAL, PARAMETER :: settlementon = .FALSE. ! settlement module on (.TRUE.) or off (.FALSE.)
!  Note: If settlement is off: set minholeid, maxholeid, minpolyid, maxpolyid, pedges, & hedges to 1
!        to avoid both wasted variable space and errors due to arrays of size 0.
!        If settlement is on and there are no holes: set minholeid, maxholeid, & hedges to 1
LOGICAL, PARAMETER :: holesExist = .TRUE.   ! Are there holes in habitat? yes(TRUE) no(FALSE)
INTEGER, PARAMETER :: minpolyid = 101001    ! Lowest habitat polygon id number
INTEGER, PARAMETER :: maxpolyid = 101004    ! Highest habitat polygon id number
INTEGER, PARAMETER :: minholeid = 100201    ! Lowest hole id number
INTEGER, PARAMETER :: maxholeid = 100401    ! Highest hole id number
INTEGER, PARAMETER :: pedges=76           ! Number of habitat polygon edge points (# of rows in habitat polygon file)
INTEGER, PARAMETER :: hedges=33           ! Number of hole edge points (number of rows in holes file)

!*** CONVERSION MODULE PARAMETERS ***
DOUBLE PRECISION, PARAMETER :: PI = 3.14159265358979     ! Pi
DOUBLE PRECISION, PARAMETER :: RCF = 180.0 / PI          ! Radian conversion factor
DOUBLE PRECISION, PARAMETER :: Earth_Radius = 6378*1000  ! Equatorial radius

!*** INPUT FILE NAME AND LOCATION PARAMETERS ***; 
!ROMS NetCDF Model Grid file
CHARACTER(LEN=46), PARAMETER :: NCgridfile = '/home/lederer/LTRANS/LTRANS_Code_wb/GRID_wb.nc'    !Filename 
  !Note: the path to the file is necessary if the file is not in the same folder as the code
  !Note: if .nc file in separate folder in Linux, then include path. For example:
  !      CHARACTER(LEN=29), PARAMETER :: NCgridfile = '/share/enorth/CPB_GRID_wUV.nc' 
  !Note: if .nc file in separate folder in Windows, then include path. For example:
  !      CHARACTER(LEN=23), PARAMETER :: NCgridfile = 'D:\ROMS\CPB_GRID_wUV.nc'

!ROMS Predictions NetCDF Input File. Filename = prefix + filenum + suffix
CHARACTER(LEN=48), PARAMETER :: prefix='/home/lederer/LTRANS/wbroms-6-200/OUT/ocean_his_'    ! NetCDF Input Filename prefix                     
INTEGER, PARAMETER :: filenum = 50                ! Number in First NetCDF Input Filename
CHARACTER(LEN=3), PARAMETER :: suffix='.nc'        ! NetCDF Input Filename suffix
  !Note: the path to the file is necessary if the file is not in the same folder as the code
  !Note: if .nc file in separate folder in Windows, then include path in prefix. For example:
  !      CHARACTER(LEN=15), PARAMETER :: prefix='D:\ROMS\y95hdr_'   
  !      if .nc file in separate folder in Linux, then include path in prefix. For example:
  !      CHARACTER(LEN=26), PARAMETER :: prefix='/share/lzhong/1995/y95hdr_'   

!Particle Location Input File 
CHARACTER(LEN=25), PARAMETER :: parfile     = 'initial_part_location.csv' !Particle locations 
  !Note: the path to the file is necessary if the file is not in the same folder as the code

!Habitat Polygon Location Input Files
CHARACTER(LEN=24), PARAMETER :: habitatfile = 'sample_habitat_edges.csv'  !Habitat polygons
CHARACTER(LEN=24), PARAMETER :: holefile    = 'sample_habitat_holes.csv'  !Holes in habitat
  !Note: the path to the file is necessary if the file is not in the same folder as the code

!*** OTHER PARAMETERS *** 
INTEGER, PARAMETER :: seed = 9               ! Seed value for random number generator (Mersenne Twister)
LOGICAL, PARAMETER :: BoundaryBLNs = .FALSE. ! Create Surfer Blanking Files of boundaries? .TRUE.=yes, .FALSE.=no
LOGICAL, PARAMETER :: SaltTempOn = .FALSE.   ! Calculate salinity and temperature at particle 
                                             ! location: yes (.TRUE.) or no (.FALSE.)