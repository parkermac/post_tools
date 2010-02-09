MODULE BEHAVIOR_MOD

! The behavior module is used to assign biological or physical characteristics to particles. 
! Currently particle movement is in the vertical direction. 
!
! Particle characteristics can include a swimming/sinking speed component and 
! a behavioral cue component that can depend upon particle age. The swimming/sinking speed 
! component controls the speed of particle motion and can be constant or set with a function. 
! The behavioral cue component regulates the direction of particle movement. For biological 
! behaviors, a random component is added to the swimming speed and direction to simulate random 
! variation in the movements of individuals (in behavior types 1 - 5, see list below). Physical 
! characteristics can also be assigned to particles, like constant sinking velocity, without 
! the additional random movements (behavior type 6). The following behavior types are currently 
! available in LTRANS and are specified using the Behavior parameter in the LTRANS.inc file:
!
!
! Passive (no behavior): Behavior = 0. In this case, the behavior module is not executed. 
!     Particle motion is based on advection, and, if turned on, horizontal and vertical turbulence.
!
! Near-surface orientation: Behavior = 1. Particles swim up if they are deeper than 1 m from the surface.  
!
! Near-bottom orientation: Behavior = 2. Particles swim down if they are shallower than 1 m from the bottom.  
!
! Diurnal vertical migration: Behavior = 3. Particles swim down if light levels at the particle location 
!     exceed a predefined threshold value.  
!
! Crassostrea virginica oyster larvae: Behavior = 4. Swimming speeds and direction of motion vary depending 
!     upon age (stage) according to field and laboratory observations (see North et al. 2008). 
!
! C. ariakensis oyster larvae: Behavior = 5. Swimming speeds and direction of motion vary depending upon age 
!     (stage) according to field and laboratory observations (see North et al. 2008).
!
! Sinking velocity: Behavior = 6. Particles move up or down with constant sinking (or floating) speeds 
!     without individual random motion. Code that calculates salinity and temperature at the particle 
!     location is included (but commented out) as a basis for calculating density-dependent sinking velocities.  
!
!
! Behavior algorithms and code created by: Elizabeth North
! Module structure created by:             Zachary Schlag
! Created on:			                   2004
! Last Modified on:	                       31 Aug 2008
!


  USE PARAM_MOD, ONLY: numpar
  IMPLICIT NONE
  PRIVATE
  SAVE

  REAL :: timer(numpar)                    !Timer for C. ariakensis downward swimming behavior
  INTEGER :: P_behave(numpar),           & !Behavior of each particle
             status(numpar)                !Status of each particle (behavior 1, behavior 2, settled 3, or dead 4) 
  DOUBLE PRECISION :: P_pediage(numpar), & !Age at which the particle will settle (become a pediveliger)
                      P_deadage(numpar), & !Age at which the particle will stop moving (die)
                      P_Sprev(numpar),   & !Salinity at particle's previous location (for calculating salt gradient)
                      P_zprev(numpar),   & !Particle's previous depth (for calculating salt gradient)
                      P_swim(numpar,3)     !Swimming speed (age-dependent, linear increase unless constant)   
                                           !(n,1)slope, (n,2)intercept, (n,3) speed at current age

  !The following procedures have been made public for the use of other program units:
  PUBLIC :: initBehave,updateStatus,Behave,getColor   

CONTAINS

  SUBROUTINE initBehave()    !Initialize the behavior module
    USE PARAM_MOD, ONLY: Behavior,swimfast,swimslow,swimstart,pediage,deadage,Sgradient,settlementon 
    USE SETTLEMENT_MOD, ONLY: initSettlement
	USE NORM_MOD,   ONLY: NORM
    IMPLICIT NONE
	INTEGER :: n

    write(*,*) 'initialize behavior'    

    do n=1,numpar
      !Set behavior to the one specified in LTRANS.inc
      P_behave(n) = Behavior  !Behavior( 0 Passive, 1 near-surface, 2 near-bottom,  
                              ! 3 DVM, 4 C.virginica, 5 C.ariakensis, 6 constant )
      P_pediage(n) = pediage  !age at which particle reaches maximum swimming
	                          !speed and can settle (becomes a pediveliger) (s)
      P_deadage(n) = deadage  !age at which particle stops moving (i.e., dies) (s)
      !Note: the following code assigns different veliger and pediveliger stage duration
	  !P_pediage(n) = (14. + norm()*0.5)*24.*3600.
	  !P_deadage(n) = P_pediage(n) + (7. + norm()*0.5)*24.*3600.

      !Calculate slope and intercept for age-dependent linear swimming speed
      P_swim(n,1) = (swimfast - swimslow)/(P_pediage(n) - swimstart) !slope
      P_swim(n,2) = swimfast - P_swim(n,1)*P_pediage(n)              !intercept
      P_swim(n,3) = 0.0                                              !swimming speed (m/s) 
      !Note: P_swim(n,3) is updated at each time step in Subroutine Behave
    enddo

    !The following variables are used by the C. virginica and C. ariakensis behavior routines
    timer = 0.0         !to count how long C. arikensis particles swim down
    status = 1          !status of the particle (1 - First Behavior (veliger), 
	                    ! 2 - Second Behavior (pediveliger), 3 - Settled, 4 - Dead)
    ! Initialize salt storage matrices 
    P_Sprev = 0.0		!Initialized to 0.0
    P_zprev = 0.0		!Initialized to 0.0

    !if Settlement is turned on then inform Settlement module of the age at which particle 
	!can settle (i.e., become pediveligers)
    if(settlementon)then
      CALL initSettlement(P_pediage)
    endif

  END SUBROUTINE initBehave


  SUBROUTINE updateStatus(P_age,n)  !Update particle status                                 
    USE SETTLEMENT_MOD, ONLY: SETTLED,DIE
    IMPLICIT NONE

	INTEGER, INTENT(IN) :: n
	DOUBLE PRECISION, INTENT(IN) :: P_age

    !If particle settled, update status    
    if(SETTLED(n)) status(n) = 3

    !Determine if particle dies, then update settle code and status    
    if (P_age .GE. P_deadage(n) .AND. .NOT. SETTLED(n)) then
      CALL DIE(n)              !sets settle(n) = 2                              
      status(n) = 4
    endif

  END SUBROUTINE updateStatus


  SUBROUTINE Behave(Xpar,Ypar,Zpar,Pwc_zb,Pwc_zc,Pwc_zf,P_zb,P_zc,P_zf,P_zetac,P_age,P_depth,n,it,ex,ix,daytime,p,Behav)
    USE PARAM_MOD, ONLY: us,dt,idt,twistart,twiend,Em,pi,daylength,Kd,thresh,Sgradient,swimfast,swimstart,constant   
    USE HYDRO_MOD, ONLY: WCTS_ITPI
	USE RANDOM_MOD, ONLY: genrand_real1
	IMPLICIT NONE

	REAL, INTENT(IN) :: daytime
	DOUBLE PRECISION, INTENT(IN) :: Xpar,Ypar,Zpar,Pwc_zb(us),Pwc_zc(us),Pwc_zf(us),P_zb,P_zc,P_zf,P_zetac,P_age,P_depth,ex(3),ix(3)
	INTEGER, INTENT(IN) :: n,it,p
	DOUBLE PRECISION, INTENT(OUT) :: Behav

	INTEGER :: btest,i,deplvl
	REAL :: negpos,dev1,devB,switch,switchslope
	DOUBLE PRECISION :: P_S,P_T,parBehav,Sslope,deltaS,deltaz
    DOUBLE PRECISION :: dtime,tst,E0,P_light


!   ************************ Update swimming speeds

	if(P_age .GE. swimstart) P_swim(n,3) = P_swim(n,1)*P_age+P_swim(n,2)
    if(P_age .GE. P_pediage(n)) P_swim(n,3) = swimfast

!   ************************ Prepare for TYPE 4 & 5 Behaviors

	!Update pediveliger behavior/stauts and timer
    IF(P_behave(n) .EQ. 4 .OR. P_behave(n) .EQ. 5) THEN

      !Set behavior code and status for pediveligers
      if (P_age .GE. P_pediage(n) .AND. P_age .LT. P_deadage(n)) then
        P_behave(n) = 2.0
        status(n) = 2
      endif

      !decrement timer
      timer(n) = max(0., timer(n)-dt) 
	ENDIF

	!obtain salinity at particle location (P_S) to cue oyster larvae behavior
	IF ((P_behave(n).EQ.4) .OR. (P_behave(n).EQ.5 .AND. timer(n).EQ.0.0)) THEN 

        do i=3,us-2
          if ((Zpar .LT. Pwc_zb(i)) .OR. (Zpar .LT. Pwc_zc(i)) .OR. (Zpar .LT. Pwc_zf(i))) exit
        enddo
        deplvl = i-2   !depth level

        !Salinity at particle location
        P_S = WCTS_ITPI("salt",Xpar,Ypar,deplvl,Pwc_zb,Pwc_zc,Pwc_zf,us,P_zb,P_zc,P_zf,ex,ix,p,4)

    ENDIF

!			*********************************************************
!			*														*
!			*					Behaviors					        *
!			*														*
!			*********************************************************

	parBehav = 0.0

    !TYPE 1. Ssurface oriented. Particle swims up if deeper than 1 m.
    IF (P_behave(n).EQ.1) THEN 
      btest = 0   !switch to control behavior

      !particle has 80% chance of swimming up if deeper than 1.0 m of bottom
      if (P_zc .LT. (P_zetac-1.0)) then
        negpos = 1.0
        dev1=genrand_real1()
        switch = 0.80 
        if (dev1.GT.switch) negpos = -1.0
        devB=genrand_real1()
        parBehav=negpos*devB*P_swim(n,3) 
        btest = 1
      end if

      !if within 1 m of surface, just swim randomly (50% chance of swimming up)
      if (btest.EQ.0) then    
        negpos = 1.0
        dev1=genrand_real1()
        switch = 0.5 
        if (dev1.GT.switch) negpos = -1.0
        devB=genrand_real1()
        parBehav=negpos*devB*P_swim(n,3)  
      end if  
   
    END IF


    !TYPE 2. Near-bottom. Particle swim down if not within 1 m of bottom.
    IF (P_behave(n).EQ.2 .OR. (P_behave(n).EQ.5 .AND. timer(n).GT.0.0)) THEN 
      btest = 0   !switch to control behavior

      !particle has 80% change of swimming down if greater than 1.0 m from bottom
      if (P_zc .GT. (P_depth+1.0)) then
        negpos = 1.0
        dev1=genrand_real1()
        switch = 0.20 
        if (dev1.GT.switch) negpos = -1.0
        devB=genrand_real1()
        parBehav=negpos*devB*P_swim(n,3)
        btest = 1
      end if

      !if within 1 m of bottom, just swim randomly
      if (btest.EQ.0) then    
        negpos = 1.0
        dev1=genrand_real1()
        switch = 0.5 
        if (dev1.GT.switch) negpos = -1.0
        devB=genrand_real1()
        parBehav=negpos*devB*P_swim(n,3)   
      end if

    END IF

    !TYPE 3: Diurnal Vertical Migration
    IF (P_behave(n).EQ.3) THEN

      !A. Find daytime in hrs since midnight (dtime)
      dtime = (daytime - aint(daytime))*24.0  !time of day 
      !This assumes that model simulations start at midnight

      !B. Calcluate irradiance at the water's surface (E0)
      tst = 0.0  !seconds since twilight start
      E0 = 0.0   !irradiance at the water's surface
      if (dtime.GT.twiStart .AND. dtime.LT.twiEnd) then
        tst=(dtime-twiStart)*3600.  
        E0= Em*SIN(PI*tst/(daylength*3600.))*SIN(PI*tst/(daylength*3600.))
      else 
        E0 = 0.0
      end if

      !C. Calcluate irradiance at depth of the particle
      P_light = E0 * exp(1*Kd*P_zc)

      !If light at particle location is less than threshold, random swimming
      if (P_light.LT.thresh ) then
        negpos = 1.0
        dev1=genrand_real1()
        switch = 0.5 
        if (dev1.GT.switch) negpos = -1.0
        devB=genrand_real1()
        parBehav=negpos*devB*P_swim(n,3)   
      end if
      !If light at particle > threshold, then have 80% chance of swimming down
      if (P_light.GT.thresh ) then
        negpos = 1.0
        dev1=genrand_real1()
        switch = 0.20 
        if (dev1.GT.switch) negpos = -1.0
        devB=genrand_real1()
        parBehav=negpos*devB*P_swim(n,3)  
      end if  	  	  

    END IF


	!TYPE 4. Crassostrea virginica -- above the halocline
	IF (P_behave(n).EQ.4) THEN 
	  if (it.EQ.1) then
		P_Sprev(n) = P_S  !for first iteration
		P_zprev(n) = P_zc
	  endif
	  btest = 0   !switch to control behavior
	  Sslope = 0.0  !salinity gradient that larvae swam through

	  !determine if larva swam through salinity gradient large enough to cue behavior.
	  !if so, then 80% chance of swimming up                                            
	  deltaS = P_Sprev(n) - P_S
	  deltaz = P_zprev(n) - P_zc
	  if (it.GT.1) Sslope = deltaS/deltaz
	  if (abs(Sslope).GT.Sgradient) then
		negpos = 1.0
		dev1=genrand_real1()
		switch = 0.80 
		if (dev1.GT.switch) negpos = -1.0
		parBehav=negpos*P_swim(n,3)
		btest = 1
	  endif

	  !if no directed swimming, just swim randomly with probabilities that result
	  !in particles moving up initially, then slowly moving toward bottom with increasing age  
	  if (btest.EQ.0) then    
		negpos = 1.0
		dev1=genrand_real1()
		if (P_age .LT. 1.5*24.*3600.) switch = 0.1
		if (P_age .GT. 1.5*24.*3600. .AND. P_age .LT. 5.*24.*3600.) switch = 0.49
		if (P_age .GT. 5.*24.*3600. .AND. P_age .LT. 8.*24.*3600.) switch = 0.50
		if (P_age .GT. 8.*24.*3600.) then 
		  switchslope =(0.50-0.517)/(8.0*24.*3600.- P_pediage(n))
		  switch =switchslope*P_age+0.50 - switchslope*8.0*24.*3600. 
		  if (P_zc .LT. P_depth+1.) switch = 0.5
		endif
		if (dev1.GT.(1-switch)) negpos = -1.0
		devB=genrand_real1()
		parBehav=negpos*devB*P_swim(n,3)  
	  endif

	  !update previous salt and depth matrix for next iteration
	  P_Sprev(n) = P_S
	  P_zprev(n) = P_zc
	ENDIF


	!TYPE 5. Crassostrea ariakensis -- below the halocline
	IF (P_behave(n).EQ.5 .AND. timer(n).EQ.0.0) THEN 
	  if (it.EQ.1) then
		P_Sprev(n) = P_S  !for first iteration
		P_zprev(n) = P_zc
	  endif
	  btest = 0   !switch to control behavior
	  Sslope = 0.0  !salinity gradient that larvae swam through

	  !determine if larva swam through salinity gradient large enough to cue behavior.   
	  !if so, then 80% chance of swimming down. Set timer that keep particle near bottom for 2 hrs
	  deltaS = P_Sprev(n) - P_S
	  deltaz = P_zprev(n) - P_zc
	  if (it.GT.1) Sslope = deltaS/deltaz
	  if (abs(Sslope).GT.Sgradient) then
		negpos = 1.0
		dev1=genrand_real1()
		switch = 0.20 
		btest = 1
		timer(n) = 2.*3600.  !2 hr times 3600 s
		if (dev1.GT.switch) negpos = -1.0
		parBehav=negpos*P_swim(n,3) 
        !keep bottom oriented behavior from starting until after particle is 3.5 days old  
		if (P_age .LT. 3.5*24.*3600.) then  
		  btest = 0
		  timer(n) = 0.
		endif         
	  endif

	  !if no directed swimming, just swim randomly with probabilities that result      
	  !in particles moving up initially, then moving toward bottom with increasing age
	  if (btest.EQ.0) then    
		negpos = 1.0
		dev1=genrand_real1()
		switch = 0.495 
		if (P_age .LT. 1.5*24.*3600.) switch = 0.9
		if (P_age .GT. 2.0*24.*3600. .AND. P_age .LT. 3.5*24.*3600.) then
		  switchslope=(0.3-0.495)/(2.0*24.*3600.- 3.5*24.*3600.)
		  switch =switchslope*P_age+0.3 - switchslope*2.0*24.*3600. 
		endif
		if (dev1.GT.switch) negpos = -1.0
		devB=genrand_real1()
		parBehav=negpos*devB*P_swim(n,3)  
	  endif

	  !update previous salt and depth matrix for next iteration        
	  P_Sprev(n) = P_S
	  P_zprev(n) = P_zc
	ENDIF

	!TYPE 6. Constant -- no random motion to vertical movement
	IF ((P_behave(n).EQ.6)) THEN        
	  
	  if(P_age .GE. swimstart) then
          parBehav = constant 
      else
	      parBehav = P_swim(n,3)
	  endif

    !Note: the code below is included if someone wants to calculate density
	!   !To calculate salinity (P_S) and temperature (P_T) at particle location 
    !    do i=3,us-2
    !      if ((Zpar .LT. Pwc_zb(i)) .OR. (Zpar .LT. Pwc_zc(i)) .OR. (Zpar .LT. Pwc_zf(i))) exit
    !    enddo
    !    deplvl = i-2   !depth level
    !
    !    !Salinity at particle location
    !    P_S = WCTS_ITPI("salt",Xpar,Ypar,deplvl,Pwc_zb,Pwc_zc,Pwc_zf,us,P_zb,P_zc,P_zf,ex,ix,p,4)
    ! 
	!	!Temperature at particle location 
    !    P_T = WCTS_ITPI("temp",Xpar,Ypar,deplvl,Pwc_zb,Pwc_zc,Pwc_zf,us,P_zb,P_zc,P_zf,ex,ix,p,4)

    ENDIF


	!Calculate movement due to behavior 
	Behav = parBehav * idt

! ******************* End Particle Behavior ******************************

  END SUBROUTINE Behave


  INTEGER FUNCTION getColor(n)
    !This function returns an identification number that describes a particle's  
	!behavior type or status for use in visualization routines. It was initially 
	!developed to contain the color code for plotting in Surfer.)                
    USE PARAM_MOD, ONLY: SETTLEMENTON
    USE SETTLEMENT_MOD, ONLY: SETTLED,DEAD
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n

    if(P_behave(n).EQ.0)getColor = 0   
    if(P_behave(n).EQ.1)getColor = 1   
    if(P_behave(n).EQ.2)getColor = 2   
    if(P_behave(n).EQ.3)getColor = 3   
	if(P_behave(n).EQ.4)getColor = 4   
    if(P_behave(n).EQ.5)getColor = 5   
    if(P_behave(n).EQ.6)getColor = 6   

    if(SETTLEMENTON)then
      if(SETTLED(n))getColor = 7         
      if(DEAD(n))getColor = 8            
    endif

  END FUNCTION getColor


END MODULE BEHAVIOR_MOD


! remove status? it is not being used to control anything here    