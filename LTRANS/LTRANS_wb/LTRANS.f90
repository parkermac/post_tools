! LTRANS - Larval TRANSport Lagrangian model v.1                                
! Date: 28 August 2008
!
! Description: The Larval TRANSport Lagrangian model (LTRANS) is an 
! off-line particle-tracking model that runs with the stored predictions 
! of a 3D hydrodynamic model, specifically the Regional Ocean Modeling System 
! (ROMS). Although LTRANS was built to simulate oyster larvae, it can  
! be adapted to simulate passive particles and other planktonic organisms. 
! LTRANS is written in Fortran 90 and is designed to track the trajectories 
! of particles in three dimensions. It includes a 4th order Runge-Kutta scheme 
! for particle advection and a random displacement model for vertical turbulent 
! particle motion. Reflective boundary conditions, larval behavior, and 
! settlement routines are also included. Components of LTRANS have been in 
! development since 2002 and are described in the following publications:
! North et al. 2004, North et al. 2006a, North et al. 2006b, and North et al. 2008.
!
! Developers:
!	Elizabeth North: enorth@hpl.umces.edu
!	Zachary Schlag: zschlag@hpl.umces.edu
!
! Mailing Address:  
!	University of Maryland
!	Center for Envir. Science
!	Horn Point Laboratory
!	Cambridge, MD 21613 USA
!
! Funding was provided by the National Science Foundation Biological Oceanography  
! Program, Maryland Department of Natural Resources, NOAA Chesapeake Bay Office, 
! NOAA Maryland Sea Grant College Program, and NOAA-funded UMCP Advanced Study 
! Institute for the Environment. 
! 
! **********************************************************************
! **********************************************************************
! **                      Copyright (c) 2008                          **
! **   The University of Maryland Center for Environmental Science    **
! **********************************************************************
! **                                                                  **
! ** This Software is open-source and licensed under the following    **
! ** conditions as stated by MIT/X License:                           **
! **                                                                  **
! **  (See http://www.opensource.org/licenses/mit-license.php ).      **
! **                                                                  **
! ** Permission is hereby granted, free of charge, to any person      **
! ** obtaining a copy of this Software and associated documentation   **
! ** files (the "Software"), to deal in the Software without          **
! ** restriction, including without limitation the rights to use,     **
! ** copy, modify, merge, publish, distribute, sublicense,            **
! ** and/or sell copies of the Software, and to permit persons        **
! ** to whom the Software is furnished to do so, subject to the       **
! ** following conditions:                                            **
! **                                                                  **
! ** The above copyright notice and this permission notice shall      **
! ** be included in all copies or substantial portions of the         **
! ** Software.                                                        **
! **                                                                  **
! ** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,  **
! ** EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE           **
! ** WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE  **
! ** AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT  **
! ** HOLDERS BE LIABLE FOR ANY CLAIMS, DAMAGES OR OTHER LIABILITIES,  **
! ** WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     **
! ** FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR    **
! ** OTHER DEALINGS IN THE SOFTWARE.                                  **
! **                                                                  **
! ** The most current official versions of this Software and          **
! ** associated tools and documentation are available at:             **
! **                                                                  **
! **  http://northweb.hpl.umces.edu/LTRANS.htm                        **
! **                                                                  **
! ** We ask that users make appropriate acknowledgement of            **
! ** The University of Maryland Center for Environmental Science,     **
! ** individual developers, participating agencies and institutions,  **
! ** and funding agencies. One way to do this is to cite one or       **
! ** more of the relevant publications listed at:                     **
! **                                                                  **
! **  http://northweb.hpl.umces.edu/LTRANS.htm#Description            **
! **                                                                  **
! **********************************************************************
! ********************************************************************** 

PROGRAM main

! LTRANS.f90 contains the main structure of the particle-tracking program. 
! It executes the external time step, internal time step, and particle loops, 
! advects particles, and writes output. It calls modules that read in 
! hydrodynamic model information, move particles due to turbulence and behavior, 
! test if particles are in habitat polygons, and apply boundary conditions to 
! keep particles in the model domain. 
!
! Program created by:   Elizabeth North
! Modified by:          Zachary Schlag
! Created on:			2004
! Last Modified on:	    3 September 2008


  USE PARAM_MOD,      ONLY: numpar,us,ws,days,dt,idt,seed,parfile,delay,HTurbOn,VTurbOn,settlementon,iprint,Behavior,SaltTempOn 

  USE BOUNDARY_MOD,   ONLY: CREATEBOUNDS,MBOUNDS,IBOUNDS,INTERSECT_REFLECT
  USE CONVERT_MOD,    ONLY: LON2X,LAT2Y,X2LON,Y2LAT
  USE RANDOM_MOD,     ONLY: INIT_GENRAND

  USE HTURB_MOD,      ONLY: HTURB
  USE VTURB_MOD,      ONLY: VTURB
  USE BEHAVIOR_MOD,   ONLY: INITBEHAVE,UPDATESTATUS,BEHAVE,GETCOLOR     
  USE SETTLEMENT_MOD, ONLY: SETTLED,DEAD,SETTLEMENT,DIE

  USE HYDRO_MOD, ONLY: initGrid,initHydro,updateHydro,setEle,setInterp,getInterp,getSlevel,getWlevel,WCTS_ITPI

IMPLICIT NONE
!	*************************************************************************
!	*																		*
!	*						Variable Declarations						    *
!	*																		*
!	*************************************************************************

!INPUT/OUTPUT FILE NAME CONSTRUCTION VARIABLES
CHARACTER(LEN=100) :: buffer2,filenm2
CHARACTER(LEN=4  ) :: prefix2,suffix2

! Initialization
REAL :: seconds,daytime
INTEGER :: stepT,stepIT,time,prcount,printdt,i,j,k,ii,deplvl,n,p,it
DOUBLE PRECISION :: idum_call_count

! Stepping
DOUBLE PRECISION :: ex(3),ix(3)

! Particle tracking
REAL :: P_latlon(numpar,3),P_nlatlon(numpar,3)
INTEGER :: startpoly(numpar),endpoly(numpar)
DOUBLE PRECISION :: P_xyz(numpar,3),newP_xyz(numpar,3),P_age(numpar,3),Xpar,Ypar,Zpar,P_depth,P_angle,  &
  P_zetab,P_zetac,P_zetaf,newXpos,newYpos,newZpos,P_zb,P_zc,P_zf,P_Salt(numpar),P_Temp(numpar)
DOUBLE PRECISION :: Pwc_zb(us),Pwc_zc(us),Pwc_zf(us),Pwc_wzb(ws),Pwc_wzc(ws),Pwc_wzf(ws)

! Behavior and Turbulence
DOUBLE PRECISION :: TurbHx,TurbHy,TurbV,Behav

! Boundaries
INTEGER :: intersectf,skipbound,in_island,inbounds,reflects
DOUBLE PRECISION :: reflect,fintersectX,fintersectY,freflectX,freflectY,Xpos,Ypos,nXpos,nYpos,island

! Settlement
INTEGER :: inpoly

! Advect
DOUBLE PRECISION :: AdvectX,AdvectY,AdvectZ,P_V,P_U,P_W,UAD,VAD,WAD,maxpartdepth,minpartdepth,          &
  kn1_u,kn1_v,kn1_w,kn2_u,kn2_v,kn2_w,kn3_u,kn3_v,kn3_w,kn4_u,kn4_v,kn4_w,x1,x2,x3,y1,y2,y3,z1,z2,z3

! Print
REAL :: color
INTEGER :: counter2,settle

!getEle Error Return Variable
INTEGER :: ele_err

! Character String for Program Termination Read Statements
CHARACTER :: anykey

! added by NCL
!logical :: particle_n_is_dead

! ***************************************************************************
! *																			*
! *							Initialize Model								*
! *																			*
! ***************************************************************************

write(*,*) ' '
write(*,*) '****** LTRANS INITIALIZATION *******'

! ***************************************************************************
! *			Initialize time steps, print counters, and constants			*
! ***************************************************************************

! THE FOLLOWING VARIABLE INITIALIZATIONS SHOULD NOT BE CHANGED:
seconds = days*24*60*60		!Total number of seconds to run model in model time
stepT = int(seconds/dt)		!number of external time steps
stepIT = int(dt/idt)		!number of internal time steps
time = 0					!initialize time to 0
prcount=0					!print counter; number of external time steps
printdt=0					!print counter
idum_call_count = 0.0		!counter for number of times random number generator is called 
CALL init_genrand(seed)     !set random number generator Seed Value

! ***************************************************************************
! *				Initialize particle tracking matrices						*
! ***************************************************************************

! Read-in lat/long of particles. If settlement module is on, read in    
! the habitat polygon on which the particle start                       
!   P_latlon(i,1) = latitude, P_latlon(i,2) = longitude,
!   P_xyz(i,3) = starting depth (negative m), startpoly(i) = starting habitat polygon

write(*,*) 'read in particle locations', numpar

OPEN(1,FILE=parfile)
1 format(F12.8,F12.8,D6.2,I6)

if(settlementon)then
  do i=1,numpar
    read (1,1) P_latlon(i,2),P_latlon(i,1),P_xyz(i,3),startpoly(i)
  enddo
else
  do i=1,numpar
    read (1,1) P_latlon(i,2),P_latlon(i,1),P_xyz(i,3)
  enddo
endif

CLOSE(1)

write(*,*) '  Particle n=5 Latitude=',P_latlon(5,1),'Longitude=',P_latlon(5,2)
write(*,*) '  Particle n=5 Depth=',P_xyz(5,3)
if(settlementon) write(*,*) '  Particle n=5 Start Polygon=',startpoly(5)

!	*******************************************************************
!	*					    Particle Attributes					      *
!	*******************************************************************

do i=1,numpar
  P_age(i,1) = Delay    !time at which particle starts
  P_age(i,2) = 0.0      !age of particle
  P_age(i,3) = 0.0      !time at which particle ends

  ! Do not specify: 
  ! Convert particle latitude and longitude to meters using equations from  
  ! sg_mercator.m and seagrid2roms.m in Seagrid.
  P_xyz(i,1) = lon2x(P_latlon(i,2))
  P_xyz(i,2) = lat2y(P_latlon(i,1)) 
enddo

write(*,*) 'particle 5 x, y, z'
write(*,*) p_xyz(5,1),p_xyz(5,2),p_xyz(5,3)

endpoly = 0                 !set end polygon location to zero

! ***************************************************************************
! *																			*
! *				Initial Read-In of Hydrodynamic Model Information			*
! *																			*
! ***************************************************************************

!Initialize Grid / Create Elements
CALL initGrid()

! ***************************************************************************
! *																			*
! *						Prepare for Particle Tracking						*
! *																			*
! ***************************************************************************

write(*,*) 'prepare boundary arrays'

!Create Boundaries
CALL createBounds()

!Initialize Behavior
CALL initBehave()

CALL initHydro()  !Read in initial hydrodynamoc model data

!	*************************************************************************
!	*																		*
!	*					Start Iterations - External Time Steps				*
!	*																		*
!	*************************************************************************

write(*,*) ' '
write(*,*) '****** BEGIN ITERATIONS *******'

! External time step
DO p=1,stepT     !External time step

  !Read in hydrodynamic model data 
  IF(p > 2) CALL updateHydro()      !do not start updating matrices until third iteration
 
  !Prepare external time step values to be used for 
  !  calculating Advection and Turbulence
  ex=0.0
  ex(1) = (p-2)*dt
  ex(2) = (p-1)*dt
  ex(3) = p*dt

!		*****************************************************************
!		*																*
!		*			Start Iterations - Internal Time Steps			    *
!		*																*
!		*****************************************************************

  DO it=1,stepIT

!		*****************************************************************
!		*					Advance Time  					            *
!		*****************************************************************

	time= time + idt  
	daytime= float(time)/(3600*24)

!		*****************************************************************
!		*		    	Adjust the Z-coordinate of Each Node			*
!		*****************************************************************

	!Prepare internal time step values to be used for 
	!  calculating Advection and Turbulence
	ix(1) = ex(2) + (it-2)*idt
	ix(2) = ex(2) + (it-1)*idt
	ix(3) = ex(2) + it*idt


!		*****************************************************************
!		*				Begin Loop to Update Each Particle			    *
!		*****************************************************************

	  DO n=1,numpar

!			*********************************************************
!			*														*
!			*	    Update Particle Age	and Characteristics		    *
!			*														*
!			*********************************************************

        !If the particle is not yet released, set new location to 
        !  current location, and cycle to next particle
        if(time <= P_age(n,1))then
          newP_xyz(n,1) = P_xyz(n,1)
          newP_xyz(n,2) = P_xyz(n,2)
          newP_xyz(n,3) = P_xyz(n,3)
          cycle
        endif

        !Update particle age
		P_age(n,2) = P_age(n,2) + float(idt)

        !Update particle status if settled or dead
		CALL updateStatus(P_age(n,2),n)

		!If particle settled or dead, skip tracking
        if(settlementon)then
          if ( SETTLED(n) .OR. DEAD(n) ) cycle
        endif
        
        ! NCL the code doesn't really work well. We want to be able to turn off particles
        ! that no longer make sense (left the domain, stuck on land, etc.) but we don't want
        ! to use the settlement module because it is overly complicated.
        if ( SETTLED(n) .OR. DEAD(n) ) then
!          write(*,*) n, ' particle is dead'
          cycle
        endif
        
        
!			*********************************************************
!			*														*
!			*			Find Element that Contains Particle		    *
!			*														*
!			*********************************************************

		!Find rho/u/v elements in which the particle is located. Store in P_element matrices  

		Xpar = P_xyz(n,1)
		Ypar = P_xyz(n,2)
		if( (p .EQ. 1) .AND. (it .EQ. 1) ) then !if the first iteration
         
		  inbounds = 0
		  !Determine if particle is within model bounadaries
		  call mbounds(Ypar,Xpar,inbounds)
		  if (inbounds.EQ.0) then 
			write(*,*) 'outside main bounds, n=',n
			call DIE(n)
!			write(*,*) ' '
!			write(*,*) 'The Program Cannot Continue and Will Terminate'
!			write(*,*) 'Press Any Key and Enter'
!			read(*,*) anykey
!			stop
!			exit
		  endif

          in_island = 0	
          call ibounds(in_island,Ypar,Xpar,island)
          if (in_island.EQ.1) then
            write(*,*) 'in island, n=',n
            call DIE(n)
!            write(*,*) ' '
!            write(*,*) 'The Program Cannot Continue and Will Terminate'
!            write(*,*) 'Press Any Key and Enter'
!            read(*,*) anykey
!            stop
          endif

          !Determine which Rho, U, & V elements the particle is in
          CALL setEle(Xpar,Ypar,n,ele_err,.TRUE.)

		else !if not the first iteration 

          !Determine which Rho, U, & V elements the particle is in
          CALL setEle(Xpar,Ypar,n,ele_err)

		endif


        !If the particle was not found to be within an element,
		!  write a message to the screen and discontinue the program
        IF(ele_err > 0)THEN
          call DIE(n)
!          SELECT CASE (ele_err)
!            CASE(1)
!              write(*,*) n,'start - particle outside model domain (rho)'
!            CASE(2)
!              write(*,*) n,'start - particle outside model domain (u)'
!            CASE(3)
!              write(*,*) n,'start - particle outside model domain (v)'
!            CASE(4)
!              write(*,*) " "
!              write(*,*) "Jumped over a rho element"
!              write(*,*) " "
!              write(*,*)" - Now Stopped -"
!            CASE(5)
!              write(*,*) " "
!              write(*,*) "Jumped over a u element"
!              write(*,*) " "
!              write(*,*)" - Now Stopped -"
!            CASE(6)
!              write(*,*) " "
!              write(*,*) "Jumped over a v element"
!              write(*,*) " "
!              write(*,*)" - Now Stopped -"
!          END SELECT
!
!          write(*,*) ' '
!          write(*,*) 'The Program Cannot Continue and Will Terminate'
!          write(*,*) 'Press Enter Key'
!          read(*,*) anykey
!          stop
!          write(*,*) 'Problem with particle ', n
          exit
        else
!          write(*,*) 'Should be good to go with particle ', n
        ENDIF

        !Set Interpolation Values for the current particle
        CALL setInterp(Xpar,Ypar,n)

!			*********************************************************
!			*														*
!			*		Ensure Particle is Within Verticle Bounds	    *
!			*														*
!			*********************************************************

		!Find depth, angle, and sea surface height at particle location

        P_depth = -1.* getInterp("depth")
        P_angle = getInterp("angle")
        P_zetab = getInterp("zetab")
        P_zetac = getInterp("zetac")
        P_zetaf = getInterp("zetaf")

		!Check if particle location above or below boundary, If so, place
		!  just within boundary (1 mm)
		if (P_xyz(n,3).LT.P_depth) P_xyz(n,3) = P_depth + 0.001
		P_zb = P_xyz(n,3)
		P_zc = P_xyz(n,3)
		P_zf = P_xyz(n,3)
		if (P_xyz(n,3).GT.P_zetab) P_zb = P_zetab - 0.001
		if (P_xyz(n,3).GT.P_zetac) P_zc = P_zetac - 0.001
		if (P_xyz(n,3).GT.P_zetaf) P_zf = P_zetaf - 0.001

		Zpar = P_xyz(n,3)

!			*********************************************************
!			*														*
!			*		 	Create Matrix of Z-coordinates			    *
!			*														*
!			*********************************************************

		!Create matrix of z-coordinates at particle and at each node for
		!  back, center, forward times
		do i=1,us

		  !Rho-coordinate depths at particle location
		  Pwc_zb(i)=getSlevel(P_zetab,P_depth,i)
		  Pwc_zc(i)=getSlevel(P_zetac,P_depth,i)
		  Pwc_zf(i)=getSlevel(P_zetaf,P_depth,i)

		  !W-coordinate depths at particle location
		  Pwc_wzb(i)= getWlevel(P_zetab,P_depth,i)
		  Pwc_wzc(i)= getWlevel(P_zetac,P_depth,i)
		  Pwc_wzf(i)= getWlevel(P_zetaf,P_depth,i)

		enddo

		!W-coordinate depths at particle location (cont.)
		Pwc_wzb(i)= getWlevel(P_zetab,P_depth,ws)
		Pwc_wzc(i)= getWlevel(P_zetac,P_depth,ws)
		Pwc_wzf(i)= getWlevel(P_zetaf,P_depth,ws)


!			*********************************************************
!			*														*
!			*			Prepare for Particle Movement			    *
!			*														*
!			*********************************************************

		AdvectX = 0.0
		AdvectY = 0.0
		AdvectZ = 0.0
		TurbHx = 0.0
		TurbHy = 0.0
		TurbV = 0.0
		Behav = 0.0

!			*********************************************************
!			*														*
!			*                 	   ADVECTION                        *
!			*														*
!			*********************************************************

		maxpartdepth = Pwc_wzb(1)
		if(Pwc_wzc(1) .GT. maxpartdepth) maxpartdepth = Pwc_wzc(1)
		if(Pwc_wzf(1) .GT. maxpartdepth) maxpartdepth = Pwc_wzf(1)

		minpartdepth = Pwc_wzb(ws)			
		if(Pwc_wzc(ws) .LT. minpartdepth) minpartdepth = Pwc_wzc(ws)
		if(Pwc_wzf(ws) .LT. minpartdepth) minpartdepth = Pwc_wzf(ws)

		!Find advection currents at original coordinates
		CALL FIND_CURRENTS(Xpar,Ypar,Zpar,Pwc_zb,Pwc_zc,Pwc_zf,      &
		  Pwc_wzb,Pwc_wzc,Pwc_wzf,P_zb,P_zc,P_zf,ex,ix,p,1,Uad,Vad,Wad)

		!Store advection currents at original coordinates
		kn1_u = Uad
		kn1_v = Vad
		kn1_w = Wad

		!Estimate new coordinates for next RK position
		x1 = Xpar + (Uad*cos(P_angle) - Vad*sin(P_angle)) * (idt/2) 
		y1 = Ypar + (Uad*sin(P_angle) + Vad*cos(P_angle)) * (idt/2) 
		z1 = Zpar +   Wad * (idt/2) 
		if(z1 .GT. minpartdepth) z1 = minpartdepth - 0.000001
		if(z1 .LT. maxpartdepth) z1 = maxpartdepth + 0.000001

		!Find advection currents at estimated next RK position
		CALL FIND_CURRENTS(x1,y1,z1,Pwc_zb,Pwc_zc,Pwc_zf,Pwc_wzb,Pwc_wzc,Pwc_wzf,   &
          P_zb,P_zc,P_zf,ex,ix,p,2,Uad,Vad,Wad)

		!Store advection currents at 2nd RK position
		kn2_u = Uad
		kn2_v = Vad
		kn2_w = Wad

		!Estimate new coordinates for next RK position
		x2 = Xpar + (Uad*cos(P_angle) - Vad*sin(P_angle)) * (idt/2) 
		y2 = Ypar + (Uad*sin(P_angle) + Vad*cos(P_angle)) * (idt/2) 
		z2 = Zpar + (Wad * (idt/2))
		if(z2 .GT. minpartdepth) z2 = minpartdepth - 0.000001
		if(z2 .LT. maxpartdepth) z2 = maxpartdepth + 0.000001

		!Find advection currents at estimated next RK position
		CALL FIND_CURRENTS(x2,y2,z2,Pwc_zb,Pwc_zc,Pwc_zf,Pwc_wzb,Pwc_wzc,Pwc_wzf,   &
          P_zb,P_zc,P_zf,ex,ix,p,2,Uad,Vad,Wad)

		!Store advection currents at 3rd RK position
		kn3_u = Uad
		kn3_v = Vad
		kn3_w = Wad

		!Calculate the coordinates at the final position
		x3 = Xpar + (Uad*cos(P_angle) - Vad*sin(P_angle)) * (idt/2) 
		y3 = Ypar + (Uad*sin(P_angle) + Vad*cos(P_angle)) * (idt/2) 
		z3 = Zpar + (Wad * idt)
		if(z3 .GT. minpartdepth) z3 = minpartdepth - 0.000001
		if(z3 .LT. maxpartdepth) z3 = maxpartdepth + 0.000001

		!Find advection currents at the final position
		CALL FIND_CURRENTS(x3,y3,z3,Pwc_zb,Pwc_zc,Pwc_zf,Pwc_wzb,Pwc_wzc,Pwc_wzf,   &
          P_zb,P_zc,P_zf,ex,ix,p,3,Uad,Vad,Wad)

		!Store advection currents at final position
		kn4_u = Uad
		kn4_v = Vad
		kn4_w = Wad

		!Use the RK formula to get the final Advection values
		P_U = (kn1_u + 2*kn2_u + 2*kn3_u + kn4_u)/6.
		P_V = (kn1_v + 2*kn2_v + 2*kn3_v + kn4_v)/6.
		P_W = (kn1_w + 2*kn2_w + 2*kn3_w + kn4_w)/6.


		AdvectX = idt*(P_U*cos(P_angle) - P_V*sin(P_angle))
		AdvectY = idt*(P_U*sin(P_angle) + P_V*cos(P_angle))
		AdvectZ = idt*P_W


!			*********************************************************
!			*														*
!			*    	      Salinity and Temperature			        *
!			*														*
!			*********************************************************

		IF (SaltTempOn) THEN
        
		  !Calculate salinity and temperture at the particle location       
	      do i=3,us-2
            if ((Zpar .LT. Pwc_zb(i)) .OR. (Zpar .LT. Pwc_zc(i)) .OR. (Zpar .LT. Pwc_zf(i))) exit
          enddo
          deplvl = i-2    
          P_Salt(n) = WCTS_ITPI("salt",Xpar,Ypar,deplvl,Pwc_zb,Pwc_zc,Pwc_zf,us,P_zb,P_zc,P_zf,ex,ix,p,4)
          P_Temp(n) = WCTS_ITPI("temp",Xpar,Ypar,deplvl,Pwc_zb,Pwc_zc,Pwc_zf,us,P_zb,P_zc,P_zf,ex,ix,p,4)

		ENDIF

!			*********************************************************
!			*														*
!			*			  Horizontal Turbulence			            *
!			*														*
!			*********************************************************

		IF (HTurbOn) THEN

		  CALL HTurb(TurbHx,TurbHy)
		  idum_call_count = idum_call_count + 2.0

		ENDIF


!			*********************************************************
!			*														*
!			*				Verticle Turbulence				        *
!			*														*
!			********************************************************* 

		IF (VTurbOn) THEN

		  CALL VTurb(P_zc,P_depth,P_zetac,p,ex,ix,Pwc_wzb,Pwc_wzc,Pwc_wzf,TurbV)

		  idum_call_count = idum_call_count + 2.0*float(idt/2)

		endif


!			*********************************************************
!			*														*
!			*					   Behavior					        *
!			*														*
!			*********************************************************

        if (Behavior.NE.0) then 
		                                                                   
          CALL Behave(Xpar,Ypar,Zpar,Pwc_zb,Pwc_zc,Pwc_zf,P_zb,P_zc,P_zf,P_zetac,P_age(n,2),P_depth,n,it,ex,ix,daytime,p,Behav)

          idum_call_count = idum_call_count + 2.0                           

        endif                                                               


!			*********************************************************
!			*														*
!			*	Update Particle Locations and Check Boundaries  	*
!			*														*
!			*********************************************************

		newXpos = 0.0
		newYpos = 0.0
		newZpos = 0.0
		
		!Update due to Advection and Turbulence		
		newXpos = P_xyz(n,1) + AdvectX + TurbHx 
		newYpos = P_xyz(n,2) + AdvectY + TurbHy
		newZpos = P_xyz(n,3) + AdvectZ + TurbV

		!Check vertical boundares and reflect
		!  if particle above surface, particle reflects off surface
		reflect=0.0
		if (newZpos.GT.P_zetac) then
		  reflect = P_zetac - newZpos
		  NewZpos = P_zetac + reflect
		endif

		!  if particle deeper than bottom, particle reflects off bottom
		if (newZpos.LT.P_depth) then
		  reflect = P_depth - newZpos
		  NewZpos = P_depth + reflect
		endif 

        !Update due to Behavior	
		newZpos = NewZpos + Behav

		!Check vertical boundares and move back in domain
		!  if particle above surface, particle moves just below surface
		if (newZpos.GT.P_zetac) NewZpos = P_zetac - 0.000001

		!  if particle deeper than bottom, particle moves just off bottom
		if (newZpos.LT.P_depth) NewZpos = P_depth + 0.000001

        !Horizontal boundary tests. Check if particle still within horizontal domain
		!If not, reflect particle off main boundary or island walls 
        Xpos = P_xyz(n,1)
        Ypos = P_xyz(n,2)
        nXpos = newXpos
        nYpos = newYpos
        skipbound = -1
        reflects = 0
        do
          if (.NOT. DEAD(n)) then 
            call intersect_reflect(Xpos,Ypos,nXpos,nYpos,fintersectX,fintersectY,         &
                   freflectX,freflectY,intersectf,skipbound)
            if(intersectf == 0)exit
            reflects = reflects + 1
            if(reflects > 3) then
              write(*,*) n,'still out after 3rd reflection.'
              call DIE(n)
!             write(*,*) ' '
!             write(*,*) 'The Program Cannot Continue and Will Terminate'
!             write(*,*) 'Press Any Key and Enter'
!             read(*,*) anykey
!             stop
             !Set particle location equal to previous location
              newXpos=P_xyz(n,1)
              newYpos=P_xyz(n,2)
		    endif
		    Xpos = fintersectX
		    Ypos = fintersectY
		    nXpos = freflectX
		    nYpos = freflectY
		  else
		    exit
		  endif
       enddo

       newXpos = nXpos
       newYpos = nYpos

        !Check to make sure new position is within model boundaries
       if (.NOT. DEAD(n)) then
        call mbounds(newYpos,newXpos,inbounds)
        if(inbounds /= 1) then
          write(*,*) 'ERROR: Particle Outside Main Boundaries After Intersect_Reflect'
!          write(*,*) 'Model Run Cannot Continue'
!          write(*,*) 'Press Any Key and Enter'
!          read(*,*) anykey
          call DIE(n)
          cycle
        endif
       else
         cycle
       endif

        !Check to make sure new position is not within an island
       if (.NOT. DEAD(n)) then
         call ibounds(in_island,newYpos,newXpos,island)
         if(in_island == 1) then
           write(*,*) 'ERROR: Particle Inside Island Boundaries After Intersect_Reflect'
!           write(*,*) 'Model Run Cannot Continue'
!           write(*,*) 'Press Any Key and Enter'
!           read(*,*) anykey
           call DIE(n)
           cycle
         endif
       else
         cycle
       endif
!      End boundary condition tests ******************* 

		!Assign new particle positions to newP_xyz matrix
		newP_xyz(n,1) = newXpos
		newP_xyz(n,2) = newYpos
		newP_xyz(n,3) = newZpos

        ! Check to make sure new position is within a rho, u and v element
        CALL setEle(nXpos,nYpos,n,ele_err)

        IF(ele_err > 0)THEN
          write(*,*) " "

		  SELECT CASE (ele_err)
            CASE(4)
              write(*,*) "Jumped over a rho element"
            CASE(5)
              write(*,*) "Jumped over a u element"
            CASE(6)
              write(*,*) "Jumped over a v element"
          END SELECT
          call DIE(n)
!
!          write(*,*) " "
!          write(*,*)" - Now Stopped -"
!          write(*,*) ' '
!          write(*,*) 'The Program Cannot Continue and Will Terminate'
!          write(*,*) 'Press Enter Key'
!          read(*,*) anykey
!          stop
          exit
        ENDIF

!			*********************************************************
!			*														*
!			*				Settlement   				            *
!			*														*
!			*********************************************************
        if(settlementon) then

		  CALL settlement(P_age(n,2),n,P_xyz(n,1),P_xyz(n,2),inpoly)
		  if (inpoly .GT. 0) then
			newP_xyz(n,3) = P_depth
			endpoly(n) = inpoly
			P_age(n,3) = P_age(n,2)
		  endif

		endif !  ********** end Settlement ******************************


!		*****************************************************************
!		*						End of Particle Loop				    *
!		*****************************************************************

	  ENDDO	!end loop for each particle

!			*********************************************************
!			*			 Update particle locations		  	        *
!			*********************************************************

	do n=1,numpar
	  P_xyz(n,1) = newP_xyz(n,1)
	  P_xyz(n,2) = newP_xyz(n,2)
	  P_xyz(n,3) = newP_xyz(n,3)
	enddo

!		*****************************************************************
!		*			     End of Internal Time Step Loop				    *
!		*****************************************************************

  ENDDO	!end internal time step

!  write(*,*) 'end track'


!		*****************************************************************
!		*																*
!		*						Print Statements					    *
!		*																*
!		*****************************************************************
!
  prcount = time/dt
  if (prcount.GE.1) printdt=printdt+dt

  IF (printdt.EQ.iprint .OR. p.EQ.1) THEN
	write(*,*) 'write output to file, day = ',daytime  

	!Convert particle position (in meters) to latitude and longitude using equations
	!from sg_mercator.m and seagrid2roms.m in Seagrid. 
	do n=1,numpar
	  P_nlatlon(n,2) = x2lon(P_xyz(n,1))
	  P_nlatlon(n,1) = y2lat(P_xyz(n,2))
	enddo

	!Create a filename and unit number for each iteration    
	counter2=prcount+10000000
	prefix2='para'
	suffix2='.csv'
	write(buffer2,'(A,I8,A)') prefix2,counter2,suffix2
	read(buffer2,'(A)') filenm2
	open(2,FILE=filenm2(1:LEN_TRIM(filenm2)),STATUS='REPLACE')
	2 format(F10.6,',',F10.6,',', F10.6,',', F10.6,',', F10.6,',', F10.6)
	5 format(F10.6,',',F10.6,',', F12.6,',', F10.6)

	!Write data to file                                      
	do ii = 1,numpar

      !Find identification number that describes a particle's behavior type or status  
	  !for use in visualization routines
      color = getColor(ii)

 	  if (SaltTempOn) then
	     !Write salinity and temperature at the particle location from the previous internal time step (idt)
	     write(2,2) P_xyz(ii,3),color,P_nlatlon(ii,2),P_nlatlon(ii,1),P_Salt(ii),P_Temp(ii)
	  else
	     write(2,5) P_xyz(ii,3),color,P_nlatlon(ii,2),P_nlatlon(ii,1)
	  endif 

	enddo

	close(2)

	printdt=0  !reset print counter	
  ENDIF


!	*************************************************************************
!	*																		*
!	*					End of External Time Step Loop					    *
!	*																		*
!	*************************************************************************


! ********************** END ITERATIONS ***************************
! *****************************************************************

ENDDO !end external time step


!	*************************************************************************
!	*		     Write final positions and status to file			        *
!	*************************************************************************

! Convert particle position (in meters) to latitude and longitude using equations
! from sg_mercator.m and seagrid2roms.m in Seagrid. 
do n=1,numpar
  P_nlatlon(n,2) = x2lon(P_xyz(n,1))
  P_nlatlon(n,1) = y2lat(P_xyz(n,2))
enddo


open(3,FILE='endfile.csv',STATUS='REPLACE')
  3 format(I8,',',I8,',', I8,',', F25.8,',',F25.15,',',F25.15,',',F25.15)
  4 format(F25.8,',',F25.15,',',F25.15,',',F25.15)
do n=1,numpar
  color = getColor(n)
  if(settlementon)then
    settle = 0
    if(SETTLED(n))settle = 1
    if(DEAD(n))settle = 2
    write(3,3) startpoly(n),endpoly(n),settle,color,P_nlatlon(n,1),P_nlatlon(n,2),P_age(n,3) 
  else
    write(3,4) color,P_nlatlon(n,1),P_nlatlon(n,2),P_age(n,3) 
  endif
enddo
close(3)

write(*,*) 'write endfile.csv  '

write(*,*) '  '
write(*,*) 'Number of times random number generator was called:'
write(*,*) '                  ', int(idum_call_count)


write(*,*) '  '
write(*,*) '****** END LTRANS *******'

!	**************************************************************
!	*		            End Main Program 			             *
!	**************************************************************


CONTAINS


  SUBROUTINE FIND_CURRENTS(Xpar,Ypar,Zpar,Pwc_zb,Pwc_zc,Pwc_zf,Pwc_wzb,Pwc_wzc,Pwc_wzf,  &
    P_zb,P_zc,P_zf,ex,ix,p,version,Uad,Vad,Wad)
    !This Subroutine calculates advection currents at the particles location in space and time

    USE PARAM_MOD,  ONLY: us,ws,z0
    USE HYDRO_MOD,  ONLY: interp,WCTS_ITPI
	USE TENSION_MOD, ONLY: TSPSI,HVAL
	USE INT_MOD,    ONLY: linint,polintd
	IMPLICIT NONE

	INTEGER, INTENT(IN) :: p,version
	DOUBLE PRECISION, INTENT(IN) :: Xpar,Ypar,Zpar,P_zb,P_zc,P_zf,ex(3),ix(3),           &
      Pwc_zb(us),Pwc_zc(us),Pwc_zf(us),Pwc_wzb(ws),Pwc_wzc(ws),Pwc_wzf(ws)
    DOUBLE PRECISION, INTENT(OUT) :: Uad,Vad,Wad

    INTEGER, PARAMETER :: nN = 4     !Number of Depth Levels to create tension spline with

	INTEGER :: i,ii,iii
    DOUBLE PRECISION :: P_Ub,P_Uc,P_Uf,P_Vb,P_Vc,P_Vf,P_Wb,P_Wc,P_Wf,ey(3),              &
	  Pwc_ub,Pwc_uc,Pwc_uf,Pwc_vb,Pwc_vc,Pwc_vf,Pwc_wb,Pwc_wc,Pwc_wf

		!version: 1 = return b, 2 = return c, 3 = return f

	! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	!		Determine the Lowest Numbered US-Level of the Closest Four
	! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    do i=3,us-2
      if ((Zpar .LT. Pwc_zb(i)) .OR. (Zpar .LT. Pwc_zc(i)) .OR. (Zpar .LT. Pwc_zf(i))) exit
    enddo
    ii = i-2

	! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	!		Determine the Lowest Numbered WS-Level of the Closest Four
	! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    do i=3,ws-2
      if ((Zpar .LT. Pwc_wzb(i)) .OR. (Zpar .LT. Pwc_wzc(i)) .OR. (Zpar .LT. Pwc_wzf(i))) exit
    enddo
    iii = i - 2

	!			*********************************************************
	!			*														*
	!			*		Calculate U,V,W in Water Column Profile			*
	!			*														*
	!			*********************************************************


	!i. Determine if particle is deep enough that velocities are affected by the bottom.
	! If so, apply log layer between deepest current velocity predicitons 
	! (deepest rho s-level for u,v and deepest w s-level for w) and bottom. 
	if ((Zpar .LT. Pwc_zb(1)) .OR. (Zpar .LT. Pwc_zc(1)) .OR. (Zpar .LT. Pwc_zf(1))) then

	  Pwc_Ub = interp(Xpar,Ypar,"uvelb",1)
	  Pwc_Uc = interp(Xpar,Ypar,"uvelc",1)
	  Pwc_Uf = interp(Xpar,Ypar,"uvelf",1)
	  Pwc_Vb = interp(Xpar,Ypar,"vvelb",1)
	  Pwc_Vc = interp(Xpar,Ypar,"vvelc",1)
	  Pwc_Vf = interp(Xpar,Ypar,"vvelf",1)
	  Pwc_Wb = interp(Xpar,Ypar,"wvelb",2)
	  Pwc_Wc = interp(Xpar,Ypar,"wvelc",2)
	  Pwc_Wf = interp(Xpar,Ypar,"wvelf",2)

	  P_Ub = Pwc_Ub*log10((Zpar-Pwc_wzb(1))/z0)/log10((Pwc_zb(1) -Pwc_wzb(1))/z0)
	  P_Uc = Pwc_Uc*log10((Zpar-Pwc_wzb(1))/z0)/log10((Pwc_zc(1) -Pwc_wzb(1))/z0)
	  P_Uf = Pwc_Uf*log10((Zpar-Pwc_wzb(1))/z0)/log10((Pwc_zf(1) -Pwc_wzb(1))/z0)
	  P_Vb = Pwc_Vb*log10((Zpar-Pwc_wzb(1))/z0)/log10((Pwc_zb(1) -Pwc_wzb(1))/z0)
	  P_Vc = Pwc_Vc*log10((Zpar-Pwc_wzb(1))/z0)/log10((Pwc_zc(1) -Pwc_wzb(1))/z0)
	  P_Vf = Pwc_Vf*log10((Zpar-Pwc_wzb(1))/z0)/log10((Pwc_zf(1) -Pwc_wzb(1))/z0)
	  P_Wb = Pwc_Wb*log10((Zpar-Pwc_wzb(1))/z0)/log10((Pwc_wzb(2)-Pwc_wzb(1))/z0)
	  P_Wc = Pwc_Wc*log10((Zpar-Pwc_wzb(1))/z0)/log10((Pwc_wzc(2)-Pwc_wzb(1))/z0)
	  P_Wf = Pwc_Wf*log10((Zpar-Pwc_wzb(1))/z0)/log10((Pwc_wzf(2)-Pwc_wzb(1))/z0)

      !		*********************************************************
      !		*		 Find Internal b,c,f and Advection Values		*
      !		*********************************************************
      !
      ! ii. fit polynomial to hydrodynamic model output and find internal b,c,f values

      !a. U velocity
      !	1. Prepare external time step values
      if (p .EQ. 1) then
        ey=0.0
        ey(1) = P_Ub
        ey(2) = P_Ub
        ey(3) = P_Uc
      else
        ey=0.0
        ey(1) = P_Ub
        ey(2) = P_Uc
        ey(3) = P_Uf
      endif

      !	2. Get Advection value
      if(version .EQ. 1) then
        Uad = polintd(ex,ey,3,ix(1))
      elseif (version .EQ. 2) then
        Uad = polintd(ex,ey,3,ix(2))
      else
        Uad = polintd(ex,ey,3,ix(3))
      endif

      !b. V velocity
      !	1. Prepare external time step values
      if (p .EQ. 1) then
        ey=0.0
        ey(1) = P_Vb
        ey(2) = P_Vb
        ey(3) = P_Vc
      else
        ey=0.0
        ey(1) = P_Vb
        ey(2) = P_Vc
        ey(3) = P_Vf
      endif

      !	2. Get Advection value
      if(version .EQ. 1) then
        Vad = polintd(ex,ey,3,ix(1))
      elseif (version .EQ. 2) then
        Vad = polintd(ex,ey,3,ix(2))
      else
        Vad = polintd(ex,ey,3,ix(3))
      endif


      !c. W velocity
      !	1. Prepare external time step values
      if (p .EQ. 1) then
        ey=0.0
        ey(1) = P_Wb
        ey(2) = P_Wb
        ey(3) = P_Wc
      else
        ey=0.0
        ey(1) = P_Wb
        ey(2) = P_Wc
        ey(3) = P_Wf
      endif

      !	2. Get Advection value
      if(version .EQ. 1) then
        Wad = polintd(ex,ey,3,ix(1))
      elseif (version .EQ. 2) then
        Wad = polintd(ex,ey,3,ix(2))
      else
        Wad = polintd(ex,ey,3,ix(3))
      endif

	else

      Uad = WCTS_ITPI("uvel",Xpar,Ypar,ii ,Pwc_zb ,Pwc_zc ,Pwc_zf ,us,P_zb,P_zc,P_zf,ex,ix,p,version)
      Vad = WCTS_ITPI("vvel",Xpar,Ypar,ii ,Pwc_zb ,Pwc_zc ,Pwc_zf ,us,P_zb,P_zc,P_zf,ex,ix,p,version)
      Wad = WCTS_ITPI("wvel",Xpar,Ypar,iii,Pwc_wzb,Pwc_wzc,Pwc_wzf,ws,P_zb,P_zc,P_zf,ex,ix,p,version)

    endif

	RETURN
  END SUBROUTINE FIND_CURRENTS


END PROGRAM
 
!	*************************************************************************
!	*																		*
!	*								The End									*
!	*																		*
!	*************************************************************************


