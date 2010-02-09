MODULE SETTLEMENT_MOD

!  The Settlement Module handles all code related to the settlement routine.  This includes 
!  reading in the habitat polygons and holes, creating variables containing the specifications 
!  of the habitat polygons and holes, keeping track of the settlement status of every particle,
!  and checking if the particle is within a habitat polygon and can settle.
!
!  Original concepts and code by:             Elizabeth North
!  Module creation and code modification by:  Zachary Schlag
!  Created on:			                      2005
!  Last modified on:	                      27 Aug 2008


  USE PARAM_MOD, ONLY: numpar,rho_elements,minholeid,maxholeid,minpolyid,maxpolyid,pedges,hedges
  IMPLICIT NONE
  SAVE
  PRIVATE

  DOUBLE PRECISION :: polys(pedges,5)
	! polys(i,1)=habitat polygon number
	! polys(i,2)=center x
	! polys(i,3)=center y
	! polys(i,4)=edge x
	! polys(i,5)=edge y

  DOUBLE PRECISION :: holes(hedges,6)
	! Holes(i,1)=hole number
	! Holes(i,2)=center longitude
	! Holes(i,3)=center latitude
	! Holes(i,4)=edge longitude
	! Holes(i,5)=edge latitude
	! Holes(i,6)=habitat polygon number

  !Maximum distance from each habitat polygon's center to its farthest edge point
  DOUBLE PRECISION :: maxbdis(minpolyid:maxpolyid)

  !Maximum distance from each hole's center to its farthest edge point
  DOUBLE PRECISION :: maxhdis(minholeid:maxholeid)

  INTEGER :: settle(numpar)                      !(0 - not settled, 1 - settled, 2 - dead)
  DOUBLE PRECISION :: settletime(numpar)         !age particles are competent to settle

  TYPE :: polyPerEle
    INTEGER :: numpoly                           !   number of polygons in each element 
	                                             !OR number of holes in each polygon
    INTEGER, ALLOCATABLE, DIMENSION(:) :: poly   !   id number of each polygon in the element
	                                             !OR id number of holes in the polygon
  END TYPE polyPerEle

  INTEGER :: polyspecs(minpolyid:maxpolyid,2)    !polyspecs(n,1) = location in polys of first point of habitat polygon n
                                                 !polyspecs(n,2) = # of edge points that make up habitat polygon n
  TYPE (polyPerEle) :: elepolys(rho_elements)    !id # of every polygon in each element
  INTEGER :: holespecs(minholeid:maxholeid,2)    !holespecs(n,1) = location in holes of first point of hole n
                                                 !holespecs(n,2) = # of edge points that make up hole n
  TYPE (polyPerEle) :: polyholes(minpolyid:maxpolyid) !id # of every hole in each polygon

  !The following procedures have been made public for the use of other program units:
  PUBLIC :: initSettlement,settlement,SETTLED,DEAD,DIE

CONTAINS

  SUBROUTINE initSettlement(P_pediage)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: P_pediage(numpar)
    INTEGER :: n

    settle = 0  !0=not settled, 1=successful settlement, 2=dead

    do n=1,numpar
      settletime(n) = P_pediage(n)
    enddo

    call readinHabitat()
    call createPolySpecs()

  END SUBROUTINE initSettlement


  SUBROUTINE readinHabitat()
    USE PARAM_MOD, ONLY: HABITATFILE,HOLEFILE,HOLESEXIST
    USE CONVERT_MOD, ONLY: LON2X,LAT2Y
    IMPLICIT NONE

	INTEGER :: i,curpoly
	DOUBLE PRECISION :: dise,P_lonlat(pedges,5),H_lonlat(hedges,6)
	! ***********************************************************************
	! *				2C.	Initialize settlement model							*
	! ***********************************************************************

    write(*,*) 'read in habitat polygon locations'

	! need to import center and edge coordinate data
    ! User specified file name
	OPEN(1,FILE=habitatfile)
    1 format(F10.0, F12.8, F11.8, F12.8, F11.8)

      do i=1,pedges
        read (1,*) P_lonlat(i,1),P_lonlat(i,2),P_lonlat(i,3),P_lonlat(i,4),P_lonlat(i,5)
      enddo
      ! P_lonlat(i,1)=habitat polygon number,
      ! P_lonlat(i,2)=center longitude, P_lonlat(i,3)=center latitude,
      ! P_lonlat(i,4)=edge longitude,   P_lonlat(i,5)=edge latitude
	CLOSE(1)

    write(*,*) '  Edge i=5 Center Lat=',P_lonlat(1,3),'Long=',P_lonlat(1,2)
    write(*,*) '  Edge i=5   Edge Lat=',P_lonlat(1,5),'Long=',P_lonlat(1,4)


	! Convert particle latitude and longitude to meters using equations from  
	! sg_mercator.m and seagrid2roms.m in Seagrid. 
	! polys(i,1)=habitat polygon number
	! polys(i,2)=center x
	! polys(i,3)=center y
	! polys(i,4)=edge x
	! polys(i,5)=edge y
	do i=1,pedges	
	  polys(i,1) = P_lonlat(i,1)

	  polys(i,2) = lon2x(P_lonlat(i,2))
	  polys(i,3) = lat2y(P_lonlat(i,3))

	  polys(i,4) = lon2x(P_lonlat(i,4))
	  polys(i,5) = lat2y(P_lonlat(i,5))
	enddo

	! Determine maximum distance between center and edge coordinates for each habitat polygon
	maxbdis = 1.0					! this is used to restrict settlement model search
	do i=1,pedges
	  curpoly = int(polys(i,1))
	  dise=sqrt((polys(i,2)-polys(i,4))**2+(polys(i,3)-polys(i,5))**2)
	  if (dise.GT.maxbdis(curpoly)) maxbdis(curpoly)=dise
	enddo


    if(holesExist)then

      ! User specified file name
      OPEN(2,FILE=holefile)
        2 format(F10.0, F12.8, F11.8, F12.8, F11.8, F10.0)

        do i=1,hedges
          read (2,*) H_lonlat(i,1),H_lonlat(i,2),H_lonlat(i,3),H_lonlat(i,4),H_lonlat(i,5),H_lonlat(i,6)
        enddo
        ! H_lonlat(i,1)=hole number,H_lonlat(i,2)=center longitude,H_lonlat(i,3)=center latitude,
        ! H_lonlat(i,4)=edge longitude,H_lonlat(i,5)=edge latitude,H_lonlat(i,6)=habitat polygon number
      CLOSE(2)

      write(*,*) '  Hole i=5 Center Lat=',H_lonlat(1,3),'Long=',H_lonlat(1,2)
      write(*,*) '  Hole i=5   Edge Lat=',H_lonlat(1,5),'Long=',H_lonlat(1,4)


      ! Convert particle latitude and longitude to meters using equations from  
      ! sg_mercator.m and seagrid2roms.m in Seagrid. 
      ! Holes(i,1)=hole number
      ! Holes(i,2)=center longitude
      ! Holes(i,3)=center latitude
      ! Holes(i,4)=edge longitude
      ! Holes(i,5)=edge latitude
      ! Holes(i,6)=habitat polygon number
      do i=1,hedges													
        holes(i,1) = H_lonlat(i,1)

        holes(i,2) = lon2x(H_lonlat(i,2))
        holes(i,3) = lat2y(H_lonlat(i,3))

        holes(i,4) = lon2x(H_lonlat(i,4))
        holes(i,5) = lat2y(H_lonlat(i,5))

        holes(i,6) = H_lonlat(i,6)
      enddo 


      ! Determine maximum distance between center and edge coordinates for each hole
      maxhdis = 1.0			! this is used to restrict settlement model hole search
      do i=1,hedges
        curpoly = int(holes(i,1))
        dise=sqrt((holes(i,2)-holes(i,4))**2+(holes(i,3)-holes(i,5))**2)
        if (dise.GT.maxhdis(curpoly)) maxhdis(curpoly)=dise
      enddo

    endif

  END SUBROUTINE readinHabitat


  SUBROUTINE createPolySpecs()
    USE PARAM_MOD, ONLY: rho_elements,holesExist
    USE HYDRO_MOD, ONLY: getR_ele
    USE GRIDCELL_MOD, ONLY: GRIDCELL
    USE PIP_MOD, ONLY: INPOLY
    IMPLICIT NONE

    LOGICAL :: check
	INTEGER :: i,j,k,count,polynums(50),triangle,checkele,P_ele
    DOUBLE PRECISION :: dis
    DOUBLE PRECISION :: r_ele_x(4,rho_elements),r_ele_y(4,rho_elements)
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: poly
  
    write(*,*) 'find polygons in elements'

    CALL getR_ele(r_ele_x,r_ele_y)

    polyspecs = 0

    !iterate through each habitat polygon and store information:
    !  polyspecs(n,1) = location in polys of first point of habitat polygon n
    !  polyspecs(n,2) = # of edge points that make up habitat polygon n
    count = 1
    polyspecs(INT(polys(1,1)),1) = 1
    do i=2,pedges
	  if(polys(i,1) /= polys(i-1,1))then
        polyspecs(INT(polys(i,1)),1) = i
        polyspecs(INT(polys(i-1,1)),2) = count
        count = 1
      elseif(i==pedges)then
        polyspecs(INT(polys(i,1)),2) = count + 1
      else
        count = count + 1
      endif
    enddo


    !iterate through every element determining which habitat polygons are inside them
	do i=1,rho_elements
      count = 0
	  polynums = 0
      do j=1,pedges
        !if the current polygon edge point has the same id as the last one added to the
        !  list skip this iteration (so same polygon isnt added multiple times)
        if(count.NE.0 .AND. polys(j,1).EQ.polynums(count)) cycle

        !check if each habitat polygon edge point is in the element
        triangle = 0
        checkele = i
        CALL gridcell(rho_elements,r_ele_y,r_ele_x,polys(j,4),polys(j,5),P_ele,triangle,checkele)
        if(triangle /= 0) then
          count = count + 1
          polynums(count) = polys(j,1)
		  cycle
        endif

        !if none of the current habitat polygon edge points are in the element
        !  check to make sure none of the element edge points are in the habitat polygon
        if(j==pedges .OR. polys(j,1)/=polys(j+1,1))then
          check = .FALSE.
          !check if any of the element edge points are in range of the habitat polygon
          do k=1,4
            dis = sqrt( (r_ele_x(k,i)-polys(j,2))**2 + (r_ele_y(k,i)-polys(j,3))**2 )
            if(dis<maxbdis(INT(polys(j,1))))check = .TRUE.
          enddo

          !if the element edge points are in range:
          if(check)then
            !allocate poly to contain the current habitat polygon edge points
            ALLOCATE(poly(polyspecs(INT(polys(j,1)),2),2))
            do k=1,polyspecs(INT(polys(j,1)),2)
              poly(k,1) = polys(polyspecs(INT(polys(j,1)),1)+k-1,4)
              poly(k,2) = polys(polyspecs(INT(polys(j,1)),1)+k-1,5)
            enddo
            !and check if any of the four element edge points are in the polygon
            do k=1,4
              if(INPOLY(r_ele_x(k,i),r_ele_y(k,i),polyspecs(INT(polys(j,1)),2),poly))then
                !if one is inside the polygon, add the polygons id to polynums
                count = count + 1
                polynums(count) = polys(j,1)
                exit
              endif
            enddo
            !make sure poly is deallocated
            DEALLOCATE(poly)
          endif
        endif          

      enddo

      !if there were any polygons inside this element, transfer that information
      !  from polynums to elepolys
      elepolys(i)%numpoly = count
      if(count>0)then
        ALLOCATE(elepolys(i)%poly(count))
        do j=1,count
          elepolys(i)%poly(j) = polynums(j)
        enddo
      endif
    enddo


    if(holesExist)then

      !iterate through each hole and store information:
      !  holespecs(n,1) = location in holes of first point of hole n
      !  holespecs(n,2) = # of edge points that make up hole n
      count = 1
      holespecs(INT(Holes(1,1)),1) = 1
      do i=2,hedges
        if(Holes(i,1) /= Holes(i-1,1))then
          holespecs(INT(Holes(i,1)),1) = i
          holespecs(INT(Holes(i-1,1)),2) = count
          count = 1
        elseif(i==hedges)then
          holespecs(INT(Holes(i,1)),2) = count + 1
        else
          count = count + 1
        endif
      enddo


      !iterate through every habitat polygon determining which holes are inside them
      do i=minpolyid,maxpolyid
        if(polyspecs(i,1)==0)cycle

        count = 0
        polynums = 0

        !iterate through all the hole edge points
        do j=1,hedges
          if(j == 1 .OR. holes(j,1) /= holes(j-1,1)) then
            if(holes(j,6)==polys(polyspecs(i,1),1))then
              !if the hole is in the current habitat polygon, add it to polynums
              count = count + 1
              polynums(count) = holes(j,1)
              cycle
            endif
          endif
        enddo

        !if there were any holes in the current habitat polygon, transfer that
        !  information from polynums to polyholes
        polyholes(i)%numpoly = count
        if(count>0)then
          ALLOCATE(polyholes(i)%poly(count))
          do j=1,count
            polyholes(i)%poly(j) = polynums(j)
          enddo
        endif
      enddo

    endif

  END SUBROUTINE createPolySpecs


  !Subroutine to determine if the particle is on any oyster polys (including holes)
  SUBROUTINE settlement(P_age,n,Px,Py,inpoly)
    USE PARAM_MOD, ONLY: holesExist
    USE HYDRO_MOD, ONLY: getP_r_element
	IMPLICIT NONE

	INTEGER, INTENT(IN) :: n
	INTEGER, INTENT(OUT) :: inpoly
	DOUBLE PRECISION, INTENT(IN) :: P_age,Px,Py

	INTEGER :: polyin,R_ele

    R_ele = getP_r_element(n)

	polyin = 0
	inpoly = 0

    if(P_age >= settletime(n))then

      !check if the particle is within the boundaries of any habitat polygon
      CALL psettle(Px,Py,R_ele,polyin)
      !if within a habitat polygon, check if particle is within the boundaries
      !	of any hole that is in that particular habitat polygon
      if (polyin .GT. 0) then
        inpoly = polyin				!if its within a habitat polygon set inpoly to the polygon id #
        if(holesExist)then
          CALL hsettle(Px,Py,polyin)
          if (polyin /= 0) inpoly = 0	!if its within a hole, reset inpoly to 0
        endif
      endif

      if (inpoly .GT. 0) then
        settle(n) = 1
      endif

    endif

  END SUBROUTINE settlement


  SUBROUTINE psettle(Px,Py,R_ele,polyin)
    USE PIP_MOD, ONLY: INPOLY
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: R_ele
    INTEGER, INTENT(OUT) :: polyin
    DOUBLE PRECISION, INTENT(IN) :: Px,Py

    INTEGER :: i,j,start,size
    DOUBLE PRECISION :: dis
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: polybnds

    polyin = 0               !initialize polyin to 0

    !if there are any habitat polygons in the element that the particle is in:
    if(elepolys(R_ele)%numpoly > 0)then

      !iterate through all the habitat polygons in that element
      do i=1,elepolys(R_ele)%numpoly
        start = polyspecs(elepolys(R_ele)%poly(i),1)
        size =  polyspecs(elepolys(R_ele)%poly(i),2)

        !if the particle is not within range of the habitat polygon, skip this habitat polygon
        dis = sqrt( (Px-polys(start,2))**2 + (Py-polys(start,3))**2 )
        if(dis>maxbdis(INT(polys(start,1))))cycle

        !allocate polybnds and fill it with the boundary point locations of the current habitat polygon
        ALLOCATE(polybnds(size,2))
        do j=1,size
          polybnds(j,1) = polys(start + j - 1, 4)
          polybnds(j,2) = polys(start + j - 1, 5)
        enddo

        if(INPOLY(Px,Py,size,polybnds))then			!  call INPOLY to see if the point is in the polygon
          polyin = polys(start,1)						!    if so set polyin to the id of the polygon it is in
          DEALLOCATE(polybnds)						!    deallocate polybnds
          exit										!    and exit the subroutine
        endif

        !make sure polybnds is deallocated
        DEALLOCATE(polybnds)
      enddo
    endif

  END SUBROUTINE psettle


  SUBROUTINE hsettle(Px,Py,holein)
    USE PIP_MOD, ONLY: INPOLY
    IMPLICIT NONE

    INTEGER, INTENT(INOUT) :: holein
    DOUBLE PRECISION, INTENT(IN) :: Px,Py

    INTEGER :: i,j,start,size,polyin
    DOUBLE PRECISION :: dis
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: polybnds

    polyin = holein          !initialize polyin to holein
    holein = 0               !initialize holein to 0

    !if there are any holes in the habitat polygon that the particle is inside:
    if(polyholes(polyin)%numpoly > 0)then

      !iterate through each hole
      do i=1,polyholes(polyin)%numpoly
        start = holespecs(polyholes(polyin)%poly(i),1)
        size =  holespecs(polyholes(polyin)%poly(i),2)

        !if the particle is not within range of the hole, skip this hole
        dis = sqrt( (Px-Holes(start,2))**2 + (Py-Holes(start,3))**2 )
        if(dis>maxhdis(INT(Holes(start,1))))cycle

        !allocate polybnds and fill it with the boundary point locations of the current hole
        ALLOCATE(polybnds(size,2))
        do j=1,size
          polybnds(j,1) = Holes(start + j - 1, 4)
          polybnds(j,2) = Holes(start + j - 1, 5)
        enddo

        if(INPOLY(Px,Py,size,polybnds,.FALSE.))then	!  call INPOLY to see if the point is in the hole
		!  NOTICE: onin is set .FALSE. meaning a 
		!  particle on the edge of a hole in a 
		!  habitat polygon is not considered to be 
		!  in the hole and thus will still settle
          holein = Holes(start,1)					!    if so set holein to the id of the hole it is in
          DEALLOCATE(polybnds)						!    deallocate polybnds
          exit										!    and exit the subroutine
        endif

        !make sure polybnds is deallocated
        DEALLOCATE(polybnds)
      enddo
    endif

  END SUBROUTINE hsettle


  LOGICAL FUNCTION SETTLED(n)
  !This function returns .TRUE. if the particle has "settled", and FALSE if not
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n

    SETTLED = .FALSE.
    IF(settle(n)==1)SETTLED = .TRUE.
  END FUNCTION SETTLED


  LOGICAL FUNCTION DEAD(n)
  !This function returns .TRUE. if the particle is "dead", and FALSE if not
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n

    DEAD = .FALSE.
    IF(settle(n)==2)DEAD = .TRUE.
  END FUNCTION DEAD


  SUBROUTINE DIE(n)
  !This subroutine sets the value of settle(n) to 2, indicating the particle is "dead"
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n

    settle(n) = 2
  END SUBROUTINE DIE


END MODULE SETTLEMENT_MOD