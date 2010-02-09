MODULE HYDRO_MOD

!  This module handles all the input from the hydrodynamic NetCDF input files.
!  It is the only module that interacts with NetCDF input files.  It contains
!  all the variables read in from the NetCDF files.  It also contains all the
!  information and variables related to the grid elements.
!
!  Created by:            Zachary Schlag			
!  Created on:			  07 Aug 2008
!  Last Modified on:	  19 Aug 2008

  USE PARAM_MOD, ONLY: numpar,ui,vi,uj,vj,us,ws,tdim,rho_nodes,u_nodes,v_nodes,max_rho_elements,max_u_elements,max_v_elements,rho_elements,u_elements,v_elements
  IMPLICIT NONE
  PRIVATE
  SAVE

  INTEGER :: iint, & !Keeps track of the input file, 0 = file 1, 1 = file 2, etc.
             stepf   !Keeps track of the forward time step

  !Used for reading in NetCDF variables one time step at a time
  INTEGER :: STARTr(4),COUNTr(4),STARTz(3),COUNTz(3)

  !Keeps track of the Rho, U, and V element that each particle is in
  INTEGER :: P_r_element(numpar),P_u_element(numpar),P_v_element(numpar)

  !The Rho, U, and V nodes that make up the Rho, U, and V element the particle is in
  INTEGER :: rnode1,rnode2,rnode3,rnode4,unode1,unode2,unode3,unode4,vnode1,vnode2,vnode3,vnode4

  !These variables keep track of the interpolation method and weights
  INTEGER :: tOK
  DOUBLE PRECISION :: t,u,Wgt1,Wgt2,Wgt3,Wgt4

  !S-Level location variables
  REAL :: SC(us),CS(us),SCW(ws),CSW(ws)

  !Depth at each rho node location
  REAL :: depth(rho_nodes)
  
  !read in zeta,salinity,temperature,vertical diffusivity, and U,V,W velocities at hydrodynamic back, center, and forward time
  REAL :: zetab(rho_nodes),zetac(rho_nodes),zetaf(rho_nodes)
  REAL :: saltb(rho_nodes,us),saltc(rho_nodes,us),saltf(rho_nodes,us),tempb(rho_nodes,us),tempc(rho_nodes,us),tempf(rho_nodes,us)
  REAL :: Wvelb(rho_nodes,ws),Wvelc(rho_nodes,ws),Wvelf(rho_nodes,ws),Khb(rho_nodes,ws),Khc(rho_nodes,ws),Khf(rho_nodes,ws)
  REAL :: Uvelb(u_nodes,us),Uvelc(u_nodes,us),Uvelf(u_nodes,us),Vvelb(v_nodes,us),Vvelc(v_nodes,us),Vvelf(v_nodes,us)

  !Rho, U, and V grid wet elements(the four node numbers that make up the element)(wet means at least 1 node is masked as water)
  INTEGER :: RE(4,rho_elements),UE(4,u_elements),VE(4,v_elements)

  !For each element, a list containing itself and all the elements that share a node with that element,
  !  used to speed up determining which element the particle has moved to, if it has moved at all
  INTEGER :: r_Adjacent(rho_elements,10),u_Adjacent(u_elements,10),v_Adjacent(v_elements,10)

  !X/Y location of all the Rho,U,V grid nodes, and the angle between x-coordinate and true east direction (radian)
  DOUBLE PRECISION :: rho_angle(rho_nodes),rx(rho_nodes),ry(rho_nodes),ux(u_nodes),uy(u_nodes),vx(v_nodes),vy(v_nodes)
  DOUBLE PRECISION :: r_ele_x(4,rho_elements),r_ele_y(4,rho_elements),u_ele_x(4,u_elements),u_ele_y(4,u_elements),v_ele_x(4,v_elements),v_ele_y(4,v_elements)

  !U, and V grid metric node locations, and Rho grid masking
  !  These variables are shared with boundary_module.f90 to create the model boundaries
  REAL :: mask_rho(vi,uj),x_u(ui,uj),y_u(ui,uj),x_v(vi,vj),y_v(vi,vj)

  !Keeps track if the grid has been read in yet or not
  !  If the grid hasn't been read in, the boundaries can't be made
  LOGICAL :: GRD_SET = .FALSE.

  !The concatenated hydrodynamic input file name
  CHARACTER(len=100) :: filenm

  !The following procedures have been made public for the use of other program units:
  PUBLIC :: initGrid,initHydro,updateHydro,setEle,setInterp,getInterp,interp,WCTS_ITPI,getSlevel,getWlevel,getMask_Rho,getUVxy,getR_ele,getP_r_element

CONTAINS

  SUBROUTINE initGrid()
    !This subroutine reads in the grid information and with it creates all the element variables
    USE PARAM_MOD, ONLY: NCgridfile,prefix,suffix,filenum
    USE CONVERT_MOD,    ONLY: LON2X,LAT2Y,X2LON,Y2LAT
    USE netcdf
    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

    INTEGER :: STATUS,NCID,VID

    INTEGER :: i,j,m,count,inele
    REAL :: romdepth(vi,uj),mask_u(ui,uj),mask_v(vi,vj),x_rho(vi,uj),y_rho(vi,uj),angle(vi,uj)
    REAL :: lat_rho(vi,uj),lon_rho(vi,uj),lat_u(ui,uj),lon_u(ui,uj),lat_v(vi,vj),lon_v(vi,vj)
    INTEGER :: r_ele(4,max_rho_elements),u_ele(4,max_u_elements),v_ele(4,max_v_elements),     &
      rho_mask(rho_nodes),u_mask(u_nodes),v_mask(v_nodes)

    WRITE(*,*) 'read-in grid information'
    WRITE(*,*) NCgridfile

    ! ******************************************* READ IN GRID INFO *******************************************

    STATUS = NF90_OPEN(NCgridfile,NF90_NOWRITE, NCID)
    if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem NF90_OPEN'

      ! Depth (m)
      STATUS = NF90_INQ_VARID(NCID,'h',VID)
      STATUS = NF90_GET_VAR(NCID,VID,romdepth)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read depth'

      ! x-coordinate at rho (m)
      STATUS = NF90_INQ_VARID(NCID,'lon_rho',VID)
      STATUS = NF90_GET_VAR(NCID,VID,lon_rho)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read lon_rho'
      do i = 1,vi
        do j = 1,uj
          x_rho(i,j) = LON2X(lon_rho(i,j))
        enddo
      enddo
      

      ! y-coordinate at rho (m)
      STATUS = NF90_INQ_VARID(NCID,'lat_rho',VID)
      STATUS = NF90_GET_VAR(NCID,VID,lat_rho)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read lat_rho'
      do i = 1,vi
        do j = 1,uj
          y_rho(i,j) = LAT2Y(lat_rho(i,j))
        enddo
      enddo

      ! x-coordinate at u (m)
      STATUS = NF90_INQ_VARID(NCID,'lon_u',VID)
      STATUS = NF90_GET_VAR(NCID,VID,lon_u)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read lon_u'
      do i = 1,ui
        do j = 1,uj
          x_u(i,j) = LON2X(lon_u(i,j))
        enddo
      enddo

      ! y-coordinate at u (m)
      STATUS = NF90_INQ_VARID(NCID,'lat_u',VID)
      STATUS = NF90_GET_VAR(NCID,VID,lat_u)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read lat_u'
      do i = 1,ui
        do j = 1,uj
          y_u(i,j) = LAT2Y(lat_u(i,j))
        enddo
      enddo

      ! x-coordinate at v (m)
      STATUS = NF90_INQ_VARID(NCID,'lon_v',VID)
      STATUS = NF90_GET_VAR(NCID,VID,lon_v)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read lon_v'
      do i = 1,vi
        do j = 1,vj
          x_v(i,j) = LON2X(lon_v(i,j))
        enddo
      enddo

      ! y-coordinate at v (m)
      STATUS = NF90_INQ_VARID(NCID,'lat_v',VID)
      STATUS = NF90_GET_VAR(NCID,VID,lat_v)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read lat_v'
      do i = 1,vi
        do j = 1,vj
          y_v(i,j) = LAT2Y(lat_v(i,j))
        enddo
      enddo
      
!     write out the lat lon of the four corners of the grid for debugging purposes
      
      write(*,*) 'ul: ', x_rho(1,uj), ' ', y_rho(1,uj)
      write(*,*) 'll: ', x_rho(1,1), ' ', y_rho(1,1)
      write(*,*) 'ur: ', x_rho(vi,uj), ' ', y_rho(vi,uj)
      write(*,*) 'lr: ', x_rho(vi,1), ' ', y_rho(vi,1)
      
      ! mask on rho grid
      STATUS = NF90_INQ_VARID(NCID,'mask_rho',VID)
      STATUS = NF90_GET_VAR(NCID,VID,mask_rho)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read mask_rho'

      ! mask on u grid
      STATUS = NF90_INQ_VARID(NCID,'mask_u',VID)
      STATUS = NF90_GET_VAR(NCID,VID,mask_u)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read mask_u'

      ! mask on v grid
      STATUS = NF90_INQ_VARID(NCID,'mask_v',VID)
      STATUS = NF90_GET_VAR(NCID,VID,mask_v)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read mask_v'

      ! angle between x-coordinate and true east direction (radian)
      STATUS = NF90_INQ_VARID(NCID,'angle',VID)
      STATUS = NF90_GET_VAR(NCID,VID,angle)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read angle'

    STATUS = NF90_CLOSE(NCID)


    ! ******************************************* READ IN S LEVEL INFO *******************************************

    WRITE(filenm,'(A,I4.4,A)') prefix,filenum,suffix
!    write(*,*) filenm(1:LEN_TRIM(filenm))

    STATUS = NF90_OPEN(filenm, NF90_NOWRITE, NCID)
    if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem NF90_OPEN'

      ! s-coordinate on rho grid (sc_r)
      STATUS = NF90_INQ_VARID(NCID,'s_rho',VID)
      STATUS = NF90_GET_VAR(NCID,VID,SC)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read SC'

      ! Cs value on rho grid (Cs_r)
      STATUS = NF90_INQ_VARID(NCID,'Cs_r',VID)
      STATUS = NF90_GET_VAR(NCID,VID,CS)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read CS'

      ! s-coordinate on w grid (sc_w)
      STATUS = NF90_INQ_VARID(NCID,'s_w',VID)
      STATUS = NF90_GET_VAR(NCID,VID,SCW)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read SCW'

      ! Cs value on w grid (Cs_w)
      STATUS = NF90_INQ_VARID(NCID,'Cs_w',VID)
      STATUS = NF90_GET_VAR(NCID,VID,CSW)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read CSW'

    !close the dataset and reassign the NCID
    STATUS = NF90_CLOSE(NCID)


    ! ******************************************* CREATE ELEMENTS *******************************************


	! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	! ~  4B. Prepare Elements (i.e., assign ID numbers to rectangular grids)  ~
	! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	write(*,*) 'create elements'

	! Assign mask values to rho nodes 
	count = 0
	do j=1,uj
	  do i=1,vi
		count = count + 1	!move to next node number
			!cycles through each variable replacing the vi,uj part with count
			!  essentially giving it node numbers
		rho_mask(count) = mask_rho(i,j)
		rho_angle(count) = angle(i,j)
	  enddo
	enddo

	! Assign mask values to u nodes
	count = 0
	do j=1,uj
	  do i=1,ui
		count = count + 1	!move to next node number
			!cycles through each variable replacing the vi,uj part with count
			!  essentially giving it node numbers
		u_mask(count) = mask_u(i,j)
	  enddo
	enddo

	! Assign mask values to v nodes
	count = 0
	do j=1,vj
	  do i=1,vi
		count = count + 1	!move to next node number
			!cycles through each variable replacing the vi,uj part with count
			!  essentially giving it node numbers
		v_mask(count) = mask_v(i,j)
	  enddo
	enddo


	! Create matrix that contains the node numbers for each rho element
	count = 0
	do j=1,uj-1							!z2v3.2
	  do i=1,vi-1
		count = count + 1
		r_ele(1,count) = i + (j-1)*vi
		r_ele(2,count) = i + 1 + (j-1)*vi
		r_ele(3,count) = i + 1 + j*vi
		r_ele(4,count) = i + j*vi
	  enddo
	enddo

	! Create matrix that contains the node numbers for each u element
	count = 0
	do j=1,uj-1							!z2v3.2
	  do i=1,ui-1
		count = count + 1
		u_ele(1,count) = i + (j-1)*ui
		u_ele(2,count) = i + 1 + (j-1)*ui
		u_ele(3,count) = i + 1 + j*ui
		u_ele(4,count) = i + j*ui
	  enddo
	enddo

	! Create matrix that contains the node numbers for each v element
	count = 0
	do j=1,vj-1							!z2v3.2
	  do i=1,vi-1
		count = count + 1
		v_ele(1,count) = i + (j-1)*vi
		v_ele(2,count) = i + 1 + (j-1)*vi
		v_ele(3,count) = i + 1 + j*vi
		v_ele(4,count) = i + j*vi
	  enddo
	enddo


	! Create matrix that contains only the rho elements that contain a node
	! whose mask value = 1 (i.e., it has at least one water point). 
	count = 0
	do i=1,max_rho_elements
	  inele = 0
	  !using the mask determine if any of the nodes for the current
	  !  element are inbounds, if so set inele to 1
	  if( rho_mask(r_ele(1,i)) .EQ. 1) inele=1
	  if( rho_mask(r_ele(2,i)) .EQ. 1) inele=1
	  if( rho_mask(r_ele(3,i)) .EQ. 1) inele=1
	  if( rho_mask(r_ele(4,i)) .EQ. 1) inele=1				!z2v3.2
	  !if inele = 1 then at least one of the three nodes for this element
	  !  are in bounds.
	  if( inele .EQ. 1 ) then
		count = count + 1
		!create array of elements that contain at least one node in bounds
		RE(1,count) = r_ele(1,i)
		RE(2,count) = r_ele(2,i)
		RE(3,count) = r_ele(3,i)
		RE(4,count) = r_ele(4,i)							!z2v3.2
	  endif
	enddo

	! Create matrix that contains only the u elements that contain a node
	! whose mask value = 1 (i.e., it has at least one water point). 
	count = 0
	do i=1,max_u_elements
	  inele = 0
	  !using the mask determine if any of the nodes for the current
	  !  element are inbounds, if so set inele to 1
	  if( u_mask(u_ele(1,i)) .EQ. 1) inele=1
	  if( u_mask(u_ele(2,i)) .EQ. 1) inele=1
	  if( u_mask(u_ele(3,i)) .EQ. 1) inele=1
	  if( u_mask(u_ele(4,i)) .EQ. 1) inele=1					!z2v3.2
	  !if inele = 1 then at least one of the three nodes for this element
	  !  are in bounds.
	  if( inele .EQ. 1 ) then
		count = count + 1
		!create array of elements that contain at least one node in bounds
		UE(1,count) = u_ele(1,i)
		UE(2,count) = u_ele(2,i)
		UE(3,count) = u_ele(3,i)
		UE(4,count) = u_ele(4,i)							!z2v3.2
	  endif 
	enddo

	! Create matrix that contains only the v elements that contain a node
	! whose mask value = 1 (i.e., it has at least one water point). 
	count = 0
	do i=1,max_v_elements
	  inele = 0
	  !using the mask determine if any of the nodes for the current
	  !  element are inbounds, if so set inele to 1
	  if( v_mask(v_ele(1,i)) .EQ. 1) inele=1
	  if( v_mask(v_ele(2,i)) .EQ. 1) inele=1
	  if( v_mask(v_ele(3,i)) .EQ. 1) inele=1
	  if( v_mask(v_ele(4,i)) .EQ. 1) inele=1					!z2v3.2
	  !if inele = 1 then at least one of the three nodes for this element
	  !  are in bounds.
	  if( inele .EQ. 1 ) then
		count = count + 1
		!create array of elements that contain at least one node in bounds
		VE(1,count) = v_ele(1,i)
		VE(2,count) = v_ele(2,i)
		VE(3,count) = v_ele(3,i)
		VE(4,count) = v_ele(4,i)							!z2v3.2
	  endif 
	enddo

	! Create matrices of  x/y for rho nodes and depth values in rho node number format 
	count = 0
	do j=1,uj
	  do i=1,vi
		count = count + 1	!move to next node number
		!cycles through each variable replacing the vi,uj part with count
		!  essentially giving it node numbers
		rx(count) = x_rho(i,j)
		ry(count) = y_rho(i,j)
		depth(count) = romdepth(i,j)
	  enddo
	enddo

	! Create matrices of x/y values for u nodes in u node number format 
	count = 0
	do j=1,uj
	  do i=1,ui
		count = count + 1	!move to next node number
		!cycles through each variable replacing the ui,uj part with count
		!  essentially giving it node numbers
		ux(count) = x_u(i,j)
		uy(count) = y_u(i,j)
	  enddo
	enddo

	! Create matrices of x/y values for v nodes in v node number format
	count = 0
	do j=1,vj
	  do i=1,vi
		count = count + 1	!move to next node number
		!cycles through each variable replacing the vi,vj part with count
		!  essentially giving it node numbers
		vx(count) = x_v(i,j)
		vy(count) = y_v(i,j)
	  enddo
	enddo


	! Create matrices of x/y node values for each rho, u, and v element
	do j=1,rho_elements
	  do i=1,4											!z2v3.3
		r_ele_x(i,j) = rx(RE(i,j))
		r_ele_y(i,j) = ry(RE(i,j))
	  enddo
	enddo

	do j=1,u_elements
	  do i=1,4											!z2v3.3
		u_ele_x(i,j) = ux(UE(i,j))
		u_ele_y(i,j) = uy(UE(i,j))
	  enddo
	enddo

	do j=1,v_elements
	  do i=1,4											!z2v3.3
		v_ele_x(i,j) = vx(VE(i,j))
		v_ele_y(i,j) = vy(VE(i,j))
	  enddo
	enddo


    ! **************************************** FIND ADJACENT ELEMENTS ****************************************

	! Create search restriction algorithms
	write(*,*) 'find adjacent elements'
! I. For each element, list all elements that are adjacent to it 
!$OMP PARALLEL
!$OMP SECTIONS
!$OMP SECTION
	do i=1,rho_elements							!z2v4.2 - start
	  r_Adjacent(i,1) = i
	  m=1
	  do j=1,rho_elements
		if(j.EQ.i) cycle
		if(  (RE(1,i) .EQ. RE(1,j)) .OR. (RE(1,i) .EQ. RE(2,j))											&
		.OR. (RE(1,i) .EQ. RE(3,j)) .OR. (RE(1,i) .EQ. RE(4,j))											&
		.OR. (RE(2,i) .EQ. RE(1,j)) .OR. (RE(2,i) .EQ. RE(2,j))											&
		.OR. (RE(2,i) .EQ. RE(3,j)) .OR. (RE(2,i) .EQ. RE(4,j))											&
		.OR. (RE(3,i) .EQ. RE(1,j)) .OR. (RE(3,i) .EQ. RE(2,j))											&
		.OR. (RE(3,i) .EQ. RE(3,j)) .OR. (RE(3,i) .EQ. RE(4,j))											&
		.OR. (RE(4,i) .EQ. RE(1,j)) .OR. (RE(4,i) .EQ. RE(2,j))											&
		.OR. (RE(4,i) .EQ. RE(3,j)) .OR. (RE(4,i) .EQ. RE(4,j))  )then
		  m=m+1
		  r_Adjacent(i,m) = j
		endif
	  enddo 
	enddo
!$OMP SECTION
	do i=1,u_elements
	  u_Adjacent(i,1) = i
	  m=1
	  do j=1,u_elements
		if(j.EQ.i) cycle
		if(  (UE(1,i) .EQ. UE(1,j)) .OR. (UE(1,i) .EQ. UE(2,j))											&
		.OR. (UE(1,i) .EQ. UE(3,j)) .OR. (UE(1,i) .EQ. UE(4,j))											&
		.OR. (UE(2,i) .EQ. UE(1,j)) .OR. (UE(2,i) .EQ. UE(2,j))											&
		.OR. (UE(2,i) .EQ. UE(3,j)) .OR. (UE(2,i) .EQ. UE(4,j))											&
		.OR. (UE(3,i) .EQ. UE(1,j)) .OR. (UE(3,i) .EQ. UE(2,j))											&
		.OR. (UE(3,i) .EQ. UE(3,j)) .OR. (UE(3,i) .EQ. UE(4,j))											&
		.OR. (UE(4,i) .EQ. UE(1,j)) .OR. (UE(4,i) .EQ. UE(2,j))											&
		.OR. (UE(4,i) .EQ. UE(3,j)) .OR. (UE(4,i) .EQ. UE(4,j))  )then
		  m=m+1
		  u_Adjacent(i,m) = j
		endif
	  enddo 
	enddo 
!$OMP SECTION
	do i=1,v_elements
	  v_Adjacent(i,1) = i
	  m=1
	  do j=1,v_elements
		if(j.EQ.i) cycle
		if(  (VE(1,i) .EQ. VE(1,j)) .OR. (VE(1,i) .EQ. VE(2,j))											&
		.OR. (VE(1,i) .EQ. VE(3,j)) .OR. (VE(1,i) .EQ. VE(4,j))											&
		.OR. (VE(2,i) .EQ. VE(1,j)) .OR. (VE(2,i) .EQ. VE(2,j))											&
		.OR. (VE(2,i) .EQ. VE(3,j)) .OR. (VE(2,i) .EQ. VE(4,j))											&
		.OR. (VE(3,i) .EQ. VE(1,j)) .OR. (VE(3,i) .EQ. VE(2,j))											&
		.OR. (VE(3,i) .EQ. VE(3,j)) .OR. (VE(3,i) .EQ. VE(4,j))											&
		.OR. (VE(4,i) .EQ. VE(1,j)) .OR. (VE(4,i) .EQ. VE(2,j))											&
		.OR. (VE(4,i) .EQ. VE(3,j)) .OR. (VE(4,i) .EQ. VE(4,j))  )then     
		  m=m+1
		  v_Adjacent(i,m) = j
		endif
	  enddo 
	enddo
!$OMP END SECTIONS
!$OMP END PARALLEL
!write(*,*) r_Adjacent
	GRD_SET = .TRUE.
    
  END SUBROUTINE initGrid



  SUBROUTINE initHydro()
    !This Subroutine reads in the hydrodynamic information for the first iteration
    USE PARAM_MOD, ONLY: prefix,suffix,filenum
    USE netcdf
    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

    INTEGER :: STATUS,NCID,VID

	INTEGER :: i,j,k,count,counter,stepb,stepc
	REAL :: romZb(vi,uj,1),romZc(vi,uj,1),romZf(vi,uj,1),romWb(vi,uj,ws,1),romWc(vi,uj,ws,1),romWf(vi,uj,ws,1),    &
      romSb(vi,uj,us,1),romSc(vi,uj,us,1),romSf(vi,uj,us,1),romTb(vi,uj,us,1),romTc(vi,uj,us,1),romTf(vi,uj,us,1), &
      romUb(ui,uj,us,1),romUc(ui,uj,us,1),romUf(ui,uj,us,1),romVb(vi,vj,us,1),romVc(vi,vj,us,1),romVf(vi,vj,us,1), &      
      romKHb(vi,uj,ws,1),romKHc(vi,uj,ws,1),romKHf(vi,uj,ws,1)


	!Open netCDF file
	iint = 0
	counter=iint+filenum	! 177 --> June 26,1995
	WRITE(filenm,'(A,I4.4,A)') prefix,counter,suffix
	write(*,*) filenm(1:LEN_TRIM(filenm))
	
	if (tdim >= 3) then
		stepb=1    !Back step is 1st time step of file
		stepc=2    !Center step is 2nd time step of file
		stepf=3    !Forward step is 3rd time step of file

		! Read in data for first three external time steps
		STATUS = NF90_OPEN(filenm, NF90_NOWRITE, NCID)
		if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem NF90_OPEN'
!		write(*,*) 'ncid = ', ncid

		startr(1)=1
		startr(2)=1
		startr(3)=1
      
		startz(1)=1
		startz(2)=1

		! **** Zeta ****
		countz(1)=vi
		countz(2)=uj
		countz(3)=1
		STATUS = NF90_INQ_VARID(NCID,'zeta',VID)

		startz(3)=stepb
		STATUS = NF90_GET_VAR(NCID,VID,romZb,STARTz,COUNTz)
		if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read zeta array b'

		startz(3)=stepc
		STATUS = NF90_GET_VAR(NCID,VID,romZc,STARTz,COUNTz)
		if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read zeta array c'

		startz(3)=stepf
		STATUS = NF90_GET_VAR(NCID,VID,romZf,STARTz,COUNTz)
		if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read zeta array f'

		! **** Salt ****
		countr(1)=vi
		countr(2)=uj
		countr(3)=us
		countr(4)=1
		STATUS = NF90_INQ_VARID(NCID,'salt',VID)

		startr(4)=stepb
		STATUS = NF90_GET_VAR(NCID,VID,romSb,STARTr,COUNTr)
		if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read salt array b'

		startr(4)=stepc
		STATUS = NF90_GET_VAR(NCID,VID,romSc,STARTr,COUNTr)
		if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read salt array c'

		startr(4)=stepf
		STATUS = NF90_GET_VAR(NCID,VID,romSf,STARTr,COUNTr)
		if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read salt array f'

		! **** Temp ****
		countr(1)=vi
		countr(2)=uj
		countr(3)=us
		countr(4)=1
		STATUS = NF90_INQ_VARID(NCID,'temp',VID)

		startr(4)=stepb
		STATUS = NF90_GET_VAR(NCID,VID,romTb,STARTr,COUNTr)
		if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read temp array b'

		startr(4)=stepc
		STATUS = NF90_GET_VAR(NCID,VID,romTc,STARTr,COUNTr)
		if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read temp array c'

		startr(4)=stepf
		STATUS = NF90_GET_VAR(NCID,VID,romTf,STARTr,COUNTr)
		if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read temp array f'

		! **** U velocity ****
		countr(1)=ui
		countr(2)=uj
		countr(3)=us
		countr(4)=1
		STATUS = NF90_INQ_VARID(NCID,'u',VID)

		startr(4)=stepb
		STATUS = NF90_GET_VAR(NCID,VID,romUb,STARTr,COUNTr)
		if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read u array b'

		startr(4)=stepc
		STATUS = NF90_GET_VAR(NCID,VID,romUc,STARTr,COUNTr)
		if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read u array c'

		startr(4)=stepf
		STATUS = NF90_GET_VAR(NCID,VID,romUf,STARTr,COUNTr)
		if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read u array f'

		! **** V velocity ****
		countr(1)=vi
		countr(2)=vj
		countr(3)=us
		countr(4)=1
		STATUS = NF90_INQ_VARID(NCID,'v',VID)

		startr(4)=stepb
		STATUS = NF90_GET_VAR(NCID,VID,romVb,STARTr,COUNTr)
		if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read v array b'

		startr(4)=stepc
		STATUS = NF90_GET_VAR(NCID,VID,romVc,STARTr,COUNTr)
		if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read v array c'

		startr(4)=stepf
		STATUS = NF90_GET_VAR(NCID,VID,romVf,STARTr,COUNTr)
		if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read v array f'
		
		! **** W velocity ****
		countr(1)=vi
		countr(2)=uj
		countr(3)=ws
		countr(4)=1
		STATUS = NF90_INQ_VARID(NCID,'w',VID)
		
		startr(4)=stepb
		STATUS = NF90_GET_VAR(NCID,VID,romWb,STARTr,COUNTr)
		if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read w array b'
		
		startr(4)=stepc
		STATUS = NF90_GET_VAR(NCID,VID,romWc,STARTr,COUNTr)
		if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read w array c'
		
		startr(4)=stepf
		STATUS = NF90_GET_VAR(NCID,VID,romWf,STARTr,COUNTr)
		if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read w array f'
		
		! **** Vertical diffusivity for salt (Aks) ****
		STATUS = NF90_INQ_VARID(NCID,'AKs',VID)
		
		startr(4)=stepb
		STATUS = NF90_GET_VAR(NCID,VID,romKHb,STARTr,COUNTr)
		if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read AKs array b'
		
		startr(4)=stepc
		STATUS = NF90_GET_VAR(NCID,VID,romKHc,STARTr,COUNTr)
		if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read AKs array c'
		
		startr(4)=stepf
		STATUS = NF90_GET_VAR(NCID,VID,romKHf,STARTr,COUNTr)
		if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read AKs array f'
		
		!close the dataset and reassign the NCID
		STATUS = NF90_CLOSE(NCID)
	else
		startr(1)=1
		startr(2)=1
		startr(3)=1
		startr(4)=1
		
		startz(1)=1
		startz(2)=1
		startz(3)=1
		
		do j = 1, 3
			if (j > tdim + iint) then
				iint = iint + 1
				counter = iint + filenum
				write(filenm,'(A,I4.4,A)') prefix,counter,suffix
				write(*,*) filenm(1:LEN_TRIM(filenm))
				startr(4) = 1
				startz(3) = 1
			endif
			STATUS = NF90_OPEN(filenm, NF90_NOWRITE, NCID)
			if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem NF90_OPEN'
			
			select case(j)
				case (1)
					write(*,*) 'initHydro: read arrays at first time step'
					
					! **** zeta ****
					countz(1)=vi
					countz(2)=uj
					countz(3)=1
					
					STATUS = NF90_INQ_VARID(NCID,'zeta',VID)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem finding zeta array b'
					STATUS = NF90_GET_VAR(NCID,VID,romZb,STARTz,COUNTz)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem reading zeta array b'
					
					! **** salt ****
					countr(1)=vi
					countr(2)=uj
					countr(3)=us
					countr(4)=1
					
					STATUS = NF90_INQ_VARID(NCID,'salt',VID)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem finding salt array b'
					STATUS = NF90_GET_VAR(NCID,VID,romSb,STARTr,COUNTr)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem reading salt array b'
					
					! **** Temp ****
					countr(1)=vi
					countr(2)=uj
					countr(3)=us
					countr(4)=1
					
					STATUS = NF90_INQ_VARID(NCID,'temp',VID)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem finding temp array b'
					STATUS = NF90_GET_VAR(NCID,VID,romTb,STARTr,COUNTr)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read temp array b'
					
					
					! **** U velocity ****
					countr(1)=ui
					countr(2)=uj
					countr(3)=us
					countr(4)=1
					STATUS = NF90_INQ_VARID(NCID,'u',VID)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem finding u array b'
					STATUS = NF90_GET_VAR(NCID,VID,romUb,STARTr,COUNTr)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read u array b'
					
					! **** V velocity ****
					countr(1)=vi
					countr(2)=vj
					countr(3)=us
					countr(4)=1
					STATUS = NF90_INQ_VARID(NCID,'v',VID)
					STATUS = NF90_GET_VAR(NCID,VID,romVb,STARTr,COUNTr)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read v array b'
					
					! **** W velocity ****
					countr(1)=vi
					countr(2)=uj
					countr(3)=ws
					countr(4)=1
					STATUS = NF90_INQ_VARID(NCID,'w',VID)
					STATUS = NF90_GET_VAR(NCID,VID,romWb,STARTr,COUNTr)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read w array b'
					
					! **** Vertical diffusivity for salt (Aks) ****
					STATUS = NF90_INQ_VARID(NCID,'AKs',VID)
					STATUS = NF90_GET_VAR(NCID,VID,romKHb,STARTr,COUNTr)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read AKs array b'
					
				case (2)
					write(*,*) 'initHydro: read arrays at second time step'
					! **** zeta ****
					countz(1)=vi
					countz(2)=uj
					countz(3)=1
					
					STATUS = NF90_INQ_VARID(NCID,'zeta',VID)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem finding zeta array c'
					STATUS = NF90_GET_VAR(NCID,VID,romZc,STARTz,COUNTz)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem reading zeta array c'
!					if (STATUS .NE. NF90_NOERR) write(*,*) STATUS, startz, countz 
					
					! **** salt ****
					countr(1)=vi
					countr(2)=uj
					countr(3)=us
					countr(4)=1
					
					STATUS = NF90_INQ_VARID(NCID,'salt',VID)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem finding salt array c'
					STATUS = NF90_GET_VAR(NCID,VID,romSc,STARTr,COUNTr)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem reading salt array c'
					
					! **** Temp ****
					countr(1)=vi
					countr(2)=uj
					countr(3)=us
					countr(4)=1
					
					STATUS = NF90_INQ_VARID(NCID,'temp',VID)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem finding temp array c'
					STATUS = NF90_GET_VAR(NCID,VID,romTc,STARTr,COUNTr)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read temp array c'
					
					! **** U velocity ****
					countr(1)=ui
					countr(2)=uj
					countr(3)=us
					countr(4)=1
					STATUS = NF90_INQ_VARID(NCID,'u',VID)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem finding u array c'
					STATUS = NF90_GET_VAR(NCID,VID,romUc,STARTr,COUNTr)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read u array c'
					
					! **** V velocity ****
					countr(1)=vi
					countr(2)=vj
					countr(3)=us
					countr(4)=1
					STATUS = NF90_INQ_VARID(NCID,'v',VID)
					STATUS = NF90_GET_VAR(NCID,VID,romVc,STARTr,COUNTr)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read v array c'
					
					! **** W velocity ****
					countr(1)=vi
					countr(2)=uj
					countr(3)=ws
					countr(4)=1
					STATUS = NF90_INQ_VARID(NCID,'w',VID)
					STATUS = NF90_GET_VAR(NCID,VID,romWc,STARTr,COUNTr)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read w array c'
					
					! **** Vertical diffusivity for salt (Aks) ****
					STATUS = NF90_INQ_VARID(NCID,'AKs',VID)
					STATUS = NF90_GET_VAR(NCID,VID,romKHc,STARTr,COUNTr)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read AKs array c'
					
				case (3)
					write(*,*) 'initHydro: read arrays at third time step'
					! **** zeta ****
					countz(1)=vi
					countz(2)=uj
					countz(3)=1
					
					STATUS = NF90_INQ_VARID(NCID,'zeta',VID)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem finding zeta array f'
					STATUS = NF90_GET_VAR(NCID,VID,romZf,STARTz,COUNTz)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem reading zeta array f'
					
					! **** salt ****
					countr(1)=vi
					countr(2)=uj
					countr(3)=us
					countr(4)=1
					
					STATUS = NF90_INQ_VARID(NCID,'salt',VID)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem finding salt array f'
					STATUS = NF90_GET_VAR(NCID,VID,romSf,STARTr,COUNTr)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem reading salt array f'
					
					! **** Temp ****
					countr(1)=vi
					countr(2)=uj
					countr(3)=us
					countr(4)=1
					
					STATUS = NF90_INQ_VARID(NCID,'temp',VID)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem finding temp array f'
					STATUS = NF90_GET_VAR(NCID,VID,romTf,STARTr,COUNTr)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read temp array f'
					
					! **** U velocity ****
					countr(1)=ui
					countr(2)=uj
					countr(3)=us
					countr(4)=1
					STATUS = NF90_INQ_VARID(NCID,'u',VID)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem finding u array f'
					STATUS = NF90_GET_VAR(NCID,VID,romUf,STARTr,COUNTr)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read u array f'
					
					! **** V velocity ****
					countr(1)=vi
					countr(2)=vj
					countr(3)=us
					countr(4)=1
					STATUS = NF90_INQ_VARID(NCID,'v',VID)
					STATUS = NF90_GET_VAR(NCID,VID,romVf,STARTr,COUNTr)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read v array f'
					
					! **** W velocity ****
					countr(1)=vi
					countr(2)=uj
					countr(3)=ws
					countr(4)=1
					STATUS = NF90_INQ_VARID(NCID,'w',VID)
					STATUS = NF90_GET_VAR(NCID,VID,romWf,STARTr,COUNTr)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read w array f'
					
					! **** Vertical diffusivity for salt (Aks) ****
					STATUS = NF90_INQ_VARID(NCID,'AKs',VID)
					STATUS = NF90_GET_VAR(NCID,VID,romKHf,STARTr,COUNTr)
					if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read AKs array f'
					
				case default
					write(*,*) 'initHydro: ERROR. This case should never happen'
			end select
			
			STATUS = NF90_CLOSE(NCID)
			if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem NF90_CLOSE'
			
			startr(4) = startr(4) + 1
			startz(3) = startz(3) + 1
			
		enddo
	endif

      !Reshape input to fit node numbers assigned to elements
      count = 0
      do j=1,uj
        do i=1,vi
          count = count + 1
          do k=1,us	
            saltb(count,k) = romSb(i,j,k,1)
            saltc(count,k) = romSc(i,j,k,1)
            saltf(count,k) = romSf(i,j,k,1)
            tempb(count,k) = romTb(i,j,k,1)
            tempc(count,k) = romTc(i,j,k,1)
            tempf(count,k) = romTf(i,j,k,1)
            Wvelb(count,k) = romWb(i,j,k,1)
            Wvelc(count,k) = romWc(i,j,k,1)
            Wvelf(count,k) = romWf(i,j,k,1)
            KHb(count,k) = romKHb(i,j,k,1)
            KHc(count,k) = romKHc(i,j,k,1)
            KHf(count,k) = romKHf(i,j,k,1)
          enddo
          Wvelb(count,ws) = romWb(i,j,ws,1)
          Wvelc(count,ws) = romWc(i,j,ws,1)
          Wvelf(count,ws) = romWf(i,j,ws,1)
          KHb(count,ws) = romKHb(i,j,ws,1)
          KHc(count,ws) = romKHc(i,j,ws,1)
          KHf(count,ws) = romKHf(i,j,ws,1)
          zetab(count) = romZb(i,j,1)
          zetac(count) = romZc(i,j,1)
          zetaf(count) = romZf(i,j,1)
        enddo
      enddo

      count = 0
      do j=1,uj
        do i=1,ui
          count = count + 1
          do k=1,us
            Uvelb(count,k) = romUb(i,j,k,1)
            Uvelc(count,k) = romUc(i,j,k,1)
            Uvelf(count,k) = romUf(i,j,k,1)
          enddo
        enddo
      enddo

      count = 0
      do j=1,vj
        do i=1,vi
          count = count + 1
          do k=1,us
            Vvelb(count,k) = romVb(i,j,k,1)
            Vvelc(count,k) = romVc(i,j,k,1)
            Vvelf(count,k) = romVf(i,j,k,1)
          enddo
        enddo		
      enddo

  END SUBROUTINE initHydro


  SUBROUTINE updateHydro()
    USE PARAM_MOD, ONLY: prefix,suffix,filenum
    USE netcdf
    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

    INTEGER :: STATUS,NCID,VID

	INTEGER :: i,j,k,count,counter
	REAL :: romZf(vi,uj,1),romSf(vi,uj,us,1),romTf(vi,uj,us,1), romUf(ui,uj,us,1),romVf(vi,vj,us,1),    &
      romWf(vi,uj,ws,1),romKHf(vi,uj,ws,1)

    !if the current input file is not yet finished, just increment stepf to the next time step
	IF (stepf .LT. tdim) THEN

	  stepf=stepf+1

	ELSE
    !if the current input file is finished, update filnm to next input file, and reset stepf to 1

	  !Open netCDF file
	  iint = iint+1
	  counter=iint+filenum	!176 + 1 = 177 --> June 26,1995
	  WRITE(filenm,'(A,I4.4,A)') prefix,counter,suffix
	  write(*,*) filenm(1:LEN_TRIM(filenm))

      stepf = 1

    ENDIF


	!Incremental loop that updates back and center data for the new timestep
	do i = 1,rho_nodes
	  do j = 1,us
		Wvelb(i,j) = Wvelc(i,j)
		Wvelc(i,j) = Wvelf(i,j)
		saltb(i,j) = saltc(i,j)
		saltc(i,j) = saltf(i,j)
		tempb(i,j) = tempc(i,j)
		tempc(i,j) = tempf(i,j)
		Khb(i,j) = Khc(i,j)
		Khc(i,j) = Khf(i,j)
	  enddo
	  Wvelb(i,ws) = Wvelc(i,ws)
	  Wvelc(i,ws) = Wvelf(i,ws)
	  Khb(i,ws) = Khc(i,ws)
	  Khc(i,ws) = Khf(i,ws)
	  zetab(i) = zetac(i)
	  zetac(i) = zetaf(i)
	enddo

	do i=1,u_nodes
	  do j=1,us
		Uvelb(i,j) = Uvelc(i,j)
		Uvelc(i,j) = Uvelf(i,j)
	  enddo
	enddo

	do i=1,v_nodes
	  do j=1,us
		Vvelb(i,j) = Vvelc(i,j)
		Vvelc(i,j) = Vvelf(i,j)
	  enddo
	enddo

	! Read in forward time step data 
	STATUS = NF90_OPEN(filenm, NF90_NOWRITE, NCID)
	if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem NF90_OPEN'
!	write(*,*) 'ncid = ', ncid

	! **** Zeta ****
	startz(3)=stepf
	STATUS = NF90_INQ_VARID(NCID,'zeta',VID)
	STATUS = NF90_GET_VAR(NCID,VID,romZf,STARTz,COUNTz)
	if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read zeta array 3'

	startr(4)=stepf
	countr(4)=1

	! **** U velocity ****
	countr(1)=ui
	countr(2)=uj
	countr(3)=us
	STATUS = NF90_INQ_VARID(NCID,'u',VID)
	STATUS = NF90_GET_VAR(NCID,VID,romUf,STARTr,COUNTr)
	if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read u array'

	! **** V velocity ****
	countr(1)=vi
	countr(2)=vj
	STATUS = NF90_INQ_VARID(NCID,'v',VID)
	STATUS = NF90_GET_VAR(NCID,VID,romVf,STARTr,COUNTr)
	if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read v array'


	! **** Salt ****
	countr(1)=vi
	countr(2)=uj
	STATUS = NF90_INQ_VARID(NCID,'salt',VID)
	STATUS = NF90_GET_VAR(NCID,VID,romSf,STARTr,COUNTr)
	if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read salt array'

	! **** Temp ****
	STATUS = NF90_INQ_VARID(NCID,'temp',VID)
	STATUS = NF90_GET_VAR(NCID,VID,romTf,STARTr,COUNTr)
	if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read temp array'


	! **** W velocity ****
	countr(3)=ws
	STATUS = NF90_INQ_VARID(NCID,'w',VID)
	STATUS = NF90_GET_VAR(NCID,VID,romWf,STARTr,COUNTr)
	if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read w array'

	! **** Vertical diffusivity for salt (Aks) ****
	STATUS = NF90_INQ_VARID(NCID,'AKs',VID)
	STATUS = NF90_GET_VAR(NCID,VID,romKHf,STARTr,COUNTr)
	if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read AKs array'

	!close the dataset and reassign the NCID
	STATUS = NF90_CLOSE(NCID)


    !Reshape input to fit node numbers assigned to elements
	count = 0
	do j=1,uj
	  do i=1,vi
		count = count + 1
		do k=1,us		
		  saltf(count,k) = romSf(i,j,k,1)
		  tempf(count,k) = romTf(i,j,k,1)
		  Wvelf(count,k) = romWf(i,j,k,1)
		  KHf(count,k) = romKHf(i,j,k,1)
		enddo
		Wvelf(count,ws) = romWf(i,j,ws,1)
		KHf(count,ws) = romKHf(i,j,ws,1)
		zetaf(count) = romZf(i,j,1)
	  enddo
	enddo

	count = 0
	do j=1,uj
	  do i=1,ui
		count = count + 1
		do k=1,us		
		  Uvelf(count,k) = romUf(i,j,k,1)
		enddo
	  enddo
	enddo

	count = 0
	do j=1,vj
	  do i=1,vi
		count = count + 1
		do k=1,us		
		  Vvelf(count,k) = romVf(i,j,k,1)
		enddo
	  enddo		
	enddo

	write(*,*) 'existing matrix,stepf=',stepf

  END SUBROUTINE updateHydro



  SUBROUTINE setEle(Xpar,Ypar,n,err,first)
    !This Subroutine determines which Rho, U, and V grid elements the particle is in
    USE GRIDCELL_MOD, ONLY: GRIDCELL
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: Xpar,Ypar
	INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(OUT), OPTIONAL :: err
    LOGICAL, INTENT(IN), OPTIONAL :: first

    LOGICAL :: fst
    INTEGER :: i,triangle,checkele,P_r_ele,P_u_ele,P_v_ele,oP_ele,P_ele,error
	CHARACTER :: anykey

    error = 0

    if( PRESENT(first) ) then
      fst = first
    else
      fst = .FALSE.
    endif

	if(fst) then !if the first iteration

      !Find rho element in which particle is located
	  P_r_ele=0
	  triangle=0
	  call gridcell(rho_elements,r_ele_y,r_ele_x,Xpar,Ypar,P_r_ele,triangle)
	  if (triangle.EQ.0) error = 1
	  P_r_element(n)=P_r_ele

      !Find u element in which particle is located
	  P_u_ele=0
	  triangle=0
	  call gridcell(u_elements,u_ele_y,u_ele_x,Xpar,Ypar,P_u_ele,triangle)
	  if (triangle.EQ.0) error = 2
	  P_u_element(n)=P_u_ele

      !Find v element in which particle is located
	  P_v_ele=0
	  triangle=0
	  call gridcell(v_elements,v_ele_y,v_ele_x,Xpar,Ypar,P_v_ele,triangle)
	  if (triangle.EQ.0) error = 3
	  P_v_element(n)=P_v_ele


	else !if not the first iteration 

	  !Find rho element in which particle is located
	  oP_ele = P_r_element(n)
	  do i=1,10
		if(r_Adjacent(oP_ele,i).EQ.0) error = 4

		triangle = 0
		checkele = r_Adjacent(oP_ele,i)
		call gridcell(rho_elements,r_ele_y,r_ele_x,Xpar,Ypar,P_ele,triangle,checkele)
	  	if(triangle .NE. 0) then
		  P_r_element(n) = P_ele
		  exit
		endif

	  enddo !r_singlecellloop


	  !Find u element in which particle is located
	  oP_ele = P_u_element(n)
	  do i=1,10
		if(u_Adjacent(oP_ele,i).EQ.0) error = 5

		triangle = 0
		checkele = u_Adjacent(oP_ele,i)
		call gridcell(u_elements,u_ele_y,u_ele_x,Xpar,Ypar,P_ele,triangle,checkele)
	  	if(triangle .NE. 0) then
		  P_u_element(n) = P_ele
		  exit
		endif

	  enddo


	  !Find v element in which particle is located
	  oP_ele = P_v_element(n)
	  do i=1,10
		if(v_Adjacent(oP_ele,i).EQ.0) error = 6

		triangle = 0
		checkele = v_Adjacent(oP_ele,i)
		call gridcell(v_elements,v_ele_y,v_ele_x,Xpar,Ypar,P_ele,triangle,checkele)
	  	if(triangle .NE. 0) then
		  P_v_element(n) = P_ele
		  exit
		endif

	  enddo

	endif

    !Assign node numbers for rho,u,v calculations
    rnode1 = RE(1,P_r_element(n))
    rnode2 = RE(2,P_r_element(n))
    rnode3 = RE(3,P_r_element(n))
    rnode4 = RE(4,P_r_element(n))

    unode1 = UE(1,P_u_element(n))
    unode2 = UE(2,P_u_element(n))
    unode3 = UE(3,P_u_element(n))
    unode4 = UE(4,P_u_element(n))

    vnode1 = VE(1,P_v_element(n))
    vnode2 = VE(2,P_v_element(n))
    vnode3 = VE(3,P_v_element(n))
    vnode4 = VE(4,P_v_element(n))

    if(PRESENT(err)) err = error

  END SUBROUTINE setEle



  SUBROUTINE setInterp(xp,yp,n)
    !This subroutine calculates and stores the interpolation method and values for the current particle
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    DOUBLE PRECISION, INTENT(IN) :: xp,yp

    double precision x1,x2,x3,x4,y1,y2,y3,y4
	double precision Dis1,Dis2,Dis3,Dis4,TDis

    x1 = rx(rnode1)
    x2 = rx(rnode2)
    x3 = rx(rnode3)
    x4 = rx(rnode4)
    y1 = ry(rnode1)
    y2 = ry(rnode2)
    y3 = ry(rnode3)
    y4 = ry(rnode4)

    Wgt1 = 0
    Wgt2 = 0
    Wgt3 = 0
    Wgt4 = 0

	tOK = 0 !to store information on the interpolation method

	! bilinear interpolation of first triangle
    t =  ( (xp-x1)*(y3-y1) + (y1-yp)*(x3-x1) ) / ( (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1) )
	u =  ( (xp-x1)*(y2-y1) + (y1-yp)*(x2-x1) ) / ( (x3-x1)*(y2-y1) - (y3-y1)*(x2-x1) )
	tOK = 1 !first triangle

	! if point outside triangle, then do bilinear interpolation of other triangle
	if( t.LT.0. .OR. u.LT.0. .OR. (t+u).GT.1.0 ) then
	  t = ( (xp-x3)*(y1-y3) + (y3-yp)*(x1-x3) ) / ( (x4-x3)*(y1-y3) - (y4-y3)*(x1-x3) )
	  u = ( (xp-x3)*(y4-y3) + (y3-yp)*(x4-x3) ) / ( (x1-x3)*(y4-y3) - (y1-y3)*(x4-x3) )
	  tOK = 2 !second triangle

	  !if bilinear techniques are undefined, then use inverse weighted distance
	  if( t.LT.0. .OR. u.LT.0. .OR. (t+u).GT.1.0 ) then
		!if particle on node, then set equal to node value
        if(  (xp.EQ.x1 .AND. yp.EQ.y1).OR.(xp.EQ.x2 .AND. yp.EQ.y2)									&
		 .OR.(xp.EQ.x3 .AND. yp.EQ.y3).OR.(xp.EQ.x4 .AND. yp.EQ.y4))then
		  if (xp.EQ.x1 .AND. yp.EQ.y1) Wgt1 = 1.0
		  if (xp.EQ.x2 .AND. yp.EQ.y2) Wgt2 = 1.0
		  if (xp.EQ.x3 .AND. yp.EQ.y3) Wgt3 = 1.0
		  if (xp.EQ.x4 .AND. yp.EQ.y4) Wgt4 = 1.0
		else !use inverse weighted distance instead  
		  Dis1=1./( SQRT( (x1-xp)**2 + (y1-yp)**2 ) ) 
		  Dis2=1./( SQRT( (x2-xp)**2 + (y2-yp)**2 ) ) 
		  Dis3=1./( SQRT( (x3-xp)**2 + (y3-yp)**2 ) ) 
		  Dis4=1./( SQRT( (x4-xp)**2 + (y4-yp)**2 ) ) 
		  TDis = Dis1+Dis2+Dis3+Dis4
		  Wgt1= Dis1/TDis
		  Wgt2= Dis2/TDis
		  Wgt3= Dis3/TDis
		  Wgt4= Dis4/TDis
		  tOK = 3 !no triangle - used inverse weighted distance
		endif
	  endif
	endif     

  END SUBROUTINE setInterp


  DOUBLE PRECISION FUNCTION getInterp(var,i)
    !This Function returns the interpolated value at the particle's location using
    !  the interpolation variables stored from function setInterp, and the hydrodynamic
    !  variables that have been read in
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: var
    INTEGER, INTENT(IN), OPTIONAL :: i

    DOUBLE PRECISION :: v1,v2,v3,v4

    !For Error State, & 'Press Any Key'
    CHARACTER :: anykey

    !Determine which data to interpolate from
    SELECT CASE(var)
      CASE("depth")
        v1 = depth(rnode1)
        v2 = depth(rnode2)
        v3 = depth(rnode3)
        v4 = depth(rnode4)
      CASE("angle")
        v1 = rho_angle(rnode1)
        v2 = rho_angle(rnode2)
        v3 = rho_angle(rnode3)
        v4 = rho_angle(rnode4)
      CASE("zetab")
        v1 = zetab(rnode1)
        v2 = zetab(rnode2)
        v3 = zetab(rnode3)
        v4 = zetab(rnode4)
      CASE("zetac")
        v1 = zetac(rnode1)
        v2 = zetac(rnode2)
        v3 = zetac(rnode3)
        v4 = zetac(rnode4)
      CASE("zetaf")
        v1 = zetaf(rnode1)
        v2 = zetaf(rnode2)
        v3 = zetaf(rnode3)
        v4 = zetaf(rnode4)
      CASE("saltb")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          write(*,*) 'Press Any Key and Enter'
          read(*,*) anykey
          stop
        endif
        v1 = saltb(rnode1,i)
        v2 = saltb(rnode2,i)
        v3 = saltb(rnode3,i)
        v4 = saltb(rnode4,i)
      CASE("saltc")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          write(*,*) 'Press Any Key and Enter'
          read(*,*) anykey
          stop
        endif
        v1 = saltc(rnode1,i)
        v2 = saltc(rnode2,i)
        v3 = saltc(rnode3,i)
        v4 = saltc(rnode4,i)
      CASE("saltf")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          write(*,*) 'Press Any Key and Enter'
          read(*,*) anykey
          stop
        endif
        v1 = saltf(rnode1,i)
        v2 = saltf(rnode2,i)
        v3 = saltf(rnode3,i)
        v4 = saltf(rnode4,i)
      CASE("tempb")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          write(*,*) 'Press Any Key and Enter'
          read(*,*) anykey
          stop
        endif
        v1 = tempb(rnode1,i)
        v2 = tempb(rnode2,i)
        v3 = tempb(rnode3,i)
        v4 = tempb(rnode4,i)
      CASE("tempc")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          write(*,*) 'Press Any Key and Enter'
          read(*,*) anykey
          stop
        endif
        v1 = tempc(rnode1,i)
        v2 = tempc(rnode2,i)
        v3 = tempc(rnode3,i)
        v4 = tempc(rnode4,i)
      CASE("tempf")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          write(*,*) 'Press Any Key and Enter'
          read(*,*) anykey
          stop
        endif
        v1 = tempf(rnode1,i)
        v2 = tempf(rnode2,i)
        v3 = tempf(rnode3,i)
        v4 = tempf(rnode4,i)
      CASE("khb")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          write(*,*) 'Press Any Key and Enter'
          read(*,*) anykey
          stop
        endif
        v1 = KHb(rnode1,i)
        v2 = KHb(rnode2,i)
        v3 = KHb(rnode3,i)
        v4 = KHb(rnode4,i)
      CASE("khc")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          write(*,*) 'Press Any Key and Enter'
          read(*,*) anykey
          stop
        endif
        v1 = KHc(rnode1,i)
        v2 = KHc(rnode2,i)
        v3 = KHc(rnode3,i)
        v4 = KHc(rnode4,i)
      CASE("khf")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          write(*,*) 'Press Any Key and Enter'
          read(*,*) anykey
          stop
        endif
        v1 = KHf(rnode1,i)
        v2 = KHf(rnode2,i)
        v3 = KHf(rnode3,i)
        v4 = KHf(rnode4,i)
      CASE DEFAULT
        write(*,*) 'Problem interpolating ',var
        write(*,*) ' '
        write(*,*) 'The Program Cannot Continue and Will Terminate'
        write(*,*) 'Press Any Key and Enter'
        read(*,*) anykey
        stop
    END SELECT

    !interpolate using the variables from setInterp
    if(tOK == 1) then 
      getInterp = v1 + (v2-v1)*t + (v3-v1)*u
    elseif(tOK == 2) then
      getInterp = v3 + (v4-v3)*t + (v1-v3)*u
    else 
      getInterp = Wgt1*v1 + Wgt2*v2 + Wgt3*v3 + Wgt4*v4 
    endif

  END FUNCTION getInterp


  DOUBLE PRECISION FUNCTION interp(xp,yp,var,i)
    !This Function determines the method of interpolation and returns the interpolated 
	!  value at the particle's location using the hydrodynamic variables that have been read in
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: xp,yp
    CHARACTER(LEN=*), INTENT(IN) :: var
    INTEGER, INTENT(IN), OPTIONAL :: i

    INTEGER :: RUV
  	DOUBLE PRECISION :: tt,uu,x1,x2,x3,x4,y1,y2,y3,y4,v1,v2,v3,v4,vp,     &
					    Dis1,Dis2,Dis3,Dis4,TDis,Wt1,Wt2,Wt3,Wt4

    !For Error State, & 'Press Any Key'
    CHARACTER :: anykey

    RUV = 1

    !determine which data to interpolate from
    SELECT CASE(var)
      CASE("depth")
        v1 = depth(rnode1)
        v2 = depth(rnode2)
        v3 = depth(rnode3)
        v4 = depth(rnode4)
      CASE("angle")
        v1 = rho_angle(rnode1)
        v2 = rho_angle(rnode2)
        v3 = rho_angle(rnode3)
        v4 = rho_angle(rnode4)
      CASE("zetab")
        v1 = zetab(rnode1)
        v2 = zetab(rnode2)
        v3 = zetab(rnode3)
        v4 = zetab(rnode4)
      CASE("zetac")
        v1 = zetac(rnode1)
        v2 = zetac(rnode2)
        v3 = zetac(rnode3)
        v4 = zetac(rnode4)
      CASE("zetaf")
        v1 = zetaf(rnode1)
        v2 = zetaf(rnode2)
        v3 = zetaf(rnode3)
        v4 = zetaf(rnode4)
      CASE("saltb")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          write(*,*) 'Press Any Key and Enter'
          read(*,*) anykey
          stop
        endif
        v1 = saltb(rnode1,i)
        v2 = saltb(rnode2,i)
        v3 = saltb(rnode3,i)
        v4 = saltb(rnode4,i)
      CASE("saltc")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          write(*,*) 'Press Any Key and Enter'
          read(*,*) anykey
          stop
        endif
        v1 = saltc(rnode1,i)
        v2 = saltc(rnode2,i)
        v3 = saltc(rnode3,i)
        v4 = saltc(rnode4,i)
      CASE("saltf")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          write(*,*) 'Press Any Key and Enter'
          read(*,*) anykey
          stop
        endif
        v1 = saltf(rnode1,i)
        v2 = saltf(rnode2,i)
        v3 = saltf(rnode3,i)
        v4 = saltf(rnode4,i)
      CASE("tempb")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          write(*,*) 'Press Any Key and Enter'
          read(*,*) anykey
          stop
        endif
        v1 = tempb(rnode1,i)
        v2 = tempb(rnode2,i)
        v3 = tempb(rnode3,i)
        v4 = tempb(rnode4,i)
      CASE("tempc")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          write(*,*) 'Press Any Key and Enter'
          read(*,*) anykey
          stop
        endif
        v1 = tempc(rnode1,i)
        v2 = tempc(rnode2,i)
        v3 = tempc(rnode3,i)
        v4 = tempc(rnode4,i)
      CASE("tempf")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          write(*,*) 'Press Any Key and Enter'
          read(*,*) anykey
          stop
        endif
        v1 = tempf(rnode1,i)
        v2 = tempf(rnode2,i)
        v3 = tempf(rnode3,i)
        v4 = tempf(rnode4,i)
      CASE("uvelb")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          write(*,*) 'Press Any Key and Enter'
          read(*,*) anykey
          stop
        endif
        v1 = Uvelb(unode1,i)
        v2 = Uvelb(unode2,i)
        v3 = Uvelb(unode3,i)
        v4 = Uvelb(unode4,i)
        RUV = 2
      CASE("uvelc")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          write(*,*) 'Press Any Key and Enter'
          read(*,*) anykey
          stop
        endif
        v1 = Uvelc(unode1,i)
        v2 = Uvelc(unode2,i)
        v3 = Uvelc(unode3,i)
        v4 = Uvelc(unode4,i)
        RUV = 2
      CASE("uvelf")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          write(*,*) 'Press Any Key and Enter'
          read(*,*) anykey
          stop
        endif
        v1 = Uvelf(unode1,i)
        v2 = Uvelf(unode2,i)
        v3 = Uvelf(unode3,i)
        v4 = Uvelf(unode4,i)
        RUV = 2
      CASE("vvelb")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          write(*,*) 'Press Any Key and Enter'
          read(*,*) anykey
          stop
        endif
        v1 = Vvelb(vnode1,i)
        v2 = Vvelb(vnode2,i)
        v3 = Vvelb(vnode3,i)
        v4 = Vvelb(vnode4,i)
        RUV = 3
      CASE("vvelc")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          write(*,*) 'Press Any Key and Enter'
          read(*,*) anykey
          stop
        endif
        v1 = Vvelc(vnode1,i)
        v2 = Vvelc(vnode2,i)
        v3 = Vvelc(vnode3,i)
        v4 = Vvelc(vnode4,i)
        RUV = 3
      CASE("vvelf")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          write(*,*) 'Press Any Key and Enter'
          read(*,*) anykey
          stop
        endif
        v1 = Vvelf(vnode1,i)
        v2 = Vvelf(vnode2,i)
        v3 = Vvelf(vnode3,i)
        v4 = Vvelf(vnode4,i)
        RUV = 3
      CASE("wvelb")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          write(*,*) 'Press Any Key and Enter'
          read(*,*) anykey
          stop
        endif
        v1 = Wvelb(rnode1,i)
        v2 = Wvelb(rnode2,i)
        v3 = Wvelb(rnode3,i)
        v4 = Wvelb(rnode4,i)
      CASE("wvelc")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          write(*,*) 'Press Any Key and Enter'
          read(*,*) anykey
          stop
        endif
        v1 = Wvelc(rnode1,i)
        v2 = Wvelc(rnode2,i)
        v3 = Wvelc(rnode3,i)
        v4 = Wvelc(rnode4,i)
      CASE("wvelf")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          write(*,*) 'Press Any Key and Enter'
          read(*,*) anykey
          stop
        endif
        v1 = Wvelf(rnode1,i)
        v2 = Wvelf(rnode2,i)
        v3 = Wvelf(rnode3,i)
        v4 = Wvelf(rnode4,i)
      CASE("khb")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          write(*,*) 'Press Any Key and Enter'
          read(*,*) anykey
          stop
        endif
        v1 = KHb(rnode1,i)
        v2 = KHb(rnode2,i)
        v3 = KHb(rnode3,i)
        v4 = KHb(rnode4,i)
      CASE("khc")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          write(*,*) 'Press Any Key and Enter'
          read(*,*) anykey
          stop
        endif
        v1 = KHc(rnode1,i)
        v2 = KHc(rnode2,i)
        v3 = KHc(rnode3,i)
        v4 = KHc(rnode4,i)
      CASE("khf")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          write(*,*) 'Press Any Key and Enter'
          read(*,*) anykey
          stop
        endif
        v1 = KHf(rnode1,i)
        v2 = KHf(rnode2,i)
        v3 = KHf(rnode3,i)
        v4 = KHf(rnode4,i)
      CASE DEFAULT
        write(*,*) 'Problem interpolating ',var
        write(*,*) ' '
        write(*,*) 'The Program Cannot Continue and Will Terminate'
        write(*,*) 'Press Any Key and Enter'
        read(*,*) anykey
        stop
    END SELECT

    !determine which node locations to interpolate from
    SELECT CASE(RUV)
      CASE(1)
        x1 = rx(rnode1)
        x2 = rx(rnode2)
        x3 = rx(rnode3)
        x4 = rx(rnode4)
        y1 = ry(rnode1)
        y2 = ry(rnode2)
        y3 = ry(rnode3)
        y4 = ry(rnode4)
      CASE(2)
        x1 = ux(unode1)
        x2 = ux(unode2)
        x3 = ux(unode3)
        x4 = ux(unode4)
        y1 = uy(unode1)
        y2 = uy(unode2)
        y3 = uy(unode3)
        y4 = uy(unode4)
      CASE(3)
        x1 = vx(vnode1)
        x2 = vx(vnode2)
        x3 = vx(vnode3)
        x4 = vx(vnode4)
        y1 = vy(vnode1)
        y2 = vy(vnode2)
        y3 = vy(vnode3)
        y4 = vy(vnode4)
    END SELECT

    !determine the method of interpolation:

	! bilinear interpolation of first triangle
    tt =  ( (xp-x1)*(y3-y1) + (y1-yp)*(x3-x1) ) / ( (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1) )
	uu =  ( (xp-x1)*(y2-y1) + (y1-yp)*(x2-x1) ) / ( (x3-x1)*(y2-y1) - (y3-y1)*(x2-x1) )
	vp = v1 + (v2-v1)*tt + (v3-v1)*uu

	! if point outside triangle, then do bilinear interpolation of other triangle
	if( tt.LT.0. .OR. uu.LT.0. .OR. (tt+uu).GT.1.0 ) then
	  tt = ( (xp-x3)*(y1-y3) + (y3-yp)*(x1-x3) ) / ( (x4-x3)*(y1-y3) - (y4-y3)*(x1-x3) )
	  uu = ( (xp-x3)*(y4-y3) + (y3-yp)*(x4-x3) ) / ( (x1-x3)*(y4-y3) - (y1-y3)*(x4-x3) )
	  vp = v3 + (v4-v3)*tt + (v1-v3)*uu

	  !if bilinear techniques are undefined, then use inverse weighted distance
	  if( tt.LT.0. .OR. uu.LT.0. .OR. (tt+uu).GT.1.0 ) then
		!if particle on node, then set equal to node value
        if(  (xp.EQ.x1 .AND. yp.EQ.y1).OR.(xp.EQ.x2 .AND. yp.EQ.y2)									&
		 .OR.(xp.EQ.x3 .AND. yp.EQ.y3).OR.(xp.EQ.x4 .AND. yp.EQ.y4))then
		  if (xp.EQ.x1 .AND. yp.EQ.y1) vp=v1
		  if (xp.EQ.x2 .AND. yp.EQ.y2) vp=v2
		  if (xp.EQ.x3 .AND. yp.EQ.y3) vp=v3
		  if (xp.EQ.x4 .AND. yp.EQ.y4) vp=v4
		else !use inverse weighted distance instead  
		  Dis1=1./( SQRT( (x1-xp)**2 + (y1-yp)**2 ) ) 
		  Dis2=1./( SQRT( (x2-xp)**2 + (y2-yp)**2 ) ) 
		  Dis3=1./( SQRT( (x3-xp)**2 + (y3-yp)**2 ) ) 
		  Dis4=1./( SQRT( (x4-xp)**2 + (y4-yp)**2 ) ) 
		  TDis = Dis1+Dis2+Dis3+Dis4
		  Wt1= Dis1/TDis
		  Wt2= Dis2/TDis
		  Wt3= Dis3/TDis
		  Wt4= Dis4/TDis
		  vp = Wt1*v1 + Wt2*v2 + Wt3*v3 + Wt4*v4 
		endif
	  endif
	endif

    interp = vp

  END FUNCTION Interp

  !This function creates a Water Column Tension Spline at back, center, and forward hydrodynamic time
  !  Then uses Polynomial Interpolation to determine Internal Time values to finally get the value
  !  of the particle in space and time.
  !Name derived from: Water Column Tension Spline, Internal Time Polynomial Interpolation
  !The final variable v is for version(ie what is to be returned): 1-back, 2-center, 3-forward, 4-(b+4c+f)/6
  DOUBLE PRECISION FUNCTION WCTS_ITPI(var,Xpos,Ypos,deplvl,Pwc_zb,Pwc_zc,Pwc_zf,slvls,P_zb,P_zc,P_zf,ex,ix,p,v)
    USE TENSION_MOD, ONLY: TSPSI,HVAL
    USE INT_MOD, ONLY: LININT,POLINTD
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: deplvl,slvls,p,v
    DOUBLE PRECISION, INTENT(IN) :: Xpos,Ypos,Pwc_zb(slvls),Pwc_zc(slvls),Pwc_zf(slvls),P_zb,P_zc,P_zf,ex(3),ix(3)
    CHARACTER(LEN=*), INTENT(IN) :: var

    INTEGER, PARAMETER :: nN = 4

    CHARACTER(LEN=LEN(var)+1) :: varb,varc,varf
    INTEGER :: i
    DOUBLE PRECISION :: abb_zb(nN),abb_zc(nN),abb_zf(nN),abb_vb(nN),abb_vc(nN),abb_vf(nN),       &
      P_vb,P_vc,P_vf,P_V,ey(3),vb,vc,vf,slope

    !TSPACK Variables
    INTEGER :: IER,SigErr
    DOUBLE PRECISION :: YP(nN),SIGM(nN)

    !For Error State, & 'Press Any Key'
    CHARACTER :: anykey
      
    varb = var//"b"
    varc = var//"c"
    varf = var//"f"

    do i=1,nN
      abb_zb(i) = Pwc_zb(i+deplvl-1)
      abb_zc(i) = Pwc_zc(i+deplvl-1)
      abb_zf(i) = Pwc_zf(i+deplvl-1)
      abb_vb(i) = interp(Xpos,Ypos,varb,i+deplvl-1)
      abb_vc(i) = interp(Xpos,Ypos,varc,i+deplvl-1)
      abb_vf(i) = interp(Xpos,Ypos,varf,i+deplvl-1)
    enddo

!		*********************************************************
!		*		5Aiic6b.  Fit Tension Spline to WC Profile		*
!		*********************************************************

    !ii. call TSPACK to fit a tension spline to water column profile 
    !  of U,V,W velocities at particle x-y location and find value at particle

    P_vb=0.0
    SigErr=0
    CALL TSPSI (nN,abb_zb,abb_vb,YP,SIGM,IER,SigErr)
    IF (SigErr.EQ.0) THEN
      P_vb = HVAL (P_zb,nN,abb_zb,abb_vb,YP,SIGM,IER)
    ELSE
      CALL linint(abb_zb,abb_vb,nN,P_zb,P_vb,slope)
    ENDIF

    P_vc=0.0
    SigErr=0
    CALL TSPSI (nN,abb_zc,abb_vc,YP,SIGM,IER,SigErr)
    IF (SigErr.EQ.0) THEN
      P_vc = HVAL (P_zc,nN,abb_zc,abb_vc,YP,SIGM,IER)
    ELSE
      CALL linint(abb_zc,abb_vc,nN,P_zc,P_vc,slope)
    ENDIF

    P_vf=0.0
    SigErr=0
    CALL TSPSI (nN,abb_zf,abb_vf,YP,SIGM,IER,SigErr)
    IF (SigErr.EQ.0) THEN
      P_vf = HVAL (P_zf,nN,abb_zf,abb_vf,YP,SIGM,IER)
    ELSE
      CALL linint(abb_zf,abb_vf,nN,P_zf,P_vf,slope)
    ENDIF


!		*********************************************************
!		*				Find Internal b,c,f Values				*
!		*********************************************************

    !iii. fit polynomial to hydrodynamic model output and find 
    !  internal b,c,f values

    !	 1. Prepare external time step values
    if (p .EQ. 1) then
      ey=0.0
      ey(1) = P_vb
      ey(2) = P_vb
      ey(3) = P_vc
    else
      ey=0.0
      ey(1) = P_vb
      ey(2) = P_vc
      ey(3) = P_vf
    endif

    !	 2. Get value
    vb = polintd(ex,ey,3,ix(1))
    vc = polintd(ex,ey,3,ix(2))
    vf = polintd(ex,ey,3,ix(3))
    P_V = (vb + vc*4 + vf) / 6.0

    SELECT CASE (v)
      CASE (1)
        WCTS_ITPI = vb
      CASE (2)
        WCTS_ITPI = vc
      CASE (3)
        WCTS_ITPI = vf
      CASE (4)
        WCTS_ITPI = P_V
      CASE DEFAULT
        write(*,*) 'ERROR: Illegal WCTS_ITPI version number'
        write(*,*) ' '
        write(*,*) 'The Program Cannot Continue and Will Terminate'
        write(*,*) 'Press Any Key and Enter'
        read(*,*) anykey
        stop
    END SELECT

  END FUNCTION WCTS_ITPI

  DOUBLE PRECISION FUNCTION getSlevel(zeta,depth,i)
    !This function returns the depth of the current s-level
    USE PARAM_MOD, ONLY: hc
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i
    DOUBLE PRECISION, INTENT(IN) :: zeta,depth
    getSlevel = zeta*(1+SC(i))+hc*SC(i)+((-1.*depth)-hc)*CS(i)
  END FUNCTION getSlevel

  DOUBLE PRECISION FUNCTION getWlevel(zeta,depth,i)
    !This function returns the depth of the current w s-level
    USE PARAM_MOD, ONLY: hc
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i
    DOUBLE PRECISION, INTENT(IN) :: zeta,depth
    getWlevel = zeta*(1+SCW(i))+hc*SCW(i)+((-1.*depth)-hc)*CSW(i)
  END FUNCTION getWlevel

  SUBROUTINE getMask_Rho(mask)
    !This subroutine returns the values in the variable mask_rho
    !This is used by createBounds() in the boundary module to make the 
	!  boundaries based on mask_rho
    IMPLICIT NONE
    REAL, INTENT(OUT) :: mask(vi,uj)

    !For Error State, & 'Press Any Key'
    CHARACTER :: anykey

    if(GRD_SET)then
      mask = mask_rho
    else
      write(*,*) 'ERROR: Cannot create boundaries, rho_mask not yet read in'
      write(*,*) ' '
      write(*,*) 'The Program Cannot Continue and Will Terminate'
      write(*,*) 'Press Any Key and Enter'
      read(*,*) anykey
      stop
    endif

  END SUBROUTINE getMask_Rho

  SUBROUTINE getUVxy(ux,uy,vx,vy)
    !This subroutine returns the values in the variables x_u,y_u,x_v,y_v
    !This is used by createBounds() in the boundary module to make the 
	!  boundaries on the U & V node locations
    IMPLICIT NONE
    REAL, INTENT(OUT) :: ux(ui,uj),uy(ui,uj),vx(vi,vj),vy(vi,vj)

    !For Error State, & 'Press Any Key'
    CHARACTER :: anykey

    if(GRD_SET)then
      ux = x_u
      uy = y_u
      vx = x_v
      vy = y_v
    else
      write(*,*) 'ERROR: Cannot create boundaries, x_u, y_u, x_v, or y_v is not yet read in'
      write(*,*) ' '
      write(*,*) 'The Program Cannot Continue and Will Terminate'
      write(*,*) 'Press Any Key and Enter'
      read(*,*) anykey
      stop
    endif

  END SUBROUTINE getUVxy

  SUBROUTINE getR_ele(ele_x,ele_y)
    !This subroutine returns the values in the variables r_ele_x, and r_ele_y
    !This is used by createPolySpecs() in the settlement module to determine 
	!  which habitat polygons are in each element
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(OUT) :: ele_x(4,rho_elements),ele_y(4,rho_elements)

    !For Error State, & 'Press Any Key'
    CHARACTER :: anykey

    if(GRD_SET)then
      ele_x = r_ele_x
      ele_y = r_ele_y
    else
      write(*,*) 'ERROR: Cannot create Poly Specs, r_ele_x or r_ele_y is not yet created'
      write(*,*) ' '
      write(*,*) 'The Program Cannot Continue and Will Terminate'
      write(*,*) 'Press Any Key and Enter'
      read(*,*) anykey
      stop
    endif

  END SUBROUTINE getR_ele

  INTEGER FUNCTION getP_r_element(n)
    !This subroutine returns the id of the rho element the particle is currently in
    !This is used by settlement() to determine which habitat polygons to check for
	!  settlement
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    getP_r_element = P_r_element(n)
  END FUNCTION getP_r_element


END MODULE HYDRO_MOD