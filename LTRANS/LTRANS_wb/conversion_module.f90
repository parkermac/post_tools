MODULE CONVERT_MOD

! The Conversion Module contains all the procedures necessary to convert locations between 
!  latitude and longitude coordinates and metric (x and y) coordinates.  The conversions 
!  are done using equations from the sg_mercator.m and seagrid2roms.m matlab scripts that 
!  are found in Seagrid and are used to generate the ROMS model grid.  The module contains
!  four interface blocks to cover the four necessary conversions:  
!
!  lon2x = longitude to x- coordinate
!  lat2y = latitude  to y- coordinate
!  x2lon = x- coordinate to longitude
!  y2lat = y- coordinate to latitude.
!
! Created by:			Zachary Schlag
! Created on:			23 Jul 2008
! Last modified on:	    11 Aug 2008

USE PARAM_MOD, ONLY: PI,RCF,EARTH_RADIUS
IMPLICIT NONE
PRIVATE
SAVE

  !NOTE: The Interface blocks are so that the functions can handle both
  !  REAL and DOUBLE PRECISION variables as input 
  !NOTE: Output is DOUBLE PRECISION regardless of input

  !The functions used to convert are as follows:
  ! x   = lon / RCF * Earth_Radius
  ! y   = log( tan( pi/4.0 + lat/(RCF*2.0) ))*Earth_Radius
  ! lon = x / Earth_Radius * RCF
  ! lat = 2.0*RCF * ( atan(exp(y/Earth_Radius)) - pi/4.0 )

  !Return x location given longitude
  INTERFACE lon2x
    MODULE PROCEDURE rlon2x  !real input
    MODULE PROCEDURE dlon2x  !double precision input
  END INTERFACE lon2x


  !Return y location given latitude
  INTERFACE lat2y
    MODULE PROCEDURE rlat2y  !real input
    MODULE PROCEDURE dlat2y  !double precision input
  END INTERFACE lat2y


  !Return longitude given x location
  INTERFACE x2lon
    MODULE PROCEDURE rx2lon  !real input
    MODULE PROCEDURE dx2lon  !double precision input
  END INTERFACE x2lon


  !Return latitude given y location
  INTERFACE y2lat
    MODULE PROCEDURE ry2lat  !real input
    MODULE PROCEDURE dy2lat  !double precision input
  END INTERFACE y2lat

  !The following procedures have been made public for the use of other program units:
  PUBLIC :: lon2x,lat2y,x2lon,y2lat

CONTAINS

  DOUBLE PRECISION FUNCTION rlon2x(lon)
    IMPLICIT NONE
    REAL, INTENT(IN) :: lon

    rlon2x = lon / RCF * Earth_Radius

  END FUNCTION rlon2x

  DOUBLE PRECISION FUNCTION dlon2x(lon)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: lon

    dlon2x = lon / RCF * Earth_Radius

  END FUNCTION dlon2x

  DOUBLE PRECISION FUNCTION rlat2y(lat)
    IMPLICIT NONE
    REAL, INTENT(IN) :: lat

    rlat2y = log( tan( pi/4.0 + lat/(RCF*2.0) ))*Earth_Radius

  END FUNCTION rlat2y

  DOUBLE PRECISION FUNCTION dlat2y(lat)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: lat

    dlat2y = log( tan( pi/4.0 + lat/(RCF*2.0) ))*Earth_Radius

  END FUNCTION dlat2y

  DOUBLE PRECISION FUNCTION rx2lon(x)
    IMPLICIT NONE
    REAL, INTENT(IN) :: x

    rx2lon = x / Earth_Radius * RCF

  END FUNCTION rx2lon

  DOUBLE PRECISION FUNCTION dx2lon(x)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: x

    dx2lon = x / Earth_Radius * RCF

  END FUNCTION dx2lon

  DOUBLE PRECISION FUNCTION ry2lat(y)
    IMPLICIT NONE
    REAL, INTENT(IN) :: y

    ry2lat = 2.0*RCF * ( atan(exp(y/Earth_Radius)) - pi/4.0 )

  END FUNCTION ry2lat

  DOUBLE PRECISION FUNCTION dy2lat(y)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: y

    dy2lat = 2.0*RCF * ( atan(exp(y/Earth_Radius)) - pi/4.0 )

  END FUNCTION dy2lat

END MODULE CONVERT_MOD