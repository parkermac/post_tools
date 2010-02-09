MODULE NORM_MOD

!  The Norm Module contains the function Norm, which returns a random number (deviate) 
!   drawn from a normal distribution with zero mean and unit variance (i.e., standard 
!   deviation = 1).
!
!  Created by:            Elizabeth North
!  Modified by:           Zachary Schlag
!  Created on:			  22 Aug 2008
!  Last Modified on:	  29 Aug 2008

IMPLICIT NONE
PRIVATE

  PUBLIC :: NORM

CONTAINS

  REAL FUNCTION norm()
    ! This function returns a normally distributed deviate with zero mean and 
    ! unit variance (i.e., standard deviation = 1). By E. North, 8/22/08
    USE PARAM_MOD,  ONLY: PI
    USE RANDOM_MOD, ONLY: genrand_real1
    IMPLICIT NONE

    REAL :: dev1,dev2

    dev1 = genrand_real1()
    dev2 = genrand_real1()
    norm = sqrt(-2.*log(dev1)) * cos(2.*PI*dev2)

  END FUNCTION norm

END MODULE NORM_MOD