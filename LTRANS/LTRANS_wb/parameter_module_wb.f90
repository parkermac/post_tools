MODULE PARAM_MOD

!  The Parameter Module reads in the two include files, LTRANS.inc and GRID.inc, 
!  making the parameters declared in the include files available to all the other 
!  modules.  
!
!  Created by:            Zachary Schlag
!  Created on:			  28 Jul 2008
!  Last Modified on:	  06 Aug 2008

IMPLICIT NONE
PUBLIC
SAVE

include 'GRID_wb.inc'
include 'LTRANS_wb.inc'

END MODULE PARAM_MOD