       module modd_binsolu

       use modd_glo, only: NBSP
       IMPLICIT NONE
       PUBLIC

!**************************************************************************
!Include file: binsolu.h  
!
!Purpose: concentration of solute (micromol solute / microgram water)
!	 in a binary solution at RH = Aw = 0, 0.1, 0.2, 0.3 etc
!
!Include dependencies: used in zsr.c
!
!Revision history: Developed by Betty Pun, December 1999, under CARB funding
!                  
!***********************************************************************/
  
!/* graduation of RH scale in molalbin definition, 
!starting with RH = 0 (first row), ending with 1 (last row) */
         REAL, PARAMETER :: RHgrad = 0.1
  
!/* binary solution molality (umol/ug water) corresponding to RH 
!   at 10% graduations --> 11 rows */   
         REAL, DIMENSION(11,NBSP)   ::  molalbin 
!Values set in ZSR_INI routine

       END MODULE modd_binsolu



