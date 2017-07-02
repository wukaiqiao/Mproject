!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR_INIT_NAMELIST                                      C
!  Purpose: initialize user_defined NAMELIST variables                 C
!                                                                      C
!  Author:                                            Date:            C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE USR_INIT_NAMELIST
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      Use usr
      USE param, only: dimension_3, dimension_m, DIMENSION_USR
      USE param1, only: zero
      USE output, only: USR_DT
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER LC
!-----------------------------------------------
!     
      Allocate(USRrdfg (DIMENSION_3, DIMENSION_M) )
      ITERMU = ZERO
      u_muprime = ZERO 
      usr_ratio = ZERO
      USRrdfg(:,1) = ZERO
      
      !DO LC=1, DIMENSION_USR
      USR_DT(:) = 0.01D0      
      !ENDDO

      RETURN
      END SUBROUTINE USR_INIT_NAMELIST
