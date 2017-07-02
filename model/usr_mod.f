      MODULE usr
!     
!     
!       Declare the user-defined namelist variables (usrnlst.inc) in this module.
!       Also Include user-defined variables in this module.  To access the
!       variables from a subroutine add the statement "Use usr".  If allocatable
!       arrays are defined in this module allocate them in usr0.  To turn on the
!       user defined subroutines (usr0, usr1, and usr2) set call_usr to true in
!       mfix.dat.
!
!                       a dummy variable listed in usrnlst.inc
        DOUBLE PRECISION DUMMY_DP
! Solids phase frictional coefficient
        DOUBLE PRECISION usr_fricoef
        DOUBLE PRECISION ITERMU, u_muprime, usr_ratio
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: USRrdfg

! Remember to initialize the usr-defined variables using usr_initialization
      CONTAINS   

!``````````````````````````````````````````````````````````````````````!
!  Module name: phi_c (par_fric)                                       !
!  Author: KWU                                         Date: 27-jun-17 !
!                                                                      !
!  Purpose: Calculate phi_c from the given mu (par fric coefficient)   !
!......................................................................!
      DOUBLE PRECISION FUNCTION PHI_C (par_fric)

! Dummy Arguments:
!---------------------------------------------------------------------//
            DOUBLE PRECISION, INTENT(IN) :: par_fric

             PHI_C = -0.1891*par_fric**3 + 0.3636*par_fric**2 -  &
             0.2285*par_fric + 0.6348

            RETURN
      END FUNCTION PHI_C


!     Allocate(  USRrdfg (DIMENSION_3, DIMENSION_M) )
END MODULE usr