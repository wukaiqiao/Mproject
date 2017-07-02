!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_USR1 (L)                                         C
!  Purpose: Write user-defined output                                  C
!                                                                      C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE WRITE_USR1(L)
!... Maximum L set to 5 as default  
!... Need to define USR_DT in usr initialize conditions    
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE cdist
      USE compar
      USE cutcell
      USE exit, only: mfix_exit
      USE fldvar
      USE funits
      USE geometry
      USE machine
      USE mpi_utility
      USE output
      USE param
      USE param1
      USE physprop
      USE run
      USE rxns
      USE scalars
      use discretelement, only: DES_CONTINUUM_COUPLED
      use discretelement, only: DISCRETE_ELEMENT
      use discretelement, only: PRINT_DES_DATA
      USE usr

      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! flag whether to write a particular SPx file
      INTEGER L
!-----------------------------------------------
!
!  Local variables
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!             flag whether to write a particular SPx file

!              offset for use in post_mfix
      INTEGER  unit_add
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
! local variables
!
!//
      double precision, allocatable :: array1(:)     !//
      double precision, allocatable :: array2(:)     !//

!             loop counters
      INTEGER LC, NN
!
!             Pointer to the next record
      INTEGER NEXT_REC
!
!              Number of records written each time step
      INTEGER  NUM_REC

      INTEGER  uspx   ! UNIT_SPX + offset from post_mfix
      CHARACTER(LEN=50), DIMENSION(1) :: LINE   !error message
      double precision, dimension(:), allocatable :: TMP_VAR

      allocate(TMP_VAR(DIMENSION_3))

!-----------------------------------------------
      uspx = UNIT_SPX + unit_add

!
      if (myPE .eq.PE_IO) then
         allocate (array1(ijkmax2))   !//
         allocate (array2(ijkmax3))   !//
      else
         allocate (array1(1))   !//
         allocate (array2(1))   !//
      end if

      SELECT CASE (L)
      CASE (1)
         if (myPE.eq.PE_IO.or.bDist_IO) then
            READ (uspx + L, REC=3) NEXT_REC, NUM_REC
            NUM_REC = NEXT_REC
            WRITE (uspx + L, REC=NEXT_REC) REAL(TIME), NSTEP
            NEXT_REC = NEXT_REC + 1
         end if

         if (bDist_IO) then
            DO LC = 1, MMAX
               IF(RE_INDEXING) THEN
                  CALL UNSHIFT_DP_ARRAY(USRrdfg,TMP_VAR)
                  call OUT_BIN_R(uspx+L,TMP_VAR,size(TMP_VAR),NEXT_REC)
               ELSE
                  call OUT_BIN_R(uspx+L,USRrdfg,size(USRrdfg),NEXT_REC)
               ENDIF
            ENDDO  
!           call OUT_BIN_R(uspx+L,P_g,size(P_g),NEXT_REC)
!           call OUT_BIN_R(uspx+L,P_star,size(P_star),NEXT_REC)
         else
         DO LC = 1, MMAX
           call gatherWriteSpx (USRrdfg,array2, array1, uspx+L, NEXT_REC)   !//
         ENDDO
         
         end if
         if (myPE.eq.PE_IO.or.bDist_IO) then
            NUM_REC = NEXT_REC - NUM_REC
            WRITE (uspx + L, REC=3) NEXT_REC, NUM_REC
            if(unit_add == 0) CALL FLUSH_bin (uspx + L)
         end if

      CASE DEFAULT
            LINE(1) = 'Unknown USR1 file index'
            CALL WRITE_ERROR ('WRITE_usr1', LINE, 1)
            CALL MFIX_EXIT(myPE)
      END SELECT

      !//      call unlock_tmp_array
!
      deallocate (array1)
      deallocate (array2)
      deallocate (TMP_VAR)

      RETURN
      END SUBROUTINE WRITE_USR1
