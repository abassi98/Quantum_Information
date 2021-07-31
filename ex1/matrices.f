      PROGRAM matrices
      
!     This program compute the product of two matrices
      
      IMPLICIT NONE

!     Read the dimension of the matrices
     
      INTEGER nn,mm,kk,ll
      INTEGER ii,jj,ss
      REAL value
      REAL, DIMENSION(:,:), ALLOCATABLE :: A,B,C
      REAL start, finish !see the cpu time
      

!     Input of the dimensions
      
      PRINT *, "Insert the dimensions of the first matrix nxm"
      READ *, nn,mm
      PRINT *, "Insert the dimension of the second matrix mxk"
      READ *, kk,ll
      
      !     Check if the multiplication can be performed
      
      IF (kk.ne.mm) THEN
         PRINT *, "The multiplication cannot be performed"
         STOP
      END IF

      ALLOCATE(A(nn,mm))
      ALLOCATE(B(kk,ll))
      ALLOCATE(C(nn,ll))
      


!     Read the matrix A from the terminal

      DO jj=1,mm
         DO ii=1,nn
            PRINT *, "A(",ii,",",jj,")"
            READ *, value
            A(ii,jj)=value
         ENDDO
      ENDDO
      
!     Read the matrix B from the terminal

     
      DO jj=1,ll
         DO ii=1,kk
            PRINT *, "B(",ii,",",jj,")"
            READ *, value
            B(ii,jj)=value
         ENDDO
      ENDDO

      CALL CPU_TIME (start)

      A=TRANSPOSE(A)
      DO jj=1,ll
         DO ii=1,nn
            C(ii,jj)=0
            DO ss=1,mm
               C(ii,jj)=C(ii,jj)+A(ss,ii)*B(ss,jj)
            ENDDO
         ENDDO
      ENDDO
      CALL CPU_TIME(finish)
      
      DEALLOCATE(A)
      DEALLOCATE(B)

      PRINT *, "The computation time is ", finish -start
      PRINT *, C
      DEALLOCATE(C)
      STOP
      END
