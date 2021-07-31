      PROGRAM mat_time
      
!     This program compute the product of two matrices
      
      IMPLICIT NONE

!     Read the dimension of the matrices
     
      INTEGER nn
      INTEGER ii,jj,ss
      REAL, DIMENSION(:,:), ALLOCATABLE :: A,B,C
      REAL start, finish !see the cpu time
      

      nn=2000
      
      ALLOCATE(A(nn,nn))
      ALLOCATE(B(nn,nn))
      ALLOCATE(C(nn,nn))
      


!     Write A

      DO jj=1,nn
         DO ii=1,nn
            A(ii,jj)=ii+(jj-1)*nn
         ENDDO
      ENDDO
      
!     Write B

     
      DO jj=1,nn
         DO ii=1,nn
            B(ii,jj)=ii+jj
         ENDDO
      ENDDO

!     First loop, standard multiplication!
      
      CALL CPU_TIME (start)
    
      DO jj=1,nn
         DO ii=1,nn
            C(ii,jj)=0
            DO ss=1,nn
               C(ii,jj)=C(ii,jj)+A(ii,ss)*B(ss,jj)
            ENDDO
         ENDDO
      ENDDO
      
      CALL CPU_TIME(finish)
      

      PRINT *, "Standard loop computation time is:  ", finish - start

!     Second loop, rowxrow multiplication throigh transposition

      CALL CPU_TIME (start)
      A=TRANSPOSE(A)
      CALL CPU_TIME(finish)

      PRINT *, "TRansposing time is:  ", finish - start
      
      CALL CPU_TIME (start)
    
      DO jj=1,nn
         DO ii=1,nn
            C(ii,jj)=0
            DO ss=1,nn
               C(ii,jj)=C(ii,jj)+A(ss,ii)*B(ss,jj)
            ENDDO
         ENDDO
      ENDDO

      
      
      CALL CPU_TIME(finish)

     

      PRINT *, "Second loop computation time is:  ",finish-start


      
!     Built_in function
      A=TRANSPOSE(A)

      CALL CPU_TIME(start)
      C=MATMUL(A,B)
      CALL CPU_TIME(finish)

      PRINT *, "BUilt-in mult. time is: ", finish-start
      
      
      DEALLOCATE(A)
      DEALLOCATE(B)
      DEALLOCATE(C)

     
   
      STOP
      END
