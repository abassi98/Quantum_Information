      MODULE MATRICES
         IMPLICIT NONE
         !Define the type matrix
         TYPE MATRIX
            DOUBLE COMPLEX, DIMENSION (:,:), ALLOCATABLE :: mat
            INTEGER, DIMENSION(2) :: dimensions
            DOUBLE COMPLEX :: trace
            DOUBLE COMPLEX :: determinant
            CHARACTER*20 :: mat_name
         END TYPE MATRIX
         
         !Define the two interface operators for calculating the trace and the  adjoint matrix
         INTERFACE OPERATOR (.adj.)
         MODULE PROCEDURE adjoint
         END INTERFACE
         INTERFACE OPERATOR (.trc.)
         MODULE PROCEDURE trace
         END INTERFACE
       
         CONTAINS
         
         !Define the function to compute the trace
         DOUBLE COMPLEX FUNCTION trace(tmp)
            IMPLICIT NONE
            TYPE(MATRIX), INTENT(IN) :: tmp
            INTEGER :: ss   
            trace=(0d0,0d0)        
            IF (tmp%dimensions(1).eq.tmp%dimensions(2)) THEN !Control if the matrix is square. If the matrix is rectangular this function does nothing, so the trace remains (0,0)
               DO ss=1,tmp%dimensions(1)
                  trace=trace +tmp%mat(ss,ss)
               ENDDO
            END IF
         END
         
         !Initialize the matrix
         SUBROUTINE init_mat (tmp)
            IMPLICIT NONE
            TYPE(MATRIX) :: tmp
            INTEGER :: nn,mm,ii,jj
            DOUBLE PRECISION :: xx,yy
            tmp%mat_name="A_matrix"
            tmp%mat_name=TRIM(tmp%mat_name)
            PRINT *, "Insert the dimension of the matrix nxm"
            READ *, nn,mm
            tmp%dimensions(1)=nn    !Initialize the dimensions
            tmp%dimensions(2)=mm
            ALLOCATE(tmp%mat(nn,mm)) !Allocate the memory
            !We fill the matrix with random numbers
            DO jj=1,mm
               DO ii=1,nn
                  CALL RANDOM_NUMBER(xx)
                  CALL RANDOM_NUMBER(yy)
                  tmp%mat(ii,jj)=COMPLEX(xx,yy)
               ENDDO
            ENDDO 
            !We initialize the trace and the determinant
            tmp%trace=.trc.tmp
    	    tmp%determinant=(1d0,1d0) !We don't compute the det. actually, we set it to (1,1) for all the matrices, even if they are rectangular
         END
         
         !Calculate and initialize the adjoint matrix 
         TYPE(MATRIX) FUNCTION adjoint(tmp)
            IMPLICIT NONE
            TYPE(MATRIX),INTENT(IN) :: tmp
            INTEGER :: ii,jj,aa,bb
            adjoint%mat_name="adjoint_"//tmp%mat_name
            adjoint%mat_name=TRIM(adjoint%mat_name)
            adjoint%dimensions(1)=tmp%dimensions(2)
            adjoint%dimensions(2)=tmp%dimensions(1)
            aa=adjoint%dimensions(1)
            bb=adjoint%dimensions(2)
            ALLOCATE(adjoint%mat(aa,bb)) !Allocate the memory for the adjoint matrix
            DO jj=1,adjoint%dimensions(2)
               DO ii=1,adjoint%dimensions(1)
                  adjoint%mat(ii,jj)=CONJG(tmp%mat(jj,ii))
               ENDDO
            ENDDO
            !We use the properties of the adjoint matrix for calculating the trace end the determinant
            adjoint%trace=CONJG(tmp%trace) 
            adjoint%determinant=CONJG(tmp%determinant)
         END 
         
         !This subroutine allows us to view the matrix and its properties defined in the type
         SUBROUTINE view_mat(tmp)
            IMPLICIT NONE
            TYPE(MATRIX) :: tmp  
            INTEGER :: ii,jj  
            PRINT *, "The dimensions are: ", tmp%dimensions
            PRINT *, "The matrix ",tmp%mat_name," is (Re,Im): "
            DO ii=1,tmp%dimensions(1)
               PRINT *, (tmp%mat(ii,jj), jj=1,tmp%dimensions(2))
            ENDDO
            IF (tmp%dimensions(1).eq.tmp%dimensions(2)) THEN
               PRINT *,"The trace is:", tmp%trace
               PRINT *,"The determinant is:", tmp%determinant
            ELSE 
               PRINT *,"The trace for this matrix is not defined "
               PRINT *,"The determinant for this matrix is not defined "
            END IF 
         END
         
         !Print the matrix in a file whose name is required as an input
         SUBROUTINE MYOPEN (tmp)
            IMPLICIT NONE
            CHARACTER*20 :: output
            TYPE(MATRIX) :: tmp
            INTEGER :: ii,jj
            PRINT *,"Insert the name of the file"
            READ *, output
            output=TRIM(output)
            OPEN(UNIT=20,FILE=output,STATUS='unknown')
            WRITE(20,*) "The dimensions are:", tmp%dimensions
            WRITE(20,*) "The matrix ",tmp%mat_name," is (Re,Im): "
            DO ii=1,tmp%dimensions(1)
               WRITE(20,*) (tmp%mat(ii,jj), jj=1,tmp%dimensions(2))
            ENDDO
            IF (tmp%dimensions(1).eq.tmp%dimensions(2)) THEN
               WRITE(20,*)"The trace is:", tmp%trace
               WRITE(20,*)"The determinant is:", tmp%determinant
            ELSE 
               WRITE(20,*) "The trace is not defined"
               WRITE(20,*) "The determinant is not defined"
            END IF 
            CLOSE(20)
            RETURN
         END
      
      END MODULE 

      PROGRAM mat2
         USE MATRICES
         IMPLICIT NONE
         TYPE(MATRIX) :: A,B
         CALL init_mat(A)
         CALL view_mat(A)
         CAll MYOPEN(A)
         B=.adj.A
         CALL view_mat(B)
         CALL MYOPEN(B)
      END 
       
