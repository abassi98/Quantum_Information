   
         module multiplication
          
         contains
         subroutine std_mult(A,B,C,nnA,mmA,nnB,mmB,time)
            implicit none
            integer :: nnA,mmA,nnB,mmB
            integer :: ii,jj,ss
            real :: start,finish,time
            real, dimension(nnA,mmA) :: A
            real, dimension(nnB,mmB) :: B
            real, dimension(nnA,mmB) :: C
            call cpu_time (start) 
            do jj=1,mmB
               do ii=1,nnA
                  C(ii,jj)=0.0
                  do ss=1,mmA
                     C(ii,jj)=C(ii,jj)+A(ii,ss)*B(ss,jj)
                  enddo
                  call cpu_time(finish)
               enddo
            enddo
            call cpu_time(finish)
            time=finish - start
            !print *,"Standard mult. runs in:", time
         end
         
         subroutine clmn_mult(A,B,C,nnA,mmA,nnB,mmB,time)
            implicit none
            integer :: nnA,mmA,nnB,mmB
            integer :: ii,jj,ss
            real :: start,finish,time
            real, dimension(nnA,mmA) :: A
            real, dimension(nnB,mmB) :: B
            real, dimension(nnA,mmB) :: C
            A=transpose(A)
            call cpu_time(start)
            do jj=1,mmB
               do ii=1,nnA
                  C(ii,jj)=0.0
                  do ss=1,mmA
                     C(ii,jj)=C(ii,jj)+A(ss,ii)*B(ss,jj)
                  enddo
                  call cpu_time(finish)
               enddo
            enddo
            call cpu_time(finish)
            time=finish-start
            !print *, "Column par column mult. runs in:", time
            A=transpose(A)
         end
        
         subroutine built_in(A,B,C,nnA,mmA,nnB,mmB,time)
            implicit none
            integer :: nnA,mmA,nnB,mmB
            real :: start,finish,time
            real, dimension(nnA,mmA) :: A
            real, dimension(nnB,mmB) :: B
            real, dimension(nnA,mmB) :: C
            call cpu_time(start)
            C=matmul(A,B)
            call cpu_time(finish)
            time=finish - start
            !print *, "The built_in mult. runs in:", time
         end
      
         end module
      
         program main
      
         use multiplication
         
         implicit none

         real, dimension(:,:),allocatable  :: A,B,C !C=AB
         integer :: dmn !dimensions of the square matrices read from the file 
         real :: time !Computation time
         
!Read the dimensions from a file

        
         open(20,file='dimensions4.txt', status='unknown')
         read(20,*) dmn
         close(20)
         
         
!     allocate the memory
         

         allocate(A(dmn,dmn))
         allocate(B(dmn,dmn))
         allocate(C(dmn,dmn))

         !Initialize with random numbers
         call random_number(A)
         call random_number(B)

         call std_mult(A,B,C,dmn,dmn,dmn,dmn,time)
         open(21,file="std_mult.txt",status="unknown",access="append")
         write(21,*) dmn, time
         close(21)

         call clmn_mult(A,B,C,dmn,dmn,dmn,dmn,time)
         open(21,file="clmn_mult.txt",status="unknown",access="append")
         write(21,*) dmn, time
         close(21)

         call  built_in(A,B,C,dmn,dmn,dmn,dmn,time)
         open(21,file="built_in.txt",status="unknown",access="append")
         write(21,*) dmn, time
         close(21)

         deallocate(A)
         deallocate(B)
         deallocate(C)

         end 
