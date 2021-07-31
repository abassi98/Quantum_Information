 
      module debugging
         
         contains
         
        
         !initialize a debug variable. If it is true, then the debu gqill be performed
         subroutine init_debug(debug)
            implicit none
            character :: yes, no, var 
            logical :: debug !our debug variable
            yes="y"
            no="n"
            print *, "Do you want to start the debugging mode?[y]/[n]"
            read *, var
            if (var.eq.yes) then
               debug=.true.
               print *, "=================="
               print *, "Debugging mode: ON"
               print *, "=================="
            else if(var.eq.no) then
               debug=.false.
               print *, "==================="
               print *, "Debugging mode: OFF"
               print *, "==================="
            else 
               print *, "ERROR: insert a valid character"
               print *, "Program aborted"
               stop
            endif
         end 

      !check if the dimensiions are positive
         subroutine check_dimensions(nn,mm,myname,debug)
            implicit none
            integer :: nn,mm !dimensions of the matrix
            logical :: debug
            character*10:: myname !name of the matrix
            myname=trim(myname)
            if(debug.eqv..true.) then
               if (nn.ge.1 .and. mm.ge.1) then
                  print *,  "=================================="
                  print *,  "Check dimensions of ",myname," :OK"
                  print *,  "=================================="
               else
                  print *, "====================================="
                  print *, "Check dimensions of ",myname," :ERROR"
                  print *, "Proggram aborted"
                  print *, "====================================="
                  stop
               endif
            endif
         end

      !Check if the two matrices can be multyplied
         subroutine match_dimensions(mm_1,nn_2,debug)
            implicit none
            integer :: mm_1,nn_2
            logical :: debug
            if (debug.eqv..true.) then
               if(nn_2.ne.mm_1) then
                  print *, "Dimensions matching: ERROR"
                  print *, "Program aborted"
                  stop
               else
                  print *, "Dimensions matching: OK"
               endif
            endif
         end

      !Check the allocation status
         subroutine check_allocation(myname,stat,debug)
            implicit none
            integer :: stat !it takes and integr value. 0 is the allocation succeeded
            character*10 :: myname
            logical :: debug
            if(debug.eqv..true.) then
               if(stat.ne.0) then
                  print *, "Allocation of ",myname,": FAILED"
                  print *, "Program aborted"
                  stop
               else
                  print *, "Allocation of ",myname,": SUCCEEDED"
               endif
            endif
         end 
         
         !This subroutine checks if the component by component difference between two matrices of the same diemensions is smaller than a fixed threshold
         subroutine check_diff(A,B,nn,mm,nameA,nameB,threshold,debug)
            implicit none
            integer :: nn,mm
            real, dimension(nn,mm) :: A,B
            logical :: debug
            character*10 :: nameA,nameB
            real :: threshold,diff
            integer :: ii,jj
            integer :: counter  !Count how  many differnces are greater than the threshold
            threshold=0.00001
            counter=0
            if(debug.eqv..true.) then
               do jj=1,mm
                  do ii=1,nn
                     diff=A(ii,jj)-B(ii,jj)
                     if(diff.gt.threshold) then
                        counter=counter+1
                     endif    
                  enddo
               enddo
               if(counter.eq.0) then
                  print *, "Check equality: OK"
                  print *, nameA, " and ",nameB,"are equal"
               else
                  print *, "Check equality: ERROR"
                  print *, nameA, " and ",nameB,"are different"
                  print *, "Program aborted"
                  stop
               endif 
            endif 
         end
         
         !This subroutine checks if the computation time exceeds a certain threshold.
         subroutine check_time(start,threshold,debug)
            implicit none
            logical debug
            real :: start, finish,threshold
            if(debug.eqv..true.) then
               call cpu_time(finish)
               if(finish-start .gt. threshold) then    
                  print *, "Computation time exceeded"
                  print *,"Ran in: ",finish-start
                  print *, "Program aborted"
                  stop
               endif
            endif
         end
         
      end module 
      
      
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
            
            call cpu_time (start) !CHECK HERE
            do jj=1,mmB
               do ii=1,nnA
                  C(ii,jj)=0.0
                  do ss=1,mmA
                     C(ii,jj)=C(ii,jj)+A(ii,ss)*B(ss,jj)
                  enddo
                  call cpu_time(finish)
                  if(finish -start .gt. 300.0) then    !If the computation time exceeds a certain threshold, the program is automatically aborted. Here 60 seconds
                     print *,"Computation time exceeded in mat-mat mult"
                     print *, "Program aborted"
                     stop
                  endif
               enddo
            enddo
      
            call cpu_time(finish)
            time=finish - start
            print *,"Standard mult. runs in:", time
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
                  if(finish -start .gt. 300.0) then    !If the computation time exceeds a certain threshold, the program is automatically aborted. Here 60 seconds
                     print *,"Computation time exceeded in mat-mat mult"
                     print *, "Program abortedd"
                     stop
                  endif
               enddo
            enddo
            call cpu_time(finish)
            time=finish-start
            print *, "Column par column mult. runs in:", time
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
            print *, "The built_in mult. runs in:", time
         end
      
      end module
      
      program main
      
         use debugging
         use multiplication
         
         implicit none
   
         logical :: debug      !Debug logical variable. If it true, the debug is performed
         real:: rr             !Random number to initialize the matrices
         integer ii,jj,ss	!Iteration variables
         real,  dimension(:,:), allocatable :: A,B,C1,C2,C3  !The three matrices
         real, dimension(:,:), allocatable :: m1,m2,out1
         integer :: nnA,mmA,nnB,mmB    !The dimensions of the matrices
         character*10::matrix1,matrix2,matrix3,matrix4,matrix5  !The names of the matrices
         integer :: mystat          !Integer variable to check the allocation
         
         character*10 :: filename
         integer ::  dmn
         real :: start, finish,time1,time2,time3,max_time
         
         character :: yes, no,var,growth_kind
         real :: threshold    !Threshold for checking the differences
         
         time1=0.0
         time2=0.0
         time3=0.0
         threshold=0.00001
         
         !Initialize the names of the matrices
         matrix1="A"
         matrix2="B"
         matrix3="C1"
         matrix4="C2"
         matrix5="C3"
         
         
         call init_debug(debug)       !Ask for debugging mode
         
         !Ask for the dimensions
         print *, "Insert the dimensions (n,m) of ",matrix1
         read *, nnA,mmA
         call check_dimensions(nnA,mmA,matrix1,debug)
         
         print *, "Insert the dimensions (n,m) of ",matrix2
         read *, nnB,mmB
         call check_dimensions(nnB,mmB,matrix2,debug)
         
         call match_dimensions(mmA,nnB,debug)
         
         ! Check if the allocation succedeed
         allocate(A(nnA,mmA),stat=mystat)
         call check_allocation(matrix1,mystat,debug)
         allocate(B(nnB,mmB),stat=mystat)
         call check_allocation(matrix2,mystat,debug)
         allocate(C1(nnA,mmB),stat=mystat)
         call check_allocation(matrix3,mystat,debug)
         allocate(C2(nnA,mmB),stat=mystat)
         call check_allocation(matrix4,mystat,debug)
         allocate(C3(nnA,mmB),stat=mystat)
         call check_allocation(matrix5,mystat,debug)
         
         ! Initialize the two matrices with random numbers taken from  a uniform distribution in [0,1]. Every time we run along columns
         
         do jj=1,mmA
            do ii=1,nnA
               call random_number(rr)
               A(ii,jj)=rr
            enddo
         enddo
         !print *, A  for debugging
         
         do jj=1,mmB
            do ii=1,nnB
               call random_number(rr)
               B(ii,jj)=rr
            enddo
         enddo
         !print *, B for debugging
         
         !Perform the matrix-matrix multiplication in three different ways
        
         call std_mult(A,B,C1,nnA,mmA,nnB,mmB,time1)
         call clmn_mult(A,B,C2,nnA,mmA,nnB,mmB,time2)
         call built_in(A,B,C3,nnA,mmA,nnB,mmB,time3)
         !Check if all the products are equal
         
         call check_diff(C1,C2,nnA,mmB,matrix3,matrix4,threshold,debug)
         call check_diff(C2,C3,nnA,mmB,matrix4,matrix5,threshold,debug)
         !We create a graph to compare the different computation times
         
         yes="y"
         no="n"
         print *,"Do you want to monitor the performance?[y]/[n]"
         read *, var
            if (var.eq.yes) then
               print *, "Insert the computation time"
               read *, max_time
               
               
               filename="output.txt"
               ss=0
               dmn=1
               debug=.true. !Otherwise the program may never stop
         
               call cpu_time(start)
               !finish=start  !Set finish equal to start
        
               open(unit=20,file=filename)
               do while(ss.le.21) 
                  call check_time(start,max_time,debug)   !Check if the computation time is less than the treshsold
                  allocate(m1(dmn,dmn))
                  allocate(m2(dmn,dmn))
                  allocate(out1(dmn,dmn))
                  !Fill the matrices every time with 1! Maybe too slow, find a way to append vectors
                  m1=1
                  m2=1
                  call std_mult(m1,m2,out1,dmn,dmn,dmn,dmn,time1)
                  call clmn_mult(m1,m2,out1,dmn,dmn,dmn,dmn,time2)
                  call built_in(m1,m2,out1,dmn,dmn,dmn,dmn,time3)
                  print *,"========================================"
                  write(20,*) dmn," ",time1," ",time2," ",time3
            
                  deallocate(m1)
                  deallocate(m2)
                  deallocate(out1)
                  ss=ss+1
                  print *, dmn
                  dmn=1.5**ss  !choose the incremental rate here modified
               enddo
               close(20)
               call cpu_time(finish)   !here modified
               print *,"Total computation ran in: ", finish-start
            else if(var.eq.no) then
               print *, "Program aborted"
               stop
            else 
               print *, "ERROR: insert a valid character"
               print *, "Program aborted"
               stop
            endif
 
         deallocate(A)
         deallocate(B)
         deallocate(C1)
         deallocate(C2)
         deallocate(C3)
    
      end
      
      
