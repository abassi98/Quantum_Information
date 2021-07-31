
module debugging
        use quantum
        use manybodyQS
         contains
         !Every subroutine has an error counter 
        
         !initialize a debug variable. If it is true, then the debu gqill be performed
         subroutine init_debug(debug,printer)
            implicit none
            character :: yes, no, var 
            logical :: debug,printer !our debugging variable and printer
            printer=.false. !Default no warmings messages
            yes="y"
            print *, "Do you want to start the debugging mode?[y]/[n]"
            read *, var
            if (var.eq.yes) then
               debug=.true.
               print *, "=================="
               print *, "Debugging mode: ON"
               print *, "=================="
               print *, "Do you want  to print checkpoints?[y]/[n]"
               read *, var
               if(var.eq.yes) printer=.true.
            else
               debug=.false.
               print *, "==================="
               print *, "Debugging mode: OFF"
               print *, "==================="
            endif
            
          end subroutine init_debug
          

      !check if the dimensiions are positive
         subroutine check_dimensions(nn,mm,debug,printer,count)
            implicit none
            integer*16 :: nn,mm !dimensions of the matrix
            logical :: debug,printer
            integer*16 :: count !error counter
            if(debug.eqv..true.) then
               if (nn.ge.1 .and. mm.ge.1) then
                  if(printer.eqv..true.) then
                  print *,  "=================================="
                  print *,  "Check dimensions: OK"
                  print *,  "=================================="
                  endif
               else
                  if(printer.eqv..true.) then
                  print *, "====================================="
                  print *, "Check dimensions of: ERROR"
                  print *, "====================================="
                  endif
                  nn=1 !set default  values for dimensions if it encounters some error
                  mm=1 !try to interpret dimensions anyway
                  count=count+1
                  stop
               endif
            endif
          end subroutine check_dimensions
          
      
      !check if the input is integer and convert it to integer
      subroutine check_rtint(nn,mm,debug,printer,count)
      implicit none
      integer*16 :: count,mm
      double precision :: nn
      logical :: debug,printer
      if(debug.eqv..true.) then
         if(floor(nn).eq.nn) then
            if(printer.eqv..true.) then
               print *, "Check integerness: SUCCEEDED"
            endif
         else
            if(printer.eqv..true.) then
               print *, "Check integerness: FAILED"
            endif
            count=count+1
         endif
      endif
      mm=int(floor(nn)) !in any case convert nn in integer
      end subroutine check_rtint
    
      
      !Check the allocation status
         subroutine check_allocation(stat,debug,printer)
            implicit none
            integer*16 :: stat     !it takes and integr value. 0 is the allocation succeeded
            logical :: debug,printer
            if(debug.eqv..true.) then
               if(stat.ne.0) then
                  if(printer.eqv..true.) then
                  print *, "Allocation:  FAILED"
                  print *, "Program aborted"
               endif
               stop !if allocation fails, the program must be stopped
               else
                  if(printer.eqv..true.) then
                     print *, "Allocation: SUCCEEDED"
                  endif
               endif
            endif
          end subroutine check_allocation
          

!check eigenvalues  computation status
      subroutine check_eigen(info,debug,printer,count)
      implicit none
      integer :: info
      logical :: debug,printer
      integer*16 :: count !error counter
      if(debug.eqv..true.) then
         if(info.eq.0) then
            if(printer.eqv..true.) then
               print *, "Eigenvalues computation: SUCCEEDED"
            endif
         else
            if(printer.eqv..true.) then
               print *, "Eigenvalues computation: FAILED"
               print *, "INFO=",info
            endif
            count=count+1
         endif
      endif
      end subroutine check_eigen
    
         !This subroutine checks if the computation time exceeds a certain threshold.
         subroutine check_time(start,threshold,debug,printer)
            implicit none
            logical :: debug,printer
            double precision :: start, finish,threshold
            if(debug.eqv..true.) then
               call cpu_time(finish)
               if(finish-start .gt. threshold) then
                  if(printer.eqv..true.) then
                  print *, "Computation time exceeded"
                  print *,"Ran in: ",finish-start
                  print *, "Program aborted"
                  endif
                  stop
               endif
            endif
          end subroutine check_time

         subroutine check_diff(mat1,mat2,nn,mm,threshold,debug,printer,count)
            implicit none
            integer*16 :: nn,mm
            double complex, dimension(nn,mm) :: mat1,mat2
            logical :: debug,printer
            double precision :: threshold,diff
            integer*16 :: ii,jj,temp,count
            temp=0
            if(debug.eqv..true.) then
               do jj=1,mm
                  do ii=1,nn
                     diff=sqrt(dble((mat1(ii,jj)-mat2(ii,jj))*conjg(mat1(ii,jj)-mat2(ii,jj))))
                     if(diff.gt.threshold) then
                        temp=temp+1
                     endif    
                  enddo
               enddo
               if(temp.eq.0) then
                  if(printer.eqv..true.) then
                     print *, "Check equality: OK"
                  endif         
               else
                  if(printer.eqv..true.) then
                     print *, "Check equality: ERROR"
                  endif
                  count=count+1
               endif 
            endif 
          end subroutine check_diff
          
          subroutine check_trace(matrix,dmn,threshold,debug,printer,count)
            implicit none
            logical ::debug,printer
            integer*16:: dmn,count
            double complex, dimension(dmn,dmn)::matrix
            double precision:: threshold
            if(debug.eqv..true.) then
               if(abs(dble(trace(matrix,dmn))-1)>=threshold.or.abs(dimag(trace(matrix,dmn)))>=threshold) then
                  count=count+1
                  if (printer.eqv..true.) then
                     print *, "Check unitariety of trace: ERROR"
                  endif
               else
                  if(printer.eqv..true.) then
                     print *, "Check unitariety of trace: OK"
                  endif
               endif
            endif                    
          end subroutine check_trace

          subroutine check_norm(vector,dmn,threshold,debug,printer,count)
            implicit none
            logical::debug,printer
            integer*16::dmn,count
            double complex,dimension(dmn)::vector
            double precision::threshold
            if(debug.eqv..true.) then
               if(printer.eqv..true.) then
                  if(abs(norm(vector,dmn,1d0)-1)>=threshold) then
                     count=count+1
                     print *, "Check normalization: ERROR"
                  else
                     print *, "Check normalization: OK"
                  endif
               else
                  if(abs(norm(vector,dmn,1d0)-1)>=threshold) then
                     count=count+1
                  endif
               endif
            endif
          end subroutine check_norm
          
         
      end module 