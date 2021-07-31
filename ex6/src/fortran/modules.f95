      module debugging
         
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
            
         end 

      !check if the dimensiions are positive
         subroutine check_dimensions(nn,mm,debug,printer,count)
            implicit none
            integer :: nn,mm !dimensions of the matrix
            logical :: debug,printer
            integer :: count !error counter
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
         end
      
      !check if the input is integer and convert it to integer
      subroutine check_rtint(nn,mm,debug,printer,count)
      implicit none
      integer :: count,mm
      real :: nn
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
      end
      
      !Check the allocation status
         subroutine check_allocation(stat,debug,printer)
            implicit none
            integer :: stat     !it takes and integr value. 0 is the allocation succeeded
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
         end 

!check eigenvalues  computation status
      subroutine check_eigen(info,debug,printer,count)
      implicit none
      integer :: info
      logical :: debug,printer
      integer :: count !error counter
      if(debug.eqv..true.) then
         if(info.eq.0) then
            if(printer.eqv..true.) then
               print *, "Eigenvalues computation: SUCCEEDED"
            endif
         else
            if(printer.eqv..true.) then
               print *, "Eigenvalues computation: FAILED"
            endif
            count=count+1
         endif
      endif
      end
         !This subroutine checks if the computation time exceeds a certain threshold.
         subroutine check_time(start,threshold,debug,printer)
            implicit none
            logical :: debug,printer
            real :: start, finish,threshold
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
         end
         
      end module 
      

      module quantum
      contains

      subroutine harm_matrix(matrix,NN,xmin,omega,hh,status)
      implicit none
      real, dimension(:,:), allocatable :: matrix
      integer :: NN,status
      real :: xmin,omega,hh,xx
      integer :: ii !iterator
      
      allocate(matrix(NN+1,NN+1),stat=status) !allocate matrix, n+1 points
      matrix=0.0                !off diagonal terms
      do ii=1,NN+1
         xx=(xmin+(ii-1)*hh)*hh*omega
         matrix(ii,ii)=1+0.5*(xx**2)
         if(ii.lt.NN+1) then
            matrix(ii,ii+1)=-0.5
              matrix(ii+1,ii)=-0.5
         endif
      enddo

      end

      real function norm(matrix,NN,num,hh)
      implicit none
      integer :: NN,num
      real, dimension(NN+1,NN+1) :: matrix
      integer :: ii
      real :: hh
      norm=0
      do ii=1,NN+1
         norm=norm+(matrix(ii,num)**2)
      enddo
      norm=norm*hh
      return
      end
      
      end module
