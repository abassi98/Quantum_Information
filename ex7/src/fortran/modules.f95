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
            
          end subroutine init_debug
          

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
          end subroutine check_dimensions
          
      
      !check if the input is integer and convert it to integer
      subroutine check_rtint(nn,mm,debug,printer,count)
      implicit none
      integer :: count,mm
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
          end subroutine check_allocation
          

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
          
         
      end module 
      

      module quantum
      contains

      subroutine harm_matrix(matrix,NN,xmin,omega,hh,status)
      implicit none
      double precision, dimension(:,:), allocatable :: matrix
      integer :: NN,status
      double precision :: xmin,omega,hh,xx
      integer :: ii !iterator
      
      allocate(matrix(NN,NN),stat=status) !allocate matrix, n+1 points
      matrix=0.0                !off diagonal terms
      do ii=1,NN
         xx=(xmin+(ii-1)*hh)*hh*omega
         matrix(ii,ii)=1+0.5*(xx**2)
         if(ii.lt.NN) then
            matrix(ii,ii+1)=-0.5
              matrix(ii+1,ii)=-0.5
         endif
      enddo

      end subroutine harm_matrix
    

      double precision function norm(vec,NN,hh) !hh=1 for norm of a vector
      implicit none
      integer :: NN,ii
      double complex, dimension(NN) :: vec
      double precision :: hh
      double complex :: temp
      temp=cmplx(0d0,0d0)
      do ii=1,NN
         temp=temp+vec(ii)*conjg(vec(ii))*hh
      enddo
      norm=dble(sqrt(temp))
      return
      end function norm

      subroutine peak_search(vec,NN,xmin,hh,position,maxim)
        implicit none
        integer :: NN,ii
        double complex, dimension (NN) :: vec
        real, dimension(NN) :: temp
        double precision :: xmin,hh,position,maxim
        do ii=1,NN
           temp(ii)=sqrt(dble(vec(ii)*conjg(vec(ii))))
        enddo
        maxim=maxval(temp)
        ii=findloc(temp,maxim,dim=(NN))
        position=xmin+(ii-1)*hh+hh/2
      end subroutine peak_search
      
           
      !This subroutine compute the forward or backward FFT according to the flag
      subroutine FFT(in,out,NN,flag,count)
        use,intrinsic :: iso_c_binding     
        implicit none
        include 'fftw3.f03'   
        integer :: NN, ii, flag,count
        double complex, dimension(NN) :: in,out
        integer*8 :: plan
        if(flag==1) then
           call dfftw_plan_dft_1d(plan,NN,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
           call dfftw_execute_dft(plan,in,out)
           call dfftw_destroy_plan(plan)
           out=out/sqrt(dble(NN))
        elseif(flag==-1) then
           call dfftw_plan_dft_1d(plan,NN,in,out,FFTW_BACKWARD,FFTW_ESTIMATE)
      	   call dfftw_execute_dft(plan,in,out)
          call dfftw_destroy_plan(plan)
          out=out/sqrt(dble(NN))
        else
           print *, "Error: Invalid Flag"
           count=count+1
        endif
      end subroutine FFT
      
        
      !This subroutine applies the potential matrix to x-sapce vector
      subroutine potential(wave,NN,omega,xmin,hh,dt,time,speed)
        implicit none
        integer ::NN !dimension of vector
        double complex, dimension(NN) :: wave
        double precision :: omega, xmin, hh,V,dt,time,speed
        integer :: ii
        do ii=1,NN
           V=0.5*((xmin+hh/2+hh*dble(ii-1)-time*speed)*omega)**2
           wave(ii)=exp(-cmplx(0d0,1d0)*0.5*V*dt)*wave(ii)
        enddo

      end subroutine potential

      !This subroutine applies the kinetic matrix to p-space vector
      subroutine kinetic(wave,NN,pmax,dt)
        implicit none
        integer ::NN !dimension of vector-1
        double complex, dimension(NN) :: wave
        double precision :: pmax,K,dt
        integer :: ii
        do ii=1,NN/2
           K=0.5*((pmax*dble(ii-1))**2)
           wave(ii)=exp(-cmplx(0d0,1d0)*K*dt)*wave(ii)
        enddo
        do ii=NN/2+1,NN
           K=0.5*((pmax*dble(ii-NN-1))**2)
           wave(ii)=exp(-cmplx(0d0,1d0)*K*dt)*wave(ii)
        enddo
      end subroutine kinetic

      
      end module
