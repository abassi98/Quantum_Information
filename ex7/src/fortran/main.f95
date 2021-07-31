      program harm_oscillator_timedep
      use debugging
      use quantum
      
      implicit none
           
      double precision :: xmin,xmax              !Interval in x-space
      double precision:: omega !pulse
      integer :: NN,dmn             !number of points of the disretization
      double precision, dimension(:,:), allocatable :: matrix !Our hermitian matrix
      double precision :: hh    !delta x
      logical :: debug, printer !logical variable for debugging
      integer :: status         !allocation status
      integer :: count          !counting errors
      integer :: ii,jj !iterators
      !p-space parameters
      double precision :: pmin,pmax, hhp !interval by pmin an ddiscretization 
      !ssyev variables
      integer :: lda,lwork,info
      double precision, dimension(:),allocatable :: w,work
      double complex, dimension(:),allocatable :: wave,temp
      double precision, dimension(:,:), allocatable :: ev_temp
      double precision :: xx,position,maxim
!     time-dependence parameters
      double precision :: speed, dt,time !max time T and discretization time dt
      integer :: N_time            !number of points of time disretization
      integer :: flag !flag for FFT, forward or backward
      
      
      count=0 
      open(10,file="temp/debugging.txt",status='unknown')
      read(10,*) debug, printer
      close(10)

      open(10,file="temp/parameters.txt",status='unknown')
      read(10,*)  dmn, xmax,speed,N_time
      close(10)
      NN=int(2**dmn)
      call check_dimensions(NN,NN,debug,printer,count) !check positiveness dimenion

      
      omega=1d0 !set omega, fixed this time
      xmin=-xmax
      hh=2*xmax/dble(NN)

      !Initialize harmonic potential matrix and check allocation
      call harm_matrix(matrix,NN,xmin,omega,hh,status)
      call check_allocation(status,debug,printer)

      !allocate and check dsyev variables
      allocate(w(NN),stat=status)
      call check_allocation(status,debug,printer)

      lda=max(1,NN)
      lwork=max(1,3*NN-1)
      allocate(work(lwork),stat=status)
      call check_allocation(status,debug,printer)

      !Compute ground state 
      call dsyev('V','U',NN,matrix,NN,w,work,lwork,info)
      call check_eigen(info,debug,printer,count)

      !Set ground state in a vector and free the rest of the memory
      allocate(wave(NN),stat=status)
      call check_allocation(status,debug,printer)
      allocate(temp(NN),stat=status)
      call check_allocation(status,debug,printer)
      matrix=matrix/sqrt(hh) !normalization as a function
      wave=cmplx(matrix(:,1),0d0)
      
      
      print *, "Ground state computed. Starting time evolution."
      
      !Define time discretization
      dt=1/(speed*dble(N_time)) !speed=1/T
      
      !Define p_space parameters
      pmax=abs(acos(-1d0))/xmax !here p-space interval
      pmin=-pmax
      print *, "Integral norm of ground state: ", norm(wave,NN,hh)
      
      !Allocate a matrix where to store time evolution
      allocate(ev_temp(NN+1,N_time+2),stat=status)
      call check_allocation(status,debug,printer)
      
      do jj=1,NN
         ev_temp(jj,1)=xmin+hh*(jj-1)+hh/2
         ev_temp(jj,2)=sqrt(dble(wave(jj)*conjg(wave(jj))))
      enddo

      open(20,file="data/statistics.txt",status='unknown') !to print statistics
      
      !Time evolution. Split operator  method
      do ii=1,N_time
         time=dt*dble(ii)
         !compute space evolution
         call potential(wave,NN,omega,xmin,hh,dt,time,speed)       
         !Forward FFT to go to p-space
         call FFT(wave,temp,NN,1,count)
         wave=temp
         temp=(0d0,0d0)
         !compute kinetic evolution
         call kinetic(wave,NN,pmax,dt)           
         !Backward FFT to get back to x-space
         call FFT(wave,temp,NN,-1,count)
         wave=temp
         temp=(0d0,0d0)
         !Last sapce evolution
         call potential(wave,NN,omega,xmin,hh,dt,time,speed)   
         do jj=1,NN
            ev_temp(jj,ii+2)=sqrt(dble(wave(jj)*conjg(wave(jj))))
         enddo
         call  peak_search(wave,NN,xmin,hh,position,maxim)
         write(20,*) time, norm(wave,NN,hh),position-speed*time,maxim
      enddo
      
      close(20)
      !write in a file
      open(10,file="data/time_ev.txt",status='unknown')
      do ii=1,NN
         write(10,*) (ev_temp(ii,jj),jj=1,N_time+2)
      enddo
      close(10)
      
      !Free memory and print number of errors encountered
      
      print *, "Error encountered:", count
      !Free memory
      deallocate(wave,temp)
      deallocate(matrix,w,work)
    end program harm_oscillator_timedep
    

      
