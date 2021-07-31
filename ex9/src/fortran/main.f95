      program ising
      use debugging
      use quantum
      use manybodyQS
      use, intrinsic :: iso_fortran_env
      implicit none
      integer :: DD,NN,count,status,ii,jj,info,lwork,lda,kk,hh,zz,N_lam
      logical :: debug,printer
      double precision :: lambda
      double complex, dimension(:), allocatable ::work
      double complex, dimension(:,:),allocatable :: hamiltonian
      double precision, dimension(:,:), allocatable :: levels,der1,der2
      double precision, dimension(:),allocatable :: w,rwork
      character*20:: file_name
      real :: start,finish1,finish2
      count=0

      open(10,file="temp/debugging.txt",status='unknown')
      read(10,*) debug, printer
      close(10)
      open(10,file="temp/parameters.txt",status='unknown')
      read(10,*)  NN,DD,N_lam,kk
      close(10)


      !print *,"Hamiltonian"
      !do ii=1,DD**NN
         !print *, (hamiltonian(ii,jj),jj=1,DD**NN)
      !enddo

      !Allocate memeory for levels (eigenvalues at varying lambda
      allocate(levels(N_lam,kk),stat=status)
      call check_allocation(status,debug,printer)
      
      do zz=1,N_lam
         call cpu_time(start)
         lambda=3/dble(N_lam)*(zz-1)
         !Initialize hamiltonian
         call init_ising_hamiltonian(hamiltonian,DD,NN,lambda,status,count)
         call check_allocation(status,debug,printer)
         call cpu_time(finish1)
         !Diagonalize hamiltonian
         lda=max(1,DD**NN)
         lwork=max(1,2*DD**NN-1)
         allocate(w(DD**NN))
         allocate(work(lwork))
         allocate(rwork(max(1,3*DD**NN-2)))
         call zheev('N','U',DD**NN,hamiltonian,lda,w,work,lwork,rwork,info)     
         do hh=1,kk
            levels(zz,hh)=w(hh)
         enddo
         if (zz==N_lam) then
            print *, hamiltonian
         endif
         deallocate(hamiltonian,w,work,rwork)
         call cpu_time(finish2)
         if (zz==1) then
            open(13,file="data/time.txt",status='unknown',access='append')
            write(13,*) NN, finish2-start,finish2-finish1
            close(13)
         endif
        
         
      enddo

      allocate(der1(N_lam-2,kk),stat=status)
      call check_allocation(status,debug,printer)
      allocate(der2(N_lam-2,kk),stat=status)
      call check_allocation(status,debug,printer)
      
      !Compute finite difference central first derivative and second derivative
      
      do hh=1,kk
         call central_der(levels(:,hh),der1(:,hh),N_lam,0d0,3d0)
         call second_der(levels(:,hh),der2(:,hh),N_lam,0d0,3d0)
      enddo
      
         
      !Print results in files
      file_name='data/levels_'//char(NN)//'.txt'
      open(12,file=file_name,status='unknown')
      do zz=1,N_lam
         lambda=3/dble(N_lam)*(zz-1)
         write(12,*) lambda,(levels(zz,hh),hh=1,kk)
      enddo
      close(12)

      file_name='data/firstder_'//char(NN)//'.txt'
      open(12,file=file_name,status='unknown')
      do zz=2,N_lam-1
         lambda=3/dble(N_lam)*(zz-1)
         write(12,*) lambda,(der1(zz-1,hh),hh=1,kk)
      enddo
      close(12)

      file_name='data/secondder_'//char(NN)//'.txt'
      open(12,file=file_name,status='unknown')
      do zz=2,N_lam-1
         lambda=3/dble(N_lam)*(zz-1)
         write(12,*) lambda,(der2(zz-1,hh),hh=1,kk)
      enddo
      close(12)

      
      print *, "Total number of errors:", count

    end program ising
    

      
