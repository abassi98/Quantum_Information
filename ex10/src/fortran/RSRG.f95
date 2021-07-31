   program RSRG
      use debugging
      use quantum
      use manybodyQS
      implicit none
      integer*16 :: DD,NN,count,status,ii,jj,lwork,lda,kk,hh,zz,N_lam,ll,N_iter
      integer::info
      logical :: debug,printer
      double precision :: lambda,gs,gs_old,threshold
      double complex, dimension(:,:),allocatable ::ham_left,ham_right,int_left,int_right, hamiltonian,temp,projector,sigmax,identity
      double precision, dimension(:),allocatable :: w

      !Set error counter to zero
      count=0

      open(10,file="temp/debugging.txt",status='unknown')
      read(10,*) debug, printer
      close(10)
      open(10,file="temp/parameters.txt",status='unknown')
      read(10,*)  NN,N_lam,kk,N_iter
      close(10)

      DD=2

      !Allocate eigenvalues vector
      allocate(w(DD**(2*NN)),stat=status)
      call check_allocation(status,debug,printer)
      
      
      !Initialize Pauli matrix sigmax
     
      call init_sigmax(sigmax)

      !Temporary matrices for real space RG
      allocate(hamiltonian(DD**(2*NN),DD**(2*NN)),stat=status) 
      call check_allocation(status,debug,printer)
      allocate(temp(DD**(2*NN),DD**(2*NN)),stat=status) 
      call check_allocation(status,debug,printer)
      allocate(int_left(DD**NN,DD**NN),stat=status)
      call check_allocation(status,debug,printer)
      allocate(int_right(DD**NN,DD**NN),stat=status)
      call check_allocation(status,debug,printer)
      allocate(projector(DD**(2*NN),DD**NN),stat=status) !projector
      call check_allocation(status,debug,printer)
      

      print *,"Starting real space RG algorithm"
      open(12,file="data/enRSRG.txt",status='unknown')
      do ll=1,N_lam
         !Set lambda
         lambda=3/dble(N_lam)*(ll-1)
         !Initialize left and right hamiltonians and check equality between them
         call init_ising_hamiltonian(ham_left,DD,NN,lambda,status,count)
         call check_allocation(status,debug,printer)
         call init_ising_hamiltonian(ham_right,DD,NN,lambda,status,count)
         call check_allocation(status,debug,printer)
         threshold=0.0001
         call check_diff(ham_left,ham_right,DD**NN,DD**NN,threshold,debug,printer,count)

         call init_identity(identity,DD,NN-1)
         !Initialize left interaction term
         call matrix_tens_product(identity,DD**(NN-1),sigmax,DD,int_left,DD**NN,count)
         !Initialize right interaction term
         call matrix_tens_product(sigmax,DD,identity,DD**(NN-1),int_right,DD**NN,count)
         deallocate(identity)
        
         !Real space RG algorithm
         call init_identity(identity,DD,NN) 
         open(15,file="data/RSRG_conv.txt",status='unknown')
         do ii=1,N_iter      
            !Double the system
            call matrix_tens_product(ham_left,DD**NN,identity,DD**NN,temp,DD**(2*NN),count)
            hamiltonian=temp
            call matrix_tens_product(identity,DD**NN,ham_right,DD**NN,temp,DD**(2*NN),count)
            hamiltonian=hamiltonian+temp
            call matrix_tens_product(int_left,DD**NN,int_right,DD**NN,temp,DD**(2*NN),count)
            hamiltonian=hamiltonian+temp
            !Diagonalize the Hamiltonian
            temp=hamiltonian
            call diagonalize(temp,DD**(2*NN),w,info)
            call check_eigen(info,debug,printer,count)
            gs=w(1)/dble(NN*(2**ii))
            !Check convergence at lambda=0
            if(ll==1) then
               write(15,*) ii,gs
            endif
            !Save first  D**N eigenvectors in projector
            do hh=1,DD**NN
               projector(:,hh)=temp(:,hh)
               call check_norm(temp(:,hh),DD**(2*NN),0.0001d0,debug,printer,count)
            enddo           
            !Project the hamiltonian in the new basis
            ham_left=matmul(transpose(conjg(projector)),matmul(hamiltonian,projector))      
            ham_right=ham_left
            !Compute the new left interaction term
            call matrix_tens_product(int_left,DD**NN,identity,DD**NN,temp,DD**(2*NN),count)
            int_left=matmul(transpose(conjg(projector)),matmul(temp,projector))
            !Compute the new right interaction term
            call matrix_tens_product(identity,DD**NN,int_right,DD**NN,temp,DD**(2*NN),count)
            int_right=matmul(transpose(conjg(projector)),matmul(temp,projector))
         enddo
         deallocate(ham_left,ham_right,identity)
         write(12,*) lambda,(w(hh)/dble(NN*2**(ii-1)),hh=1,kk)
      enddo

      close(12)
      close(15)
      
      !Free memory
      deallocate(w,sigmax,temp,int_left,int_right,projector,hamiltonian)
            
      print *, "Total number of errors:", count

    end program RSRG
    
    
    

      
