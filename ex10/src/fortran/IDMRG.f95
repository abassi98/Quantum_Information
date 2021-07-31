   program IDMRG
      use debugging
      use quantum
      use manybodyQS
      implicit none
      integer*16:: DD,NN,count,status,ii,jj,kk,hh,ll,N_lam,zz,N_iter
      integer:: info
      logical :: debug,printer
      double precision :: lambda,gs,gs_old,threshold
      double complex,dimension(:), allocatable::groundstate
      double complex,dimension(:,:),allocatable::hamiltonian,ham_left,ham_right,reduced_ham,gsdensity,sigmax,sigmaz,density_right
      double complex,dimension(:,:),allocatable::prod,identity,small_temp,density_left,least_temp,projector,temp
      double precision, dimension(:),allocatable :: w,w2

      !Counting errors set to zero
      count=0

      open(10,file="temp/debugging.txt",status='unknown')
      read(10,*) debug, printer
      close(10)
      open(10,file="temp/parameters.txt",status='unknown')
      read(10,*)  NN,N_lam,kk,N_iter
      close(10)

      !Set dimension of spins
      DD=2
      
      !Check if N is too large. Showuld be 2*N+1<12
      if (NN>5) then
         NN=5
         print *, "Warning: dimension too large. Set N=5"
      endif
      
     

      !Temporary matrices for infinite density matrix RG: temp,groundstate,gsdensity,density_left,density_right,projector
      allocate(temp(DD**(2*NN+2),DD**(2*NN+2)),stat=status) 
      call check_allocation(status,debug,printer)
      allocate(groundstate(DD**(2*NN+2)),stat=status)
      call check_allocation(status,debug,printer)
      allocate(gsdensity(DD**(2*NN+2),DD**(2*NN+2)),stat=status) 
      call check_allocation(status,debug,printer)
      allocate(density_left(DD**(NN+1),DD**(NN+1)),stat=status)
      call check_allocation(status,debug,printer)
      allocate(density_right(DD**(NN+1),DD**(NN+1)),stat=status)
      call check_allocation(status,debug,printer)
      allocate(projector(DD**(NN+1),DD**NN),stat=status)
      call check_allocation(status,debug,printer)

      !Temporary matrices to construct the hamiltonian: hamiltonian,sigmax,sigmaz,sigmax tp sigmax,small_temp,red_ham,least_temp
      allocate(hamiltonian(DD**(2*NN+2),DD**(2*NN+2)),stat=status)
      call check_allocation(status,debug,printer)
      call init_sigmax(sigmax)
      allocate(prod(DD**2,DD**2),stat=status)
      call check_allocation(status,debug,printer)
      call matrix_tens_product(sigmax,DD,sigmax,DD,prod,DD**2,count)
      allocate(small_temp(DD**(NN+2),DD**(NN+2)),stat=status)
      call check_allocation(status,debug,printer)
      allocate(reduced_ham(DD**NN,DD**NN),stat=status)
      call check_allocation(status,debug,printer)
      call init_sigmaz(sigmaz)
      allocate(least_temp(DD**(NN+1),DD**(NN+1)),stat=status)
      call check_allocation(status,debug,printer)
      
      
      !Eigenvalues vectors
      allocate(w(DD**(2*NN+2)),stat=status)
      call check_allocation(status,debug,printer)
      allocate(w2(DD**(NN+1)),stat=status)
      call check_allocation(status,debug,printer)

      print *, "Starting infinite density matrix RG algorithm"
      open(12,file="data/enIDMRG.txt",status='unknown')
      do ll=1,N_lam        
         lambda=3/dble(N_lam)*(ll-1)
         !Initialize hamiltonian of left and right part
         call init_ising_hamiltonian(ham_left,DD,NN+1,lambda,status,count)
         call check_allocation(status,debug,printer)
         call init_ising_hamiltonian(ham_right,DD,NN+1,lambda,status,count)
         call check_allocation(status,debug,printer)

         !check difference
         call check_diff(ham_left,ham_right,DD**(NN+1),DD**(NN+1),threshold,debug,printer,count)
             
         !Infinite density matrix RG algorithm
         open(15,file="data/IDMRG_conv.txt",status='unknown')
         do ii=1,N_iter
            !Construct the total 2*NN+2 hamiltonian
            call init_identity(identity,DD,NN)
            call matrix_tens_product(identity,DD**(NN),prod,DD**2,small_temp,DD**(NN+2),count)
            call matrix_tens_product(small_temp,DD**(NN+2),identity,DD**NN,temp,DD**(2*NN+2),count)
            hamiltonian=temp
            deallocate(identity)
            call init_identity(identity,DD,NN+1)
            call matrix_tens_product(ham_left,DD**(NN+1),identity,DD**(NN+1),temp,DD**(2*NN+2),count)
            hamiltonian=hamiltonian+temp
            call matrix_tens_product(identity,DD**(NN+1),ham_right,DD**(NN+1),temp,DD**(2*NN+2),count)
            hamiltonian=hamiltonian+temp
            deallocate(identity)
            !Diagonalize the Hamiltonian
            call diagonalize(hamiltonian,DD**(2*NN+2),w,info)
            call check_eigen(info,debug,printer,count)
            gs=w(1)/dble(2*NN+2*ii)
            !Save the ground state
            groundstate=hamiltonian(:,1)
            !Check convergence at lambda=0
            if(ll==1) then
               write(15,*) ii,gs
            endif
            !Compute its density matrix
            call density_matrix(groundstate,DD**(2*NN+2),gsdensity)
            !Trace out half of the system
            call right_trace(gsdensity,DD**(2*NN+2),density_left,DD**(NN+1),DD**(NN+1),count)
            call left_trace(gsdensity,DD**(2*NN+2),density_right,DD**(NN+1),DD**(NN+1),count)
            !check if the two matrices are equal, as expected due to symmetry
            threshold=0.00001
            call check_diff(density_left,density_right,DD**(NN+1),DD**(NN+1),threshold,debug,printer,count)
            !Diagonalize density_left
            call diagonalize(density_left,DD**(NN+1),w2,info)
            call check_eigen(info,debug,printer,count)
            !Define the projector. We take the greatest m=DD**NN eigenvalues
            do jj=1,DD**NN
               projector(:,jj)=density_left(:,DD**(NN+1)-jj+1)
               call check_norm(density_left(:,jj),DD**(NN+1),0.0001d0,debug,printer,count)
            enddo
            !Project the hamiltonian of left part
            reduced_ham=matmul(transpose(conjg(projector)),matmul(ham_left,projector))
            !Add the reduced hamiltonian
            call init_identity(identity,DD,int(1,16))
            call matrix_tens_product(reduced_ham,DD**NN,identity,DD,least_temp,DD**(NN+1),count)
            ham_left=least_temp
            deallocate(identity)           
            !Enlarge the left system by adding 1 spin to the right
            call init_identity(identity,DD,NN)
            call matrix_tens_product(identity,DD**NN,sigmaz,DD,least_temp,DD**(NN+1),count)
            ham_left=ham_left+lambda*least_temp
            !Finally add the interaction term
            call matrix_tens_product(identity,DD**NN,sigmax,DD,least_temp,DD**(NN+1),count)
            reduced_ham=matmul(transpose(conjg(projector)),matmul(least_temp,projector)) !use reduced_ham as temporary
            call matrix_tens_product(reduced_ham,DD**NN,sigmax,DD,least_temp,DD**(NN+1),count)
            deallocate(identity)
            ham_left=ham_left+least_temp
            ham_right=ham_left
            call check_diff(ham_left,ham_right,DD**(NN+1),DD**(NN+1),threshold,debug,printer,count)            
         enddo
         close(15)
         deallocate(ham_left,ham_right)
         write(12,*) lambda,(w(hh)/dble(2*NN+2*(ii-1)),hh=1,kk)     
      enddo

      close(12)
      
      !Free memory
      deallocate(w,w2,hamiltonian,reduced_ham,temp,small_temp,least_temp)
      deallocate(groundstate,gsdensity,density_left,density_right,projector,sigmax,sigmaz,prod)
            
      print *, "Total number of errors:", count

    end program IDMRG
    
    
    
    

      
