      program many_body
      use debugging
      use quantum
      use manybodyQS
      
      implicit none
      !Varible declaration
      double complex, dimension(:),allocatable ::bigstate, state1,state2,state3,product12,product13,product23,product123
      double complex, dimension(:,:),allocatable ::density13, density12,density1,density2,densitypar1,densitypar13,density123,bigmat
      logical :: debug,printer
      integer :: count,status,NN,DD,ii,jj !NN is the number of systems, DD their dimension
      double precision :: entropy
      double precision, parameter :: threshold=0.0000001
      !Parameters reading
      count=0 !counting errors variable set to zero
      open(10,file="temp/debugging.txt",status='unknown')
      read(10,*) debug, printer
      close(10)
      open(10,file="temp/parameters.txt",status='unknown')
      read(10,*)  NN,DD
      close(10)
      
      call check_dimensions(DD,DD,debug,printer,count) !check positiveness dimensions
      call check_dimensions(NN,NN,debug,printer,count) !check positiveness dimensions

      !Initialize pure states of 1 particle, diemension D and normalise
      call genpure_state(state1,DD,1,status)
      call check_allocation(status,debug,printer)
      state1=state1/norm(state1,DD,1d0)
     
      call genpure_state(state2,DD,1,status)
      call check_allocation(status,debug,printer)
      state2=state2/norm(state2,DD,1d0)

      call genpure_state(state3,DD,1,status)
      call check_allocation(status,debug,printer)
      state3=state3/norm(state3,DD,1d0)
     
      !Compute the tensor products 12,23,13,123
      allocate(product12(DD**2),stat=status)
      call check_allocation(status,debug,printer)
      call tensor_product(state1,DD,state2,DD,product12,DD**2,count)     
      print *, "Norm of tensor product 12:",norm(product12,DD**2,1d0)

      allocate(product13(DD**2),stat=status)
      call check_allocation(status,debug,printer)
      call tensor_product(state1,DD,state3,DD,product13,DD**2,count)     
      print *, "Norm of tensor product 13:",norm(product13,DD**2,1d0)

      allocate(product23(DD**2),stat=status)
      call check_allocation(status,debug,printer)
      call tensor_product(state2,DD,state3,DD,product23,DD**2,count)     
      print *, "Norm of tensor product 23:",norm(product23,DD**2,1d0)

      allocate(product123(DD**3),stat=status)
      call check_allocation(status,debug,printer)
      call tensor_product(product12,DD**2,state3,DD,product123,DD**3,count)     
      print *, "Norm of tensor product 123:",norm(product123,DD**3,1d0)

      
      !Compute density matrices for checking 1,2,12,123
      allocate(density1(DD,DD),stat=status)
      call check_allocation(status,debug,printer)
      call density_matrix(state1,DD,density1)
     
      
      allocate(density2(DD,DD),stat=status)
      call check_allocation(status,debug,printer)
      call density_matrix(state2,DD,density2)
     
      allocate(density13(DD**2,DD**2),stat=status)
      call check_allocation(status,debug,printer)
      call density_matrix(product13,DD**2,density13)

      allocate(density12(DD**2,DD**2),stat=status)
      call check_allocation(status,debug,printer)
      call density_matrix(product12,DD**2,density12)

      allocate(density123(DD**3,DD**3),stat=status)
      call check_allocation(status,debug,printer)
      call density_matrix(product123,DD**3,density123)
      
      !Partial trace 123 over 2
      allocate(densitypar13(DD**2,DD**2),stat=status)
      call check_allocation(status,debug,printer)
      call partial_trace_gen(density123,DD,3,densitypar13,2)
      call check_diff(density13,densitypar13,DD**2,DD**2,threshold,debug,printer,count)


      !partial trace 12 over 2
      allocate(densitypar1(DD,DD),stat=status)
      call check_allocation(status,debug,printer)
      call partial_trace_gen(density12,DD,2,densitypar1,2)
      call check_diff(density1,densitypar1,DD,DD,threshold,debug,printer,count)

      print *, "Density matrix of state 1:"
      do ii=1,DD
         print *, (density1(ii,jj), jj=1,DD)
      enddo

      print *, "Reduced denisty matrix of state 12 over state2:"
      do ii=1,DD
         print *, (densitypar1(ii,jj), jj=1,DD)
      enddo

      !Compute  a generic pure state
      call genpure_state(bigstate,DD,NN,status)
      call check_allocation(status,debug,printer)
      bigstate=bigstate/norm(bigstate,DD**NN,1d0)

      allocate(bigmat(DD**NN,DD**NN),stat=status)
      call check_allocation(status,debug,printer)
      call density_matrix(bigstate,DD**NN,bigmat)

      call entropy_vn(entropy,bigmat,DD**NN)

      print *, "Von Neumann entropy of D^N pure state", entropy
      
      print *, "Total number of errors:",count

    end program many_body
    

      
