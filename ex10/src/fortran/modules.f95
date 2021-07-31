      
      

      module quantum
      contains

      subroutine harm_matrix(matrix,NN,xmin,omega,hh,status)
      implicit none
      double precision, dimension(:,:), allocatable :: matrix
      integer*16 :: NN,status
      double precision :: xmin,omega,hh,xx
      integer*16 :: ii !iterator
      
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
      integer*16 :: NN,ii
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
        integer*16 :: NN,ii
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
        integer*16 :: NN, ii, flag,count
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
        integer*16 ::NN !dimension of vector
        double complex, dimension(NN) :: wave
        double precision :: omega, xmin, hh,V,dt,time,speed
        integer*16 :: ii
        do ii=1,NN
           V=0.5*((xmin+hh/2+hh*dble(ii-1)-time*speed)*omega)**2
           wave(ii)=exp(-cmplx(0d0,1d0)*0.5*V*dt)*wave(ii)
        enddo

      end subroutine potential

      !This subroutine applies the kinetic matrix to p-space vector
      subroutine kinetic(wave,NN,pmax,dt)
        implicit none
        integer*16 ::NN !dimension of vector-1
        double complex, dimension(NN) :: wave
        double precision :: pmax,K,dt
        integer*16 :: ii
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

      module manybodyQS
      contains

        
        !This subroutine diagonalize a hermitian matrix
        subroutine diagonalize(matrix,dmn,w,info)
          implicit none
          integer*16::dmn
          integer:: info,lwork,lda
          double complex, dimension(dmn,dmn):: matrix
          double precision, dimension(dmn) :: w
          double precision,dimension(:),allocatable:: rwork
          double complex, dimension(:), allocatable :: work
          lda=max(1,dmn)
          lwork=max(1,2*dmn-1)
          allocate(work(max(1,lwork)))
          allocate(rwork(max(1,3*dmn-2)))
          call zheev('V','U',dmn,matrix,lda,w,work,lwork,rwork,info)
          deallocate(work,rwork)
        end subroutine diagonalize
        

        !THis subroutine initializes a generic pure state of dimension D^N
        subroutine genpure_state(state,DD,NN,status)
          implicit none
          double complex, dimension(:),allocatable :: state
          integer*16 :: status,DD,NN,ii
          double precision :: re,im
          allocate(state(DD**NN),stat=status)
          do ii=1,DD**NN
             call random_number(re)
             call random_number(im)
             state(ii)=cmplx(2*re-1,2*im-1)
          enddo
        end subroutine genpure_state

       

        !This subroutine computes the trace of a generic double complex matrix
        double complex function trace(mat,dim)
          implicit none
          integer*16 :: dim, ii
          double complex, dimension(dim,dim) :: mat
          trace=cmplx(0d0,0d0)
          do ii=1,dim
             trace=trace+mat(ii,ii)
          enddo
          return
        end function trace

        !THis  subroutine comnpute Von Neumann entropy of a generic state
        !We take the absolute value inside the log because the eigenvalues can become negative due to numerical instabilities
        subroutine entropy_vn(entropy,rho,dim)
          implicit none
          integer*16 :: dim,lda,info,lwork,status,ii
          double complex, dimension(dim,dim):: rho
          double complex, dimension(:),allocatable :: work
          double precision, dimension(:),allocatable :: w,rwork
          double precision :: entropy
          lda=max(1,dim)
          lwork=max(1,2*dim-1)
          allocate(w(dim))
          allocate(work(lwork))
          allocate(rwork(max(1,3*dim-2)))
          call zheev('N','U',dim,rho,lda,w,work,lwork,rwork,info)
          entropy=0d0
          do ii=1,dim
             entropy=entropy-w(ii)*log(abs(w(ii)))
          enddo
          deallocate(w,work,rwork)
        end subroutine entropy_vn
        
        
        !This subroutine computes the tensor product of two pure states
        subroutine state_tens_product(stateA,dimA,stateB,dimB,prod,dim,count)
          implicit none
          integer*16 :: dimA,dimB,dim,count,ii,jj
          double complex, dimension(dimA) :: stateA
          double complex, dimension(dimB) :: stateB
          double complex, dimension(dim) :: prod
          if(dimA*dimB.ne.dim) then
             print *, "Match dimensions failed in tensor product: ERROR"
             count=count+1
          endif
          do ii=0,dimA-1
             do jj=0,dimB-1
                prod(ii*dimB+jj+1)=stateA(ii+1)*stateB(jj+1)
             enddo
          enddo
        end subroutine state_tens_product

        !This suborutine computes the tensor product between two generic denisty matrices
        subroutine matrix_tens_product(mat1,dim1,mat2,dim2,out,dim,count)
          implicit none
          integer*16 :: dim1,dim2,dim,count,jj,ii,ll,kk
          double complex, dimension(dim1,dim1) :: mat1
          double complex, dimension(dim2,dim2) :: mat2
          double complex, dimension(dim,dim):: out
          if(dim1*dim2.ne.dim) then
             print *, "Subsystems dimensions product different: ERROR"
             count=count+1
          endif
          do jj=0,dim1-1
             do  ii=0,dim1-1
                do kk=0,dim2-1
                   do ll=0,dim2-1
                      out(ii*dim2+ll+1,jj*dim2+kk+1)=mat1(ii+1,jj+1)*mat2(ll+1,kk+1)
                   enddo
                enddo
             enddo
          enddo
        end subroutine matrix_tens_product
        
        !This subroutine computes the denisty matrix of a pure state
        subroutine density_matrix(state,dim,density)
          implicit none
          integer*16 :: dim,ii,jj
          double complex, dimension(dim) :: state
          double complex, dimension(dim,dim) :: density
          do jj=0,dim-1
             do ii=0,dim-1
                density(ii+1,jj+1)=state(ii+1)*conjg(state(jj+1))
             enddo
          enddo
        end subroutine density_matrix

        !This subroutine computes the reduced density matrix over the right part of the system
        subroutine right_trace(density,dim,densityA,dimA,dimB,count)
          implicit none
          integer*16 :: dim, dimA,dimB,ii,jj,kk,count,flag
          double complex, dimension(dim,dim) :: density
          double complex, dimension(dimA,dimA) :: densityA
          if(dimA*dimB.ne.dim) then
             print *, "Subsystems dimensions product different: ERROR"
             count=count+1
          endif
          do jj=0,dimA-1
             do ii=0,dimA-1
                densityA(ii+1,jj+1)=cmplx(0d0,0d0)
                do kk=0,dimB-1
                   densityA(ii+1,jj+1)=densityA(ii+1,jj+1)+density(ii*dimB+kk+1,jj*dimB+kk+1)
                enddo
             enddo
          enddo          
        end subroutine right_trace

         !This subroutine computes the reduced density matrix over the left part of the system
        subroutine left_trace(density,dim,densityB,dimA,dimB,count)
          implicit none
          integer*16 :: dim, dimA,dimB,ii,jj,kk,count,flag
          double complex, dimension(dim,dim) :: density
          double complex, dimension(dimB,dimB) :: densityB
          if(dimA*dimB.ne.dim) then
             print *, "Subsystems dimensions product different: ERROR"
             count=count+1
          endif
          do jj=0,dimB-1
             do ii=0,dimB-1
                densityB(ii+1,jj+1)=cmplx(0d0,0d0)
                do kk=0,dimB-1
                   densityB(ii+1,jj+1)=densityB(ii+1,jj+1)+density(ii*dimB+kk+1,jj*dimB+kk+1)!TO CHECK.LIKE THIS,  IT IS NOT CORRECT
                enddo
             enddo
          enddo          
        end subroutine left_trace
        
        

        !This subroutine computes the reduced density matrix over a generic subsystem w.r.t. index k=1,NN
        subroutine partial_trace_gen(density,DD,NN,reduced,kk)
          implicit none
          integer*16 :: DD,NN,kk,ii,jj, ss,alpha,beta
          double complex, dimension(DD**NN,DD**NN) :: density
          double complex, dimension(DD**(NN-1),DD**(NN-1)):: reduced
          do jj=0,DD**(NN-1)-1
             do ii=0,DD**(NN-1)-1
                reduced(ii+1,jj+1)=cmplx(0d0,0d0)
                do ss=0,DD-1
                   alpha=DD**(NN-kk)*ss+(DD**(NN-kk+1)-DD**(NN-kk))*int(ii/DD**(NN-kk))+ii
                   beta=DD**(NN-kk)*ss+(DD**(NN-kk+1)-DD**(NN-kk))*int(jj/DD**(NN-kk))+jj
                   reduced(ii+1,jj+1)=reduced(ii+1,jj+1)+density(alpha+1,beta+1)
                enddo
             enddo
          enddo
        end subroutine partial_trace_gen

        subroutine init_sigmax(sigmax)
          implicit none
          double complex, dimension(:,:),allocatable :: sigmax
          allocate(sigmax(2,2))
          sigmax=cmplx(0d0,0d0)
          sigmax(1,2)=cmplx(1d0,0d0)
          sigmax(2,1)=cmplx(1d0,0d0)
        end subroutine init_sigmax

        subroutine init_sigmay(sigmay)
          implicit none
          double complex, dimension(:,:),allocatable :: sigmay
          allocate(sigmay(2,2))
          sigmay=cmplx(0d0,0d0)
          sigmay(1,2)=cmplx(0d0,-1d0)
          sigmay(2,1)=cmplx(0d0,1d0)
        end subroutine init_sigmay
        
        subroutine init_sigmaz(sigmaz)
          implicit none
          double complex, dimension(:,:),allocatable :: sigmaz
          allocate(sigmaz(2,2))
          sigmaz=cmplx(0d0,0d0)
          sigmaz(1,1)=cmplx(1d0,0d0)
          sigmaz(2,2)=cmplx(-1d0,0d0)
        end subroutine init_sigmaz

        subroutine init_identity(identity,DD,kk)
          implicit none
          integer*16 :: DD,kk,ii
          double complex, dimension(:,:), allocatable::identity
          allocate(identity(DD**kk,DD**kk))
          identity=cmplx(0d0,0d0)
          do ii=1,DD**kk
             identity(ii,ii)=cmplx(1d0,0d0)
          enddo
        end subroutine init_identity
        

        subroutine init_ising_hamiltonian(hamiltonian,DD,NN,lambda,status,count)
          implicit none
          integer*16 :: ii,jj,kk,status,DD,NN,count
          double precision ::lambda
          double complex,dimension (:,:),allocatable::hamiltonian,temp,temp_up,sigmax,sigmaz,identity,prod
          if(NN<2) then
             print *, "Invalid number of particles. Set NN=3"
             NN=2
          endif          
          call init_sigmax(sigmax)
          call init_sigmaz(sigmaz)
          allocate(hamiltonian(DD**NN,DD**NN),stat=status) !allocate hamiltonian
          allocate(temp(DD**NN,DD**NN)) !allocate temp for swapping
          allocate(prod(DD**2,DD**2)) !compute sigmax tensprod(tp) sigmax
          call matrix_tens_product(sigmax,DD,sigmax,DD,prod,DD**2,count)
          hamiltonian=cmplx(0d0,0d0) !initialize to 0 hamiltonian  
          !Quadratic part of hamiltonian
          do kk=0,NN-2
             !compute tensor product  of first k identities with sigmx tp sigmx
             allocate(temp_up(DD**(kk+2),DD**(kk+2)))
             call init_identity(identity,DD,kk)
             call matrix_tens_product(identity,DD**kk,prod,DD**2,temp_up,DD**(kk+2),count)
             deallocate(identity)
             !Compute the second tensor product
             call init_identity(identity,DD,NN-kk-2)
             call matrix_tens_product(temp_up,DD**(kk+2),identity,DD**(NN-kk-2),temp,DD**NN,count)
             deallocate(identity,temp_up)
             !Sum to get the hamiltonian
             hamiltonian=hamiltonian+temp
          enddo
          
          !Linear part of hamiltonian
          do kk=0,NN-1
             !compute tensor product of firts part. Now sigmaz
             allocate(temp_up(DD**(kk+1),DD**(kk+1)))
             call init_identity(identity,DD,kk)
             call matrix_tens_product(identity,DD**kk,sigmaz,DD,temp_up,DD**(kk+1),count)
             deallocate(identity)
             !Compute  second part
             call init_identity(identity,DD,NN-kk-1)
             call matrix_tens_product(temp_up,DD**(kk+1),identity,DD**(NN-kk-1),temp,DD**NN,count)
             deallocate(identity,temp_up)
             !Sum to get hailtonian
             hamiltonian=hamiltonian+temp*lambda
          enddo
          
          deallocate(temp,sigmax,sigmaz,prod)
       
        end subroutine init_ising_hamiltonian

        !This subroutine computes the central finite difference derivative of a vector
        subroutine central_der(vec,der,dim,xmin,xmax)
          implicit none
          integer*16 :: dim,ii
          double  precision :: step,xmin,xmax
          double precision, dimension(dim) :: vec
          double precision, dimension(dim-2):: der
          step=(xmax-xmin)/dble(dim)
          do ii=2,dim-1
             der(ii-1)=(vec(ii+1)-vec(ii-1))/(2*step)
          enddo
        end subroutine central_der

        !This subroutine computes the second derivative of a vector
        subroutine second_der(vec,der,dim,xmin,xmax)
          implicit none
          integer*16 :: dim, ii
          double  precision :: step,xmin,xmax
          double precision, dimension(dim) :: vec
          double precision, dimension(dim-2):: der
          step=(xmax-xmin)/dble(dim)
          do ii=2,dim-1
             der(ii-1)=(vec(ii+1)+vec(ii-1)-2*vec(ii))/(step**2)
          enddo
        end subroutine second_der
        
          
      end module manybodyQS
      
           
      
