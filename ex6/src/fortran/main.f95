      program harm_oscillator
      use debugging
      use quantum
      
      implicit none
      real :: xmin,xmax              !Interval
      real :: omega !pulse
      integer :: NN             !number of points of the disretization
      real, dimension(:,:), allocatable :: matrix !Our hermitian matrix
      real :: hh    !delta x
      logical :: debug, printer !logical variable for debugging
      integer :: status         !allocation status
      integer :: count          !counting errors
      integer :: kk !how many eigenvectors
      integer :: ii,jj !iterators
      !dsyev variables
      integer :: lda,lwork,info
      real, dimension(:),allocatable :: w,work
      real :: xx
      real :: norma
      
      count=0 
      open(10,file="temp/debug.txt",status='unknown')
      read(10,*) debug
      close(10)
      open(10,file="temp/printer.txt",status='unknown')
      read(10,*) printer
      close(10)
     
      
      !Read number of intervals
      open(10,file="temp/N.txt",status='unknown')
      read(10,*) NN
      close(10)
      call check_dimensions(NN,NN,debug,printer,count) !check dimension
     
 
      !Read xmax
      open(12,file="temp/xmax.txt",status='unknown')
      read(12,*) xmax
      close(12)

      xmin=-xmax !symmetric interval
      !Read omega
      open(13,file="temp/omega.txt",status='unknown')
      read(13,*) omega
      close(13)

      !Read number of eigenvectors
      open(13,file="temp/k.txt",status='unknown')
      read(13,*) kk
      close(13)
      hh=(xmax-xmin)/NN*1.0

      !Initialize harmonic potential matrix and check allocation
      call harm_matrix(matrix,NN,xmin,omega,hh,status)
      call check_allocation(status,debug,printer)

      !allocate and check dsyev variables
      allocate(w(NN+1),stat=status)
      call check_allocation(status,debug,printer)
      lda=max(1,NN+1)
      lwork=max(1,3*NN+2)
      allocate(work(lwork),stat=status)
      call check_allocation(status,debug,printer)

      !Compute eigenvalues and eigenvectors
      call ssyev('V','U',NN+1,matrix,NN+1,w,work,lwork,info)
      call check_eigen(info,debug,printer,count)

      !normalize the eigenfunctions because dsyev returns orthonormalized eigenvectors
 
      !Print and save eigenvalues and eigenvectors in a file
      open(15,file="data/eigenvalues.txt",status='unknown')
      open(16,file="data/eigenvectors.txt",status='unknown')

      matrix=matrix/sqrt(hh)
      do  ii=1,NN+1
         xx=xmin+(ii-1)*1.0*hh
        write(16,*)xx,matrix(ii,1),matrix(ii,2),(matrix(ii,jj), jj=3,kk)
         write(15,*) ii-1,w(ii)/(hh**2)
      enddo

      
      close(15)
      close(16)

      !check the normalizations
      do ii=1,kk
      print *,"Norm of state",ii-1,":",norm(matrix,NN,kk,hh)
      enddo
      
      !Free memory and print number of errors encountered
      deallocate(matrix,w,work)
      print *, "Error encountered:", count
      
      end

      
