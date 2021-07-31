     
      program eigenvalues
      use random_matrix
      use debugging
      implicit none
      complex, dimension(:,:),allocatable :: H  !Complex matrix
      real, dimension(:,:), allocatable :: D   !Real matrix
      real, dimension(:), allocatable :: w, rwork, n_spac,r
      complex, dimension(:), allocatable :: work
      integer :: NN,info,div,lda,lwork !Matrix size, info, division and others for  cheev
      integer :: ii,jj
      real :: xx,yy
      logical :: debug          !debugging variable
      logical :: printer !if true, it prints debugging mesagges
      integer :: status !status of allocation
      integer :: count          !counting errors
      integer :: nbins          !number of bins of the histogram
      real :: xmin, xmax,binwidth !intervals of prob and binwidth
      count=0                   !INITIALIZE COUNTING ERRORS

      
      
      ! Read debug and printer from files
      open(30,file="fortran/temp/debug.txt",status='unknown')
      read(30,*) debug
      close(30)
      open(30, file="fortran/temp/printer.txt",status="unknown")
      read(30,*) printer
      close(30)
      
      
      !Read dimension from a file
      open(30, file="fortran/temp/dimension.txt",status='unknown')
      read(30,*) NN
      close(30)


      call  check_dimensions(NN,NN,debug,printer,count)
      
      ! Read division from a file
      open(30,file="fortran/temp/division.txt",status='unknown')
      read(30,*) div
      close(30)
      
      !Iinitialize the hermitian radnom matrix
      call init_Hmatrix(H,NN,status)
      call check_allocation(status,debug,printer)
      lda=NN
      lwork=max(1,2*NN-1)

      !allocate the memory and check. Normally no printed messages
      allocate(w(NN),stat=status)
      call check_allocation(status,debug,printer)
      allocate(n_spac(NN-1),stat=status)
      call check_allocation(status,debug,printer)
      allocate(work(max(1,lwork)),stat=status)
      call check_allocation(status,debug,printer)
      allocate(rwork(max(1,3*NN-2)),stat=status)
      call check_allocation(status,debug,printer)

      !Compute the vector r
      allocate(r(NN-2),stat=status)
      call check_allocation(status,debug,printer)

      !Compute the eigenvalues and check the corretcness of computation
      call cheev('N','U',NN,H,lda,w,work,lwork,rwork,info)
      call check_eigen(info,debug,printer,count)

      !Compute the normalized spacings for the hermitian matrix.
      call norm_spacings(n_spac,w,NN,div)
      open(30,file="data/H_spac.txt",status='unknown',access='append')
      open(31,file="data/H_r.txt",status='unknown',access='append')
      do ii=1,NN-1
         write(30,*) n_spac(ii)
         if(ii<=NN-3) then
         r(ii)=min(n_spac(ii+1),n_spac(ii))/max(n_spac(ii+1),n_spac(ii))
         write(31,*) r(ii)
         endif
      enddo
      close(30)
      close(31)

      !initialize a diagonal real random matrix and check allocation
      call init_Dmatrix(D,NN,status)
      call check_allocation(status,debug,printer)

      ! Reallocate the necessaru memory
      deallocate(work,w)
      lwork=max(1,3*NN-1) !now lwork must be it
      allocate(w(NN),stat=status)
      call check_allocation(status,debug,printer)
      allocate(work(max(1,lwork)),stat=status)
      call check_allocation(status,debug,printer)

      !Compute and check eigenvaluesfor a  diagonal real matrix
      call ssyev('N','U',NN,D,lda,w,work,lwork,info)
      call check_eigen(info,debug,printer,count)

      !Compute and save in a file the normalized spacings
      call norm_spacings(n_spac,w,NN,div)
      open(30,file="data/D_spac.txt",status='unknown',access='append')
      open(31,file="data/D_r.txt",status='unknown',access='append')
      do ii=1,NN-1
         write(30,*) n_spac(ii)
         if(ii<=NN-3) then
         r(ii)=min(n_spac(ii+1),n_spac(ii))/max(n_spac(ii+1),n_spac(ii))
         write(31,*) r(ii)
         endif
      enddo
      close(30)
      close(31)
      
      print *, "Error encountered: ", count
      
      !Free memory
      deallocate(H,D,w,n_spac,work,rwork)
      end

      
