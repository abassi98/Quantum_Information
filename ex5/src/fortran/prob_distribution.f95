      program prob_distribution
      use debugging
      use random_matrix
      implicit none
      logical :: debug, printer !debugging variables
      integer :: nbins          !number of bins
      real :: xmin,xmax,binwidth !interval and width of bins
      real, dimension(:),allocatable :: prob
      integer :: ii             !iterator
      integer :: count !error counter
      integer :: inter
      integer :: status !allocation status
      real :: average
      
      open(30,file="fortran/temp/debug.txt",status='unknown')
      read(30,*) debug
      close(30)
      open(30,file="fortran/temp/printer.txt",status='unknown')
      read(30,*) printer
      close(30)
      
      !Allocate probability density function vector prob
      allocate(prob(nbins),stat=status)
      call check_allocation(status,debug,printer)

      !Read nbins,xmin,xmax from files
      open(30,file="fortran/temp/nbins.txt",status='unknown')
      read(30,*) nbins
      close(30)
      open(30,file="fortran/temp/xmin.txt",status='unknown')
      read(30,*) xmin
      close(30)
      open(30,file="fortran/temp/xmax.txt",status='unknown')
      read(30,*) xmax
      close(30)

      binwidth=(xmax-xmin)/nbins*1.0 !set binwidth

      inter=0 !file "H_spac.txt"
      !Compute the probability density function for hermitian matrix
      call get_prob(inter,prob,nbins,xmin,xmax)
      open(30,file="data/H_prob.txt",status='unknown')
      do ii=1,nbins
         write(30,*) (ii*1.0-0.5)*binwidth,prob(ii)
      enddo
      close(30)

      inter=1 !file "D_spac.txt"
      !Compute probability density function for a real diagonal matrix
      call get_prob(inter,prob,nbins,xmin,xmax)
      open(30,file="data/D_prob.txt",status='unknown')
      do ii=1,nbins
         write(30,*) (ii*1.0-0.5)*binwidth,prob(ii)
      enddo
      close(30)

      !Compute the distribution and mean of
      xmin=0.
      xmax=1.
      binwidth=(xmax-xmin)/nbins*1.0 !set binwidth
      inter=0                   !"H_r.txt"
      call mean(average,inter)
      print *, "Average of r(H):",average
      call get_probr(inter,prob,nbins,xmin,xmax)
      open(30,file="data/Hr_prob.txt",status='unknown')
      do ii=1,nbins
         write(30,*) (ii*1.0-0.5)*binwidth,prob(ii)
      enddo
      close(30)
      
      inter=1
      call mean(average,inter)
      print *, "Average of r(D):",average
      call get_probr(inter,prob,nbins,xmin,xmax)
      open(30,file="data/Dr_prob.txt",status='unknown')
      do ii=1,nbins
         write(30,*) (ii*1.0-0.5)*binwidth,prob(ii)
      enddo
      close(30)
      end
