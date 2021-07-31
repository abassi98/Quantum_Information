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
            
         end 

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
         end


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
         end 

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
      end
         !This subroutine checks if the computation time exceeds a certain threshold.
         subroutine check_time(start,threshold,debug,printer)
            implicit none
            logical :: debug,printer
            real :: start, finish,threshold
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
         end
         
      end module 
      
      module random_matrix

      contains

      subroutine init_Hmatrix(A,dmn,status)
      !initialize a single precision random hermitian matrix
      implicit none
      complex, dimension(:,:),allocatable :: A !hermitian matrix
      integer :: dmn !matrix's dimension
      integer :: ii,jj !iterators
      real :: xx,yy             !random numbers
      integer :: status
      !allocate the memory
      allocate(A(dmn,dmn),stat=status)
      !Inititalize the hermitian matrix with uniform random numbers in [-1,1]
      do jj=1,dmn
         do ii=1,jj
            call random_number(xx)
            call random_number(yy)
            A(ii,jj)=complex(2*xx-1,2*yy-1)
            A(jj,ii)=conjg(A(ii,jj))
            if(ii==jj) A(ii,jj)=complex(2*xx-1,0)
         enddo
      enddo
      end

      !Inititalize a double precision daigonal random matrix
      subroutine init_Dmatrix(A,dmn,status)
      implicit none
      real, dimension(:,:), allocatable :: A
      integer :: dmn
      integer :: ii,jj          !iterator
      real :: xx
      integer :: status

      allocate(A(dmn,dmn),stat=status)
      do ii=1,dmn
         do jj=1,dmn
            if(ii==jj) then
                call random_number(xx)
                A(ii,ii)=2*xx-1 !uniform distribition in [-1,1]
             else
                A(ii,jj)=0.0
             endif
         enddo
      enddo
      end

      subroutine get_prob(inter,prob,nbins,xmin,xmax) !compute the probability densiy function of given a vector
      implicit none
      integer :: nbins !number of bins
      real, dimension(nbins) :: prob !Prob distribution discretized
      real :: xmin,xmax, bwidth !interval and width of bins
      real :: val !value to be read from file
      integer :: ii,jj             !iterator
      integer :: count          !frequency counting
      integer,intent(in) :: inter
      bwidth=(xmax-xmin)/nbins*1.0
      !Print and end flag in every file.Here -1.
      if(inter==0) then
         open(10,file="data/H_spac.txt",status='old',access='append')
         write(10,*) -1.0
      elseif(inter==1) then
         open(10,file="data/D_spac.txt",status='old',access='append')
         write(10,*) -1.0
      else
         stop
      endif
      close(10)
      !Compute the probability density function  by counting the frequency in 
      do jj=1,nbins
         if(inter==0) then
            open(10,file="data/H_spac.txt",status='old')
         elseif(inter==1) then
            open(10,file="data/D_spac.txt",status='old')
         else
            stop
         endif
         count=0
         ii=1
         do 
            read(10,*) val
         if(xmin+bwidth*(jj-1)<=val.and.val<xmin+bwidth*jj) then
            count=count+1
         endif
         ii=ii+1
         if(val<0) exit
         enddo
         close(10)
         prob(jj)=1.0*count/(bwidth*ii)
      enddo
     
      end

       subroutine get_probr(inter,prob,nbins,xmin,xmax) !compute the probability densiy function of given a vector
      implicit none
      integer :: nbins !number of bins
      real, dimension(nbins) :: prob !Prob distribution discretized
      real :: xmin,xmax, bwidth !interval and width of bins
      real :: val !value to be read from file
      integer :: ii,jj             !iterator
      integer :: count          !frequency counting
      integer,intent(in) :: inter
      bwidth=(xmax-xmin)/nbins*1.0
      !Print an ending flag (-1.) in each file
      if(inter==0) then
         open(10,file="data/H_r.txt",status='old',access='append')
         write(10,*) -1.0
      elseif(inter==1) then
         open(10,file="data/D_r.txt",status='old',access='append')
         write(10,*) -1.0
      else
         stop
      endif
      close(10)
      !Compute the probability density function
      do jj=1,nbins
         if(inter==0) then
            open(10,file="data/H_r.txt",status='old')
         elseif(inter==1) then
            open(10,file="data/D_r.txt",status='old')
         else
            stop
         endif
         count=0
         ii=1
         do 
            read(10,*) val
         if(xmin+bwidth*(jj-1)<=val.and.val<xmin+bwidth*jj) then
            count=count+1
         endif
         ii=ii+1
         if(val<0) exit
         enddo
         close(10)
         prob(jj)=1.0*count/(bwidth*ii)
      enddo
     
      end

      subroutine mean(average,inter) !compute the average of a file
      implicit none
      real :: val !value to be read from file
      integer :: ii,jj             !iterator
      real :: sum,average
      integer,intent(in) :: inter
      !Print an ending flag
      if(inter==0) then
         open(10,file="data/H_r.txt",status='old',access='append')
         write(10,*) -1.0
      elseif(inter==1) then
         open(10,file="data/D_r.txt",status='old',access='append')
         write(10,*) -1.0
      else
         stop
      endif
      close(10)
     !Compute the sum of the elements of a file
      if(inter==0) then
         open(10,file="data/H_r.txt",status='old')
      elseif(inter==1) then
         open(10,file="data/D_r.txt",status='old')
      else
         stop
      endif
      sum=0.0
      ii=1
      do 
         read(10,*) val
         sum=sum+val
         ii=ii+1
         if(val<0) exit
      enddo
      close(10)
      average=sum/ii*1.0
     
      end
      

      subroutine norm_spacings(n_spac,eig_val,dmn,div)
      implicit none
      integer :: dmn !dimension of the matrix
      integer :: div !how large is the interval around lambda_i
      integer :: ii    !iterators integers
      real :: delta             !average spacing local or not
      integer :: ii_min, ii_max
      real,dimension(dmn) :: eig_val !eigenvalues vector
      real, dimension(dmn-1) :: n_spac !output of normalized spacings
      !Compute norm_spacings. Here div=k
      if(div>dmn) div=dmn
      if(div<2) div=2
      do ii=1,dmn-1
         if(mod(div,2)==1) then
            ii_max=ii+(div-1)/2 
            ii_min=ii-(div-1)/2
         else
            ii_max=ii+div/2-1
            ii_min=ii-div/2
         endif
         if(ii_max>dmn) then
            ii_max=dmn
            ii_min=dmn-div+1
         else if(ii_min<1) then
            ii_min=1
            ii_max=div
            
         endif
         delta=(eig_val(ii_max)-eig_val(ii_min))/(div-1)
         n_spac(ii)=(eig_val(ii+1)-eig_val(ii))/delta
      enddo
      end

      end module
