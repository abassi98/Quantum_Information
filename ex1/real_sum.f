      PROGRAM real_sum
      IMPLICIT NONE
      REAL xx,yy,sum
      xx=acos(-1.0)*(10**32.0)*1.0
      yy=sqrt(2.0)*(10**21.0)*1.0
      sum=xx+yy
      PRINT *, "The sum is: ", sum
      END
