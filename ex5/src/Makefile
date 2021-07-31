
FC = gfortran
FCLAGS = -g -O3 -fbounds-check #-ffpe-trap=zero,overflow,underflow 
FFLAGS = -llapack -lfftw3 -lm

SOURCE1 = fortran/modules.f95 fortran/debugging.f95 

SOURCE2 = fortran/modules.f95 fortran/norm_spacings.f95

SOURCE3 = fortran/modules.f95 fortran/prob_distribution.f95



all: debugging norm_spacings prob_distribution



debugging: $(SOURCE1)
	$(FC) $(FCLAGS) -o $@ $(SOURCE1) $(FFLAGS)

norm_spacings: $(SOURCE2)
	$(FC) $(FCLAGS) -o $@ $(SOURCE2) $(FFLAGS) 

prob_distribution: $(SOURCE3)
	$(FC) $(FCLAGS) -o $@ $(SOURCE3) $(FFLAGS) 



clean:
	-rm -f *.mod *.txt debugging norm_spacings prob_distribution

