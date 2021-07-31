
FC = gfortran
FCLAGS = -g -O3 -fbounds-check #-ffpe-trap=zero,overflow,underflow 
FFLAGS = -llapack -lfftw3 -lm

SOURCE1 = fortran/modules.f95 fortran/debugging.f95 

SOURCE2 = fortran/modules.f95 fortran/main.f95





all: debugging main



debugging: $(SOURCE1)
	$(FC) $(FCLAGS) -o $@ $(SOURCE1) $(FFLAGS)

main: $(SOURCE2)
	$(FC) $(FCLAGS) -o $@ $(SOURCE2) $(FFLAGS) 



clean:
	-rm -f *.mod *.txt debugging main

