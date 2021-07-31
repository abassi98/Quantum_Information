
FC = gfortran
FCLAGS = -g -O3 -fbounds-check #-ffpe-trap=zero,overflow,underflow 
FFLAGS = -llapack -lfftw3 -lm

SOURCE1 = fortran/modules.f95 fortran/module_debugging.f95 fortran/debugging.f95  

SOURCE2 = fortran/modules.f95 fortran/module_debugging.f95 fortran/RSRG.f95

SOURCE3 = fortran/modules.f95 fortran/module_debugging.f95 fortran/IDMRG.f95



all: debugging RSRG IDMRG



debugging: $(SOURCE1)
	$(FC) $(FCLAGS) -o $@ $(SOURCE1) $(FFLAGS)

RSRG: $(SOURCE2)
	$(FC) $(FCLAGS) -o $@ $(SOURCE2) $(FFLAGS) 

IDMRG :  $(SOURCE3)
	$(FC) $(FCLAGS) -o $@ $(SOURCE3) $(FFLAGS)



clean:
	-rm -f *.mod *.txt debugging RSRG IDMRG

