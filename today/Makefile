# Compiler Configuration
FC = gfortran
FFLAGS = -O3 -fdefault-real-8 
#LDFLAGS = 

# Dependencies
all: MatrixSolve 
ortho: MatrixSolve.f90

# Generic Compilation Rules
%: %.f90
	$(FC) $(FFLAGS) $^ -o $@ $(LDFLAGS)
