# Example Makefile for gfort

# Compiler Configuration
FC = gfortran
FFLAGS = -O3 -fdefault-real-8 
#LDFLAGS = 

# Dependencies
all: fibonacci fibonacci.dat fibonacci.pdf

fibonacci.dat: fibonacci
	./fibonacci > $@


fibonacci.pdf: fibonacci.dat plot.gpl
	gnuplot plot.gpl	


# Generic Compilation Rules
%: %.f90
	$(FC) $(FFLAGS) $^ -o $@ $(LDFLAGS) 
