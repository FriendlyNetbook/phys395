# Compiler Configuration$
FC = gfortran$
FFLAGS = -O3 -fdefault-real-8$
LDFLAGS = -llapack

# Dependencies$
All: plotFit.pdf

plotFit.pdf: FitEvaluation
	gnuplot plotFit	

FitEvaluation: fit
	(./fit < "Assignment #2.dat") > Results

# Generic Compilation Rules$
%: %.f90
	$(FC) $(FFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm plotFit.pdf FitEvaluation Results fit
