# Compiler Configuration$
FC = gfortran$
FFLAGS = -O3 -fdefault-real-8$
LDFLAGS = -llapack$

All: RayleighEigenvalues.txt RayleighConvergence.pdf q3.normalizedProbabilityDistributions.pdf q2.normalizedProbabilityDistributions.pdf q1.typicalSolutions.pdf

RayleighEigenvalues.txt: q5
	./q5

RayleighConvergence.pdf: RConv.dat
	gnuplot plotRConv

RConv.dat: q4
	./q4

q3.normalizedProbabilityDistributions.pdf: q3
	./q3
	gnuplot plotQ3

q2.normalizedProbabilityDistributions.pdf: q2
	./q2
	gnuplot plotQ2


q1.typicalSolutions.pdf: q1
	./q1
	gnuplot plotQ1

# Generic Compilation Rules$
%: %.f90
	$(FC) $(FFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm *.dat *.pdf R*.txt

