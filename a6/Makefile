# Compiler Configuration$
FC = gfortran$
FFLAGS = -O3 -fdefault-real-8$
LDFLAGS = -llapack$

All: q2.gif q3.gif

q2.gif: q2data.dat
	gnuplot plotq2
q3.gif: q3data.dat
	gnuplot plotq3

q2data.dat: q1n2
	./q1n2

q3data.dat: q3
	./q3

# Generic Compilation Rules$
%: %.f90
	$(FC) $(FFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm *.dat *.gif q1n2 q3
