# Compiler Configuration$
FC = gfortran$
FFLAGS = -O3 -fdefault-real-8$
LDFLAGS =

# Dependencies$
All: ConvergencePath.pdf results.txt

ConvergencePath.pdf: data.dat
	gnuplot plotConvergence

data.dat: q4and5
	./q4and5 > data.dat

results.txt: q4and5 q2and3 q1
	./q4and5
	echo "Question 3" >> results.txt
	./q2and3 >> results.txt
	echo "Question 1" >> results.txt
	./q1 >> results.txt

# Generic Compilation Rules$
%: %.f90
	$(FC) $(FFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm results.txt ConvergencePath.pdf q1 q2and3 q4and5 data.dat
