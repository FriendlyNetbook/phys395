# Compiler Configuration$
FC = gfortran$
FFLAGS = -O3 -fdefault-real-8$
LDFLAGS = -lcfitsio$

# Dependencies$
All:fractal1.fit tracePlot.pdf energyErrorPlot.pdf

traceData: q12
	./q12

energyError: q12
	./q12

fractal1.fit: q3
	./q3

tracePlot.pdf: traceData
	gnuplot plotTrace

energyErrorPlot.pdf: energyError
	gnuplot plotEnergyErr

# Generic Compilation Rules$
%: %.f90
	$(FC) $(FFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm  fractal1.fit fractal2.fit fractal3.fit fractal4.fit q12 q3 traceData energyError animateData *.pdf
