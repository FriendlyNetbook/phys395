# Compiler Configuration$
FC = gfortran$
FFLAGS = -O3 -fdefault-real-8$
#LDFLAGS =

# Dependencies$
All: Pt1 CoeffPt2 Pt3.pdf Pt4Results

CoeffPt2: Pt2Coeff_N10_lobatto Pt2Coeff_N100_lobatto Pt2Coeff_N10_uniform Pt2Coeff_N100_uniform
Pt2Coeff_N10_lobatto: lobattoPt2_N10
	./lobattoPt2_N10 > Pt2Coeff_N10_lobatto
Pt2Coeff_N100_lobatto: lobattoPt2_N100
	./lobattoPt2_N100 > Pt2Coeff_N100_lobatto
Pt2Coeff_N10_uniform: uniformPt2_N10
	./uniformPt2_N10 > Pt2Coeff_N10_uniform
Pt2Coeff_N100_uniform: uniformPt2_N100
	./uniformPt2_N100 > Pt2Coeff_N100_uniform

data_Pt3_lobatto_N100: lobattoPt3_N100
	./lobattoPt3_N100 > data_Pt3_lobatto_N100
data_Pt3_uniform_N100: uniformPt3_N100
	./uniformPt3_N100 > data_Pt3_uniform_N100
data_Pt3_lobatto_N10: lobattoPt3_N10
	./lobattoPt3_N10 > data_Pt3_lobatto_N10
data_Pt3_uniform_N10: uniformPt3_N10
	./uniformPt3_N10 > data_Pt3_uniform_N10


Pt3.pdf: data_Pt3_uniform_N100 data_Pt3_lobatto_N100 data_Pt3_uniform_N10 data_Pt3_lobatto_N10
	gnuplot plotP3

Pt4Results: lobattoPt4 uniformPt4
	./lobattoPt4 > Pt4Results
	printf "\n" >> Pt4Results
	./uniformPt4 >> Pt4Results


# Generic Compilation Rules$
%: %.f90
	$(FC) $(FFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm  Pt1 uniformPt2_N10 uniformPt2_N100 uniformPt4 lobattoPt2_N10 lobattoPt2_N100 lobattoPt3_N10 uniformPt3_N10 lobattoPt3_N100 uniformPt3_N100 lobattoPt4 data* Pt2Coeff* *.pdf Pt4Results
