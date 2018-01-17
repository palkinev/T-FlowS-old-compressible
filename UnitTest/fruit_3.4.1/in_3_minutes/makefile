FC = gfortran
LIBS = ../../libfruit.a
INCLS = -I../../mod/
FCFLAGS = -O3 -Wall -ffree-line-length-none

EXE = calculator_test_driver

$(EXE): calculator.o calculator_test.o

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LIBS)
%.o: %.f90
	$(FC) $(FCFLAGS) $(INCLS) -c $<

.PHONY: clean
clean:
	rm -f *.o *.mod $(EXE)

