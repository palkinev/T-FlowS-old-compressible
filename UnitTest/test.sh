#!/bin/bash

# Runs FRUITPy unit tests

cat << EOF > test.py
from FRUIT import *

#test_modules = ["Library/allp_mod_test.f90 Library/all_mod_test.f90"]
test_modules = ["Library/allp_mod_test.f90", "Library/all_mod_test.f90", \
				"Process/Calc_Sgs_Coefficients_Dynamic_Smagorinsky_test.f90"]
driver = "test_driver.f90"
#build_command="make test_driver"

# run_command="./allp_mod"

suite = test_suite(test_modules)

# build test
#suite.build(driver, build_command)

suite.build_run(driver)

# run test
#suite.run(driver, run_command)

suite.summary()

EOF

cat << EOF > makefile
FC = mpif90
LIBS = libfruit.a
INCLS = -Imod
#FCFLAGS = -O3 -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -fbacktrace -Wall -fmax-errors=0 -Wno-array-temporaries -Warray-bounds -Wcharacter-truncation -Wline-truncation -Wconversion-extra -Wimplicit-interface -Wimplicit-procedure -Wunderflow -Wextra -Wuninitialized -ffpe-trap=invalid,overflow -fdump-core -finit-real=nan -fall-intrinsics -fcheck=all
FCFLAGS = -O3 -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -fbacktrace -Wall -fmax-errors=0 -Wno-array-temporaries -Warray-bounds -Wcharacter-truncation -Wline-truncation -Wconversion-extra -Wimplicit-interface -Wimplicit-procedure -Wunderflow -Wextra -Wuninitialized -ffpe-trap=invalid,overflow -fdump-core -finit-real=nan -fall-intrinsics
vpath %.f90 ../Library/\
			   Library\
	 	    ../Process/\
	 	       Process/\
	 	    ../Parallel/Double/ # search path

EXE = test_driver

\$(EXE): test_tools.o\
 		 allp_mod.o\
 		 allp_mod_test.o\
 		 all_mod.o\
 		 all_mod_test.o\
 		 pro_mod.o\
 		 les_mod.o\
 		 rans_mod.o\
 		 par_mod.o\
 		 vd1d_mms_mod.o\
 		 two_column_file_with_x_and_y.o\
 		 vd1d_mms_mod_table.o\
 		 vd2d_mms_mod_2.o\
 		 vd2d_mms_mod_3.o\
 		 vd2d_mms_mod_4_RT.o\
 		 moin_problem_mod.o\
 		 GloSum.o\
 		 GloMin.o\
 		 Exchng.o\
 		 GraPhi.o\
 		 CorBad.o\
 		 ReadC.o\
 		 Calc_Sgs_Coefficients_Dynamic_Smagorinsky.o\
 		 Calc_Sgs_Coefficients_Dynamic_Smagorinsky_test.o
%: %.o
	\$(FC) \$(FCFLAGS) -o \$@ \$^ \$(LIBS)
%.o: %.f90
	\$(FC) \$(FCFLAGS) \$(INCLS) -c \$<

.PHONY: clean
clean:
	rm -f *.o *.mod \$(EXE) *_driver.f90 *.pyc

EOF

python test.py
make clean
rm makefile
