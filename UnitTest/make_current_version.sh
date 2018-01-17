#!/bin/bash

# this bash-script will prepare FRUIT and FRUITpy
# Source files can be downloaded from github https://github.com/acroucher/FRUITPy
# Dependency: FRUIT lib in http://sourceforge.net/projects/fortranxunit/

# remove old files
rm -rf mod/ *.a *.so *.py *.pyc

# FRUITpy interface
git clone https://github.com/acroucher/FRUITPy

# make by python
cd FRUITPy
python setup.py build
cp FRUIT.py ..
cd ..

# FRUIT lib
#wget https://downloads.sourceforge.net/project/fortranxunit/fruit_3.4.1/fruit_3.4.1.zip; unzip -qo fruit_3.4.1.zip
mkdir mod/
tar -zxf fruit_3.4.1.tar.gz

# shared lib
# cp FRUITPy/fruit_makefile fruit_3.4.1/makefile

cd fruit_3.4.1/src/

#  old version compilation
#mv fruit_util.f90_2015mar fruit_util.f90
#mv fruit.f90_2012jun fruit.f90

#sed -i -e 's%$(HOME)%.%g' makefile

cd ..



# static lib

cat << EOF > makefile
F90=.f90
OBJ=.o
LIB = libfruit.a
SRC = src
BUILD = \$(SRC)

# compiler options:
FC = mpif90
FLAGS = -c -O3

SOURCES = \$(wildcard \$(SRC)/*\$(F90))
OBJS = \$(patsubst \$(SRC)/%\$(F90), \$(BUILD)/%\$(OBJ), \$(SOURCES))

\$(LIB) : \$(OBJS)
		ar cr \$(LIB) \$(BUILD)/*.o
		ranlib \$(LIB)

\$(BUILD)/%\$(OBJ): \$(SRC)/%\$(F90)
		\$(FC) \$(FLAGS) -c \$< -o \$@
EOF

make

cp libfruit.* ../
cp *.mod  ../mod/

cd ..

# DONE




# calculator test
# Correct output
# 
# All tests passed.
# Hit rate:
#   asserts:  1 / 1 (100%)
#   cases  :  1 / 1 (100%)

cd fruit_3.4.1/in_3_minutes
cp ../../FRUIT.py .

cat << EOF > test.py
from FRUIT import *

test_modules = ['calculator_test.f90']
driver = "calculator_test_driver.f90"

suite = test_suite(test_modules)
suite.build_run(driver)
suite.summary()
EOF

 
cat << EOF > makefile
FC = gfortran
LIBS = ../../libfruit.a
INCLS = -I../../mod/
FCFLAGS = -O3 -Wall -ffree-line-length-none

EXE = calculator_test_driver

\$(EXE): calculator.o calculator_test.o

%: %.o
	\$(FC) \$(FCFLAGS) -o \$@ \$^ \$(LIBS)
%.o: %.f90
	\$(FC) \$(FCFLAGS) \$(INCLS) -c \$<

.PHONY: clean
clean:
	rm -f *.o *.mod \$(EXE)

EOF

python test.py

cd ../../

# delete temp files
#rm -rf FRUITPy fruit_3.4.1