#!/bin/bash

#this bash-script will prepare VTK-Fortran for compiling
# Source files can be downloaded from github https://github.com/szaghi?tab=repositories

# remove old files
rm -r mod *.a
rm -rf BeFoR64 FoXy PENF StringiFor Lib_VTK_IO

mkdir mod

# Lib_VTK_IO
git clone https://github.com/szaghi/Lib_VTK_IO

# Dependencies:

# BeFoR64
git clone https://github.com/szaghi/BeFoR64
mkdir Lib_VTK_IO/src/third_party/BeFoR64/src
mkdir Lib_VTK_IO/src/third_party/BeFoR64/src/lib

# FoXy
git clone https://github.com/Fortran-FOSS-Programmers/FoXy
mkdir Lib_VTK_IO/src/third_party/FoXy/src
mkdir Lib_VTK_IO/src/third_party/FoXy/src/lib

# PENF
git clone https://github.com/szaghi/PENF
mkdir Lib_VTK_IO/src/third_party/PENF/src
mkdir Lib_VTK_IO/src/third_party/PENF/src/lib

# StringiFor (can be only compiled with ifort or gfortran-6+)
git clone https://github.com/szaghi/StringiFor
mkdir Lib_VTK_IO/src/third_party/StringiFor/src
mkdir Lib_VTK_IO/src/third_party/StringiFor/src/lib

# copy sources
cp BeFoR64/src/lib/befor64.F90                       Lib_VTK_IO/src/third_party/BeFoR64/src/lib
cp PENF/src/lib/penf.F90                             Lib_VTK_IO/src/third_party/BeFoR64/src/lib
cp PENF/src/lib/penf_global_parameters_variables.F90 Lib_VTK_IO/src/third_party/PENF/src/lib
cp PENF/src/lib/penf_b_size.F90                      Lib_VTK_IO/src/third_party/PENF/src/lib
cp PENF/src/lib/penf_stringify.F90                   Lib_VTK_IO/src/third_party/PENF/src/lib
cp BeFoR64/src/lib/befor64_pack_data_m.F90           Lib_VTK_IO/src/third_party/BeFoR64/src/lib
cp FoXy/src/lib/foxy*                                Lib_VTK_IO/src/third_party/FoXy/src/lib
cp StringiFor/src/lib/stringifor*                    Lib_VTK_IO/src/third_party/StringiFor/src/lib

# compile
cd Lib_VTK_IO

# this method does not require root access
#export MPICH_F90=gfortran-6
#sed -i -e 's%gfortran%gfortran-6%g' makefile

# release version
make COMPILER=gnu DEBUG=no OPTIMIZE=yes MPI=yes SHARED=no STATIC=yes
cp static/libvtkfortran.a ../libvtkfortran.a
cp -r static/mod/ ../

make clean

# debug version
make COMPILER=gnu DEBUG=yes OPTIMIZE=no MPI=yes SHARED=no STATIC=yes
cp -r static/libvtkfortran.a ../libvtkfortrand.a

#PENF lib
cd ../PENF

make COMPILER=gnu
cp lib/libpenf.a ../libpenf.a
cp -r lib/mod/ ../

cd ..

# another way to change default compiler:
# 
# sudo update-alternatives --remove-all gfortran #remove all alternatives
# 
# sudo update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/gfortran-4.4 50 # set alternative 1, 50 % priority
# sudo update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/gfortran-5   30 # set alternative 2, 30 % priority
# sudo update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/gfortran-6   20 # set alternative 3, 20 % priority
# 
# sudo update-alternatives --config gfortran # interactively change gfortran compiler
# 
# sudo update-alternatives --set gfortran /usr/bin/gfortran-4.4 # manually change gfortran compiler to gfortran-4.4
# sudo update-alternatives --set gfortran /usr/bin/gfortran-5   # manually change gfortran compiler to gfortran-5
# sudo update-alternatives --set gfortran /usr/bin/gfortran-6   # manually change gfortran compiler to gfortran-6
