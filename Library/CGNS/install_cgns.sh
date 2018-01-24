#!/bin/bash

# current dir which will contain CNGS/ HDF5/ and optionally MPICH/ folders
CGNS_DIR=$PWD
INSTALL_DIR=$CGNS_DIR/install_dir
SRC_DIR=$CGNS_DIR/src_dir

# decide if you wish to build mpich by yourself
BUILD_MPI=true
# decide if you wish to build gui cgns tools (cgnsview)
CGNS_TOOLS=true

# put your compilers here (gcc, gfortran are allowed if $MPI is built anyway)
export CC="gcc";
export FC="gfortran"; # or ifort, mpif90, mpifort, mpiifiort

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

mkdir -p $SRC_DIR/
mkdir -p $INSTALL_DIR/

#------ MPICH 3.2.1
if [ $BUILD_MPI == true ]; then

	cd $SRC_DIR/

	# download
	wget -N http://www.mpich.org/static/downloads/3.2.1/mpich-3.2.1.tar.gz; tar -zxvf mpich-3.2.1.tar.gz; mv mpich-3.2.1/ $SRC_DIR/MPICH/;
		
	# configure
	cd $SRC_DIR/MPICH/
	./configure \
	--prefix=$INSTALL_DIR/MPICH \
	--enable-fast=all,O3 \
	--enable-fortran=all \
	--disable-shared \
	--disable-dependency-tracking

	# build
	make

	# install
	make install

	# and now your mpicc and mpif90 compilers are:
	export CC=$INSTALL_DIR/MPICH/bin/mpicc
	export FC=$INSTALL_DIR/MPICH/bin/mpif90

	# return
	cd $CGNS_DIR
fi

#------ HDF5 5.1.8 (Paraview 5.4.1 & Visit 2.12.3 work with HDF5 5.1.8 and not with 5.1.10)

cd $SRC_DIR/

# download
git clone --depth=1  https://bitbucket.hdfgroup.org/scm/hdffv/hdf5.git --branch hdf5_1_8 ./HDF5;

# configure
cd HDF5/; rm -rf .git
FCFLAGS=-O3 \
./configure \
--prefix=$INSTALL_DIR/HDF5 \
--enable-fortran \
--enable-parallel \
--disable-shared \
--enable-production

# build
make

# install
make install

cd $CGNS_DIR

if [ $CGNS_TOOLS == true ]; then
#------ TCL

	cd $SRC_DIR/

	tar -zxvf tcl8.6.8-src.tar.gz; mv tcl8.6.8 TCL/; cd TCL/unix/

	# configure
	./configure \
	--prefix=$INSTALL_DIR/TCL

	# build
	make
    # no need to make install
    make install

	cd $CGNS_DIR

#------ TK (requires libx11-dev from repo)

	cd $SRC_DIR/

	tar -zxvf tk8.6.8-src.tar.gz; mv tk8.6.8/ TK/; cd TK/unix/

	# configure
	./configure \
	--prefix=$INSTALL_DIR/TK \
	--with-tcl=$INSTALL_DIR/TCL/lib/

    # build
	make
    # no need to make install
    make install

	cd $CGNS_DIR

	cp -r $INSTALL_DIR/TK/lib $INSTALL_DIR/TK/unix/

fi

#------ CGNS (latest version)

	cd $SRC_DIR/

	# download
	git clone --depth=1  https://github.com/CGNS/CGNS.git CGNS/; cd CGNS/; rm -rf .git; cd src/

if [ $CGNS_TOOLS == true ]; then

	# configure
	FLIBS=-Wl,--no-as-needed\ -ldl\ -lz \
	LIBS=-Wl,--no-as-needed\ -ldl\ -lz \
	./configure \
	--prefix=$INSTALL_DIR/CGNS \
	--with-hdf5=$INSTALL_DIR/HDF5 \
	--with-fortran \
	--enable-lfs \
	--enable-64bit \
	--disable-shared \
	--disable-debug \
	--with-zlib \
	--enable-64bit \
	--enable-parallel \
	--enable-cgnstools \
	--with-tcl=$SRC_DIR/TCL \
	--with-tk=$SRC_DIR/TK \
	--datarootdir=$INSTALL_DIR/CGNS/tcl_scripts

else # no cgns GUI tools

	FLIBS=-Wl,--no-as-needed\ -ldl\ -lz \
	LIBS=-Wl,--no-as-needed\ -ldl\ -lz \
	./configure \
	--prefix=$CGNS_DIR/CGNS/ \
	--with-hdf5=$CGNS_DIR/HDF5/ \
	--with-fortran \
	--enable-lfs \
	--enable-64bit \
	--disable-shared \
	--disable-debug \
	--with-zlib \
	--enable-64bit \
	--enable-parallel \
	--disable-cgnstools
fi

# build
make
# install
make install

# run tests
cd tests                                   # for self-confidence
make                                       # for self-confidence
make test                                  # for self-confidence
cd ../examples/fortran                     # for self-confidence
make                                       # for self-confidence
make test                                  # for self-confidence
cd ../../Test_UserGuideCode/Fortran_code   # for self-confidence
make                                       # for self-confidence
make test                                  # for self-confidence
cd ../C_code                               # for self-confidence
make                                       # for self-confidence
make test                                  # for self-confidence


cd $CGNS_DIR

#---------

echo CNGS is now built in "$CGNS_DIR/CGNS/install_dir"