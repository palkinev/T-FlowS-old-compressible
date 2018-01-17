#//////////////////////////////////////////////////////////////////////#
#              ___________ __________ _______________  ___ 
#              \__    ___/ \______   \\_   _____/\   \/  /  
#                |    |     |       _/ |    __)_  \     /   
#                |    |     |    |   \ |        \ /     \   
#                |____|     |____|_  //_______  //___/\  \  
#                                  \/         \/       \_/  
#                            Processor Makefile
#----------------------------------------------------------------------#
#  You must have enviroment variable COMPILER set to specific platform.
#  OR
#  compile like this:
#  make 
#  make DEBUG=yes
#----------------------------------------------------------------------#

#-----------------
# User options
#-----------------

COMPILER      = gnu
DEBUG         = no
TESTS         = no

#-----------------
# External Libs
#-----------------

#-----------------
# Fortran VTK format IO lib (supports ASCII and Base64)
VTKFortranLib = no
vtk_lib = ../VTKFortran/VTKFortran-current
# https://github.com/szaghi/VTKFortran
# gcc 5+ (Fortran 2003+) only
# or
# ifort ????

# C lib with fortran interfaces for IO CGNS format
CGNSLib = yes
# http://cgns.github.io/
cgns_lib = ../Library/CGNS/CGNS/install_dir
# depends on:
# HDF5 lib
hdf5_lib = ../Library/CGNS/HDF5/install_dir
# https://support.hdfgroup.org/HDF5/examples/intro.html#fortran

#----------------------------------------------------------------------#

ifeq ($(VTKFortranLib),yes)
	# static libs dirs (for linker)
	ifeq ($(DEBUG),yes)
		STATIC_LIBS = $(vtk_lib)/libvtkfortrand.a
	else
		STATIC_LIBS = $(vtk_lib)/libvtkfortran.a 
	endif
	STATIC_LIBS += $(vtk_lib)/libpenf.a
endif

ifeq ($(CGNSLib),yes)
	STATIC_LIBS += $(cgns_lib)/lib/libcgns.a
	STATIC_LIBS += $(hdf5_lib)/lib/libhdf5.a
	STATIC_LIBS += -Wl,--no-as-needed -ldl -lz
endif


# where to store compiled .o and .mod files
DOBJ  = ./obj/
DMOD  = ./mod/
MKDIRS := $(shell mkdir -p $(DOBJ) $(DMOD))

# mod files dirs (for compiler)
HEADER_DIRS = $(DMOD) 

ifeq ($(VTKFortranLib),yes)
	HEADER_DIRS += $(vtk_lib)/mod/
endif
ifeq ($(CGNSLib),yes)
	HEADER_DIRS += $(cgns_lib)/include/
	HEADER_DIRS += $(hdf5_lib)/include/
endif

ifeq ($(PRECISION),single)
	INCLUDE_DIRS = ../Parallel/Single
else
	INCLUDE_DIRS = ../Parallel/Double
endif

ifeq ($(FORCE),seq)
	INCLUDE_DIRS = ../Sequential
endif

vpath %.f90 ../Library/ ../User ./  $(INCLUDE_DIRS) # source .f90 files search pattern
vpath %.c   ../Library/ ../User ./  $(INCLUDE_DIRS) # source .c   files search pattern

#-------------------------
PROGRAM = TFLowS
#-------------------------

#=============================
#....Fortran compiler & linker
#-----------------------------

ifeq ($(COMPILER),gnu)
	FC     = mpif90
	#FC = /home/l_palkin_e/eclipse/compressible/eclipse-project/CGNS/MPICH/install_dir/bin/mpif90
	FL     = $(FC)
	ifeq ($(DEBUG),yes)
		FCOPTS = -g\
		-O0\
		-J $(DMOD)\
		-fdefault-real-8 -fdefault-double-8\
		-ffree-line-length-none\
		-fbacktrace\
		-fmax-errors=0\
		-Wall\
		-Wno-array-temporaries\
		-Warray-bounds\
		-Wcharacter-truncation\
		-Wline-truncation\
		-Wimplicit-interface\
		-Wunderflow\
		-Wextra\
		-Wuninitialized\
		-fmodule-private\
		-fimplicit-none\
		-ffpe-trap=invalid,overflow -fdump-core -finit-real=nan
		#-Wconversion-extra -Wimplicit-procedure -fcheck=all
		#-fno-range-check  ,precision,denormal,underflow -std=f2003 -fall-intrinsics
	else
		FCOPTS =    -O3 -J $(DMOD) -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none
	endif 
	FLOPTS =
endif

ifeq ($(COMPILER),intel)
	FC     = mpiifort
	FL     = $(FC)
	ifeq ($(DEBUG),yes)
		FCOPTS = -g -O0 -module $(DMOD) -r8 -traceback -warn all -check all -debug all -extend-source 132 -fpe-all=0 -fp-stack-check -fstack-protector-all -ftrapuv -no-ftz -traceback -gen-interfaces -std03 -ipo -inline all -ipo-jobs4 -vec-report1
	else
		FCOPTS =    -O3 -module $(DMOD) -r8
	endif
	FLOPTS =    
endif

#//////////////////////////////////////////////////////////////////////#

#....Fortran sources
#-------------------
FSRCS :=	\
			allp_mod.f90\
			all_mod.f90\
			pro_mod.f90\
			les_mod.f90\
			par_mod.f90\
			sol_mod.f90\
			rans_mod.f90\
			two_column_file_with_x_and_y.f90\
			vd1d_mms_mod.f90\
			vd1d_mms_mod_table.f90\
			vd2d_mms_mod_2.f90\
			vd2d_mms_mod_3.f90\
			vd2d_mms_mod_4_RT.f90\
			vd2d_mms_mod_5_burner.f90\
			moin_problem_mod.f90\
			Base64_Mod.f90\
			ReadC.f90\
			CnsLoa.f90\
			UnkAloc.f90\
			IniPar.f90\
			IniVar.f90\
			StaPar.f90\
			EndPar.f90\
			Timex.f90\
			BiCG.f90\
			CorUVW.f90\
			CGS.f90\
			BouLoa.f90\
			Resid.f90\
			Calc3.f90\
			Processor.f90\
			CalBou.f90\
			GraPhi.f90\
			ReaCom.f90\
			TopolM.f90\
			CalcMn.f90\
			CG.f90\
			GradP.f90\
			SavRes.f90\
			CalcPS.f90\
			LoaRes.f90\
			CalcG.f90\
			NewUVW.f90\
			CalcShear.f90\
			Dist.f90\
			GeoAloc.f90\
			ToUppr.f90\
			GeoLoa.f90\
			ISort.f90\
			NamFil.f90\
			Prec1.f90\
			GradP3.f90\
			BufLoa.f90\
			Exchng.f90\
			GloMax.f90\
			GloMin.f90\
			GloSum.f90\
			ModOut.f90\
			IGlSum.f90\
			FindBad.f90\
			CorBad.f90\
			CalcFlux.f90\
			CalcConvect.f90\
			Wait.f90\
			Approx.f90\
			CalMinMax.f90\
			ConvScheme.f90\
			CorBou.f90\
			Calc_Sgs_Coefficients_Dynamic_Smagorinsky.f90\
			UserPerturb2.f90\
			UserForce.f90\
			UserGraPhi.f90\
			GraCorNew.f90\
			CalcVort.f90\
			NearWallCell.f90\
			LoaIni.f90\
			SavIni.f90\
			UserCutLines_channel.f90\
			SSort.f90\
			SavParView.f90\
			UserCutLines_annulus.f90\
			EOS.f90\
			CalcZmix.f90\
			CalcDens.f90\
			Mass_Fraction_Equation_Source.f90\
			Service_output_Processor.f90\
			SaveVTK_ascii_base64.f90\
			#Rnd.f90\
			#init_random_seed.f90

ifeq ($(VTKFortranLib),yes)
	FSRCS += SaveVTK_VTKFortran_Lib.f90
endif

ifeq ($(CGNSLib),yes)
	FSRCS += Save_Mesh_CGNS_Lib.f90
	#FSRCS += Add_BC_To_Mesh_CGNS_Lib.f90
	FSRCS += Add_Fields_To_Mesh_CGNS_Lib.f90
endif

#....C sources
#-------------
CSRCS := logo_pro.c

#-----
FOBJS := $(FSRCS:%.f90=$(DOBJ)%.o)
COBJS := $(CSRCS:%.c=$(DOBJ)%.o)
  
  

#//////////////////////////////////////////////////////////////////////#
all      : $(PROGRAM)

#----------------------------------------------------------------------#
$(PROGRAM): $(FOBJS) $(COBJS)
	@echo "   Linking $(PROGRAM) ..."
	$(FL) -o $(PROGRAM) $(FOBJS) $(COBJS) $(FLOPTS) $(STATIC_LIBS)

$(FOBJS): $(DOBJ)%.o: %.f90
	$(FC) -c $(FCOPTS)  $(addprefix -I, $(HEADER_DIRS)) -o $@ $<

$(COBJS): $(DOBJ)%.o: %.c
	$(CC) -c $(CCOPTS) $(addprefix -I, $(HEADER_DIRS)) -o $@ $<

#======================================================================#
#...explicit targets
#-----------------------------
clean:
	rm -r $(DOBJ) $(DMOD) $(PROGRAM)
