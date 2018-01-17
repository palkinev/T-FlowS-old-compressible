#//////////////////////////////////////////////////////////////////////#
#              ___________ __________ _______________  ___ 
#              \__    ___/ \______   \\_   _____/\   \/  /  
#                |    |     |       _/ |    __)_  \     /   
#                |    |     |    |   \ |        \ /     \   
#                |____|     |____|_  //_______  //___/\  \  
#                                  \/         \/       \_/  
#                             Divisor Makefile
#----------------------------------------------------------------------#
#  You must have enviroment variable COMPILER set to specific platform.
#  OR
#  compile like this:
#  make 
#  make DEBUG=yes
#----------------------------------------------------------------------#
# User options
COMPILER = gnu
DEBUG    = no
TESTS    = no

# where to store compiled .o and .mod files
DOBJ  = ./obj/
DMOD  = ./mod/
MKDIRS := $(shell mkdir -p $(DOBJ) $(DMOD))

# mod files dirs (for compiler)
HEADER_DIRS = $(DMOD)

vpath %.f90 ../Library/ ./ # source .f90 files search pattern
vpath %.c   ../Library/ ./ # source .c   files search pattern
#-------------------------
  PROGRAM = divide
#-------------------------

#=============================
#....Fortran compiler & linker
#-----------------------------

  ifeq ($(COMPILER),gnu)
    FC     = mpif90 
    FL     = $(FC)
    ifeq ($(DEBUG),yes)
    	FCOPTS = -g -O0 -J$(DMOD) -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -fbacktrace -fmax-errors=0 -Wall -Wno-array-temporaries -Warray-bounds -Wcharacter-truncation -Wline-truncation -Wconversion-extra -Wimplicit-interface -Wimplicit-procedure -Wunderflow -Wextra -Wuninitialized -fcheck=all -fmodule-private -fimplicit-none -ffpe-trap=invalid,overflow -fbacktrace -fdump-core -finit-real=nan #-fno-range-check  ,precision,denormal,underflow -std=f2003 -fall-intrinsics
    else
        FCOPTS =    -O3 -J$(DMOD) -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none
    endif 
    FLOPTS =    
  endif
  
  ifeq ($(COMPILER),intel)
    FC     = mpiifort 
    FL     = $(FC)
    ifeq ($(DEBUG),yes)
    	FCOPTS = -g -O0 -module $(DMOD) -r8 -warn all -check all -debug all -extend-source 132 -fpe-all=0 -fp-stack-check -fstack-protector-all -ftrapuv -no-ftz -traceback -gen-interfaces -std03 -ipo -inline all -ipo-jobs4 -vec-report1
    else
        FCOPTS =    -O3 -module $(DMOD) -r8 -traceback
    endif 
    FLOPTS =    
  endif

#//////////////////////////////////////////////////////////////////////#

#....Fortran sources
#-------------------
  FSRCS := allp_mod.f90 \
            all_mod.f90 \
            par_mod.f90 \
            div_mod.f90 \
            pro_mod.f90 \
            gen_mod.f90 \
            Divisor.f90 \
            DivLoa.f90   Factor.f90   Number.f90   CasSav.f90 \
            IniDiv.f90   Split.f90    RISort.f90   Dist.f90 \
            ReadC.f90    EpsPar.f90   GenSav.f90   ComSav.f90 \
            GeoLoa.f90   GeoSav.f90   TestLn.f90   NamFil.f90 \
            INSort.f90   BCelLoa.f90  ToUppr.f90 \
            Inertia.f90  Jacobi.f90   CouMat.f90   RNSort.f90

#....C sources
#-------------
  CSRCS := logo_div.c 

#-----
  FOBJS := $(FSRCS:%.f90=$(DOBJ)%.o)
  COBJS := $(CSRCS:%.c=$(DOBJ)%.o)
  
  

#//////////////////////////////////////////////////////////////////////#
all      : $(PROGRAM)
#----------------------------------------------------------------------#
$(PROGRAM): $(FOBJS) $(COBJS)
	@echo "   Linking $(PROGRAM) ..."
	$(FL) -o $(PROGRAM) $(FOBJS) $(COBJS) $(FLOPTS)

$(FOBJS): $(DOBJ)%.o: %.f90
	$(FC) -c $(FCOPTS) $(addprefix -I, $(HEADER_DIRS)) -o $@ $<

$(COBJS): $(DOBJ)%.o: %.c
	$(CC) -c $(CCOPTS) $(addprefix -I, $(HEADER_DIRS)) -o $@ $<

#======================================================================#
#...explicit COMPILERs
#-----------------------------
clean:
	rm -r $(DOBJ) $(DMOD) $(PROGRAM)