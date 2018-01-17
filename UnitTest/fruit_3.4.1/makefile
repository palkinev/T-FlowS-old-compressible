F90=.f90
OBJ=.o
LIB = libfruit.a
SRC = src
BUILD = $(SRC)

# compiler options:
FC = mpif90
FLAGS = -c -O3

SOURCES = $(wildcard $(SRC)/*$(F90))
OBJS = $(patsubst $(SRC)/%$(F90), $(BUILD)/%$(OBJ), $(SOURCES))

$(LIB) : $(OBJS)
		ar cr $(LIB) $(BUILD)/*.o
		ranlib $(LIB)

$(BUILD)/%$(OBJ): $(SRC)/%$(F90)
		$(FC) $(FLAGS) -c $< -o $@
