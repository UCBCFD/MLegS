# Compiler (ifx or gfortran)
OMPI_FC = gfortran

# Compiler flags
ifeq ($(OMPI_FC), ifx)
  FC = mpiifx
  FFLAGS = -O3 -r8 -mcmodel=medium -g -qopenmp -module $(MOD_DIR) -L$(LIB_DIR) $(EXTLIBS)
else ifeq ($(OMPI_FC), gfortran)
  FC = mpifort
  FFLAGS = -O3 -fdefault-real-8 -mcmodel=medium -g -fopenmp -J$(MOD_DIR) -L$(LIB_DIR) $(EXTLIBS)
else
  $(error Untested compiler: $(FC). Tested compilers are gfortran <v11.2 or later> and ifx <v2024.1.0 or later>)
endif

EXTLIBS = -lffte -lfm -lblas -llapack

# Source and build directories
MODULE_DIR = src/modules
SUBMODULE_DIR = src/submodules
MAIN_DIR = test
BUILD_DIR = build
MOD_DIR = $(BUILD_DIR)/mod
OBJ_DIR = $(BUILD_DIR)/obj
BIN_DIR = $(BUILD_DIR)/bin
LIB_DIR = external/lib

# Find all Fortran source files
MODULE_SOURCES = $(wildcard $(MODULE_DIR)/*.f90)
SUBMODULE_SOURCES = $(wildcard $(SUBMODULE_DIR)/*.f90)
MAIN_SOURCES = $(wildcard $(MAIN_DIR)/*.f90)

# Create object file names
MODULE_OBJECTS = $(patsubst $(MODULE_DIR)/%.f90,$(OBJ_DIR)/%.o,$(MODULE_SOURCES))
SUBMODULE_OBJECTS = $(patsubst $(SUBMODULE_DIR)/%.f90,$(OBJ_DIR)/%.o,$(SUBMODULE_SOURCES))

# Create binary file names
MAIN_EXES = $(patsubst $(MAIN_DIR)/%.f90,%,$(MAIN_SOURCES))

# Default target
all: $(MAIN_EXES)

-include Makefile.dep

# Link object files to create the executable
$(MAIN_EXES): $(MODULE_OBJECTS) $(SUBMODULE_OBJECTS)
	@mkdir -p $(MOD_DIR) $(OBJ_DIR) $(BIN_DIR)
	$(FC) $(FFLAGS) $(MAIN_DIR)/$@.f90 -o $(BIN_DIR)/$@ $^

# Compile each module source file into an object file
$(OBJ_DIR)/%.o: $(MODULE_DIR)/%.f90
	@mkdir -p $(MOD_DIR) $(OBJ_DIR)
	$(FC) $(FFLAGS) -c -o $@ $<

$(OBJ_DIR)/%.o: $(SUBMODULE_DIR)/%.f90
	@mkdir -p $(MOD_DIR) $(OBJ_DIR)
	$(FC) $(FFLAGS) -c -o $@ $<

# Compile only modules
mods: $(MODULE_OBJECTS) $(SUBMODULE_OBJECTS)

# Clean up build files
clean:
	rm -rf $(BUILD_DIR)

.PHONY: all clean mods
