# Compiler (ifx or gfortran)
FC = gfortran

# Compiler flags
ifeq ($(FC), ifx)
  FFLAGS = -O3 -g -module $(MOD_DIR)
else ifeq ($(FC), gfortran)
  FFLAGS = -O3 -fdefault-double-8 -g -J$(MOD_DIR)
else
  $(error Untested compiler: $(FC). Tested compilers are gfortran <v11.2 or later> and ifx <v2024.1.0 or later>)
endif

# Source and build directories
BUILD_DIR = build
MOD_DIR = $(BUILD_DIR)/mod
OBJ_DIR = $(BUILD_DIR)/obj

# Find all Fortran source files
SOURCES = $(wildcard *.f90)

# Create object file names
OBJECTS = $(patsubst %.f90,$(OBJ_DIR)/%.o,$(SOURCES))

# Target static library
LIB_DIR = $(BUILD_DIR)/lib
LIB = $(LIB_DIR)/libfm.a

# Configure build directory
ifneq ($(OBJ_DIR),)
  $(shell test -d $(OBJ_DIR) || mkdir -p $(OBJ_DIR))
endif
ifneq ($(MOD_DIR),)
  $(shell test -d $(MOD_DIR) || mkdir -p $(MOD_DIR))
endif

# Default target
all: $(OBJECTS) $(LIB)

-include Makefile.dep

# Create static library from object files
$(LIB) : $(OBJECTS)
	@mkdir -p $(LIB_DIR)
	ar rcs $@ $^

# Compile objects
$(OBJ_DIR)/%.o: %.f90
	@mkdir -p $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -c -o $@ $<

# Clean up build files
clean:
	rm -rf $(BUILD_DIR)

.PHONY: all clean
