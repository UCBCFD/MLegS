# Compiler
FC = gfortran
# Compiler flags
FFLAGS = -O3 -fdefault-double-8 -g -J$(MOD_DIR)

# Source and build directories
BUILD_DIR = build
MOD_DIR = $(BUILD_DIR)/mod
OBJ_DIR = $(BUILD_DIR)/obj

# Find all Fortran source files
SOURCES = $(wildcard *.f90)

# Create object file names
OBJECTS = $(patsubst %.f90,$(OBJ_DIR)/%.o,$(SOURCES))

# Target static library
LIB = libfm.a

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
	ar rcs $@ $^

# Compile objects
$(OBJ_DIR)/%.o: %.f90
	@mkdir -p $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -c -o $@ $<

# Clean up build files
clean:
	rm -rf $(LIB) $(BUILD_DIR)

.PHONY: all clean
