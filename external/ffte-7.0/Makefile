# Compiler (ifx or gfortran)
FC = ifx

# Compiler flags
ifeq ($(FC), ifx)
  FFLAGS = -O3 -g -qopenmp -fomit-frame-pointer -f77rtl
else ifeq ($(FC), gfortran)
  FFLAGS = -O3 -g -fopenmp -fomit-frame-pointer -std=legacy
else
  $(error Untested compiler: $(FC). Tested compilers are gfortran <v11.2 or later> and ifx <v2024.1.0 or later>)
endif

# Source and build directories
BUILD_DIR = build
OBJ_DIR = $(BUILD_DIR)/obj

# Find all Fortran source files
SOURCES = $(wildcard *.f)
HEADERS = $(wildcard *.h)

# Create object file names
OBJECTS = $(patsubst %.f,$(OBJ_DIR)/%.o,$(SOURCES))

# Target static library
LIB_DIR = $(BUILD_DIR)/lib
LIB = $(LIB_DIR)/libffte.a

# Default target
all: $(OBJECTS) $(LIB)

# Create static library from object files
$(LIB) : $(OBJECTS)
	@mkdir -p $(LIB_DIR)
	ar rcs $@ $^

# Compile objects
$(OBJ_DIR)/%.o: %.f $(HEADERS)
	@mkdir -p $(OBJ_DIR)
	$(FC) $(FFLAGS) -c -o $@ $<

# Clean up build files
clean:
	rm -rf $(BUILD_DIR)

.PHONY: all clean
