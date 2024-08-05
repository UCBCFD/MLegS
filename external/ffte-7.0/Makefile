# Compiler
FC = gfortran
# Compiler flags
FFLAGS = -O3 -g -fopenmp -fomit-frame-pointer -std=legacy

# Source and build directories
BUILD_DIR = build
OBJ_DIR = $(BUILD_DIR)/obj

# Find all Fortran source files
SOURCES = $(wildcard *.f)
HEADERS = $(wildcard *.h)

# Create object file names
OBJECTS = $(patsubst %.f,$(OBJ_DIR)/%.o,$(SOURCES))

# Target static library
LIB = libffte.a

# Default target
all: $(OBJECTS) $(LIB)

# Create static library from object files
$(LIB) : $(OBJECTS)
	ar rcs $@ $^

# Compile objects
$(OBJ_DIR)/%.o: %.f $(HEADERS)
	@mkdir -p $(OBJ_DIR)
	$(FC) $(FFLAGS) -c -o $@ $<

# Clean up build files
clean:
	rm -rf $(LIB) $(BUILD_DIR)

.PHONY: all clean
