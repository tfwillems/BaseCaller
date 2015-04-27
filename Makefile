##
## Makefile for StutterTrainer
##

## Default compilation flags.
## Override with:
##   make CXXFLAGS=XXXXX
CXXFLAGS= -O3 -g -D_FILE_OFFSET_BITS=64 -std=c++0x -DMACOSX

## Source code files, add new files to this list
SRC_COMMON  = error.cpp io.cpp base_caller.cpp

# For each CPP file, generate an object file
OBJ_COMMON  := $(SRC_COMMON:.cpp=.o)

.PHONY: all
all: BaseCaller

# Clean the generated files of the main project only (leave Bamtools/vcflib alone)
.PHONY: clean
clean:
	rm -f *.o *.d BaseCaller *~

# Clean all compiled files, including bamtools/vcflib
.PHONY: clean-all
clean-all: clean

# The GNU Make trick to include the ".d" (dependencies) files.
# If the files don't exist, they will be re-generated, then included.
# If this causes problems with non-gnu make (e.g. on MacOS/FreeBSD), remove it.
include $(subst .cpp,.d,$(SRC))

BaseCaller: $(OBJ_COMMON)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

# Build each object file independently
%.o: %.cpp $(BAMTOOLS_LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

# Auto-Generate header dependencies for each CPP file.
%.d: %.cpp $(BAMTOOLS_LIB)
	$(CXX) -c -MP -MD $(CXXFLAGS) $(INCLUDE) $< > $@

