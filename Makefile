###################################################################
# This Makefile was created using the bat-project script
# for project Cygno_PMTs_fit
# bat-project is part of Bayesian Analysis Toolkit (BAT).
# BAT can be downloaded from http://mpp.mpg.de/bat
###################################################################
#
# Run 'make' to compile the program and 'make clean' to remove
# all compiled parts and 'clean' the directory.
#
# You might need to adjust the CXXFLAGS and LIBS based on
# the BAT installation on your system. Consult the gmake manual
# for details.
#
###################################################################

# Directories
SRC_DIR = src
BUILD_DIR = build
BIN_DIR = bin

# List of all class (model) sources used in the program,
# separated by spaces. A backslash indicates continuation
# on the next line
CXXSRCS = $(SRC_DIR)/PMT_association.cpp \
          $(SRC_DIR)/helper_lib.cpp \
          $(SRC_DIR)/PMT_calibration.cpp \
          $(SRC_DIR)/PMT_FindAlpha.cpp

# List of all program sources used in the program,
# separated by spaces. A backslash indicates continuation
# on the next line
PRGSRCS = $(SRC_DIR)/runfit.cpp

TARGET = $(BIN_DIR)/fitter

# compiler and flags
CXX       = g++
CXXFLAGS  = -g -O2 -Wall -fPIC -Wno-deprecated
LD        = /usr/bin/ld -m elf_x86_64
LDFLAGS   = -g -O2  -fopenmp


# ----------------------------------------------------------------------
# The following definitions rely on the script bat-config being
# available in $PATH. If BAT is not installed in the standard system
# directories, update $PATH accordingly.

CXXFLAGS += $(shell bat-config --cflags)
LIBS := $(shell bat-config --libs)

#--------------------------------------------------

# ----------------------------------------------------------------------
# don't change lines below unless you know what you're doing
#

CXXOBJS  = $(addprefix $(BUILD_DIR)/,$(addsuffix .o,$(basename $(notdir $(CXXSRCS)))))
PRGOBJS  = $(addprefix $(BUILD_DIR)/,$(addsuffix .o,$(basename $(notdir $(PRGSRCS)))))
MYPROGS  = $(BIN_DIR)/$(basename $(notdir $(PRGSRCS)))

GARBAGE  = $(CXXOBJS) $(PRGOBJS) link.d $(MYPROGS)

# targets
all : $(MYPROGS)

.PHONY : all clean print

link.d: $(CXXSRCS) $(PRGSRCS)
	$(CXX) -MM $(CXXFLAGS) $(CXXSRCS) $(PRGSRCS) > link.d
	@$(foreach prog,$(MYPROGS), echo $(prog) : $(BUILD_DIR)/$(notdir $(prog)).o >> link.d;)
-include link.d


# compilation rules
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# linking rules
$(MYPROGS) : $(CXXOBJS) $(PRGOBJS)
	$(CXX) $(LDFLAGS) $^ $(LIBS) -o $@

clean :
	rm -f $(GARBAGE)

print :
	@echo "compiler  : $(CXX)"
	@echo "c++ srcs  : $(CXXSRCS) $(PRGSRCS)"
	@echo "c++ objs  : $(CXXOBJS) $(PRGOBJS)"
	@echo "c++ flags : $(CXXFLAGS)"
	@echo "ld flags  : $(LDFLAGS)"
	@echo "libs      : $(LIBS)"
