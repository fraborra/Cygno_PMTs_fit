# List of all class (model) sources used in the program,
# separated by spaces. A backslash indicates continuation
# on the next line
CXXSRCS = PMT_association.cpp helper_lib.cpp PMT_calibration.cpp

# List of all program sources used in the program,
# separated by spaces. A backslash indicates continuation
# on the next line
PRGSRCS = runfit.cpp

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

CXXOBJS = $(addsuffix .o,$(basename $(CXXSRCS)))
MYPROGS = $(basename $(PRGSRCS))
PRGOBJS = $(addsuffix .o,$(basename $(PRGSRCS)))

GARBAGE = $(CXXOBJS) $(PRGOBJS) link.d $(MYPROGS)

# targets
all : $(MYPROGS)

.PHONY : all clean print

link.d : $(addsuffix .hpp,$(basename $(CXXSRCS))) $(CXXSRCS) $(PRGSRCS)
	$(CXX) -MM $(CXXFLAGS) $(filter-out %.h,$^) > link.d;
	@$(foreach prog,$(MYPROGS), echo $(prog) : $(prog).o >> link.d;)

-include link.d


cose.o :
	$(CXX) $(CXXFLAGS) -c $(addsuffix .cpp,$(basename $@)) -o $@


$(MYPROGS) : $(CXXOBJS) $(PRGOBJS)
	$(CXX) $(LDFLAGS) $^ $(LIBS) -o $@

clean :
	rm -f $(GARBAGE)

print :
	@echo compiler  : $(CXX)
	@echo c++ srcs  : $(CXXSRCS) $(PRGSRCS)
	@echo c++ objs  : $(CXXOBJS) $(PRGOBJS)
	@echo c++ flags : $(CXXFLAGS)
	@echo ld flags  : $(LDFLAGS)
	@echo libs      : $(LIBS)
