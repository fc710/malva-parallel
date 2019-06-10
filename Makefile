#CXX = $(shell icc -v 2>&1 >/dev/null | grep -o "icc")
CXXFLAGS = -Wall -g -std=c++14 -O3 -fopenmp 

ifeq ($(CXX), icc)
 CXXFLAGS += -xCORE_AVX2
else
 #CXX = clang++
CXX = g++
endif

INCLUDES = -I. -I./sdsl-lite/build/include -I./htslib/htslib -I./KMC
LDFLAGS = -L./sdsl-lite/build/lib -L./sdsl-lite/build/external/libdivsufsort/lib -L./htslib/
LDLIBS = -lhts -lz -lsdsl -ldivsufsort -ldivsufsort64 -ltbb
OBJS = main_exp.o MurmurHash3.o ./KMC/kmc_api/kmc_file.o ./KMC/kmc_api/kmer_api.o ./KMC/kmc_api/mmer.o
PROG = malva-geno

.PHONY: all

all: $(PROG)

$(PROG): $(OBJS) 
	@echo "* Linking $@"
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@ $(LDFLAGS) $(LDLIBS)

%.o: %.cpp
	@echo '* Compiling $<'
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ -c $<

clean:
	rm -rf *.o $(PROG)
	rm -rf KMC/kmc_api/*.o
