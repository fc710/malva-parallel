CXX = $(shell icpc -v 2>&1 >/dev/null | grep -o "icpc")
CXXFLAGS = -Wall -g -std=c++17 -O3 

ifeq ($(CXX), icpc)
 CXXFLAGS += -qopenmp -axCORE_AVX2 -parallel-source-info=2
else
 CXX = g++
CXXFLAGS += -fopenmp -msse4.2
#CXX = g++
endif

INCLUDES = -I. -I./sdsl-lite/build/include -I./htslib/htslib -I./KMC -I/opt/intel/advisor_2019.3.0.591490/include
#INCLUDES = -I/home/fabio/include -I./htslib/htslib -I./KMC
LDFLAGS = -L./sdsl-lite/build/lib -L./sdsl-lite/build/external/libdivsufsort/lib -L./htslib/
#LDFLAGS = -L/home/fabio/lib -L./htslib/
LDLIBS = -lhts -lz -lsdsl -ldivsufsort -ldivsufsort64 -ltbb -ldl
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
