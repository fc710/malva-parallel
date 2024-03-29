COMPILER1 = #$(shell which icpc | grep -o "icpc")
COMPILER2 = #$(shell which clang++ | grep -o "clang++")
CXXFLAGS = -Wall -g -std=c++17 -O3
ifeq ($(COMPILER1), icpc)
	CXX = icpc
	CXXFLAGS += -qopenmp -axcore-avx2 -parallel-source-info=2
else ifeq ($(COMPILER2), clang++)
	CXX = clang++
	CXXFLAGS += -fopenmp -axcore-avx2
else
	CXX = g++
	CXXFLAGS += -fopenmp -msse4.2

endif

INCLUDES = -I. -I./sdsl-lite/build/include -I./htslib/htslib -I./KMC
LDFLAGS = -L./sdsl-lite/build/lib -L./sdsl-lite/build/external/libdivsufsort/lib -L./htslib/
LDLIBS = -lhts -lz -lsdsl -ldivsufsort -ldivsufsort64 
OBJS = main.o MurmurHash3.o ./KMC/kmc_api/kmc_file.o ./KMC/kmc_api/kmer_api.o ./KMC/kmc_api/mmer.o
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
