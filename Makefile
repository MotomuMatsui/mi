CXX   := g++
TAR   := tar xzf
MKDIR := mkdir -p
CD    := cd
CP    := cp
CMAKE := cmake

VPATH := src
INC   := -Ilib
LIB   := -Llib -llapacke -llapack -lcblas -lrefblas -lgfortran -lm

CXXFLAGS := -O3
CXXFLAGS += -std=c++11
CXXFLAGS += -march=native
CXXFLAGS += -fno-exceptions
CXXFLAGS += -Wall

OBJECTS := main.o
OBJECTS += distance.o
OBJECTS += format.o
OBJECTS += mmseqs.o
OBJECTS += mafft.o
OBJECTS += nj.o
OBJECTS += transitivity.o

FILE  := lapack-3.7.1/make.inc
EXIST := $(shell ls | grep ${FILE})

.PHONY: all
all: mafft mmseqs lapack mi clean

.PHONY: mafft
mafft:
	$(TAR) mafft-7.407-without-extensions-src.tgz
	$(MAKE) -C mafft-7.407-without-extensions/core clean
	$(MAKE) -C mafft-7.407-without-extensions/core
	$(MAKE) -C mafft-7.407-without-extensions/core install

.PHONY: mmseqs
mmseqs:
	$(TAR) MMseqs2.tar.gz
	$(MKDIR) MMseqs2/build
	$(CD) MMseqs2/build; $(CMAKE) -DHAVE_SSE4_1=1 -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
	$(MAKE) -C MMseqs2/build
	$(MAKE) -C MMseqs2/build install 

.PHONY: lapack
lapack:
	$(TAR) lapack-3.7.1.tar.gz 
	$(MKDIR) lib
    ifneq (${EXIST}, ${FILE})
	$(CP) lapack-3.7.1/make.inc.example lapack-3.7.1/make.inc
    endif
	$(MAKE) -C lapack-3.7.1 blaslib
	$(MAKE) -C lapack-3.7.1 cblaslib
	$(MAKE) -C lapack-3.7.1 lapacklib
	$(MAKE) -C lapack-3.7.1 lapackelib
	$(CP) -f lapack-3.7.1/*.a lib
	$(CP) -f lapack-3.7.1/CBLAS/include/*.h lib
	$(CP) -f lapack-3.7.1/LAPACKE/include/*.h lib

mi: $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIB)

$(OBJECTS): %.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INC) -c $<

.PHONY: clean
clean:
	$(RM) $(OBJECTS)
