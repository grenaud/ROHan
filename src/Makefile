
CXX      = g++ -std=c++11 #-g #-pg #-fpermissive #-g # -pg
LIBGAB   = ../lib/libgab/

CXXFLAGS = -Wall -lm -O3 -lz -I../lib/libgab/ -I../lib/libgab//gzstream/ -I../lib/samtools/ -I../lib/htslib/  -I../lib/libharu/include/   -c
LDFLAGS  = -lgsl -lcurl -lgslcblas -llzma -lbz2 -lz  -lpthread 
LDLIBS   = ../lib/samtools/bedidx.o ../lib/htslib/libhts.a ../lib/samtools/libst.a ../lib/samtools/libst.a 


all:  ../lib/libgab/libgab.a ../lib/libharu/src/.libs/libhpdf.a ../lib/htslib/libhts.a ../lib/samtools/bedidx.o ../bin/rohan refPlusSegToGeno sitesPlusGeno ref2DNAprof #testComp testCompCont

%.o: %.cpp #../lib/libgab/libgab.a  
	${CXX} ${CXXFLAGS} $^ -o $@

#testComp:	testComp.o 
#	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)
#testCompCont:	testCompCont.o ../lib/libgab/libgab.a  
#	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)


refPlusSegToGeno:	refPlusSegToGeno.o ../lib/libgab/libgab.a   ../lib/libgab/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

sitesPlusGeno:	sitesPlusGeno.o ../lib/libgab/libgab.a   ../lib/libgab/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

../bin/rohan:       rohan.o PdfWriter.o Hmm.o HmmState.o DataChunk.o DataToWrite.o GenoResults.o PositionResult.o miscfunc.o  GenomicWindows.o GenomicRange.o  ../lib/libgab/FastQObj.o ../lib/libgab/FastQParser.o ../lib/libharu/src/.libs/libhpdf.a ../lib/libgab/libgab.a  ../lib/libgab/gzstream/libgzstream.a 
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

static:       rohan.o PdfWriter.o Hmm.o HmmState.o DataChunk.o DataToWrite.o GenoResults.o PositionResult.o miscfunc.o  GenomicWindows.o GenomicRange.o  ../lib/libgab/FastQObj.o ../lib/libgab/FastQParser.o ../lib/libharu/src/.libs/libhpdf.a ../lib/libgab/libgab.a  ../lib/libgab/gzstream/libgzstream.a 
	${CXX} -static -o ../bin/rohan $^ $(LDLIBS) $(LDFLAGS)

ref2DNAprof.o:	ref2DNAprof.cpp
	${CXX} ${CXXFLAGS} ref2DNAprof.cpp

ref2DNAprof:	ref2DNAprof.o ${LIBGAB}libgab.a  ${LIBGAB}gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

../lib/htslib/libhts.a:  ../lib/htslib/hts_internal.h
	cd ../lib/htslib/ && git submodule update --init --recursive && make && cd ../..
#	cd ../lib/htslib/ && git submodule update --init --recursive && autoreconf -i && ./configure --disable-libcurl  && make && cd ../..

../lib/htslib/hts_internal.h:
	rm -rfv ../lib/htslib/
	git clone --recursive https://github.com/samtools/htslib.git ../lib/htslib/

../lib/samtools/bedidx.o:  ../lib/samtools/samtools.h
	cd ../lib/samtools/  &&   make && cd ../..
#	cd ../lib/samtools/  && autoreconf -i && ./configure --disable-libcurl &&  make && cd ../..

../lib/samtools/samtools.h:
	rm -rf ../lib/samtools/
	git clone --recursive https://github.com/samtools/samtools.git ../lib/samtools/

../lib/libgab/libgab.a:  ../lib/libgab/libgab.h
	cd ../lib/libgab/ &&  make libgab.a && make -C gzstream/ && cd ../..

../lib/libgab/libgab.h:
	rm -rfv ../lib/libgab/
	mkdir ../lib/libgab/
	git clone --depth 1 https://github.com/grenaud/libgab.git ../lib/libgab/

../lib/libharu/include/hpdf.h:
	rm -rf ../lib/libharu/
	mkdir ../lib/libharu/
	git clone  https://github.com/libharu/libharu.git ../lib/libharu/
	cd ../lib/libharu/ && 	git checkout 4ae1d5f4c84459f130bf1b1ef4c5c330af8eca5d && cd ../..

../lib/libharu/src/.libs/libhpdf.a: ../lib/libharu/include/hpdf.h
	cd ../lib/libharu/ && ./buildconf.sh --force && ./configure  --without-zlib --without-png && make

clean :
	rm -f *.o testComp  testCompCont  refPlusSegToGeno sitesPlusGeno  testPrior  ../bin/rohan  ref2DNAprof 


