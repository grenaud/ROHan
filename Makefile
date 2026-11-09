all:	bin/rohan bam2prof/bam2prof
	@echo "Done"


#git clone instead of wget on the zipball: the zipball URL goes through GitHub's archive
#service, which returned 502/504 for minutes at a time and wget's default retries made the
#build hang ~5 minutes before failing. Retry a few times, then fail clearly.
bam2prof/bam2prof.cpp:
	rm -rfv bam2prof/ bam2prof.tmp/
	for i in 1 2 3; do git clone --depth 1 -b master https://github.com/grenaud/bam2prof.git bam2prof.tmp && break; rm -rf bam2prof.tmp; sleep 15; done; test -d bam2prof.tmp
	mkdir -p bam2prof/
	mv -v bam2prof.tmp/src/*  bam2prof/
	rm -rf bam2prof.tmp/


#bam2prof's own Makefile hardcodes its link libraries, with no -ldeflate, but it links against
#the htslib that src/Makefile configured. Where libdeflate-dev is installed that htslib needs
#-ldeflate, so ask htslib_static.mk for the real list (same as src/Makefile) and override
#LDFLAGS. bin/rohan is a prerequisite so htslib is built and configured before this is read.
HTSLIB_STATIC_LIBS = $(shell sed -n 's/^HTSLIB_static_LIBS[ 	]*=[ 	]*//p' lib/htslib/htslib_static.mk 2>/dev/null)
HTSLIB_STATIC_LDFLAGS = $(shell sed -n 's/^HTSLIB_static_LDFLAGS[ 	]*=[ 	]*//p' lib/htslib/htslib_static.mk 2>/dev/null)

bam2prof/bam2prof: bam2prof/bam2prof.cpp bin/rohan
	make -C bam2prof/ LDFLAGS="-lgsl -lgslcblas $(if $(HTSLIB_STATIC_LIBS),$(HTSLIB_STATIC_LIBS),-llzma -lbz2 -lz -lpthread) $(HTSLIB_STATIC_LDFLAGS)"
	rm -rfv bam2prof/grenaud*/
	rm -f bam2prof.zip

#bin/rohan is declared phony on purpose: src/Makefile is the one that knows what rohan is
#built from, so always recurse into it and let it decide what needs recompiling. Testing
#whether bin/rohan merely exists would report "Done" without rebuilding after a source edit.
bin/rohan:
	make -C  src/

#end-to-end test on simulated data, needs no network access, see testData/Makefile
test:	bin/rohan
	make -C testData/ test

clean:
	make -C src/ clean
	make -C bam2prof/ clean
	make -C testData/ test-clean

.PHONY: all test clean bin/rohan

