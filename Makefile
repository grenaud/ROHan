all:	bin/rohan bam2prof/bam2prof
	@echo "Done"


bam2prof/bam2prof.cpp:
	rm -rfv bam2prof/
	mkdir -p bam2prof/
	rm -f bam2prof.zip
	wget --no-check-certificate -O bam2prof.zip https://github.com/grenaud/bam2prof/zipball/master
	unzip bam2prof.zip -d bam2prof
	mv -v bam2prof/*/src/*  bam2prof/


bam2prof/bam2prof: bam2prof/bam2prof.cpp
	make -C bam2prof/
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

