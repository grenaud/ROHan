
#"make all" runs the two demonstration data sets and needs to download them
#first.  "make test" is the self contained one: it simulates its own data and
#needs no network access, see the bottom of this file.
all: HG02891.summary.txt Loschbour.hg19_1000g.summary.txt

HG02891.21.bam:
	samtools view -b ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG02891/alignment/HG02891.mapped.ILLUMINA.bwa.GWD.low_coverage.20121211.bam 21 > HG02891.21.bam 

HG02891.21.bam.bai: HG02891.21.bam 
	samtools index HG02891.21.bam 

human_g1k_v37.masked.fa.gz:
	wget --no-check-certificate https://www.dropbox.com/s/4ayoefse8a5ri1n/human_g1k_v37.masked.fa.gz?dl=0 -O human_g1k_v37.masked.fa.gz

human_g1k_v37.masked.fa: human_g1k_v37.masked.fa.gz
	gunzip -f human_g1k_v37.masked.fa.gz

human_g1k_v37.masked.fa.fai: human_g1k_v37.masked.fa
	samtools faidx human_g1k_v37.masked.fa


HG02891.summary.txt: human_g1k_v37.masked.fa.fai HG02891.21.bam.bai
	../bin/rohan   --auto autosomes  --err  ../illuminaProf/error.prof  -t  4  -o  HG02891  human_g1k_v37.masked.fa HG02891.21.bam	


Loschbour.hg19_1000g.5p.prof: Loschbour.hg19_1000g.chr21.bam.bai
	samtools view -F1 -h -q 30 Loschbour.hg19_1000g.chr21.bam| awk '{if( substr($$1,1,1)=="@"  || length($$10)>45){ print $$0; }}'  | head -n 1000000 | samtools view -bS - | ../bam2prof/bam2prof -5p Loschbour.hg19_1000g.5p.prof -3p Loschbour.hg19_1000g.3p.prof -minq 30  -length 40  - 


Loschbour.hg19_1000g.chr21.bam.bai:
	wget --no-check-certificate -O Loschbour.hg19_1000g.chr21.bam     https://www.dropbox.com/s/1n24l8jijh7a08m/Loschbour.hg19_1000g.chr21.bam?dl=0
	wget --no-check-certificate -O Loschbour.hg19_1000g.chr21.bam.bai https://www.dropbox.com/s/nhdvjcr1no71re3/Loschbour.hg19_1000g.chr21.bam.bai?dl=0


Loschbour.hg19_1000g.summary.txt: Loschbour.hg19_1000g.chr21.bam.bai Loschbour.hg19_1000g.5p.prof Loschbour.hg19_1000g.chr21.bam.bai  human_g1k_v37.masked.fa.fai
	../bin/rohan --rohmu 5e-5 --deam5p Loschbour.hg19_1000g.5p.prof --deam3p Loschbour.hg19_1000g.3p.prof --out Loschbour.hg19_1000g --auto autosomes  --err ../illuminaProf/error.prof human_g1k_v37.masked.fa Loschbour.hg19_1000g.chr21.bam



############################################################################
#  make test: a self contained end-to-end test, run it from the top level  #
#  directory with "make test".  Takes a few minutes and needs no network   #
#  access: simulateTestData.py simulates a diploid chromosome with a known #
#  heterozygosity rate and one known run of homozygosity, rohan is run on  #
#  it and checkTestResults.py compares the results against what was        #
#  simulated.                                                              #
############################################################################

TESTPREFIX   = simulated
TESTOUT      = simulated.rohan
TESTCHR      = testchr
TESTLENGTH   = 3000000
#the ROH is a whole number of windows so that rohan can recover it exactly
TESTROHSTART = 1000000
TESTROHEND   = 2000000
TESTTHETA    = 1e-3
TESTCOV      = 8
TESTWINDOW   = 500000
#fewer chains than the 50,000 default: enough for a 6 window genome and it
#keeps the test down to a few minutes
TESTCHAINS   = 10000
TESTTHREADS  = 4
#seed for the simulation, and the one handed to rohan --seed.  Note that the
#local het. estimates still vary slightly between runs, so checkTestResults.py
#compares against the truth with a tolerance rather than byte for byte.
TESTSEED     = 42

SAMTOOLS     = ../lib/samtools/samtools
ROHAN        = ../bin/rohan

#one rule for the whole data set: the reference, the truth file and the BAM are
#produced together and must stay in sync, so they are all made by the recipe of
#the single target the rest of the test depends on
$(TESTPREFIX).bam: simulateTestData.py
	@if [ ! -x $(SAMTOOLS) ]; then echo "$(SAMTOOLS) is missing, run make in the top level directory first"; exit 1; fi
	python3 simulateTestData.py --prefix $(TESTPREFIX) --chr $(TESTCHR) --length $(TESTLENGTH) --roh-start $(TESTROHSTART) --roh-end $(TESTROHEND) --theta $(TESTTHETA) --coverage $(TESTCOV) --seed $(TESTSEED)
	$(SAMTOOLS) view -b $(TESTPREFIX).sam > $@
	$(SAMTOOLS) index $@
	$(SAMTOOLS) faidx $(TESTPREFIX).fa
	echo $(TESTCHR) > $(TESTPREFIX).autosomes
	rm -f $(TESTPREFIX).sam

test: $(TESTPREFIX).bam
	@if [ ! -x $(ROHAN) ]; then echo "$(ROHAN) is missing, run make in the top level directory first"; exit 1; fi
	$(ROHAN) --size $(TESTWINDOW) --chains $(TESTCHAINS) --seed $(TESTSEED) --auto $(TESTPREFIX).autosomes --err ../illuminaProf/error.prof -t $(TESTTHREADS) -o $(TESTOUT) $(TESTPREFIX).fa $(TESTPREFIX).bam
	@echo ""
	python3 checkTestResults.py --prefix $(TESTOUT) --truth $(TESTPREFIX).truth.txt --roh-start $(TESTROHSTART) --roh-end $(TESTROHEND) --theta $(TESTTHETA)

test-clean:
	rm -f $(TESTPREFIX).fa $(TESTPREFIX).fa.fai $(TESTPREFIX).sam $(TESTPREFIX).bam $(TESTPREFIX).bam.bai $(TESTPREFIX).truth.txt $(TESTPREFIX).autosomes
	rm -f $(TESTOUT).*

.PHONY: all test test-clean
