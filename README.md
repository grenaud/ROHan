# ROHan: inference of heterozygosity rates and runs of homozygosity for modern and ancient samples
==========================================================

QUESTIONS :
   gabriel [dot] reno [ at sign ] gmail [dot] com


About
----------------------

ROHan is a Bayesian framework to estimate local rates of heterozygosity, infer runs of homozygosity (ROH) and compute global rates of heterozygosity. ROHan can work on modern and ancient samples with signs of ancient DNA damage.

Description
----------------------

Diploid organisms have minor differences between corresponding pairs of autosomes. Differences of a single nucleotide create heterozygous sites. Very recent inbreeding can lead to large stretches of genomic deserts of heterozygosous sites. Such deserts are called runs of homozygosity (ROH). 

Due to the lack of heterozygosous sites, ROHs can cause an underestimate of global estimates of heterozygosity. Such global estimates of rates of heterozygous sites is useful to infer population genetics parameters such as effective population size (Ne). 

Analysis of ancient samples obtained using ancient DNA (aDNA), in addition to being plagued by potential ROHs due to inbreeding, often have presence of post-mortem molecular damage and can often have low coverage due to a paucity of endogenous material.

ROHan aims at:
1) providing an estimate of local rates of heterozygosity that is robust to low coverage and post-mortem damage. 
2) identifying ROHs using the results from 1). This is done via an hidden markov model whose parameters are optimized using a Markov Chain Monte Carlo
3) using both 1) and 2), provide a global rate of heterozygosity

Please refer to the LICENSE file for license information.

Downloading:
----------------------

Go to https://github.com/grenaud/rohan and either:

1) Download the ZIP 

or

2) Do a "git clone --depth 1 https://github.com/grenaud/rohan.git"

Installation
----------------------

1) make sure you have "aclocal", "cmake", "libtool", "libpng" and "git" installed, check for it by typing " git --version" and "cmake --version".  

For Ubuntu:

     sudo apt-get install -y autotools-dev git  cmake libtool libpng-dev liblzma-dev  libbz2-dev  libcurl4-openssl-dev  libgsl-dev

For MacOS, if you have Homebrew (https://brew.sh/) installed: 

    brew install autoconf
    brew install automake
    brew install git
    brew install cmake
    brew install libtool
    brew install libpng

2) make sure you have gcc that supports -std=c++11, gcc version 4.7 or above. Type "gcc -v" and check the version. For both Ubuntu and MacOS, 

3) As the makefile uses "git clone" to download subpackages, please make sure that the computer on which you are installing ROHan has access to the internet. Once done, simply type:
   
```bash
cd ROHan
make
```

For MacOS, if you get the problem: fatal error: 'lzma.h' file not found, this is a problem building htslib with homebrew, please refer to the following htslib page: https://github.com/samtools/htslib/issues/493


4) (optional) Either put the executable in the overall path or add the path to your $PATH environment or add an alias to be able to run "rohan" from any directory.


Testing your installation
----------------------

To check that your build works, type:

```bash
make test
```

This takes a few minutes and needs no internet access. It simulates a diploid chromosome of 3 Mbp at 8X coverage with a heterozygosity rate of 1e-3 and a single 1 Mbp run of homozygosity in the middle, runs ROHan on it and compares the results against what was simulated. The test passes if ROHan recovers the local and genome-wide heterozygosity rates and the position of the ROH. The files it produces live in testData/ and are removed by "make clean".

The parameters of the simulation are at the bottom of testData/Makefile if you want to try harder cases (lower coverage, a smaller ROH, a lower heterozygosity rate). To re-generate the data after changing them, run "make -C testData/ test-clean" first.


Quick start
----------------------

For test data, make sure you are connected to the internet and type:
     cd testData/
     make

This will run ROHan on chromosome 21 for: 
* a West African individual from Phase 3 of the 1000G project
* the ancient Loschbour individual (Lazardis et al.). Please be aware that the command for the profile uses no quality score filter: "-minq 0" because the quality scores were artificially decreased to account for damage during calling. We generally do not recommend this as this will inflate substitution rate. 

For these cases, the theta estimate have very large confidence due to the use of a single chromosome. 


Preparing the BAM file
-----------------

1) Make sure you know the transition/transversion ratio for the species/population you are working with. This ratio will be specified via:

    rohan --tstv [TSTV ratio] 

2) Do not apply any filters as mapping quality and base quality are all informative in the model. Duplicate removal is very recommended. Simply use the fail QC flag to remove reads/fragments that have failed basic quality control (e.g. for duplicates). Simply sort and index and provide ROHan the same reference used for mapping. 


3) Is your sample ancient DNA? if not skip to 5). If the sample has aDNA damage, it should be quantified in incorporated in the calculation.

**Important: if you are going to correct for damage, your reads must be single-end.** The deamination rates are indexed on the distance to the 5' and 3' ends of the original molecule. For a read that is one half of a pair, the end of the read is not the end of the molecule, so the rates cannot be placed correctly. ROHan therefore **skips every paired-end read** as soon as "--deam5p"/"--deam3p" are used. If all of your reads are paired-end, nothing at all will be retained and ROHan will tell you so and stop; if only some of them are, the paired-end portion is silently dropped and ROHan will print a warning.

The usual way to deal with this, and the standard practice for ancient DNA, is to merge (collapse) overlapping read pairs into single reads before mapping, using e.g. [leeHom](https://github.com/grenaud/leeHom), "AdapterRemoval --collapse", "fastp -m" or SeqPrep, and then align the collapsed reads as single-end. Ancient and museum samples usually have short fragments, so the large majority of pairs will overlap and collapse. If you do not need the damage correction, you can run on paired-end data without "--deam5p"/"--deam3p".

There are 2 ways to quantify the damage:

i) If the sample is not too divergent from the reference used you can simply use:

    bam2prof/bam2prof

 To get the best results, take a substantial subsample of your original BAM file and use it in bam2prof. 
 a) Find the ideal length for the -length parameter. Try increasing until substition rates level off.
 b) Keep increasing  -minq from 0 until the damage levels off.
 c) If you have a high rate of substitutions outside of expected ones (e.g. C->T, G->A) consider using the second option (see ii) )
 d) Have a look at the command line in the test data directory.

ii) Use the script to estimate damage by masking polymorphic positions, this is especially important if the divergence between the reference genome and sample is great and there are a lot of polymorphic positions:

    src/estimateDamage.pl    [reference fasta] [input .bam file]

 For example:

    src/estimateDamage.pl    Equcab20_nucl.fasta ancientHorse.bam

 type the command without any arguments for a description of the different arguments.
 
4) Do you have extensive ancient DNA damage (e.g. 20% at the ends of greater) which could potentially affect the mapping? If so, we recommend filtering for reads in highly mappable regions, see : http://lh3lh3.users.sourceforge.net/snpable.shtml  to create mappability tracks.

5) Create a file with the name of the autosomes, 1 chromosome per line ex:
    chr1
    chr2
    ...

6) Run ROHan, for modern samples run:

     rohan --rohmu 2e-5   -o  [output prefix]  [reference genome]  [BAM file]

For ancient samples run:

     rohan --rohmu 2e-5 --deam5p  deamfile.5p.prof  --deam3p  deamfile.3p.prof    -o  [output prefix]  [reference genome]  [BAM file]

where deamfile.5p.prof and deamfile.3p.prof are the substitution rates due to aDNA damage and not due to sequencing errors. You can add more threads in the calculation using "-t". Be aware that as this is a full likelihood model, the data needs to live in memory and therefore, using multiple threads can result in higher memory usage.

If a mappability track was used, you can add the option: --map to filter wrt to mappable regions.

7) Inspect your results. Refer to [output prefix].het_1_X.pdf and [output prefix].hmm_1_X.pdf, if there are regions labeled as non ROH despite the fact that they clearly show a depression in heterozygosity, re-run using --hmm and a higher value of --rohmu (e.g. for the example above, one could run using 5e-5).

8) For ancient samples, if the result if really unexpected, you can re-run using --tvonly which will limit the estimate to transversions, the heterozygosity needs to be multiplied by (Ts/Tv+1) (ex: if Ts/Tv  = 2.1, multiply the results it by 3.1). This estimate can be an underestimate due to the lack of data.


Description of output files
-----------------


|                  File                |                           Contents                               | 
| ------------------------------------ |------------------------------------------------------------------| 
[output prefix].hEst.gz                | local estimates of heterozygosity                                | 
[output prefix].het_1_X.pdf            | Plot of the local estimates of heterozygosity 1 of X             | 
[output prefix].hmm_1_X.pdf            | Plot of the posterior decoding for HMM 1 of X                    | 
[output prefix].max.hmmmcmc.gz         | Log of the MCMC using maximum estimates of heterozygosity        | 
[output prefix].max.hmmp.gz            | HMM posterior decoding using maximum estimates of heterozygosity | 
[output prefix].mid.hmmmcmc.gz         | Log of the MCMC using mid. estimates of heterozygosity           | 
[output prefix].mid.hmmp.gz            | HMM posterior decoding using mid. estimates of heterozygosity    | 
[output prefix].min.hmmmcmc.gz         | Log of the MCMC using minimum estimates of heterozygosity        | 
[output prefix].min.hmmp.gz            | HMM posterior decoding using minimum estimates of heterozygosity | 
[output prefix].rginfo.gz              | List of read groups and average coverage                         | 
[output prefix].summary.txt            | Fraction of the genome in ROH and heterozygosity outside of ROH  | 
[output prefix].vcf.gz                 | Posterior probabilities for the genotypes                        | 


In [output prefix].summary.txt, there are 2 lines:

     Genome-wide theta outside ROH:	thetamid thetamin thetamax
     Genome-wide theta inc. ROH:	thetamid thetamin thetamax

The first line is for Watterson's theta for genomic regions excluding runs of homozygosity. The second is for Watterson's theta for the entire genome, including runs of homozygosity. 

FAQ
----------------------



### I get "No data was found for the entire BAM file for the regions you defined", why?

- ROHan estimates the average coverage on a subsample of genomic windows and found no read at all in any of them. Either the BAM genuinely has no data there, which you can check with "samtools view [your bam] [one of the regions listed]", or every read was filtered out. When the BAM does contain reads in those regions, by far the most common cause is using "--deam5p"/"--deam3p" on paired-end data: paired-end reads are all skipped when deamination profiles are given, see point 3) of "Preparing the BAM file" above. ROHan detects this case and reports it explicitly.

### The build fails with "fatal error: curses.h: No such file or directory", why?

- This came from the bundled samtools, which builds its interactive alignment viewer ("samtools tview") unless it is told not to. ROHan does not use the viewer, or any other part of the samtools binary; it only links four object files. As of v1.0.5 the makefile asks samtools for exactly those, so ncurses is no longer needed at all and this error should not occur any more. If you are on an older version, either install the ncurses development headers ("sudo apt-get install libncurses-dev") or pull a current ROHan.

### ROHan builds but will not start: "libgsl.so.25: cannot open shared object file", why?

- ROHan uses the GNU Scientific Library, so libgsl has to be present at run time and not only at build time. On a cluster this usually means the library exists on the login node where you compiled but not on the compute node where the job runs; loading the same module in your job script ("module load gsl/...") fixes it. If you would rather not depend on the environment at all, "make -C src/ static" links everything into the binary, and the "bin/rohan" that ships in the repository is already a static binary.

### The coverage step prints "Results bp=0 sites=0 lambda=...", is something wrong?

- No. Prior to v1.0.5 that line was printed with zeroes whenever ROHan reused a cached "[out prefix].rginfo.gz" from an earlier run rather than scanning the BAM again: the lambda is real, but the bp and sites counters belong to the run that wrote the cache and were never filled in. Current versions say "reused from ..." instead. Note that the cache is reused whenever the file exists, so if you have changed the BAM, the reference or the read filters since, pass "-f" to force the coverage and read groups to be recomputed.

### The heterozygosity step sits at 0% and the .hEst.gz file does not grow, has it hung?

- Almost certainly not, it is just slow, and neither of the two things you are watching can show you otherwise. ROHan estimates heterozygosity in every window in turn, and on a 3 Gbp genome the default 1 Mbp windows come to a few thousand of them. Measured on such an assembly (2.9 Gbp, 4,235 windows) with four threads on four cores, the progress bar took 42 minutes to move from 0% to 1%; on a single thread it is a few hours before the first percent appears. The file is worse: ".hEst.gz" is written through bgzip, which does not put anything on disk until 64 kbytes of text have accumulated, and a window is only about 70 bytes, so the file stays at exactly zero bytes for the first ~936 windows however well the run is going. From v1.0.5 onwards ROHan prints the number of windows and counts them off individually, which is the number to watch. This step is the bulk of the run time and it parallelizes well, so use "-t" with as many cores as you have, and note that memory scales with it (1.3 GB resident at "-t 4" on the genome above), so ask your scheduler for enough. "--auto" with a list of your assembled chromosomes also avoids spending time on unplaced scaffolds.

### Why do I have such huge confidence intervals for potential ROH regions?

This is normal, as we use the second derivative for the confidence intervals, the shape of the likelihood function is somewhat flat around 0 which leads to overestimated confidence intervals. 


### I have regions that have clear signs of being runs of homozygosity however they do not get flagged as such, why?

It is possible that you have genuine ROH but if you have a slight overestimate of the rate of local heterozygosity do to him properly calibrated weapon qualities or  base qualities. Therefore such regions will not be flagged as being ROHs. It is advisable to rerun with --rohmu XXX where XXX is the "heterozygosity rate" in ROH regions you are willing to tolerate

### What are the 3 different lines on the HMM plot?

- The green line is the output using the point estimates. Margenta are the estimates using upper bounds for local het. rates whereas red were computed using lower bounds of het. rates.


### What is the difference between --bed and --map?

- By default, ROHan generates windows of "--size"bp along autosomes, you can override this and use the regions in a bed file. With either option, you can specify individual sites to consider using --map, this is recommended for ancient samples with short fragment size and with some damage.


### I get slightly different results every time I run on the same data, why?

- The HMM parameters are optimized using a Markov Chain Monte Carlo, which starts from a random point. The seed is taken from the wall clock and reported both on stderr and in [output prefix].summary.txt, re-running with "--seed [that seed]" reproduces the MCMC. Note that the local heterozygosity estimates still vary in the last few digits from one run to the next, so runs are not bit for bit identical even with the same seed.


### The percentage of classified segments does not sum up to one (100%), why?

- The % of segments should not necessarily sum up to one. Segments can be:

* in ROH
* not in ROH
* not sure

Therefore, because of the segments in the "not sure" category, this percentages will not necessarily sum up to 100%.

### How can I cite ROHan?

Please cite the following:
Gabriel Renaud, Kristian Hanghøj, Thorfinn Sand Korneliussen, Eske Willerslev, Ludovic Orlando, Joint Estimates of Heterozygosity and Runs of Homozygosity for Modern and Ancient Samples, Genetics, Volume 212, Issue 3, 1 July 2019, Pages 587–614, https://doi.org/10.1534/genetics.119.302057

or in bibtex:

         @article{10.1534/genetics.119.302057,
        author = {Renaud, Gabriel and Hanghøj, Kristian and Korneliussen, Thorfinn Sand and Willerslev, Eske and Orlando, Ludovic},
        title = "{{Joint Estimates of Heterozygosity and Runs of Homozygosity for Modern and Ancient Samples}}",
        journal = {Genetics},
        volume = {212},
        number = {3},
        pages = {587-614},
        year = {2019},
        month = {05},
        abstract = "{Both the total amount and the distribution of heterozygous sites within individual genomes are informative about the genetic diversity of the population they belong to. Detecting true heterozygous sites in ancient genomes is complicated by the generally limited coverage achieved and the presence of post-mortem damage inflating sequencing errors. Additionally, large runs of homozygosity found in the genomes of particularly inbred individuals and of domestic animals can skew estimates of genome-wide heterozygosity rates. Current computational tools aimed at estimating runs of homozygosity and genome-wide heterozygosity levels are generally sensitive to such limitations. Here, we introduce ROHan, a probabilistic method which substantially improves the estimate of heterozygosity rates both genome-wide and for genomic local windows. It combines a local Bayesian model and a Hidden Markov Model at the genome-wide level and can work both on modern and ancient samples. We show that our algorithm outperforms currently available methods for predicting heterozygosity rates for ancient samples. Specifically, ROHan can delineate large runs of homozygosity (at megabase scales) and produce a reliable confidence interval for the genome-wide rate of heterozygosity outside of such regions from modern genomes with a depth of coverage as low as 5–6× and down to 7–8× for ancient samples showing moderate DNA damage. We apply ROHan to a series of modern and ancient genomes previously published and revise available estimates of heterozygosity for humans, chimpanzees and horses.}",
        issn = {1943-2631},
        doi = {10.1534/genetics.119.302057},
        url = {https://doi.org/10.1534/genetics.119.302057},
        eprint = {https://academic.oup.com/genetics/article-pdf/212/3/587/42104856/genetics0587.pdf},
        }


`

