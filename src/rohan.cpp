#include <iostream>
#include <fstream>
#include <queue>
#include <algorithm>   
#include <cfloat>   
#include <random>

//TODO
// maybe skip mpq that do not appear to save RAM?
// fix rginfo
// TODO: GC bias for coverage?

//#include "api/internal/io/BgzfStream_p.h"
//#include <api/BamConstants.h>
//#include <api/BamMultiReader.h>
//#include <utils/bamtools_fasta.h>
//#include <utils/bamtools_options.h>
//#include <utils/bamtools_pileup_engine.h>
//#include <utils/bamtools_utilities.h>


extern "C" {
    //#include "tabix.h"
    //#include "bam.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/tbx.h"
#include "htslib/bgzf.h"
#include "bam.h"

#include "samtools.h"
#include "sam_opts.h"
#include "bedidx.h"
}

#define bam_is_properpair(b)    (((b)->core.flag&BAM_FPROPER_PAIR) != 0)
#define bam_is_paired(b)    (((b)->core.flag&BAM_FPAIRED) != 0)
#define bam_is_pair1(b)     (((b)->core.flag&BAM_FREAD1)  != 0)
#define bam_is_pair2(b)     (((b)->core.flag&BAM_FREAD2)  != 0)
#define bam_is_qcfailed(b)  (((b)->core.flag&BAM_FQCFAIL)     != 0)
#define bam_is_rmdup(b)     (((b)->core.flag&BAM_FDUP)        != 0)
#define bam_is_failed(b)    ( bam_is_qcfailed(b) || bam_is_rmdup(b) )
#define bam_mqual(b)        ((b)->core.qual)
#define bam1_qname(b)       (bam_get_qname((b)))

typedef struct {     // auxiliary data structure
    samFile *fp;     // the file handle
    bam_hdr_t *hdr;  // the file header
    hts_itr_t *iter; // NULL if a region not specified
    int min_mapQ, min_len; // mapQ filter; length filter
} aux_t;


#include "PdfWriter.h"
#include "GenomicWindows.h"
#include "Hmm.h"

#define MAXCOV             50     // maximal coverage
//#define TESTHMMSIMS


#include "miscfunc.h"
#include "libgab.h"

using namespace std;
// using namespace BamTools;


//#define MIN(a,b) (((a)<(b))?(a):(b))
//#define MAX(a,b) (((a)>(b))?(a):(b))

//#define DEBUGFIRSTWINDOWS 3
//#define CORRECTCOV
//#define ONLYUSECOV 12
// #define ONLYUSECOVMIN 5
// #define ONLYUSECOVMAX 49


//#define PRECOMPUTELOG
//#define DEBUGCOV

//#define DEBUGILLUMINAFREQ
//#define DEBUGINITSCORES

// #define DEBUGINITLIKELIHOODSCORES
// #define DEBUGINITLIKELIHOODSCORES2 55

//#define DEBUGDEFAULTFREQ//print the default base frequency 
//#define DEBUGDEAM //to print deamination scores
//#define DEBUGHCOMPUTE
//#define DEBUGCOMPUTELLGENO
//#define DEBUGCOMPUTELL
//#define DEBUGCOMPUTELLEACHBASE
//#define DEBUGPRECOMPUTEBABD 2610

// #define DEBUGQUEUEDATATOWRITE

// #define HETVERBOSE
// #define DEBUGBOXINGALG

//#define COVERAGETVERBOSE
//#define DUMPTRIALLELIC //hack to remove tri-allelic, we need to account for them

// #define MINLENGTHFRAGMENT     35      // mininam length for fragment
// #define MAXLENGTHFRAGMENT     250     // maximal length for fragment


//#define MAXMAPPINGQUAL        257     // maximal mapping quality, should be sufficient as mapping qualities are encoded using 8 bits
#define MAXMAPPINGQUAL          37     // maximal mapping quality
#define MAXBASEQUAL             64      // maximal base quality score, greater qual score do not make a lot of sense

//#define MAXMAPPINGBASEQUAL    64      // maximal base quality score, should be sufficient as mapping qualities are encoded using 8 bits




//0 N
//1 A
//2 C
//3 X
//4 G
//5 X
//6 X
//7 X
//8 T
//9-15 X
string alphabetHTSLIB = "NACNGNNNTNNNNNNN";

int alphabetHTSLIB2idx [16] = {-1, 0, 1,-1,  // N A C N
			       2,-1,-1,-1,   // G N N N
			       3,-1,-1,-1,   // T N N N
			       -1,-1,-1,-1}; // N N N N



unsigned int MINLENGTHFRAGMENT  =    35;      // minimal length for fragment
unsigned int MAXLENGTHFRAGMENT  =    MINLENGTHFRAGMENT;     //  maximal length for fragment

long double stepHMM = 1000;
char offsetQual=33;
//                                         4 u       Ne
int minSegSitesPer1M =   80; // 0.00008	 = 4*2e-8* 1000
int maxSegSitesPer1M = 5000; // 0.00500  = 4*2e-8*62500

double  rho1M=20; 

long double likeMatch           [MAXBASEQUAL+1];
long double likeMismatch        [MAXBASEQUAL+1];

long double likeMatchProb       [MAXBASEQUAL+1];
long double likeMismatchProb    [MAXBASEQUAL+1];

// long double likeMatchMap        [MAXMAPPINGQUAL+1];
// long double likeMismatchMap     [MAXMAPPINGQUAL+1];

long double likeMatchProbMap    [MAXMAPPINGQUAL+1];
long double likeMismatchProbMap [MAXMAPPINGQUAL+1];


long double overestimateFactor [MAXCOV+1];

long double TStoTVratio;

vector< vector<long double> > binomVec (MAXCOV+1,vector<long double>(MAXCOV+1,0)) ;
unsigned int totalBasesSum;
unsigned int totalSitesSum;

//Substitution rates due to deamination
vector<probSubstition> sub5p;
vector<probSubstition> sub3p;

vector<diNucleotideProb> sub5pDiNuc;
vector<diNucleotideProb> sub3pDiNuc;


vector< vector<probSubstition> > subDeam;        //first dimension is fragment length, second is position
vector <vector<diNucleotideProb> > subDeamDiNuc; //first dimension is fragment length, second is position


probSubstition   defaultSubMatch;
diNucleotideProb defaultSubMatchMatrix;

probSubstition   illuminaErrorsProb;
diNucleotideProb illuminaErrorsProbMatrix;

//default basepairs frequencies
alleleFrequency         dnaDefaultBases;
vector<alleleFrequency> defaultDNA5p;
vector<alleleFrequency> defaultDNA3p;
vector< vector<alleleFrequency> > defaultDNAFragmentLength;

string fastaFile;
string fastaIndex;

map<string,rgInfo> rg2info;
bool specifiedDeam=false;
bool verbose=false;
bool outputgenol = false;
bool tvonly=false;

bool usesuggestH=false;
double sumhsuggest=0.0;
int    nhsuggest  =0;

bool suggestH=false;
double suggestHval;

bool covcutoffs   =false;
int  covcutoffsmin=0;
int  covcutoffsmax=MAXCOV;

// 1D: mapping quality 
// 2D: length of fragment
// 3D: base qual
// 4D: pos fragment
// 5D: 4X4 matrix
//vector< vector< vector< vector<diNucleotideProb> > > > mpq2Length2BaseQual2Pos2SubMatrix;

// // 1D: mapping quality 
// // 2D: pos away from 5'/3' end
// // 3D: base qual
// vector< vector< vector<diNucleotideProb> > >  mpq2Pos2BaseQual2SubMatrix5p;
// vector< vector< vector<diNucleotideProb> > >  mpq2Pos2BaseQual2SubMatrix3p;



typedef vector< vector<diNucleotideProb> > mpq2bsq2submatrix;

 // 1D: pos away from 5'/3' end
 // 2D: mapping quality 
 // 3D: base qual
// vector< mpq2bsq2submatrix >  pos2mpq2BaseQual2SubMatrix5p;
// vector< mpq2bsq2submatrix >  pos2mpq2BaseQual2SubMatrix3p;

// 1D: length of fragment
// 2D: pos fragment from the 5' end
// 3D: mapping quality 
// 4D: base qual
// 5D: 4X4 matrix

vector< vector< mpq2bsq2submatrix * > >  length2pos2mpq2bsq2submatrix;
//vector< vector< mpq2bsq2submatrix  > >  length2pos2mpq2bsq2submatrix;


//vector< vector< mpq2bsq2submatrix * > > length2pos2mpq2bsq2submatrix;
// pointer to data structure

//vector< 



//vector< vector< vector<diNucleotideProb> > > mqp;//FOR EACH FRAGMENT LENGTH



// //dummy values
// for(int L=0;L<MINLENGTHFRAGMENT;L++){ //for each fragment length
// 	vector<diNucleotideProb> baseQualEffect;//FOR EACH QUAL SCORE
// 	fragmentLengthEffect.push_back(baseQualEffect);
// }

// //int L=35;//fragment length    
// for(int L=MINLENGTHFRAGMENT;L<MAXLENGTHFRAGMENT;L++){ //for each fragment length
// 	cerr<<L<<endl;

// 	vector<diNucleotideProb> baseQualEffect;//FOR EACH QUAL SCORE




vector<long double> * cov2ProbSite;

long double contrate=0.0;
long double rateForPoissonCov;
long double pdfRateForPoissonCov;

// // Returns logl( expl(x)+expl(y) )
// inline long double oplusl(long double x, long double y ){
//     return x > y 
//         ? x + log1pl( expl( y-x ) )
//         : y + log1pl( expl( x-y ) )  ;
// }

//for gradient descent
long double alpha=1.0e-9;
long double beta =0.5;

long double randomPMismatch4Bases =  ( (long double)(3) ) /  ( (long double)(4) );  // 3/4
long double randomPMatch4Bases    =  1.0-randomPMismatch4Bases;                     // 1/4

long double randomPMmToABase      =  ( (long double)(1) ) /  ( (long double)(3) );  // 1/3, prob of a mismatch to given base

string genoIdx2Code10 [10] = {"AA","AC","AG","AT","CC","CG","CT","GG","GT","TT"};
string babdIdx2Code16 [16] = {"AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"};
//                              0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
char  geno16_2_geno10 [16] = {  0 ,  1,   2,   3,   1,   4,   5,   6,   2,   5,   7 ,  8 ,  3 ,  6,   8,   9 };

 // genoPriority[0] = 0;//AA
 //    genoPriority[1] = 4; //CC
 //    genoPriority[2] = 7; //GG
 //    genoPriority[3] = 9; //TT

 //    //hetero
 //    genoPriority[4] = 1; //AC
 //    genoPriority[5] = 2; //AG
 //    genoPriority[6] = 3; //AT
 //    genoPriority[7] = 5; //CA
 //    genoPriority[8] = 6; //CG
 //    genoPriority[9] = 8; //CT

// 0	0	0
// 4	1	1
// 7	2	2
// 9	3	3
vector<int> genoPriority (10,0);


#include "DataChunk.h"
#include "DataToWrite.h"
#include "GenoResults.h"






















int    timeThreadSleep =    10;
int    timeSleepWrite  =    1;

bool      readDataDone = false;
unsigned int sizeChunk =  1000000;

string                                                             bamFileToOpen;
string                                                             bedfilename;//only for mappability 
queue< DataChunk * >                                               queueDataToprocess;
queue< DataChunk * >                                               queueDataAllGenomic;//in case we cheat and use the first windows of data

queue< DataChunk * >                                               queueDataForCoverage;

priority_queue<DataToWrite *, vector<DataToWrite *>, CompareDataToWrite> queueDataTowrite;


pthread_mutex_t  mutexQueueToRead           = PTHREAD_MUTEX_INITIALIZER;

// also used for coverage counter
pthread_mutex_t  mutexQueueToWrite          = PTHREAD_MUTEX_INITIALIZER; 
pthread_mutex_t  mutexRank                  = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t  mutexCERR                  = PTHREAD_MUTEX_INITIALIZER;


pthread_mutex_t  mutexCOUNTCOV              = PTHREAD_MUTEX_INITIALIZER;
unsigned int queueDataForCoverageOrigsizeSum=0;

// pthread_mutex_t  mutexCoverageCounter = PTHREAD_MUTEX_INITIALIZER;


//GLOBALLY accessed
map<unsigned int, int>       threadID2Rank;


//! Chunk of code to check if a certain thread call failed
/*!
  This block is calls by the pthread

*/				
#define checkResults(string, val) {             \
 if (val) {                                     \
     cerr<<"Failed with "<<val<<" at "<<string<<endl;	\
   exit(1);                                     \
 }                                              \
}
 

//#define debugInitscores


// typedef struct{
//     long double ll;
//     long double expal1; //expectation of # of allele 1
//     long double expal2; //expectation of # of allele 2

// } computeLLRes;


// //! A method to compute the pdf of poisson dist.
// /*!
//   This method computes the probability density function for a poisson dist.
// */
// long double pdfPoisson(const long double l,const long double k ) {
//     return expl(k*logl(l)-lgammal(k+1.0)-l);
// }



//! A method to initialize various probability scores to avoid recomputation
/*!
  This method is called by the main after capturing the arguments
*/
void initScores(){
    initMiscFuncVariables();

    totalBasesSum=0;
    totalSitesSum=0;


    //homoz
    genoPriority[0] = 0;
    genoPriority[1] = 4;
    genoPriority[2] = 7;
    genoPriority[3] = 9;

    //hetero
    genoPriority[4] = 1;
    genoPriority[5] = 2;
    genoPriority[6] = 3;
    genoPriority[7] = 5;
    genoPriority[8] = 6;
    genoPriority[9] = 8;


    //BASE QUALs, cannot go down to qual zero
    for(int i=0;i<2;i++){
        likeMatch[i]          = log1pl(    -randomPMatch4Bases );          
        likeMismatch[i]       = logl  (     randomPMismatch4Bases );

        likeMatchProb[i]              =    randomPMatch4Bases;    // 1/4
        likeMismatchProb[i]           =    randomPMismatch4Bases; //3/4
    }

    for(int i=2;i<=MAXBASEQUAL;i++){
        likeMatch[i]          = log1pl(    -powl(10.0,i/-10.0) );          
        likeMismatch[i]       = logl  (     powl(10.0,i/-10.0) );

        likeMatchProb[i]              = 1.0-powl(10.0,i/-10.0);
        likeMismatchProb[i]           =     powl(10.0,i/-10.0);
    }

    //mapping QUALs
    for(int i=0;i<=MAXMAPPINGQUAL;i++){
        // likeMatchMap[i]        = log1pl(    -powl(10.0,i/-10.0) );          
        // likeMismatchMap[i]     = logl  (     powl(10.0,i/-10.0) );

        likeMatchProbMap[i]           = 1.0-powl(10.0,i/-10.0);
        likeMismatchProbMap[i]        =     powl(10.0,i/-10.0);
    }

    for(int i=1;i<=MAXCOV;i++){
	//cout<<i<<endl;

	for(int j=0;j<=i;j++){	    
	    binomVec[i][j] = ( logl(nChoosek(i,j))+logl(powl(0.5,i)) );	     
	}

    }


    overestimateFactor[0] = 0.0;
    overestimateFactor[1] = 1.0;

    for(int i=2;i<=MAXCOV;i++){

	// 2*(1/2^i)=1/2^{i-1} will appear as homoz. So the fraction of het. will be 1-1/2^{i-1}
	// the underestimate will be by a factor of 1.0/ (  1.0- (1/powl(2.0,i-1))  );
	
	overestimateFactor[i] = 1.0/ (  1.0- (1/powl(2.0,i-1))  );

    }


#ifdef DEBUGINITSCORES
    cerr<<"qQUAL"<<"\t"<<"likeMatch"<<"\t"<<"likeMismatch"<<"\t"<<"likeMatchProb"<<"\t"<<"likeMismatchProb"<<endl;
    for(int i=0;i<=MAXBASEQUAL;i++){
	cerr<<"q= "<<i<<"\t"<<likeMatch[i]<<"\t"<<likeMismatch[i]<<"\t"<<likeMatchProb[i]<<"\t"<<likeMismatchProb[i]<<endl;
    }

    cerr<<"mQUAL"<<"\t"<<"likeMatchProb"<<"\t"<<"likeMismatchProb"<<endl;
    for(int i=0;i<=MAXMAPPINGQUAL;i++){
	//cerr<<"m= "<<i<<"\t"<<likeMatchMap[i]<<"\t"<<likeMismatchMap[i]<<"\t"<<likeMatchProbMap[i]<<"\t"<<likeMismatchProbMap[i]<<endl;
	cerr<<"m= "<<i<<"\t"<<likeMatchProbMap[i]<<"\t"<<likeMismatchProbMap[i]<<endl;
    }
#endif

}//END initScores()


//! A method to combined 5' deam rates and 3' deam rates
/*!
  This method is called by the initDeamProbabilities.
  It uses the "worse" deamination
*/
void combineDeamRates(long double f1[4],long double f2[4],long double f[4],int b){
    
    long double minFreq = MIN2( f1[b] , f2[b] );
    //cout<<b<<"\t"<<f1[b]<<"\t"<<f2[b]<<"\t"<<f[b]<<"\t"<<minFreq<<endl;

    if(f1[b] == minFreq){//use f1
	for(int i=0;i<4;i++)
	    f[i] = f1[i];

    }else{
	if(f2[b] == minFreq){//use f2
	    for(int i=0;i<4;i++)
		f[i] = f2[i];
	}else{    
	    cerr<<"ERROR in combineDeamRates(), wrong state"<<endl;
	    exit(1);
	}
    }
}

//! A method to initialize the deamination probabilities
/*!
  This method is called by the main 
*/
void initDeamProbabilities(const string & deam5pfreqE,const string & deam3pfreqE){

    
    vector<substitutionRates> sub5pT;
    vector<substitutionRates> sub3pT;

    readNucSubstitionRatesFreq(deam5pfreqE,sub5pT);
    readNucSubstitionRatesFreq(deam3pfreqE,sub3pT);;

    //5'
    for(unsigned int i=0;i<sub5pT.size();i++){
	probSubstition toadd;
	for(int nuc1=0;nuc1<4;nuc1++){
	    double probIdentical=1.0;

	    for(int nuc2=0;nuc2<4;nuc2++){
		int nuc = nuc1*4+nuc2;
		if(nuc1==nuc2) continue;
		int ind2 = dimer2indexInt( nuc1,nuc2  );
		probIdentical = probIdentical-sub5pT[i].s[ ind2 ];
		toadd.s[nuc]  =               sub5pT[i].s[ ind2 ];
	    }
	    
	    if(probIdentical<0){
		cerr<<"Error with deamination profile, identity probability is less than 0"<<endl;
		exit(1);
	    }

	    int nuc = nuc1*4+nuc1;
	    toadd.s[nuc] = probIdentical;	    	    	    
	}
	sub5p.push_back(toadd);
    }

    //copying to the rest of the fragment the last position
    for(unsigned int i=(sub5p.size()-1);i<MAXLENGTHFRAGMENT;i++){     
	sub5p.push_back( sub5p[sub5p.size()-1] );
    }

    //copying to di-nucleotides
    for(unsigned int i=0;i<MAXLENGTHFRAGMENT;i++){     
	diNucleotideProb diNuctoadd;
	for(int nuc1=0;nuc1<4;nuc1++){
	    for(int nuc2=0;nuc2<4;nuc2++){
		int nuc = nuc1*4+nuc2;
		diNuctoadd.p[nuc1][nuc2]=sub5p[i].s[nuc];		
	    }	
	}
	sub5pDiNuc.push_back(diNuctoadd);
    }


    //3'
    for(unsigned int i=0;i<sub3pT.size();i++){
	probSubstition toadd;
	for(int nuc1=0;nuc1<4;nuc1++){
	    double probIdentical=1.0;

	    for(int nuc2=0;nuc2<4;nuc2++){
		int nuc = nuc1*4+nuc2;
		if(nuc1==nuc2) continue;
		int ind2 = dimer2indexInt( nuc1,nuc2  );
		probIdentical = probIdentical-sub3pT[i].s[ ind2 ];
		toadd.s[nuc]  =               sub3pT[i].s[ ind2 ];
	    }
	    
	    if(probIdentical<0){
		cerr<<"Error with deamination profile, identity probability is less than 0"<<endl;
		exit(1);
	    }

	    int nuc = nuc1*4+nuc1;
	    toadd.s[nuc] = probIdentical;	    	    	    
	}
	sub3p.push_back(toadd);
    }

    //copying to the rest of the fragment the last position
    for(unsigned int i=(sub3p.size()-1);i<MAXLENGTHFRAGMENT;i++){
	sub3p.push_back( sub3p[sub3p.size()-1] );
    }

    //copying to di-nucleotides
    for(unsigned int i=0;i<MAXLENGTHFRAGMENT;i++){     
	diNucleotideProb diNuctoadd;
	for(int nuc1=0;nuc1<4;nuc1++){
	    for(int nuc2=0;nuc2<4;nuc2++){
		int nuc = nuc1*4+nuc2;
		diNuctoadd.p[nuc1][nuc2]=sub3p[i].s[nuc];		
	    }	
	}
	sub3pDiNuc.push_back(diNuctoadd);
    }


#ifdef DEBUGDEAM

    cerr<<endl<<"-- 5' deamination rates --"<<endl;
    for(unsigned int i=0;i<sub5p.size();i++){
    	cerr<<"i="<<i<<" - ";
    	for(int nuc1=0;nuc1<4;nuc1++){
    	    for(int nuc2=0;nuc2<4;nuc2++){
    		int nuc = nuc1*4+nuc2;
    		//cout<<nuc1<<"\t"<<nuc2<<"\t"
    		cerr<<sub5p[i].s[nuc]<<" ";
    	    }
    	    cerr<<" - ";
    	}
    	cerr<<endl;
    }

    cerr<<"------------------"<<endl;

    for(unsigned int i=0;i<sub5p.size();i++){
    	cerr<<"i="<<i<<endl<<"\t";
    	for(int nuc1=0;nuc1<4;nuc1++)
	    cerr<<"ACGT"[nuc1]<<"\t";
	cerr<<endl;
    	for(int nuc1=0;nuc1<4;nuc1++){
	    cerr<<"ACGT"[nuc1]<<"\t";
    	    for(int nuc2=0;nuc2<4;nuc2++){
    		//int nuc = nuc1*4+nuc2;
    		//cout<<nuc1<<"\t"<<nuc2<<"\t"
    		cerr<<sub5pDiNuc[i].p[nuc1][nuc2]<<"\t";
    	    }
    	    cerr<<endl;
    	}
    	cerr<<endl;
    }

    cerr<<"-- 3' deamination rates --"<<endl;
    for(unsigned int i=0;i<sub3p.size();i++){
    	cerr<<"i="<<i<<" - ";
    	for(int nuc1=0;nuc1<4;nuc1++){
    	    for(int nuc2=0;nuc2<4;nuc2++){
    		int nuc = nuc1*4+nuc2;
    		//cout<<nuc1<<"\t"<<nuc2<<"\t"
    		cerr<<sub3p[i].s[nuc]<<" ";
    	    }
    	    cerr<<" - ";
    	}
    	cerr<<endl;
    }


    cerr<<"------------------"<<endl;

    for(unsigned int i=0;i<sub3p.size();i++){
    	cerr<<"i="<<i<<endl<<"\t";
    	for(int nuc1=0;nuc1<4;nuc1++)
	    cerr<<"ACGT"[nuc1]<<"\t";
	cerr<<endl;
    	for(int nuc1=0;nuc1<4;nuc1++){
	    cerr<<"ACGT"[nuc1]<<"\t";
    	    for(int nuc2=0;nuc2<4;nuc2++){
    		//int nuc = nuc1*4+nuc2;
    		//cout<<nuc1<<"\t"<<nuc2<<"\t"
    		cerr<<sub3pDiNuc[i].p[nuc1][nuc2]<<"\t";
    	    }
    	    cerr<<endl;
    	}
    	cerr<<endl;
    }

#endif

    //dummy values
    for(unsigned int L=0;L<MINLENGTHFRAGMENT;L++){     //for each fragment length
	vector<probSubstition> subDeam_;
	vector<diNucleotideProb> subDeamDiNuc_;

	subDeam.push_back(      subDeam_      );
	subDeamDiNuc.push_back( subDeamDiNuc_ );       
    }


    for(unsigned int L=MINLENGTHFRAGMENT;L<=MAXLENGTHFRAGMENT;L++){     //for each fragment length
	printprogressBarCerr( float(L-MINLENGTHFRAGMENT)/float(MAXLENGTHFRAGMENT-MINLENGTHFRAGMENT) );
	// if( (L%8)==0){
	//     cerr<<".";
	// }
	//vector<alleleFrequency > defaultDNA;
	vector<probSubstition> subDeam_;
	vector<diNucleotideProb> subDeamDiNuc_;
	for(unsigned int l=0;l<L;l++){     //position
	    probSubstition    subDeam__;
	    diNucleotideProb  subDeamDiNuc__;

	    for(int b1=0;b1<4;b1++){//original base
		//pick the highest deam rates for that position
		combineDeamRates(sub5pDiNuc[l].p[b1] ,  sub3pDiNuc[L-l-1].p[b1] , subDeamDiNuc__.p[b1], b1);
		
		for(int b2=0;b2<4;b2++){//post deam base
		    int b = b1*4+b2;
		    subDeam__.s[b]  = subDeamDiNuc__.p[b1][b2];
		}

	    }
	    subDeam_.push_back(      subDeam__      );
	    subDeamDiNuc_.push_back( subDeamDiNuc__ );
	}

	subDeam.push_back(      subDeam_      );
	subDeamDiNuc.push_back( subDeamDiNuc_ );       
    }
    cerr<<endl; //flushing progress bar
    


#ifdef DEBUGDEAM

    cerr<<"-- per length deamination rates --"<<endl;

    for(unsigned int L=MINLENGTHFRAGMENT;L<=MAXLENGTHFRAGMENT;L++){     //for each fragment length

	cerr<<endl<<"L="<<L<<endl;
	for(unsigned int l=0;l<L;l++){     //position

	    cerr<<"l="<<l<<" - ";
	    for(int nuc1=0;nuc1<4;nuc1++){
		for(int nuc2=0;nuc2<4;nuc2++){
		    int nuc = nuc1*4+nuc2;
		    cerr<<subDeam[L][l].s[nuc]<<" ";
		}
		cerr<<" - ";
	    }
	    cerr<<endl;
	}
    }



    for(unsigned int L=MINLENGTHFRAGMENT;L<=MAXLENGTHFRAGMENT;L++){     //for each fragment length

	cerr<<endl<<"L="<<L<<endl;
	for(unsigned int l=0;l<L;l++){     //position

	    cerr<<"l="<<l<<" - "<<endl;
	    
	    for(int nuc1=0;nuc1<4;nuc1++){
		cerr<<"ACGT"[nuc1]<<"\t";
		for(int nuc2=0;nuc2<4;nuc2++){
		    cerr<<subDeamDiNuc[L][l].p[nuc1][nuc2]<<"\t";
		}
		cerr<<endl;
	    }
	    cerr<<endl;
	}
    }
    


#endif

    //exit(1);
    
    //if no deamination, cannot have a mismatch
    for(int b1=0;b1<4;b1++){
	for(int b2=0;b2<4;b2++){
	    int b = b1*4+b2;
	    if(b1==b2){
		defaultSubMatch.s[ b ]           = 1.0;	
		defaultSubMatchMatrix.p[b1][b2]  = 1.0;	
	    }else{
		defaultSubMatch.s[ b ]           = 0.0;	
		defaultSubMatchMatrix.p[b1][b2]  = 0.0;	
	    }
	}
    }



}//end initDeamProbabilities







// //! A method to initialize the probability of seeing a base if a fragment is mismapped
// /*!
//   This method is called by the main after reading the deamination profiles
//   as deamination will bias this.
// */
// void initDefaultBaseFreq(const string & dnafreqFile){

//     readDNABaseFreq( dnafreqFile  ,  dnaDefaultBases );

//     for(unsigned int i=0;i<MAXLENGTHFRAGMENT;i++){     
// 	//defaultDNA5p
// 	alleleFrequency dnaFreq5p_;
// 	alleleFrequency dnaFreq3p_;

// 	for(int b=0;b<4;b++){//original base
// 	    dnaFreq5p_.f[b]=0;
// 	    dnaFreq3p_.f[b]=0;
// 	}

// 	for(int b1=0;b1<4;b1++){//    observed base

// 	    for(int b2=0;b2<4;b2++){ //original base base
// 		int b = b2*4+b1;
// 		//cerr<<"pos="<<i<<"\t"<<"ACGT"[b1]<<"\t"<<"ACGT"[b2]<<"\t"<<b<<"\t"<<dnaDefaultBases.f[b1]<<"\t"<<sub5p[i].s[b]<<"\t"<<sub3p[i].s[b]<<endl;
// 		dnaFreq5p_.f[b1] += dnaDefaultBases.f[b2]*sub5p[i].s[b];
// 		dnaFreq3p_.f[b1] += dnaDefaultBases.f[b2]*sub3p[i].s[b];
// 		//cerr<<"\t"<<dnaFreq5p_.f[b1]<<"\t"<<dnaFreq3p_.f[b1]<<endl;
// 	    }
// 	}
	
// 	defaultDNA5p.push_back(dnaFreq5p_);
// 	defaultDNA3p.push_back(dnaFreq3p_);	
//     }

// #ifdef DEBUGDEFAULTFREQ

//     cerr<<"-- Default base frequencies --"<<endl;
//     cerr<<dnaDefaultBases.f[0]<<"\t"<<dnaDefaultBases.f[1]<<"\t"<<dnaDefaultBases.f[2]<<"\t"<<dnaDefaultBases.f[3]<<endl;

//     cerr<<"-- 5' base frequencies --"<<endl;
    
//     for(unsigned int i=0;i<MAXLENGTHFRAGMENT;i++){     
//     	cerr<<"i="<<i<<" - ";
// 	cerr<<defaultDNA5p[i].f[0]<<"\t"<<defaultDNA5p[i].f[1]<<"\t"<<defaultDNA5p[i].f[2]<<"\t"<<defaultDNA5p[i].f[3]<<"\t"<<(defaultDNA5p[i].f[0]+defaultDNA5p[i].f[1]+defaultDNA5p[i].f[2]+defaultDNA5p[i].f[3])<<endl;
//     }

//     cerr<<"-- 3' base frequencies --"<<endl;

//     for(unsigned int i=0;i<MAXLENGTHFRAGMENT;i++){     
//     	cerr<<"i="<<i<<" - ";
// 	cerr<<defaultDNA3p[i].f[0]<<"\t"<<defaultDNA3p[i].f[1]<<"\t"<<defaultDNA3p[i].f[2]<<"\t"<<defaultDNA3p[i].f[3]<<"\t"<<(defaultDNA3p[i].f[0]+defaultDNA3p[i].f[1]+defaultDNA3p[i].f[2]+defaultDNA3p[i].f[3])<<endl;
//     }

    
// #endif


// }//end initDeamProbabilities


void initDefaultBaseFreq(const string & dnafreqFile){


    readDNABaseFreq( dnafreqFile  ,  dnaDefaultBases );
    //vector< vector<alleleFrequency> > defaultDNAFragmentLength;
    for(unsigned int L=0;L<MINLENGTHFRAGMENT;L++){     //for each fragment length
	vector<alleleFrequency > dummy;
	defaultDNAFragmentLength.push_back(dummy);
    }

    for(unsigned int L=MINLENGTHFRAGMENT;L<=MAXLENGTHFRAGMENT;L++){     //for each fragment length
	// if( (L%8)==0){
	//     cerr<<".";
	// }
	printprogressBarCerr( float(L-MINLENGTHFRAGMENT)/float(MAXLENGTHFRAGMENT-MINLENGTHFRAGMENT) );
	
	vector<alleleFrequency > defaultDNA;
	for(unsigned int l=0;l<L;l++){     //position

	    alleleFrequency dnaFreq_;

	    for(int b=0;b<4;b++){//original base
		dnaFreq_.f[b]=0;		
	    }
	   

	    for(int b1=0;b1<4;b1++){     // observed base		
		for(int b2=0;b2<4;b2++){ // original base base
		    int b = b2*4+b1;
		    dnaFreq_.f[b1] += dnaDefaultBases.f[b2]*subDeam[L][l].s[b];
		}
	    }
	    defaultDNA.push_back(dnaFreq_);
	}
    
#ifdef DEBUGDEFAULTFREQ
	cerr<<"-- Defautl base frequencies --"<<endl;
	cerr<<dnaDefaultBases.f[0]<<"\t"<<dnaDefaultBases.f[1]<<"\t"<<dnaDefaultBases.f[2]<<"\t"<<dnaDefaultBases.f[3]<<endl;
    
	//    for(unsigned int i=0;i<MAXLENGTHFRAGMENT;i++){
	//for(unsigned int l=MINLENGTHFRAGMENT;L<MAXLENGTHFRAGMENT;L++){     //for each fragment length
	cerr<<"--  base frequencies at length L="<<L<<" --"<<endl;
	
	for(unsigned int l=0;l<L;l++){     
	    cerr<<"pos="<<l<<" - ";
	    cerr<<defaultDNA[l].f[0]<<"\t"<<defaultDNA[l].f[1]<<"\t"<<defaultDNA[l].f[2]<<"\t"<<defaultDNA[l].f[3]<<"\t"<<(defaultDNA[l].f[0]+defaultDNA[l].f[1]+defaultDNA[l].f[2]+defaultDNA[l].f[3])<<endl;
	}

#endif
    
	defaultDNAFragmentLength.push_back(defaultDNA);
    }
    cerr<<endl; //flushing progress bar
    //exit(1);
}//end initDefaultBaseFreq 







//! A method to initialize the likelihood of observing data to avoid recomputation
/*!
  This method is called by the main after calling initDefaultBaseFreq
*/
void initLikelihoodScores(){


    //vector< vector< mpq2bsq2submatrix * > >  length2pos2mpq2bsq2submatrix;
    //insert dummy values
    for(int L=0;L<int(MINLENGTHFRAGMENT);L++){//for each fragment length
	vector< mpq2bsq2submatrix * >  vectorToAdd ;	
	length2pos2mpq2bsq2submatrix.push_back(vectorToAdd);
    }

    mpq2bsq2submatrix  * mpq2BaseQualSubMatrix;    
    bool firstIteration=true;
    for(int L=int(MINLENGTHFRAGMENT);L<=int(MAXLENGTHFRAGMENT);L++){//for each fragment length
	// if( (L%8)==0){
	//     cerr<<".";
	// }
	printprogressBarCerr( float(L-MINLENGTHFRAGMENT)/float(MAXLENGTHFRAGMENT-MINLENGTHFRAGMENT) );

	vector< mpq2bsq2submatrix * > pos2mpq2BaseQual2SubMatrix;

	for(int l=0;l<L;l++){//for each pos

	    if(specifiedDeam){//specied deamination rates
		mpq2BaseQualSubMatrix = new mpq2bsq2submatrix();
	    }else{
		if(!firstIteration){//if it's not the first iteration, since every position is the same we reuse the one from the first iteration
		    pos2mpq2BaseQual2SubMatrix.push_back(mpq2BaseQualSubMatrix);
		    continue;//skip to next iteration
		}
		mpq2BaseQualSubMatrix = new mpq2bsq2submatrix();//we will compute for a single iteration
		firstIteration=false;
	    }
	    
	    for(int mq=0;mq<=MAXMAPPINGQUAL;mq++){ //for each mapping quality 
		vector<diNucleotideProb>  baseQual2SubMatrix;//FOR EACH QUAL SCORE
		
		for(int q=0;q<=MAXBASEQUAL;q++){ //for each base quality
		    diNucleotideProb toAddForBaseQual;
		    diNucleotideProb toAddForBaseQual_;

		    for(int bTheo=0;bTheo<4;bTheo++){
			for(int bObs=0;bObs<4;bObs++){
			    toAddForBaseQual_.p[bTheo][bObs]=0.0;			 
			}
		    }
		    
		    //marginalize over each potential base post deamination bpstDeam
		    //This computes P[bObs | bTheo] in the absence of mismapping
		    //P[bObs | bTheo]  =  sum_{bpstDeam = A,C,G,T} P[bpstDeam|bTheo] * P[bObs | bpstDeam] 		   
		    for(int bTheo=0;bTheo<4;bTheo++){                   // each possible theoritical base
			for(int bObs=0;bObs<4;bObs++){	                // each possible observed base
			    for(int bpstDeam=0;bpstDeam<4;bpstDeam++){	// each possible deaminated base					    
				toAddForBaseQual_.p[bTheo][bObs] +=       
				    
				    subDeamDiNuc[L][l].p[bTheo][bpstDeam] * //P[ bTheo->bpstDeam ]

				    // P[ bObs|bpstDeam ] = 
				    // P[ bObs|bpstDeam , !E ] * (1-e) + P[ bObs,bpstDeam , E ] * (e) + 
				    (
				     likeMatchProb[q]    *   defaultSubMatchMatrix.p[bpstDeam][bObs]  // P[ bpstDeam->bObs | !E ] * (1-e) 
				     +
				     likeMismatchProb[q] * illuminaErrorsProbMatrix.p[bpstDeam][bObs] // P[ bpstDeam->bObs |  E ] *   (e) 
				     );
			    }//end each possible deaminated base
			}//end each possible observed base
		    }//end each possible theoritical base
		    


		    //This computes P[bObs|bTheo] and includes the mapping quality
		    
		    for(int bTheo=0;bTheo<4;bTheo++){               // each possible theoretical base
			for(int bObs=0;bObs<4;bObs++){	        // each possible observed base
			    
			    //0.5 since this is the prior prob of having sampled a given chromosome
			    // m=correctly mapped
			    // P[bObs|bTheo] =  P[bObs|bTheo,m]  P[m] + P[bObs|bTheo,!m]  P[!m] 
			    toAddForBaseQual.p[bTheo][bObs] = logl(0.5 * 
								   (likeMatchProbMap[mq]    * toAddForBaseQual_.p[bTheo][bObs] + 
								    likeMismatchProbMap[mq] * defaultDNAFragmentLength[L][l].f[bObs])
								   );  
			    //toAddForBaseQual3p.p[bTheo][bObs] = logl(0.5 * (likeMatchProbMap[mq]*toAddForBaseQual3p_.p[bTheo][bObs] + likeMismatchProbMap[mq] * defaultDNAFragmentLength[L][l].f[bObs]defaultDNA3p[l].f[bObs]) );
						
			}//for each bObs
		    }//for each bTheo, bpstDeam and bObs
		    


		    baseQual2SubMatrix.push_back(toAddForBaseQual);
		    		
		    
		}//end for each base quality

		mpq2BaseQualSubMatrix->push_back( baseQual2SubMatrix );
		
	    }//end for each mapping quality mq

	    
	    
	    pos2mpq2BaseQual2SubMatrix.push_back(mpq2BaseQualSubMatrix);
	}//for each position l
	

	length2pos2mpq2bsq2submatrix.push_back(pos2mpq2BaseQual2SubMatrix);	

    }//for each fragment length L
    cerr<<endl; //flushing progress bar

    // cerr<<setprecision(20)<<"TEST2 AT 38 "<<pos2mpq2BaseQual2SubMatrix5p[  12  ][37][38].p[0][3]<<"\tL2P="<<length2pos2mpq2bsq2submatrix[113][12]->at(37)[38].p[0][3]<<endl ;
    // cerr<<setprecision(20)<<"TEST2 AT 18 "<<pos2mpq2BaseQual2SubMatrix5p[  12  ][37][18].p[0][3]<<"\tL2P="<<length2pos2mpq2bsq2submatrix[113][12]->at(37)[18].p[0][3]<<endl ;

#ifdef DEBUGINITLIKELIHOODSCORES2

    for(unsigned int L=MINLENGTHFRAGMENT;L<=MAXLENGTHFRAGMENT;L++){//for each fragment length

	//vector< vector< mpq2bsq2submatrix  > >  length2pos2mpq2bsq2submatrix_;
	for(unsigned int l=0;l<L;l++){//for each pos
	    cerr<<"pos "<<l<<"/"<<L<<"\t"<<length2pos2mpq2bsq2submatrix[L][l]->size()<<endl;
	    
	    if( l == DEBUGINITLIKELIHOODSCORES2){

		for(int mq=0;mq<=MAXMAPPINGQUAL;mq++){ //for each mapping quality 
		    //cerr<<"mq="<<mq<<"\t"<<endl;//length2pos2mpq2bsq2submatrix[L][l]->at(mq).size()<<endl;

		    for(int q=0;q<=MAXBASEQUAL;q++){ //for each base quality
			cerr<<"pos "<<l<<"/"<<L<<" mq="<<mq<<" bq="<<q<<endl;

			//cerr<<" mq = "<<mq<<" ("<<likeMatchProbMap[mq]<<" "<<likeMismatchProbMap[mq]<<" )  q = "<<q<<endl;
			for(int nuc1=0;nuc1<4;nuc1++){
			    cerr<<"ACGT"[nuc1]<<"\t";		    
			    long double s=0;//sum of probs

			    for(int nuc2=0;nuc2<4;nuc2++){
				cerr<<2*exp(length2pos2mpq2bsq2submatrix[L][l]->at(mq)[q].p[nuc1][nuc2])<<"\t";
				//cerr<<2*exp(length2pos2mpq2bsq2submatrix[L][l][mq][q].p[nuc1][nuc2])<<"\t";
				//length2pos2mpq2bsq2submatrix[L][l]->at(mq)
				//toAddForBaseQual5p.p[nuc1][nuc2]<<"\t";
				s+=exp( length2pos2mpq2bsq2submatrix[L][l]->at(mq)[q].p[nuc1][nuc2] );
			    }

			    cerr<<"\tsum="<<s<<endl;
			}
			cerr<<endl;
		    
		    }//for each base qual	   		    

		}// for each mapping quality	   
	    }
	}//for each pos
    }//for each fragment length

#endif    

    //cerr<<"done "<<endl;

    //exit(1);


}// END initLikelihoodScores()



void deleteLikelihoodScores(){


    //vector< vector< mpq2bsq2submatrix * > >  length2pos2mpq2bsq2submatrix;
    //insert dummy values
    // for(int L=0;L<int(MINLENGTHFRAGMENT);L++){//for each fragment length
    // 	vector< mpq2bsq2submatrix * >  vectorToAdd ;	
    // 	length2pos2mpq2bsq2submatrix.push_back(vectorToAdd);
    // }

    bool firstIteration=true;
    for(int L=int(MINLENGTHFRAGMENT);L<=int(MAXLENGTHFRAGMENT);L++){//for each fragment length
	for(int l=0;l<L;l++){//for each pos
	    if(specifiedDeam){//specied deamination rates
		delete(length2pos2mpq2bsq2submatrix[L][l]);	   	  
	    }else{
		if(firstIteration){//if it's not the first iteration, since every position is the same we reuse the one from the first iteration
		    delete(length2pos2mpq2bsq2submatrix[L][l]);	   	  
		}       	
		firstIteration=false;
	    }
	}//for each position l
    }//for each fragment length L

}


// //! A method to compute the prob. of seeing an observed base given an original allele
// /*!
//   This method is called by the computeLL to compute the prob. of seeing observed base ob 
//   with error rate q, a of deamination probDeam and original allele al. This method is called 
//   for each base.
// */
// inline long double computeBaseAl2Obs(const int al,
// 				     const int ob,
// 				     const int q,
// 				     const probSubstition * probDeam      , //rate of deamination
// 				     const bool isRev,
// 				     const long double mismappingProb){
//     int dinucal2al=-1;

//     if( isRev ){		    
// 	dinucal2al =     complementInt(al)*4+complementInt(al);     //genotype is 1, observed is obsBase
//     }else{
// 	dinucal2al =                   al*4+               al;       //genotype is 1, observed is obsBase
//     }

//     if(probDeam->s[dinucal2al] == 1.0){ //simple case, no deamination
// 	int dinucIndexal2ob;


// 	if( isRev ){
// 	    dinucIndexal2ob =     complementInt(al)*4+complementInt(ob);     //genotype is al, observed is ob
// 	}else{
// 	    dinucIndexal2ob =                   al *4+              ob;      //genotype is al, observed is ob
// 	}
                                                            
	
// 	return ( (1.0-mismappingProb)*(
// 				       likeMatchProb[    q ]*defaultSubMatch.s[   dinucIndexal2ob]
// 				       + 
// 				       likeMismatchProb[ q ]*illuminaErrorsProb.s[dinucIndexal2ob] 
// 				       )
// 		 +
// 		 mismappingProb*randomPMatch4Bases);

//     }else{
	
// 	long double sumProbToReturn = 0.0;
// 	for(int alpostdeam=0;alpostdeam<4;alpostdeam++){
// 	    int dinucal2ald = -1;
// 	    if( isRev ){		    
// 		dinucal2ald =     complementInt(al)*4+complementInt(alpostdeam);     //genotype is 1, observed is obsBase
// 	    }else{
// 		dinucal2ald =                   al*4+               alpostdeam;       //genotype is 1, observed is obsBase
// 	    }
	    

// 	    if(probDeam->s[dinucal2ald] > 0.0){ //has deamination       	
// 		int dinucald2ob = -1;
// 		if( isRev ){		    
// 		    dinucald2ob =     complementInt(alpostdeam)*4+complementInt(ob);     //genotype is 1, observed is obsBase
// 		}else{
// 		    dinucald2ob =                   alpostdeam*4+               ob;       //genotype is 1, observed is obsBase
// 		}

// 		sumProbToReturn += probDeam->s[dinucal2ald]*( (1.0-mismappingProb)*(
// 										    likeMatchProb[    q ]*defaultSubMatch.s[   dinucald2ob]
// 										    + 
// 										    likeMismatchProb[ q ]*illuminaErrorsProb.s[dinucald2ob] 
// 										    )
// 							      +
// 							      mismappingProb*randomPMatch4Bases);
// 	    }
// 	}

// 	return sumProbToReturn;
//     }

//     return -1;
// }//end computeBaseAl2Obs



//! A method to compute the prior prob.  for each genotype
/*!
   This method is called by the computeLL to compute the prior prob. for each genotype
   
   \param h : The heterozygosity rate

*/

inline void computePriorMatrix(long double h,
			       diNucleotideProb * priorGenotype,
			       // diNucleotideProb * priorGenotypeProb,
			       // diNucleotideProb * priorGenotypeProbD,
			       diNucleotideProb * priorGenotypeD1,
			       diNucleotideProb * priorGenotypeD2){

    for(int ba=0;ba<4;ba++){//ancestral base

	for(int bd=0;bd<4;bd++){//derived base
	    if(ba == bd){
		//priorGenotype.p[ba][bd]     = dnaDefaultBases.f[ba]  *    (1.0-h);
		priorGenotype->p[ba][bd]       = logl(dnaDefaultBases.f[ba])  +    logl(1.0-h);
		
		//                                     DER( k*(1-h) ) = -k
		// priorGenotypeProb->p[ba][bd]   =      dnaDefaultBases.f[ba]   *        (1.0-h);
		// priorGenotypeProbD->p[ba][bd]  =      dnaDefaultBases.f[ba]   *        (   -1.0);		    
		//priorGenotypeD1->p[ba][bd]     = 1/( (h-1.0)*( (dnaDefaultBases.f[ba])  +  logl(1.0-h) ));
		priorGenotypeD1->p[ba][bd]     = 1/( (h-1) );               // d/dh(log(1 - h)) = 1/(h - 1)
		priorGenotypeD2->p[ba][bd]     = -1.0/( powl( (1-h),2) ); // d^2/dh^2(log(1 - h)) = -1/(1 - h)^2

		
		//cerr<<"ACGT"[ba]<<"\t"<<"ACGT"[bd]<<"\t"<<dnaDefaultBases.f[ba]<<" X "<<(1.0-h)<<endl;
	    }else{//mutation
		if( (ba%2)==(bd%2) ){//transition
		    //priorGenotype->p[ba][bd] = dnaDefaultBases.f[ba]  *  ( (h) * (TStoTVratio/(TStoTVratio+1.0)) );
		    
		    //                                     DER( k*(h) ) = k
		    // priorGenotypeProb->p[ba][bd]     =      dnaDefaultBases.f[ba]   *         (h) * (TStoTVratio/(TStoTVratio+1.0))     ;
		    // priorGenotypeProbD->p[ba][bd]    =      dnaDefaultBases.f[ba]   *         1.0 * (TStoTVratio/(TStoTVratio+1.0))     ;
		    if(tvonly){
			priorGenotype->p[ba][bd]         = logl(dnaDefaultBases.f[ba])  +   logl(  0   );
			priorGenotypeD1->p[ba][bd]       = 0;
			priorGenotypeD2->p[ba][bd]       = 0;
		    }else{
			priorGenotype->p[ba][bd]         = logl(dnaDefaultBases.f[ba])  +   logl( (h) * (TStoTVratio/(TStoTVratio+1.0))     );
			priorGenotypeD1->p[ba][bd]       = 1/(  h );// (d)/(dh)( log(c h) ) = 1/h
			priorGenotypeD2->p[ba][bd]       = -1.0/(  powl(h,2) ); //d^2/dh^2(log(h)) = -1/h^2			
		    }

		    
		    //cerr<<"ACGT"[ba]<<"\t"<<"ACGT"[bd]<<"\t"<<dnaDefaultBases.f[ba]<<" X "<< ( (h) * (TStoTVratio/(TStoTVratio+1.0)) )<<endl; 
		}else{//transversion


		    //                                     DER( k*(h) ) = k
		    // priorGenotypeProb->p[ba][bd]  =      dnaDefaultBases.f[ba]   *       ( (h) * ((        1.0/(TStoTVratio+1.0)) /2.0));
		    // priorGenotypeProbD->p[ba][bd] =      dnaDefaultBases.f[ba]   *       ( 1.0 * ((        1.0/(TStoTVratio+1.0)) /2.0));
		    //priorGenotype->p[ba][bd] = dnaDefaultBases.f[ba]  * (( (h) * (        1.0/(TStoTVratio+1.0)) )/2.0);
		    if(tvonly){
			priorGenotype->p[ba][bd]      = logl(dnaDefaultBases.f[ba])  +   logl( (h) * (         1.0/2.0  ));
			priorGenotypeD1->p[ba][bd]   = 1/h; // (d)/(dh)( log(c h) ) = 1/h
			priorGenotypeD2->p[ba][bd]   = -1.0/( powl(h,2) ); // d^2/dh^2(log(h)) = -1/h^2
		    }else{
			priorGenotype->p[ba][bd]      = logl(dnaDefaultBases.f[ba])  +   logl( (h) * ((        1.0/(TStoTVratio+1.0)) /2.0));
			priorGenotypeD1->p[ba][bd]   = 1/h; // (d)/(dh)( log(c h) ) = 1/h
			priorGenotypeD2->p[ba][bd]   = -1.0/( powl(h,2) ); // d^2/dh^2(log(h)) = -1/h^2
		    }
		    //cerr<<"ACGT"[ba]<<"\t"<<"ACGT"[bd]<<"\t"<<dnaDefaultBases.f[ba]<<" X "<< ( (h) * (        1.0/(TStoTVratio+1.0)) )/2.0<<endl; 
		}
	    }
	}
    }

}//END computePriorMatrix





inline void preComputeBaBdLikelihood(const vector<positionInformation> * piForGenomicWindow, 
				     vector<babdlikelihood> * vectorBaBdLikelihood,
				     const unsigned int numberOfSitesInWindow){

    for(unsigned int p=0;p<numberOfSitesInWindow;p++){
	
	babdlikelihood bblikeForGivenPos;

	int         babdIdx          =0;

#ifdef DEBUGPRECOMPUTEBABD
	bool printDEBUG = (piForGenomicWindow->at(p).posAlign==DEBUGPRECOMPUTEBABD);
#endif
       	
	for(uint8_t ba=0;ba<4;ba++){//ancestral base
	    uint8_t ba_c = 3-ba;

	    for(uint8_t bd=0;bd<4;bd++){//derived base
		uint8_t bd_c = 3-bd;

				
		long double loglikelihoodForGivenBaBd          =0.0;

		
		for(unsigned int i=0;i<piForGenomicWindow->at(p).readsVec.size();i++){ //for each fragment at pos p
		    //Likelihood it comes from A
		    //char bObs=piForGenomicWindow->at(p).readsVec[i].base;
		    
		    //#ifdef PRECOMPUTELOG
		    long double llA; //Likelihood it comes from A
		    long double llD; //Likelihood it comes from D
		    // cerr<<p<<" "<<i<<" "<<int(piForGenomicWindow->at(p).readsVec[i].lengthF)<<" "<<int(piForGenomicWindow->at(p).readsVec[i].pos5p)<<" "<<int(piForGenomicWindow->at(p).readsVec[i].mapq)<<" "<<int(piForGenomicWindow->at(p).readsVec[i].qual)<<" "<<int(ba_c)<<" "<<int(piForGenomicWindow->at(p).readsVec[i].base)<<endl;

#ifdef DEBUGSINGLEREAD
		    cerr<<"name: "<<piForGenomicWindow->at(p).readsVec[i].name<<endl;
#endif
		    
		    // cerr<<"1 "<<length2pos2mpq2bsq2submatrix[piForGenomicWindow->at(p).readsVec[i].lengthF].size()<<endl;
		    // cerr<<"2 "<<length2pos2mpq2bsq2submatrix
		    // 	    [piForGenomicWindow->at(p).readsVec[i].lengthF]
		    // 	[piForGenomicWindow->at(p).readsVec[i].pos5p].size()<<endl;
		    
		    // cerr<<"3 "<<length2pos2mpq2bsq2submatrix
		    // 	    [piForGenomicWindow->at(p).readsVec[i].lengthF]
		    // 	    [piForGenomicWindow->at(p).readsVec[i].pos5p]
		    // 	[piForGenomicWindow->at(p).readsVec[i].mapq].size()<<endl;


		    // cerr<<"4 "<<length2pos2mpq2bsq2submatrix
		    // 	    [piForGenomicWindow->at(p).readsVec[i].lengthF]
        	    // 	    [piForGenomicWindow->at(p).readsVec[i].pos5p]
		    // 	    [piForGenomicWindow->at(p).readsVec[i].mapq]
         	    // 	    [piForGenomicWindow->at(p).readsVec[i].qual].p[0][0]<<endl;


		    // cerr<<"1 "<<length2pos2mpq2bsq2submatrix
		    // 	    [piForGenomicWindow->at(p).readsVec[i].lengthF]
		    // 	    [piForGenomicWindow->at(p).readsVec[i].pos5p]
		    // 	    [piForGenomicWindow->at(p).readsVec[i].mapq]
		    // 	[piForGenomicWindow->at(p).readsVec[i].qual].p[ba_c][piForGenomicWindow->at(p).readsVec[i].base];

		    
		    if(piForGenomicWindow->at(p).readsVec[i].isrv){
			llA = length2pos2mpq2bsq2submatrix
			    [piForGenomicWindow->at(p).readsVec[i].lengthF]
			    [piForGenomicWindow->at(p).readsVec[i].pos5p]->at
			    (piForGenomicWindow->at(p).readsVec[i].mapq)
			    [piForGenomicWindow->at(p).readsVec[i].qual].p[ba_c][piForGenomicWindow->at(p).readsVec[i].base];

			llD = length2pos2mpq2bsq2submatrix
			    [piForGenomicWindow->at(p).readsVec[i].lengthF]
			    [piForGenomicWindow->at(p).readsVec[i].pos5p]->at
			    (piForGenomicWindow->at(p).readsVec[i].mapq)
			    [piForGenomicWindow->at(p).readsVec[i].qual].p[bd_c][piForGenomicWindow->at(p).readsVec[i].base];

		    }else{

			llA = length2pos2mpq2bsq2submatrix
			    [piForGenomicWindow->at(p).readsVec[i].lengthF]
			    [piForGenomicWindow->at(p).readsVec[i].pos5p]->at
			    (piForGenomicWindow->at(p).readsVec[i].mapq)
			    [piForGenomicWindow->at(p).readsVec[i].qual].p[ba][piForGenomicWindow->at(p).readsVec[i].base];

			llD = length2pos2mpq2bsq2submatrix
			    [piForGenomicWindow->at(p).readsVec[i].lengthF]
			    [piForGenomicWindow->at(p).readsVec[i].pos5p]->at
			    (piForGenomicWindow->at(p).readsVec[i].mapq)
			    [piForGenomicWindow->at(p).readsVec[i].qual].p[bd][piForGenomicWindow->at(p).readsVec[i].base];
		    }
		    //cerr<<"done"<<endl;

#ifdef DEBUGPRECOMPUTEBABD
		    if(printDEBUG){
		 
			cout<<"ba "<<"ACGT"[ba]<<" r"<<"ACGT"[ba_c]<<" bc "<<"ACGT"[bd]<<" r"<<"ACGT"[bd_c]<<" i="<<i<<" l="<<
			    int(piForGenomicWindow->at(p).readsVec[i].lengthF)<<" 5p="<<
			    int(piForGenomicWindow->at(p).readsVec[i].pos5p)<<" mp="<<			   
			    int(piForGenomicWindow->at(p).readsVec[i].mapq)<<" q="<<
			    int(piForGenomicWindow->at(p).readsVec[i].qual)<<" b="<<"ACGT"[int(piForGenomicWindow->at(p).readsVec[i].base)]<<" "<< piForGenomicWindow->at(p).readsVec[i].isrv  <<" "<<llA<<" "<<2*exp(llA)<<" "<<llD<<" "<<2*exp(llD)<<endl;
		

		    }
#endif


		    // #else
		    
		    // 		    // //without precompute log
		    // 		    long double llA = logl(0.5* length2pos2mpq2bsq2submatrix
		    // 		    			   [piForGenomicWindow->at(p).readsVec[i].lengthF]
		    // 		    			   [piForGenomicWindow->at(p).readsVec[i].pos5p]->at( piForGenomicWindow->at(p).readsVec[i].mapq)
		    // 		    			   [piForGenomicWindow->at(p).readsVec[i].qual].p[ba][piForGenomicWindow->at(p).readsVec[i].base]);

		    // 		    //Likelihood it comes from D
		    // 		    long double llD = logl(0.5* length2pos2mpq2bsq2submatrix
		    // 		    			   [piForGenomicWindow->at(p).readsVec[i].lengthF]
		    // 		    			   [piForGenomicWindow->at(p).readsVec[i].pos5p]->at( piForGenomicWindow->at(p).readsVec[i].mapq)
		    // 		    			   [piForGenomicWindow->at(p).readsVec[i].qual].p[bd][piForGenomicWindow->at(p).readsVec[i].base]);
// #endif
		    
    
		       
		    
		    loglikelihoodForGivenBaBd += oplusnatl( llA, llD); //, adding probs of llA and llD, multiplying probabilities for each site, assuming independence, (\prod_{fragment} P(D|G))

		}//END  for each fragment at pos p

		//  product of (\prod_{fragment} P(D|G)) times the prior P(G) for the genotype
		
		bblikeForGivenPos.gl[babdIdx] = loglikelihoodForGivenBaBd;
		
		
		

		babdIdx++;
	    }//END for each derived base
	}//END for each ancestral base

#ifdef DEBUGPRECOMPUTEBABD
	if(printDEBUG){ exit(1); }
#endif

	vectorBaBdLikelihood->push_back(  bblikeForGivenPos );
    }//END for each genomic position

}//preComputeBaBdLikelihood





inline void genotypePositions(const int                   mostLikelyBaBdIdx,
			      vector<long double>       * vectorOfloglikelihoodForGivenBaBd,
			      vector<long double>       * vectorOfloglikelihoodForGivenGeno,
			      PositionResult * pr){


    //Compute likelihood of all minus best BABD
    long double loglikelihoodForEveryBaBd_minusBest          =0.0;
	
    for(int g=0;g<16;g++){
	if(g!= mostLikelyBaBdIdx)
	    loglikelihoodForEveryBaBd_minusBest = oplusInitnatl( loglikelihoodForEveryBaBd_minusBest , vectorOfloglikelihoodForGivenBaBd->at(g) ); // \sum_{genotype} (\prod_{fragment} P(D|G))*P(G)
	//cerr<<babdIdx2Code16[g]<<"\t"<<vectorOfloglikelihoodForGivenBaBd[g]<<endl;
	//cerr<<vectorToString(vectorOfloglikelihoodForGivenBaBd,"\t")<<endl;
    }
	
	
    // if( (mostLikelyBaBdIdx !=  0) && 
    //     (mostLikelyBaBdIdx !=  5) && 
    //     (mostLikelyBaBdIdx != 10) && 
    //     (mostLikelyBaBdIdx != 15) ){
    //     //exit(1);
    //     cerr<<endl<<"---------------------------------------------------------------------------------"<<endl;
    //     cerr<<p<<"\tpos="<<piForGenomicWindow->at(p).posAlign<<"\t"<<piForGenomicWindow->at(p).readsVec.size()<<endl;

    //     for(unsigned int g=0;g<16;g++){
    // 	cerr<<babdIdx2Code16[g]<<"\t"<<vectorOfloglikelihoodForGivenBaBd[g]<<endl;
    // 	//cerr<<vectorToString(vectorOfloglikelihoodForGivenBaBd,"\t")<<endl;
    //     }

    //     cerr<<mostLikelyBaBd<<"\t"<<babdIdx2Code16[mostLikelyBaBdIdx]<<"\t"<<loglikelihoodForEveryBaBd_minusBest <<"\t"<< loglikelihoodForEveryBaBd<<"\t"<<( loglikelihoodForEveryBaBd_minusBest - loglikelihoodForEveryBaBd)<<endl;

    // } 

    //add to 10 genotypes
    // string babdIdx2Code10 [10] = {"AA","AC","AG","AT","CC","CG","CT","GG","GT","TT"};
    // string babdIdx2Code16 [16] = {"AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"};
    //                                 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15

    //---------------------------------------------
    // BEGIN convert to 10 genotypes
    //---------------------------------------------
    //homo
    vectorOfloglikelihoodForGivenGeno->at(0) = vectorOfloglikelihoodForGivenBaBd->at( 0);
    vectorOfloglikelihoodForGivenGeno->at(4) = vectorOfloglikelihoodForGivenBaBd->at( 5);
    vectorOfloglikelihoodForGivenGeno->at(7) = vectorOfloglikelihoodForGivenBaBd->at(10);
    vectorOfloglikelihoodForGivenGeno->at(9) = vectorOfloglikelihoodForGivenBaBd->at(15);

    //hetero
    vectorOfloglikelihoodForGivenGeno->at( 1) = oplusnatl(vectorOfloglikelihoodForGivenBaBd->at( 1), vectorOfloglikelihoodForGivenBaBd->at( 4));
    vectorOfloglikelihoodForGivenGeno->at( 2) = oplusnatl(vectorOfloglikelihoodForGivenBaBd->at( 2), vectorOfloglikelihoodForGivenBaBd->at( 8));
    vectorOfloglikelihoodForGivenGeno->at( 3) = oplusnatl(vectorOfloglikelihoodForGivenBaBd->at( 3), vectorOfloglikelihoodForGivenBaBd->at(12));
    vectorOfloglikelihoodForGivenGeno->at( 5) = oplusnatl(vectorOfloglikelihoodForGivenBaBd->at( 6), vectorOfloglikelihoodForGivenBaBd->at( 9));
    vectorOfloglikelihoodForGivenGeno->at( 6) = oplusnatl(vectorOfloglikelihoodForGivenBaBd->at( 7), vectorOfloglikelihoodForGivenBaBd->at(13));
    vectorOfloglikelihoodForGivenGeno->at( 8) = oplusnatl(vectorOfloglikelihoodForGivenBaBd->at(11), vectorOfloglikelihoodForGivenBaBd->at(14));

    long double mostLikelyGeno                     =-1.0*numeric_limits<long double>::infinity();
    int         mostLikelyGenoIdx                  =-1;
    long double mostLikelyGenoHet                  =-1.0*numeric_limits<long double>::infinity();
    int         mostLikelyGenoHetIdx               =-1;

    long double loglikelihoodForEveryGeno          =0.0;
    long double loglikelihoodForEveryGeno_minusBest=0.0;

    for(unsigned int g=0;g<10;g++){	   
	//
	pr->ll[g] = vectorOfloglikelihoodForGivenGeno->at(g); 
	if(vectorOfloglikelihoodForGivenGeno->at(g) > mostLikelyGeno){
	    mostLikelyGeno    = vectorOfloglikelihoodForGivenGeno->at(g);
	    mostLikelyGenoIdx = g;
	}

	//if het
	if( (g !=  0) && 
	    (g !=  4) && 
	    (g !=  7) && 
	    (g !=  9) ){
	    
	    if(vectorOfloglikelihoodForGivenGeno->at(g) > mostLikelyGenoHet){
		mostLikelyGenoHet    = vectorOfloglikelihoodForGivenGeno->at(g);
		mostLikelyGenoHetIdx = g;
	    }

	}

    }
    pr->mostLikelyGenoIdx    = mostLikelyGenoIdx;
    pr->mostLikelyGenoHetIdx = mostLikelyGenoHetIdx;

    
    for(int g=0;g<10;g++){
	if(g != mostLikelyGenoIdx)
	    loglikelihoodForEveryGeno_minusBest = oplusInitnatl( loglikelihoodForEveryGeno_minusBest , vectorOfloglikelihoodForGivenGeno->at(g) ); // 
	loglikelihoodForEveryGeno               = oplusInitnatl( loglikelihoodForEveryGeno           , vectorOfloglikelihoodForGivenGeno->at(g) ); // 
    }
    pr->gq = ( loglikelihoodForEveryGeno_minusBest - loglikelihoodForEveryGeno);



    //---------------------------------------------
    // END convert to 10 genotypes
    //---------------------------------------------

    //	if(true){
    if( (mostLikelyGenoIdx !=  0) && 
	(mostLikelyGenoIdx !=  4) && 
	(mostLikelyGenoIdx !=  7) && 
	(mostLikelyGenoIdx !=  9) ){

#ifdef DEBUGCOMPUTELLGENO
	cerr<<endl<<"---------------------------------------------------------------------------------"<<endl;
	cerr<<p<<"\tpos="<<piForGenomicWindow->at(p).posAlign<<"\t"<<piForGenomicWindow->at(p).readsVec.size()<<endl;
	for(unsigned int i=0;i<piForGenomicWindow->at(p).readsVec.size();i++){ //for each fragment at pos p
	    
	    cerr<<"ACGT"[piForGenomicWindow->at(p).readsVec[i].base]<<"\t"<<
		"ACGT"[piForGenomicWindow->at(p).readsVec[i].isrv?(3-piForGenomicWindow->at(p).readsVec[i].base):piForGenomicWindow->at(p).readsVec[i].base]<<"\t"
		<<"Q="<<int(piForGenomicWindow->at(p).readsVec[i].qual)<<"\t"
		<<"M="<<int(piForGenomicWindow->at(p).readsVec[i].mapq)<<"\t"
		<<"5="<<int(piForGenomicWindow->at(p).readsVec[i].pos5p)<<"\t"
		<<"L="<<int(piForGenomicWindow->at(p).readsVec[i].lengthF)<<"\t"
		<<"R="<<piForGenomicWindow->at(p).readsVec[i].isrv<<"\t"
		//		<<piForGenomicWindow->at(p).readsVec[i].name<<"\t"
		<<endl;		
	}
	for(unsigned int g=0;g<16;g++){
	    cerr<<babdIdx2Code16[g]<<"\t"<<vectorOfloglikelihoodForGivenBaBd[g]<<endl;
	}
	cerr<<endl<<"---------------------------------------------------------------------------------"<<endl;
	for(unsigned int g=0;g<10;g++){	       
	    cerr<<genoIdx2Code10[g]<<"\t"<<vectorOfloglikelihoodForGivenGeno[g]<<endl;
	}
	cerr<<endl<<"---------------------------------------------------------------------------------"<<endl;
	    
	cerr<<p<<"\tpos="<<piForGenomicWindow->at(p).posAlign<<"\tmostLikelyBaBd\t"<<mostLikelyBaBd<<"\t"<<babdIdx2Code16[mostLikelyBaBdIdx]<<"\t"<<loglikelihoodForEveryBaBd_minusBest <<"\t"<< loglikelihoodForEveryBaBd<<"\t"<<( loglikelihoodForEveryBaBd_minusBest - loglikelihoodForEveryBaBd)<<"\t"<<piForGenomicWindow->at(p).readsVec.size()<<endl;
	cerr<<p<<"\tpos="<<piForGenomicWindow->at(p).posAlign<<"\tmostLikelyGeno\t"<<mostLikelyGeno<<"\t"<<genoIdx2Code10[mostLikelyGenoIdx]<<"\t"<<loglikelihoodForEveryGeno_minusBest <<"\t"<< loglikelihoodForEveryGeno<<"\t"<<( loglikelihoodForEveryGeno_minusBest - loglikelihoodForEveryGeno)<<"\t"<<piForGenomicWindow->at(p).readsVec.size()<<endl;
	cerr<<endl<<"---------------------------------------------------------------------------------"<<endl;

#endif
    }//end for each geno


} //genotypePositions


inline void computeLLinternal(const diNucleotideProb priorGenotype,
			      const diNucleotideProb priorGenotypeD1,
			      const diNucleotideProb priorGenotypeD2,
			      const vector<babdlikelihood>      & vectorBaBdLikelihood,
			      const vector<positionInformation> * piForGenomicWindow,
			      const bool lastIteration,
			      vector<PositionResult *>    * vecPositionResults,
			      unsigned int covDist [MAXCOV+1],
			      long double & loglikelihoodForEveryPositionForEveryBaBd,
			      long double & loglikelihoodForEveryPositionForEveryBaBdD1,
			      long double & loglikelihoodForEveryPositionForEveryBaBdD2){

    for(unsigned int p=0;p<piForGenomicWindow->size();p++){//every genomic position
	

#ifdef DEBUGCOMPUTELLEACHBASE
	//if(p>10000 && p<11000){	    
	cerr<<p<<"\tpos="<<piForGenomicWindow->at(p).posAlign<<"\t"<<piForGenomicWindow->at(p).readsVec.size()<<endl;

	cerr<<"B\tQ\tMQ\t5p\tL"<<endl;

	for(unsigned int i=0;i<piForGenomicWindow->at(p).readsVec.size();i++){
	    cerr<<"ACGT"[piForGenomicWindow->at(p).readsVec[i].base]<<"\t"
		<<int(piForGenomicWindow->at(p).readsVec[i].qual)<<"\t"
		<<int(piForGenomicWindow->at(p).readsVec[i].mapq)<<"\t"
		<<int(piForGenomicWindow->at(p).readsVec[i].pos5p)<<"\t"
		<<int(piForGenomicWindow->at(p).readsVec[i].lengthF)<<"\t"
		<<piForGenomicWindow->at(p).readsVec[i].isrv<<"\t"
		//		<<piForGenomicWindow->at(p).readsVec[i].name<<"\t"
		<<endl;	    
	}
	//}    
#endif




	// typedef vector< vector<diNucleotideProb> > mpq2bsq2submatrix;
	
	
	//  // 2D: mapping quality 
	//  // 3D: base qual
	// vector< mpq2bsq2submatrix >  pos2mpq2BaseQual2SubMatrix5p;
	// vector< mpq2bsq2submatrix >  pos2mpq2BaseQual2SubMatrix3p;
	
	// // 1D: length of fragment
	// // 2D: pos fragment from the 5' end
	// vector< vector< mpq2bsq2submatrix * > > length2pos2mpq2bsq2submatrix;
	
	//for this position
	long double loglikelihoodForEveryBaBd          =0.0;
	long double loglikelihoodForEveryBaBdD1        =0.0;
	long double loglikelihoodForEveryBaBdD2        =0.0;


	vector<long double> vectorOfloglikelihoodForGivenBaBd     (16,0.0) ;
	vector<long double> vectorOfloglikelihoodForGivenGeno     (10,0.0) ;

	long double mostLikelyBaBd    = -1.0*numeric_limits<long double>::infinity();
	int         mostLikelyBaBdIdx = -1;
	int         babdIdx           =  0;





	// long double sumProbForPriors             =0.0; //the derivative of the GL should be 0, just the prior
	// long double sumDerProbForPriors          =0.0; //the derivative of the GL should be 0, just the prior
	//if the function is log( f(x) )

	long double internalFunctionD2   = 0.0; // f''(x)
	long double internalFunctionD1   = 0.0; // f'(x)
	long double internalFunction     = 0.0; // f(x)

	for(uint8_t ba=0;ba<4;ba++){//ancestral base
#ifdef DEBUGCOMPUTELLEACHBASE
	    uint8_t ba_c = 3-ba;
#endif
	    for(uint8_t bd=0;bd<4;bd++){//derived base
		    
#ifdef DEBUGCOMPUTELLEACHBASE
		uint8_t bd_c = 3-bd;
#endif

#ifdef DEBUGCOMPUTELLEACHBASE
		//if(p>10000 && p<11000)
		cerr<<endl<<"GENO:A="<<"ACGT"[ba]<<" ("<<"ACGT"[ba_c]<<") \tD="<<"ACGT"[bd]<<" ("<<"ACGT"[bd_c]<<")\tprior "<<expl(priorGenotype.p[ba][bd])<<endl;
		//cerr<<int(ba)<<"\t"<<int(ba_c)<<endl;
#endif
		
		long double loglikelihoodForGivenBaBdTimesPrior  =0.0;
		long double loglikelihoodForGivenBaBdTimesPriorD1=0.0;
		long double loglikelihoodForGivenBaBdTimesPriorD2=0.0;
		    
		// long double loglikelihoodForGivenBaBd          =0.0;

			       	      
		//  product of (\prod_{fragment} P(D|G)) times the prior P(G) for the genotype
		//loglikelihoodForGivenBaBdTimesPrior = loglikelihoodForGivenBaBd + priorGenotype.p[ba][bd]; // (\prod_{fragment} P(D|G))*P(G)
		loglikelihoodForGivenBaBdTimesPrior   = vectorBaBdLikelihood[p].gl[babdIdx] + priorGenotype.p[ba][bd]; // (\prod_{fragment} P(D|G))*P(G)
		loglikelihoodForGivenBaBdTimesPriorD1 = 0                                   + priorGenotypeD1.p[ba][bd]; // (\prod_{fragment} P(D|G))*P(G)
		loglikelihoodForGivenBaBdTimesPriorD2 = 0                                   + priorGenotypeD2.p[ba][bd]; // (\prod_{fragment} P(D|G))*P(G)
		
#ifdef DEBUGCOMPUTELLEACHBASE
		//if(p>10000 && p<11000)
		cerr<<"GENO:A="<<"ACGT"[ba]<<"\tD="<<"ACGT"[bd]<<"\tprior "<<expl(priorGenotype.p[ba][bd])<<"\tllForBaBD "<<loglikelihoodForGivenBaBd<<"\tllForBaBD*Prior "<<loglikelihoodForGivenBaBdTimesPrior<<" ("<<expl(loglikelihoodForGivenBaBdTimesPrior) <<")"<<" priorGenotypeProb "<<priorGenotypeProb.p[ba][bd]<<" D "<<priorGenotypeProbD.p[ba][bd]<<"\tllForEveryBaBD "<<loglikelihoodForEveryBaBd<<"\tgenoLike "<<mostLikelyBaBd<<"\tmostLikeG "<<mostLikelyBaBdIdx<<endl;
#endif


		//adding probabilities for each genotype
		loglikelihoodForEveryBaBd  = oplusInitnatl( loglikelihoodForEveryBaBd , loglikelihoodForGivenBaBdTimesPrior); // \sum_{genotype} (\prod_{fragment} P(D|G))*P(G)
		    
		    
		//(d^2)/(d^2h) (exp(f(h))) =  e^(f(h))f'(h) * f'(h) + e^(f(h)) f''(h) 
		internalFunctionD2  += 
		    expl( loglikelihoodForGivenBaBdTimesPrior )*loglikelihoodForGivenBaBdTimesPriorD1*loglikelihoodForGivenBaBdTimesPriorD1 
		    + 
		    expl( loglikelihoodForGivenBaBdTimesPrior )*loglikelihoodForGivenBaBdTimesPriorD2;  

		//(d)/(dh)(exp(f(h))) =  e^(f(h)) *f'(h)
		internalFunctionD1  += expl( loglikelihoodForGivenBaBdTimesPrior )*loglikelihoodForGivenBaBdTimesPriorD1;  

		internalFunction    += expl( loglikelihoodForGivenBaBdTimesPrior );

		    
		// sumProbForPriors          += priorGenotypeProb.p[ba][bd];
		// sumDerProbForPriors       += priorGenotypeProbD.p[ba][bd];
		    
		vectorOfloglikelihoodForGivenBaBd[babdIdx] = loglikelihoodForGivenBaBdTimesPrior ;

		if(loglikelihoodForGivenBaBdTimesPrior>mostLikelyBaBd){ //finding most likely Ba Bd pair
		    mostLikelyBaBd    = loglikelihoodForGivenBaBdTimesPrior;
		    mostLikelyBaBdIdx = babdIdx;
		}
		
#ifdef DEBUGCOMPUTELL
		//if(p>10000 && p<11000)
		//cerr<<endl<<"---------------------------------------------------------------------------------"<<endl;
#endif
		babdIdx++;
	    }//END for each derived base
	}//END for each ancestral base

	//The derivative for the log( exp(p1) + exp(p2) +...exp(p10) ) is
	//
	//             der( exp(p1) + exp(p2) +...exp(p10) )
	//             -----------------------------------
	//                ( exp(p1) + exp(p2) +...exp(p10) )
	loglikelihoodForEveryBaBdD1 +=  ( internalFunctionD1 / internalFunction );  
	//loglikelihoodForEveryBaBdD1 +=  ( sumNumerator / sumDenominator);  
	//loglikelihoodForEveryBaBdD1 +=  ( sumDerProbForPriors/sumProbForPriors);  
	    
	//If the function is log( f(x) ) where f(x) (called internalFunction) is the likelihood function then
	//  d/dx log( f(x) ) = f'(x) / f(x)
	//  and 
	//  d^2/d^2x log( f(x) ) = [f''(x) f(x) + f'(x)^2 ] / [ f(x) ]^2
	loglikelihoodForEveryBaBdD2 +=  ( internalFunctionD2 * internalFunction - powl(internalFunctionD1,2) ) / powl( internalFunction , 2 );

#ifdef DEBUGCOMPUTELL	
	//if(p>10000 && p<11000)
	//cerr<<endl<<"---------------------------------------------------------------------------------"<<endl;
#endif


	//////////////////////////////////////////////////////
	// BEGIN GENOTYPING
	//////////////////////////////////////////////////////
	if(lastIteration){//do it at the last iteration		
	    if(outputgenol){
		PositionResult * pr=new PositionResult();


		pr->refB   =     piForGenomicWindow->at(p).refBase;
		pr->pos    =     piForGenomicWindow->at(p).posAlign;
		pr->avgMQ  =     piForGenomicWindow->at(p).avgMQ;

		for(int n=0;n<4;n++)
		    pr->baseC[n] = piForGenomicWindow->at(p).baseC[n];

		pr->dp     = int(piForGenomicWindow->at(p).readsVec.size());
	    
		genotypePositions( mostLikelyBaBdIdx                 ,
				   &vectorOfloglikelihoodForGivenBaBd,
				   &vectorOfloglikelihoodForGivenGeno,
				   pr);
		vecPositionResults->push_back(pr);
	    }
	    //vecPositionResults


		
	    //count cov dist
	    covDist[ piForGenomicWindow->at(p).readsVec.size() ]++;
		
	}
	////////////////////////////////////////////////////
	// END GENOTYPING
	////////////////////////////////////////////////////
	    





	//product for each genomic position
                                                        
	// loglikelihoodForEveryPositionForEveryBaBd   += loglikelihoodForEveryBaBd;   // \prod_{site} \sum_{genotype} (\prod_{fragment} P(D|G))*P(G)
	// loglikelihoodForEveryPositionForEveryBaBdD1 += loglikelihoodForEveryBaBdD1; // \prod_{site} \sum_{genotype} (\prod_{fragment} P(D|G))*P(G)
	// loglikelihoodForEveryPositionForEveryBaBdD2 += loglikelihoodForEveryBaBdD2; // \prod_{site} \sum_{genotype} (\prod_{fragment} P(D|G))*P(G)

	// \prod_{site} \sum_{genotype} (\prod_{fragment} P(D|G))*P(G)
	//#if defined(ONLYUSECOV) || defined(ONLYUSECOVMIN) || defined(ONLYUSECOVMAX)
	// if(piForGenomicWindow->at(p).readsVec.size() < ONLYUSECOVMIN){
	//     continue;
	// }

	// if(piForGenomicWindow->at(p).readsVec.size() > ONLYUSECOVMAX){
	//     continue;
	// }

	if(covcutoffs){//do not rely on coverage prior, use hard cutoffs instead
	    
	    if( int(piForGenomicWindow->at(p).readsVec.size()) < covcutoffsmin){
		continue;
	    }

	    if( int(piForGenomicWindow->at(p).readsVec.size()) > covcutoffsmax){
		continue;
	    }

	    loglikelihoodForEveryPositionForEveryBaBd   +=                                                                 loglikelihoodForEveryBaBd;
	    loglikelihoodForEveryPositionForEveryBaBdD1 +=                                                                 loglikelihoodForEveryBaBdD1;
	    loglikelihoodForEveryPositionForEveryBaBdD2 +=                                                                 loglikelihoodForEveryBaBdD2; 
	}else{// do not use hard cutoffs, use coverage prior
	    // cerr<<loglikelihoodForEveryPositionForEveryBaBd<<"\t"<<cov2ProbSite->at( piForGenomicWindow->at(p).readsVec.size() ) <<"\t"<< loglikelihoodForEveryBaBd<<endl;
	    // cerr<<loglikelihoodForEveryPositionForEveryBaBdD1<<"\t"<<cov2ProbSite->at( piForGenomicWindow->at(p).readsVec.size() ) <<"\t"<< loglikelihoodForEveryBaBdD1<<endl;
	    // cerr<<loglikelihoodForEveryPositionForEveryBaBdD2<<"\t"<<cov2ProbSite->at( piForGenomicWindow->at(p).readsVec.size() ) <<"\t"<< loglikelihoodForEveryBaBdD2<<endl;

	    loglikelihoodForEveryPositionForEveryBaBd   += cov2ProbSite->at( piForGenomicWindow->at(p).readsVec.size() ) * loglikelihoodForEveryBaBd;
	    loglikelihoodForEveryPositionForEveryBaBdD1 += cov2ProbSite->at( piForGenomicWindow->at(p).readsVec.size() ) * loglikelihoodForEveryBaBdD1; 
	    loglikelihoodForEveryPositionForEveryBaBdD2 += cov2ProbSite->at( piForGenomicWindow->at(p).readsVec.size() ) * loglikelihoodForEveryBaBdD2; 
	}
	
	//#endif
	// loglikelihoodForEveryPositionForEveryBaBd   +=
	//     //#ifdef CORRECTCOV
	//     cov2ProbSite->at( piForGenomicWindow->at(p).readsVec.size() ) *
	//     //#endif
	//     loglikelihoodForEveryBaBd;
	// #ifdef DEBUGCOV

	// 	    cerr<<loglikelihoodForEveryPositionForEveryBaBd<<" "<<cov2ProbSite->at( piForGenomicWindow->at(p).readsVec.size() )<<" "<<loglikelihoodForEveryBaBd<<" l="<<piForGenomicWindow->at(p).readsVec.size()<<"#"<<endl;
	// 	    if( isnan( cov2ProbSite->at( piForGenomicWindow->at(p).readsVec.size() ) ) ){
	// 		cerr<<piForGenomicWindow->at(p).readsVec.size()<<endl;
	// 		exit(1);
	// 	    }
	// #endif
	    
	// 	loglikelihoodForEveryPositionForEveryBaBdD1 +=
	// #ifdef CORRECTCOV
	// 	    cov2ProbSite->at( piForGenomicWindow->at(p).readsVec.size() ) *
	// #endif
	// 	    loglikelihoodForEveryBaBdD1; 
	// 	loglikelihoodForEveryPositionForEveryBaBdD2 +=
	// #ifdef CORRECTCOV
	// 	    cov2ProbSite->at( piForGenomicWindow->at(p).readsVec.size() ) *
	// #endif
	// 	    loglikelihoodForEveryBaBdD2; 

    }//END for each genomic position

}





inline hResults computeLL(const vector<positionInformation> * piForGenomicWindow,
			  vector<PositionResult *>    * vecPositionResults,
			  const int threadID){


    int rc;

    if(verbose){
	rc = pthread_mutex_lock(&mutexCERR);
	checkResults("pthread_mutex_lock()\n", rc);
	
	cerr<<"Thread#"<<threadID<<" starting pre-computations size of window: "<<thousandSeparator(piForGenomicWindow->size())<<endl;
	
	rc = pthread_mutex_unlock(&mutexCERR);
	checkResults("pthread_mutex_unlock()\n", rc);
    }
    
    hResults hresToReturn;
    hresToReturn.sites = piForGenomicWindow->size();

    if(piForGenomicWindow->size() == 0){
	hresToReturn.h            = 0;
	hresToReturn.hLow         = 0;
	hresToReturn.hHigh        = 0;
	hresToReturn.hasConverged = false;
	
	return hresToReturn;
    }

    
    /////////////////////////////////////////
    //BEGIN pre-computing babdlikelihood
    /////////////////////////////////////////
    vector<babdlikelihood> vectorBaBdLikelihood;
    
    preComputeBaBdLikelihood(piForGenomicWindow,&vectorBaBdLikelihood,piForGenomicWindow->size());
    
    /////////////////////////////////////////
    //END pre-computing babdlikelihood
    /////////////////////////////////////////

    if(verbose){
	rc = pthread_mutex_lock(&mutexCERR);
	checkResults("pthread_mutex_lock()\n", rc);

	cerr<<"Thread#"<<threadID<<" done pre-computing computing "<<endl;

	rc = pthread_mutex_unlock(&mutexCERR);
	checkResults("pthread_mutex_unlock()\n", rc);
    }
    //max iterations


    int sitesPer1M=randomInt(minSegSitesPer1M,maxSegSitesPer1M);// between minSegSitesPer1M and maxSegSitesPer1M  and 1000 sites per million
    // long double theta     = double(sitesPer1M)/double(1000000); //theta
    // long double theta_t_1 = theta;
    long double hp;
    long double h=double(sitesPer1M)/double(1000000);

    long double hPosD1 = double(minSegSitesPer1M)/double(1000000);//last h with lowest positive D1
    //long double hPosD1LLD1=-1;//last h with lowest positive D1
    bool        hPosD1init=false;//last h with lowest positive D1
    
    long double hNegD1 = double(maxSegSitesPer1M)/double(1000000);//last h with lowest positive D1 
    //    long double hNegD1LLD1=-1;//last h with lowest positive D1 
    bool        hNegD1init=false;//last h with lowest positive D1 

    if(suggestH){
	h=suggestHval;	
    }


    long double z=h;

    long double alpha_ = alpha;
    long double beta_  = beta;

    //long double hLastIteration=LDBL_MAX;
    //long double hPrecision=0.0000001;
    long double errb;
    long double loglikelihoodForEveryPositionForEveryBaBdP=LDBL_MAX;

    bool hasConverged    = false;
    bool lastIteration   = false;
    bool bottomIteration = false;
    long double deltaBottom=1.0e-11;

    // long double mu ;

    // //for(long double h=0.000001;h<0.001000;h+=0.0000250){
    // //while(1){
    // long double lambda   = 0.0000000001;
    //long double lambdaHW  = 0.0000000001;

    int iterationsMax     = 200;
    long double loglikelihoodForEveryPositionForEveryBaBdD1ForConvergence = 100.0;//target to assess convergence
    int iterationsGrad    = 1;
    int iterationWithMinHRate    =  0;
    int maxiterationWithMinHRate = 10;
    //check jumping around min
    int iterationWithPositiveD1  = 0;//iterations with positive first derivative
    int iterationWithNegativeD1  = 0;//iterations with negative first derivative
    
    double minHRateNonROH = double(minSegSitesPer1M)/double(1000000);

    unsigned int covDist [MAXCOV+1];
    for(int i=0;i<=MAXCOV;i++){
	covDist[i] = 0.0;
    }
    
    while( iterationsGrad<iterationsMax){

	//mu = 1.0 - double(3.0)/(5+double(iterationsGrad));


	//h = (1.0+mu) * theta -  mu*theta_t_1;
	//for(long double h=0.000001;h<0.001000;h+=0.000100){
	//long double h=0.00000001;
	//long double h=0.000824;
	//long double h=0.000735;
	//long double h=0.010000;
	diNucleotideProb priorGenotype;
	// diNucleotideProb priorGenotypeProb;
	// diNucleotideProb priorGenotypeProbD;

	diNucleotideProb priorGenotypeD1;
	diNucleotideProb priorGenotypeD2;

	/////////////////////////////////////////
	//compute prior genotype matrix given h
	/////////////////////////////////////////

	computePriorMatrix(h,
			   &priorGenotype,
			   // &priorGenotypeProb,
			   // &priorGenotypeProbD,
			   &priorGenotypeD1,
			   &priorGenotypeD2);

	
	// cerr<<setprecision(20)<<"TEST3 "<<"\tL2P="<<length2pos2mpq2bsq2submatrix[113][12]->at(37)[38].p[0][3]<<"\t"<<length2pos2mpq2bsq2submatrix[113][12]->at(37)[18].p[0][3]<<"\t"<<expl(length2pos2mpq2bsq2submatrix[113][12]->at(37)[38].p[0][3])<<"\t"<<expl(length2pos2mpq2bsq2submatrix[113][12]->at(37)[18].p[0][3])<<"\t"<<logl(0.5*length2pos2mpq2bsq2submatrix[113][12]->at(37)[38].p[0][3])<<"\t"<<logl(0.5*length2pos2mpq2bsq2submatrix[113][12]->at(37)[18].p[0][3])<<endl ;

#ifdef DEBUGCOMPUTELL
	cerr<<"h="<<h<<endl;
	cerr<<"\t";
	for(int ba=0;ba<4;ba++)
	    cerr<<"ACGT"[ba]<<"\t";
	cerr<<endl;
	long double sumProb_=0.0;
	for(int ba=0;ba<4;ba++){
	    cerr<<"ACGT"[ba]<<"\t";
	    for(int bd=0;bd<4;bd++){
		cerr<<expl(priorGenotype.p[ba][bd])<<"\t";
		sumProb_+=expl(priorGenotype.p[ba][bd]);
	    }
	    cerr<<endl;
	}
	cerr<<endl;
	cerr<<"sum prior for each geno = "<<sumProb_<<endl;	   
#endif


	//for every position
	long double loglikelihoodForEveryPositionForEveryBaBd          =0.0;
	long double loglikelihoodForEveryPositionForEveryBaBdD1        =0.0;
	long double loglikelihoodForEveryPositionForEveryBaBdD2        =0.0;


	computeLLinternal(priorGenotype,
			  priorGenotypeD1,
			  priorGenotypeD2,
			  vectorBaBdLikelihood,
			  piForGenomicWindow,
			  lastIteration,
			  vecPositionResults,
			  covDist ,
			  loglikelihoodForEveryPositionForEveryBaBd,
			  loglikelihoodForEveryPositionForEveryBaBdD1,
			  loglikelihoodForEveryPositionForEveryBaBdD2	);
        errb = 1.96/sqrt(-1.0*loglikelihoodForEveryPositionForEveryBaBdD2);
	
	
	//cout<<setprecision(14)<<h<<"\t"<<loglikelihoodForEveryPositionForEveryBaBd<<"\t"<<loglikelihoodForEveryPositionForEveryBaBdD1<<"\t"<<loglikelihoodForEveryPositionForEveryBaBdD2<<"\t"<<errb<<"\t"<<(h-errb)<<"\t"<<(h+errb)<<"\t"<<hnew<<endl;

	// long double hnew = h + lambda*loglikelihoodForEveryPositionForEveryBaBdD1;
	// h=hnew;
	// theta_t_1 = theta;
	// theta  = h + lambda*loglikelihoodForEveryPositionForEveryBaBdD1;
	
	z = beta_*z + loglikelihoodForEveryPositionForEveryBaBdD1;
	long double hnew = h + alpha_*z;


	//Begin avoiding overshooting
	if(loglikelihoodForEveryPositionForEveryBaBdD1>0){//in positive derivative, need to increase h

	    if(hPosD1init==false){//init the positive boundary

		hPosD1      =h;
		//hPosD1LLD1  =loglikelihoodForEveryPositionForEveryBaBdD1;
		hPosD1init  =true;
		
#ifdef DEBUGBOXINGALG
		cerr<<"Thread#"<<threadID<<" found a + lower bound for h "<<h<<"  range  "<<hPosD1<<" - "<<hNegD1<<endl;
#endif
		
	    }else{//was initialized, should be better, otherwise it would be blocked
	    
		hPosD1      =h;
		//hPosD1LLD1  =loglikelihoodForEveryPositionForEveryBaBdD1;
#ifdef DEBUGBOXINGALG
		cerr<<"Thread#"<<threadID<<" found a + lower bound for h "<<h<<"  range  "<<hPosD1<<" - "<<hNegD1<<endl;
#endif
			    	  
	    }
	    
	    if(hPosD1init && hNegD1init){//check if hnew will overshoot
		// previous      new one
		if(hNegD1<hnew || hnew<hPosD1  ){//wrong direction		

		    long double hnewCorrect = (hPosD1 + hNegD1)/2.0;
#ifdef DEBUGBOXINGALG
		    cerr<<"Thread#"<<threadID<<" overshot in - with new h "<<hnew<<" not in  "<<hPosD1<<" - "<<hNegD1<<" moving instead to "<<hnewCorrect<<endl;
#endif
		    
		    hnew = hnewCorrect;
		    
		}else{
		    //ok
		}		    
	    }
	    
	}else{ //	if(loglikelihoodForEveryPositionForEveryBaBdD1<0){ //in negatvive derivative, need to decrease h

	    if(hNegD1init==false){
	    
		hNegD1      =h;
		///hNegD1LLD1  =loglikelihoodForEveryPositionForEveryBaBdD1;
		hNegD1init  =true;
		
#ifdef DEBUGBOXINGALG
		cerr<<"Thread#"<<threadID<<" found a - lower bound for h "<<h<<"  range  "<<hPosD1<<" - "<<hNegD1<<endl;
#endif    
			    	  	    
	    }else{//was initialized, should be better, otherwise it would be blocked

		hNegD1      =h;
		//hNegD1LLD1  =loglikelihoodForEveryPositionForEveryBaBdD1;
			
#ifdef DEBUGBOXINGALG
		cerr<<"Thread#"<<threadID<<" found a - lower bound for h "<<h<<"  range  "<<hPosD1<<" - "<<hNegD1<<endl;
#endif
		
	   	    
	    }

	    if(hPosD1init && hNegD1init){//check if hnew will overshoot
		// previous      new one
		if(hNegD1<hnew || hnew<hPosD1  ){//wrong direction

		    long double hnewCorrect = (hPosD1 + hNegD1)/2.0;
#ifdef DEBUGBOXINGALG
		    cerr<<"Thread#"<<threadID<<" overshot in + with new h "<<hnew<<" not in  "<<hPosD1<<" - "<<hNegD1<<" moving instead to "<<hnewCorrect<<endl;
#endif
		    
		    hnew = hnewCorrect;
		    
		}else{
		    //ok
		}		    
	    }
	    
	}
	

		    
	
	if(loglikelihoodForEveryPositionForEveryBaBdD1 > 0){
	    iterationWithPositiveD1++;
	}else{
	    iterationWithNegativeD1++;
	}



	
	
	if(verbose){
	    rc = pthread_mutex_lock(&mutexCERR);
	    checkResults("pthread_mutex_lock()\n", rc);
	    
	    cerr<<"Thread#"<<threadID<<setprecision(14)<<"\t"<<h<<"\t"<<loglikelihoodForEveryPositionForEveryBaBd<<"\t"<<loglikelihoodForEveryPositionForEveryBaBdD1<<"\t"<<loglikelihoodForEveryPositionForEveryBaBdD2<<"\te="<<errb<<"\t"<<(h-errb)<<"\t"<<(h+errb)<<"\thnew="<<hnew<<"\tz="<<z<<"\t"<<loglikelihoodForEveryPositionForEveryBaBdP<<"\t"<<hp<<"\t"<<alpha_<<"\t"<<beta_<<endl;
	    
	
	    rc = pthread_mutex_unlock(&mutexCERR);
	    checkResults("pthread_mutex_unlock()\n", rc);
	}
	
	if(lastIteration){//was set at the previous loop
	    break;
	}
	
	//if(abs(hnew-h)<hPrecision){
	if(fabsl(loglikelihoodForEveryPositionForEveryBaBdD1) < loglikelihoodForEveryPositionForEveryBaBdD1ForConvergence){

	    if(verbose){
		rc = pthread_mutex_lock(&mutexCERR);
		checkResults("pthread_mutex_lock()\n", rc);
		
		cerr<<"Thread#"<<threadID<<" converged "<<endl;
		
		rc = pthread_mutex_unlock(&mutexCERR);
		checkResults("pthread_mutex_unlock()\n", rc);
	    }
	    
	    lastIteration=true;
	    hasConverged =true;
	    //break;
	}

	//setting h to hnew
	if(loglikelihoodForEveryPositionForEveryBaBdP == LDBL_MAX){//if previous hadn't been set
	    loglikelihoodForEveryPositionForEveryBaBdP = loglikelihoodForEveryPositionForEveryBaBd;
	    hp = h;
	    h  = hnew;
	}else{
	
	    // if( loglikelihoodForEveryPositionForEveryBaBdP > loglikelihoodForEveryPositionForEveryBaBd){
	    // 	//reject move, this is dumb because we end up re-calculating the same iteration as before

	    // 	h =  hp;
	    // 	//is 0.5 for large values of the first derivative
	    // 	//is 0.7 at 1000
	    // 	long double factorAlpha = 0.5+1/(expl(fabsl(loglikelihoodForEveryPositionForEveryBaBdD1)/1000)+1);
	    // 	cout<<"reject move , reverting h to "<<hp<<"\t"<<alpha<<"\t"<<factorAlpha<<endl;
	    // 	alpha =  alpha * factorAlpha;
	    // 	continue;
	    // }else{
	    loglikelihoodForEveryPositionForEveryBaBdP = loglikelihoodForEveryPositionForEveryBaBd;
	    hp = h;
	    h  = hnew;
	    //}
	}


	//adjusting alpha/beta
	if(h<0){
	    if(verbose){
		rc = pthread_mutex_lock(&mutexCERR);
		checkResults("pthread_mutex_lock()\n", rc);
		
		cerr<<"Thread#"<<threadID<<setprecision(14)<<" new proposed h is less than 0 "<<h<<"\talpha:"<<alpha_<<"\tbeta:"<<beta_<<endl;
				
		rc = pthread_mutex_unlock(&mutexCERR);
		checkResults("pthread_mutex_unlock()\n", rc);
	    }


	    //h     = 1.0e-11;
	    h = hp;//resetting to previous h
	    alpha_ = alpha_*0.5;
	    beta_  = beta_ *0.9;

	    //resetting 
	    iterationWithPositiveD1=0;
	    iterationWithNegativeD1=0;

	}else{

	    if(h>0.1){
		h     = 0.1-1.0e-11;
		alpha_ = alpha_*0.5;
		beta_  = beta_ *0.9;

		//resetting 
		iterationWithPositiveD1=0;
		iterationWithNegativeD1=0;
	    }else{

		//if we spent too much on either side, reduce alpha
		if( ((iterationWithPositiveD1+iterationWithNegativeD1 )>12 ) &&
		    ((double(MIN2(iterationWithPositiveD1,iterationWithNegativeD1))/double(iterationWithPositiveD1+iterationWithNegativeD1))>0.25) ){ //if more than 25% of iterations are on the other side, adjust alpha

		    if(verbose){
			rc = pthread_mutex_lock(&mutexCERR);
			checkResults("pthread_mutex_lock()\n", rc);
			
			cerr<<"Thread#"<<threadID<<setprecision(14)<<"\t found signs of oscillation +:"<<iterationWithPositiveD1<<"  -:"<<iterationWithNegativeD1<<" "<<h<<"\talpha:"<<alpha_<<"\tbeta:"<<beta_<<endl;
			
			rc = pthread_mutex_unlock(&mutexCERR);
			checkResults("pthread_mutex_unlock()\n", rc);
		    }

		    alpha_ = alpha_*0.5;
		    beta_  = beta_ *0.9;
		    //resetting 
		    iterationWithPositiveD1=0;
		    iterationWithNegativeD1=0;
		}else{

		    //if spend too much time on one side, increase alpha

		    if( ((iterationWithPositiveD1+iterationWithNegativeD1 )>12 ) &&
			((double(MIN2(iterationWithPositiveD1,iterationWithNegativeD1))/double(iterationWithPositiveD1+iterationWithNegativeD1))<0.05) ){ //if more than 25% of iterations are on the other side, adjust alpha
			
			if(verbose){
			    rc = pthread_mutex_lock(&mutexCERR);
			    checkResults("pthread_mutex_lock()\n", rc);
			    
			    cerr<<"Thread#"<<threadID<<setprecision(14)<<"\t spent too much time on one side +:"<<iterationWithPositiveD1<<"  -:"<<iterationWithNegativeD1<<" "<<h<<"\talpha:"<<alpha_<<"\tbeta:"<<beta_<<endl;
			    
			    rc = pthread_mutex_unlock(&mutexCERR);
			    checkResults("pthread_mutex_unlock()\n", rc);
			}
			
			alpha_ = alpha_*2.0;
			//beta_  = beta_ *0.9;
			//resetting 
			iterationWithPositiveD1=0;
			iterationWithNegativeD1=0;
		    }

		    
		}

	    }
	}
	
	if( h < minHRateNonROH){
	    iterationWithMinHRate++;
	}else{
	    iterationWithMinHRate=0;
	}

	if(iterationWithMinHRate == maxiterationWithMinHRate ){
	    hasConverged =true;
	    lastIteration=true;

	    if(verbose){
		rc = pthread_mutex_lock(&mutexCERR);
		checkResults("pthread_mutex_lock()\n", rc);
		
		cerr<<"Thread#"<<threadID<<setprecision(14)<<" has reached the bottom limit ("<<minHRateNonROH<<") for "<<iterationWithMinHRate<<" iterations now"<<endl;		
		
		rc = pthread_mutex_unlock(&mutexCERR);
		checkResults("pthread_mutex_unlock()\n", rc);
	    }

	    //has not converged in terms of the first derivative, need to use another mean to asses the errb
	    if(fabsl(loglikelihoodForEveryPositionForEveryBaBdD1) > loglikelihoodForEveryPositionForEveryBaBdD1ForConvergence){

		if(verbose){
		    rc = pthread_mutex_lock(&mutexCERR);
		    checkResults("pthread_mutex_lock()\n", rc);
		    
		    cerr<<"Thread#"<<threadID<<setprecision(14)<<" reached bottom but first derivative did not reach close to 0"<<endl;
		    
		    rc = pthread_mutex_unlock(&mutexCERR);
		    checkResults("pthread_mutex_unlock()\n", rc);
		}
		
		long double h_        ;
		
		for(int iterationmininter=0;iterationmininter<200;iterationmininter++){
		    h_ = h+deltaBottom;
		    
		    diNucleotideProb priorGenotype;
		    diNucleotideProb priorGenotypeD1;
		    diNucleotideProb priorGenotypeD2;

		
		    computePriorMatrix(h_,
				       &priorGenotype,
				       // &priorGenotypeProb,
				       // &priorGenotypeProbD,
				       &priorGenotypeD1,
				       &priorGenotypeD2);
		    
	

		
		    long double loglikelihoodForEveryPositionForEveryBaBd_          =0.0;
		    long double loglikelihoodForEveryPositionForEveryBaBdD1_        =0.0;
		    long double loglikelihoodForEveryPositionForEveryBaBdD2_        =0.0;
		    

		    computeLLinternal(
			priorGenotype,
			priorGenotypeD1,
			priorGenotypeD2,
			vectorBaBdLikelihood,
			piForGenomicWindow,
			false,
			vecPositionResults,
			covDist ,
			loglikelihoodForEveryPositionForEveryBaBd_,
			loglikelihoodForEveryPositionForEveryBaBdD1_,
			loglikelihoodForEveryPositionForEveryBaBdD2_	);

		    if(verbose){
			rc = pthread_mutex_lock(&mutexCERR);
			checkResults("pthread_mutex_lock()\n", rc);
		    
			cerr<<"Thread#"<<threadID<<setprecision(14)<<" bottom "<<h_<<" "<<deltaBottom<<" "<<loglikelihoodForEveryPositionForEveryBaBd_<<" "<<(2.0*fabsl(loglikelihoodForEveryPositionForEveryBaBd_-loglikelihoodForEveryPositionForEveryBaBd))<<endl;
		    
			rc = pthread_mutex_unlock(&mutexCERR);
			checkResults("pthread_mutex_unlock()\n", rc);
		    }
		    //cerr<<h_<<" "<<deltaBottom<<" "<<loglikelihoodForEveryPositionForEveryBaBd_<<" "<<(2.0*fabsl(loglikelihoodForEveryPositionForEveryBaBd_-loglikelihoodForEveryPositionForEveryBaBd))<<endl;
		    if( (2.0*fabsl(loglikelihoodForEveryPositionForEveryBaBd_-loglikelihoodForEveryPositionForEveryBaBd)) > 3.84){
			bottomIteration=true;
			break;
		    }else{
			deltaBottom = deltaBottom*2.0;
		    }
		    
		}		
		//cout<<"h new h_ "<<h_<<endl;

	    }//end if first derivative has worked

	    

	}
	
	if(iterationsGrad == (iterationsMax-2)){
	    if(verbose){
		rc = pthread_mutex_lock(&mutexCERR);
		checkResults("pthread_mutex_lock()\n", rc);


		cerr<<"Thread#"<<threadID<<" did not converge"<<endl;

	    
		rc = pthread_mutex_unlock(&mutexCERR);
		checkResults("pthread_mutex_unlock()\n", rc);
	    }
	    hasConverged =false;
	    lastIteration=true;
	    //break;
	}

	//cout<<setprecision(14)<<h<<"\t"<<loglikelihoodForEveryPositionForEveryBaBd<<"\t"<<loglikelihoodForEveryPositionForEveryBaBdD1<<"\t"<<loglikelihoodForEveryPositionForEveryBaBdD2<<"\t"<<errb<<"\t"<<(h-errb)<<"\t"<<(h+errb)<<"\t"<<theta<<"\t"<<theta_t_1<<"\t"<<mu<<"\t"<<hLastIteration<<endl;

	iterationsGrad++;	
    }//end iteration

    if(bottomIteration){
	errb=deltaBottom;
    }
    
    unsigned int covDistSum=0.0;
    for(int i=0;i<=MAXCOV;i++){	
	covDistSum += covDist[ i ];
    }
    
    long double correctionFactor=0.0;
    
    for(int i=0;i<=MAXCOV;i++){
	correctionFactor += overestimateFactor[i] * (double(covDist[ i ])/double(covDistSum));
	//cout<<"covDist["<<i<<"] "<<covDist[ i ]<<" overestimateFactor["<<i<<"] "<<overestimateFactor[ i ]<<" "<<correctionFactor<<endl;
    }
    // cout<<setprecision(14)<<h<<"\t"<<correctionFactor<<"\t";

    //h = h*correctionFactor;
    //cout<<setprecision(14)<<h<<"\t"<<(h-errb)<<"\t"<<(h+errb)<<endl;

    hresToReturn.h            = (h)     *correctionFactor;
    hresToReturn.hLow         = (h-errb)*correctionFactor;
    hresToReturn.hHigh        = (h+errb)*correctionFactor;
    hresToReturn.errb         = errb    *correctionFactor;
    // if(tvonly){
    //   hresToReturn.h            =       hresToReturn.h     * (1.0+TStoTVratio); 
    //   hresToReturn.hLow         =       hresToReturn.hLow  * (1.0+TStoTVratio); 
    //   hresToReturn.hHigh        =       hresToReturn.hHigh * (1.0+TStoTVratio); 
    //   hresToReturn.errb         =       hresToReturn.errb  * (1.0+TStoTVratio);       
    // }
    
    hresToReturn.hasConverged = hasConverged;


    return hresToReturn;
    
}//end computeLL()





// class heteroComputerVisitor : public PileupVisitor {
  
// public:
//     heteroComputerVisitor(const RefVector& references, 
// 			  const int refID,
// 			  const unsigned int leftCoord,
// 			  const unsigned int rightCoord,
// 			  //vector<PositionResult *> * dataToWriteOut,
// 			  const int threadID,
// 			  vector<positionInformation> * piForGenomicWindow,
// 			  Fasta * fastaReference)
// 	: PileupVisitor()
// 	, m_references(references)
// 	, m_refID(refID)
// 	, m_leftCoord(leftCoord)
// 	, m_rightCoord(rightCoord)
// 	  //, m_dataToWriteOut( dataToWriteOut)
// 	, m_threadID( threadID )
// 	, m_numberOfSites( 0 )
// 	, totalBases(0)
// 	, totalSites(0)   
// 	, m_piForGenomicWindow( piForGenomicWindow )
// 	, m_fastaReference(fastaReference)
//     { 
// 	//cerr<<"heteroComputerVisitor constructor"<<endl;
//     }
//     ~heteroComputerVisitor(void) { }
  
//     // PileupVisitor interface implementation

    
//     void Visit(const PileupPosition& pileupData) {   
	


// 	if(pileupData.Position < int(m_leftCoord)   || 
// 	   pileupData.Position > int(m_rightCoord) ){
// 	    return ;
// 	}

// 	if( (m_numberOfSites%50000)==0 &&
// 	    m_numberOfSites != 0 ){
// 	    int rc = pthread_mutex_lock(&mutexCERR);
// 	    checkResults("pthread_mutex_lock()\n", rc);

// 	    if(verbose)
// 		cerr<<"Thread#"<<m_threadID<<" reading: "<<m_references[m_refID].RefName<<":"<<pileupData.Position<<" valid sites:\t"<<thousandSeparator(totalSites)<<endl;

// 	    rc = pthread_mutex_unlock(&mutexCERR);
// 	    checkResults("pthread_mutex_unlock()\n", rc);

// 	}

// 	m_numberOfSites++;


// 	// int                 totalBases=0 ;
// 	//int                 counterB  [4];
// 	//long double         llBaseDeam[4];


// 	// vector<int>              obsBase      ;
// 	// vector<int>              obsQual      ;
// 	// vector<int>              mmQual       ; //mismapping probability
// 	// vector<bool>             isRevVec;
// 	//vector<singleRead> singleReadToAdd;



// 	unsigned int                posAlign = pileupData.Position+1;

// 	char referenceBase = 'N';
// 	if ( !m_fastaReference->GetBase(pileupData.RefId, posAlign-1, referenceBase ) ) {
// 	    cerr << "rohan convert ERROR: pileup conversion - could not read reference base from FASTA file" << endl;
// 	    return;
// 	}
// 	referenceBase = toupper(referenceBase);
	
// 	if(!isResolvedDNA(referenceBase)){ //avoid Ns
// 	    return; 
// 	}

// 	positionInformation piToAdd;

// 	piToAdd.skipPosition = false;

// 	piToAdd.posAlign                     = posAlign;
// 	//piToAdd.refID                        = posAlign;
// 	piToAdd.refBase                      = referenceBase;
// 	//
// 	double probMM=0;
// 	int    basesRetained=0;
// 	bool foundSites=false;

// 	for(int n=0;n<4;n++){
// 	    piToAdd.baseC[n]=0;
// 	}

// 	for(unsigned int i=0;i<pileupData.PileupAlignments.size();i++){
// 	    //cerr<<i<<" "<<pileupData.PileupAlignments[i].Alignment.Name<<endl;
	    
// 	    if( pileupData.PileupAlignments[i].IsCurrentDeletion   ||
// 	    	pileupData.PileupAlignments[i].IsNextInsertion     ||
// 	    	pileupData.PileupAlignments[i].IsNextDeletion      ||
// 		(pileupData.PileupAlignments[i].DeletionLength>0)  ||
// 		(pileupData.PileupAlignments[i].InsertionLength>0) ){		
// 		//includeFragment was initialized as false
// 	    	continue;
// 	    }

// 	    //skip reads that were QC failed
// 	    if(  pileupData.PileupAlignments[i].Alignment.IsFailedQC() ){
// 		continue;
// 	    }

// 	    //skip fragments below the minimum length
// 	    if(  pileupData.PileupAlignments[i].Alignment.Length < int(MINLENGTHFRAGMENT) ){ 
// 		cerr<<"skipped minlength "<<pileupData.PileupAlignments[i].Alignment.Name<<endl;
// 		continue;
// 	    }

// 	    if(  pileupData.PileupAlignments[i].Alignment.Length > int(MAXLENGTHFRAGMENT) ){ 
// 		cerr<<"skipped maxlength "<<pileupData.PileupAlignments[i].Alignment.Name<<endl;
// 		continue;
// 	    }
	    
// 	    if(i>=MAXCOV){
// 		break;
// 	    }

// 	    char  b   =     pileupData.PileupAlignments[i].Alignment.QueryBases[ pileupData.PileupAlignments[i].PositionInAlignment ];
// 	    if(!isResolvedDNA(b)){ //avoid Ns
// 		continue; 
// 	    }

// 	    int bIndex = baseResolved2int(b);
// 	    int   q    = MIN2( int(pileupData.PileupAlignments[i].Alignment.Qualities[  pileupData.PileupAlignments[i].PositionInAlignment ]-offsetQual), MAXBASEQUAL);
// 	    int   m    = MIN2( int(pileupData.PileupAlignments[i].Alignment.MapQuality), MAXMAPPINGQUAL );
// 	    bool isRev = pileupData.PileupAlignments[i].Alignment.IsReverseStrand();
// 	    // if(posAlign == 74310){
// 	    // 	cout<<m<<" "<<probMM<<" "<<likeMismatchProbMap[m]<<endl;
// 	    // }
// 	    piToAdd.baseC[bIndex]++;
// 	    probMM += likeMismatchProbMap[m]; 
// 	    basesRetained++;

// 	    totalBases++;
// 	    foundSites=true;

// 	    singleRead sr_;
// 	    sr_.base    = uint8_t(bIndex);
// 	    sr_.qual    = uint8_t(q);	    
// 	    sr_.mapq    = uint8_t(m);
// 	    sr_.lengthF = uint8_t(pileupData.PileupAlignments[i].Alignment.Length);

// 	    if(isRev){
// 		sr_.pos5p = uint8_t(  pileupData.PileupAlignments[i].Alignment.Length-pileupData.PileupAlignments[i].PositionInAlignment-1 ); 
// 		sr_.base = 3 - sr_.base;//complement
// 	    }else{
// 		sr_.pos5p = uint8_t(  pileupData.PileupAlignments[i].PositionInAlignment ); 
// 	    }
// 	    sr_.isrv=isRev;

// #ifdef DEBUGSINGLEREAD
// 	    sr_.name=pileupData.PileupAlignments[i].Alignment.Name;//to remove
// #endif
	    
// 	    piToAdd.readsVec.push_back(sr_);
// 	    // obsBase.push_back( bIndex  );
// 	    // obsQual.push_back( q        );
// 	    // mmQual.push_back(  m        );
// 	    // isRevVec.push_back(isRev);	    
// 	    //mmProb.push_back(  likeMismatchProbMap[m]  );
// 	    // substitutionRatesPerRead.push_back( probSubMatchToUseEndo );
// 	    //	    includeFragment[i]=true;

// 	}//END FOR EACH READ

// 	if( foundSites ){
// 	    piToAdd.avgMQ =  round(-10*log10(probMM/double(basesRetained)));
// 	    // if(posAlign == 74310){
// 	    // 	cout<<piToAdd.avgMQ<<" "<<probMM<<" "<<basesRetained<<endl;
// 	    // }

// 	    totalSites++;
// 	}
	
// 	m_piForGenomicWindow->push_back(piToAdd);

//     }//end Visit()
    
//     unsigned int getTotalSites() const{
// 	return totalSites;
//     }

//     unsigned int getTotalBases() const{
// 	return totalBases;
//     }

// private:
//     RefVector m_references;
//     int          m_refID;
//     unsigned int m_leftCoord;
//     unsigned int m_rightCoord;
//     int          m_threadID;
//     unsigned int m_numberOfSites;

//     unsigned int totalBases;
//     unsigned int totalSites;

//     vector<positionInformation> * m_piForGenomicWindow;
//     Fasta * m_fastaReference;

//     //vector<PositionResult *> * m_dataToWriteOut;
// };//heteroComputerVisitor


	// , m_references(references)
	// , m_refID(refID)
	// , m_leftCoord(leftCoord)
	// , m_rightCoord(rightCoord)
	//   //, m_dataToWriteOut( dataToWriteOut)
	// , m_threadID( threadID )
	// , m_numberOfSites( 0 )
	// , totalBases(0)
	// , totalSites(0)   
	// , m_piForGenomicWindow( piForGenomicWindow )
	// , m_fastaReference(fastaReference)



//#define DEBUGHTS

static int read_bamHET(void *data, bam1_t *b){ // read level filters better go here to avoid pileup
    //cerr<<"read_bamHET"<<endl;
    aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
    int ret;
    while (1){
        ret = aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b);
        if ( ret<0 ) break;
	// int32_t qlen   = bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b));
	int32_t qlen   = b->core.l_qseq;

	//int32_t isize  = b->core.n_cigar;
	// cerr<<"read_bamHET1: "<<bam1_qname(b)<<endl;	  

	//skip fragments marked as:
	//1: unmapped
	//2: secondary alignment, only retain primary
	//3: QC failed
	//4: PCR duplicates
        if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
	// cerr<<"read_bamHET2: "<<bam1_qname(b)<<endl;	  
	//if a reads is paired and not properly paired, skip
	if( bam_is_paired( b) &&    
	    !bam_is_properpair(b) ){ 
	    continue; 
	}
		


        if ( (int)b->core.qual < aux->min_mapQ ) continue;
	// cerr<<"read_bamHET3: "<<bam1_qname(b)<<endl;	  
	// cerr<<aux->min_len<<" "<<qlen<<" "<<int(MINLENGTHFRAGMENT)<<" "<<int(MAXLENGTHFRAGMENT)<<endl;
        if ( aux->min_len && ( (qlen <  int(MINLENGTHFRAGMENT) ) || (qlen > int(MAXLENGTHFRAGMENT) ) )) continue;
	//cerr<<"read_bamHET4: "<<bam1_qname(b)<<endl;	  
	if(specifiedDeam){
	    if( bam_is_paired( b) ) continue; //skip paired-end reads because we cannot get proper deamination

	    uint8_t *rgptr = bam_aux_get(b, "RG");
	    // cerr<<"rg1 "<<rgptr<<endl;
	    // cout<<"isize "<<isize<<endl;
	    string rg = "UNKNOWN";
	    
	    if(rgptr){
		rg = string( (const char*)(rgptr+1));
	    }
	    //cerr<<"rg2 "<<rg<<endl;
	    //TODO
	    if(rg2info.find(rg) == rg2info.end()){
		cerr << "Heterozygosity computation: found an unknown RG tag from  "<<bam1_qname(b)<<" for BAM file:" << bamFileToOpen <<endl;
		//exit(1);	
	    }else{
		if(!rg2info[rg].isPe){ //if single end
		    if(qlen == rg2info[rg].maxReadLength){//probably reached the end of the read length, we will keep maxlength-1 and under		
			//cerr<<"name "<<bam1_qname(b)<<"\t"<<rg<<endl;
			continue;
		    }
		}else{
		    //since we skipped the paired reads, if we have reached here
		    //we can add single-end fragments 		    
		}
	    }	    
	    // }else{
	    // 	//cerr<<"WARNING: Could not retrieve RG tag from : "<<bam1_qname(b)<<" discarding"<<endl;	  	    
	    // }
	}
	// if(specifiedDeam){
	//     if( bam_is_paired( b) ) continue; //skip paired-end reads because we cannot get proper deamination
	// }
	// cerr<<"read_bamHET5: "<<bam1_qname(b)<<endl;	  

        break;
    }
    return ret;
}

#define AROUNDINDELS //if enabled, will remove the base after an indel at the cost of a bit of time


//! Method called for each thread
/*!
  

*/				
void *mainHeteroComputationThread(void * argc){
    if(verbose)
	cerr<<"mainHeteroComputationThread started"<<endl;

    int   rc;

    //#ifdef HETVERBOSE    
    int rankThread=0;
    //#endif

    // if(verbose)    
    // 	cerr<<"mainHeteroComputationThread mutex1"<<endl;

    rc = pthread_mutex_lock(&mutexRank);
    checkResults("pthread_mutex_lock()\n", rc);
    //     if(verbose)
    //     cerr<<"mainHeteroComputationThread mutex21 ID="<<(*(int *)pthread_self())<<endl;    
    //     cerr<<"mainHeteroComputationThread mutex22 ID="<<(threadID2Rank.size()+1)<<endl;
    // }
    //cerr<<"mainHeteroComputationThread mutex2"<<endl;    

    threadID2Rank[*(int *)pthread_self()]  = threadID2Rank.size()+1;

    //cerr<<"mainHeteroComputationThread mutex3"<<endl;
    //#ifdef HETVERBOSE    
    rankThread = threadID2Rank[*(int *)pthread_self()];
    //#endif
    
    rc = pthread_mutex_unlock(&mutexRank);
    checkResults("pthread_mutex_unlock()\n", rc);

 checkqueue:    
    // stackIndex=-1;
    //check stack

    
    rc = pthread_mutex_lock(&mutexQueueToRead);
    checkResults("pthread_mutex_lock()\n", rc);


    bool foundData=false;
    
#ifdef HETVERBOSE
    if(verbose)    {
	rc = pthread_mutex_lock(&mutexCERR);
	checkResults("pthread_mutex_lock()\n", rc);


	cerr<<"Thread #"<<rankThread <<" started and is requesting data"<<endl;

	rc = pthread_mutex_unlock(&mutexCERR);
	checkResults("pthread_mutex_unlock()\n", rc);
    }
#endif

    DataChunk    * currentChunk;


    if(!queueDataToprocess.empty()){    
 	foundData=true;
 	currentChunk = queueDataToprocess.front();
 	queueDataToprocess.pop();

#ifdef HETVERBOSE
	if(verbose){
	    rc = pthread_mutex_lock(&mutexCERR);
	    checkResults("pthread_mutex_lock()\n", rc);

	    cerr<<"Thread #"<<rankThread<<" is reading chunk rank#"<<currentChunk->rank<<endl;

	    rc = pthread_mutex_unlock(&mutexCERR);
	    checkResults("pthread_mutex_unlock()\n", rc);
	}
#endif
	//cout<<"rank "<< &(currentChunk->dataToProcess) <<endl;
    }

    
  

    if(!foundData){
	rc = pthread_mutex_unlock(&mutexQueueToRead);
	checkResults("pthread_mutex_unlock()\n", rc);


	if(readDataDone){
#ifdef HETVERBOSE
	    cerr<<"Thread #"<<rankThread<<" is done"<<endl;
#endif
	    return NULL;	
	}else{

#ifdef HETVERBOSE
	    cerr<<"Thread #"<<rankThread<<" sleeping for "<<timeThreadSleep<<endl;
#endif
	    sleep(timeThreadSleep);
	    goto checkqueue;
	}

    }else{
	//release stack
	rc = pthread_mutex_unlock(&mutexQueueToRead);
	checkResults("pthread_mutex_unlock()\n", rc);
    }

    //////////////////////////////////////////////////////////////////
    //                BEGIN HET COMPUTATION                         //
    //////////////////////////////////////////////////////////////////

    //cout<<currentChunk->rangeGen<<endl;

#ifdef HETVERBOSE
    cerr<<"Thread #"<<rankThread<<" is reading region "<<currentChunk->rangeGen<<endl;
#endif

    //sleep(10);


    


#ifdef HETVERBOSE
    cerr<<"Thread #"<<rankThread<<" is reading BAM file: "<<bamFileToOpen<<endl;
#endif

    // BamReader reader;
    // if ( !reader.Open(bamFileToOpen) ) {
    // 	cerr << "Could not open input BAM file:" << bamFileToOpen <<endl;
    // 	exit(1);
    // }

    // reader.LocateIndex();

    // if(!reader.HasIndex()){
    // 	cerr << "The BAM file: " << bamFileToOpen <<" does not have an index"<<endl;
    // 	exit(1);
    // }

    // // retrieve reference data
    // const RefVector  references = reader.GetReferenceData();
    // const int        refID      = reader.GetReferenceID( currentChunk->rangeGen.getChrName() );

    // #ifdef HETVERBOSE
    //     cerr<<"Thread #"<<rankThread<<" refID   "<<refID<<endl;    
    //      cerr<<"Thread #"<<rankThread<<" refName "<<references[refID].RefName<<endl;
    // #endif
    
    // BamRegion bregion (refID, 
    // 		       currentChunk->rangeGen.getStartCoord(), 
    // 		       refID, 
    // 		       currentChunk->rangeGen.getEndCoord()   );

    // bool setRegionRes=reader.SetRegion( bregion   );

    //if(verbose){
    // if(verbose){
    // 	rc = pthread_mutex_lock(&mutexCERR);
    // 	checkResults("pthread_mutex_lock()\n", rc);
	
    // 	cerr<<"Thread #"<<rankThread<<" setting region: "<<references[refID].RefName<<":"<<currentChunk->rangeGen.getStartCoord()<<"-"<<currentChunk->rangeGen.getEndCoord()<<"\t"<<(setRegionRes?"success!":"failed")<<endl;
	
    // 	rc = pthread_mutex_unlock(&mutexCERR);
    // 	checkResults("pthread_mutex_unlock()\n", rc);
    // }
    
    //}
#ifdef HETVERBOSE

#endif

    // if( refID==-1 ||
    //    !setRegionRes){
    // 	cerr << "Heterozygous computation: could not set region "<<currentChunk->rangeGen<<" for BAM file:" << bamFileToOpen <<" "<<currentChunk->rangeGen.getStartCoord()<<" "<< currentChunk->rangeGen.getEndCoord()<< endl;
    // 	exit(1);
    // }

    DataToWrite  * dataToWrite = new DataToWrite();
    //cout<<"range1 "<<currentChunk->rangeGen<<endl;
    dataToWrite->rangeGen      =  currentChunk->rangeGen;
    //cout<<"range2 "<<currentChunk->rangeGen<<endl;
    dataToWrite->rank          =  currentChunk->rank;
    //dataToWrite->refID         =  refID;

    string region      = currentChunk->rangeGen.getChrName()+":"+stringify(currentChunk->rangeGen.getStartCoord())+"-"+stringify(currentChunk->rangeGen.getEndCoord());

    faidx_t *fai = fai_load(fastaFile.c_str());
   
    if ( !fai ) {
         cerr << "ERROR: failed to open fasta file " <<fastaFile<<" and fasta index " << fastaIndex<<""<<endl;
         exit(1);
    }
    int seq_len;
    char *seq = fai_fetch(fai, region.c_str(), &seq_len);
    if ( seq_len < 0 ) {
	cerr << "ERROR: failed to retrieve sequence from fasta file " <<fastaFile<<" and fasta index " << fastaIndex<<" and region "<<region<<endl;
	exit(1);
    }

    if(seq_len != int(   currentChunk->rangeGen.getEndCoord() - currentChunk->rangeGen.getStartCoord() +1 ) ){
	cerr << "ERROR: failed to retrieve sequence from  file " <<fastaFile<<" and fasta index " << fastaIndex<<" asked for "<<int(   currentChunk->rangeGen.getEndCoord() - currentChunk->rangeGen.getStartCoord() )<<" "<<seq_len<<endl;
	exit(1);	
    }


    // Fasta fastaReference;
    
    // if ( !fastaReference.Open(fastaFile , fastaIndex) ){        //from bamtools
    //      cerr << "ERROR: failed to open fasta file " <<fastaFile<<" and fasta index " << fastaIndex<<""<<endl;
    //      exit(1);
    //  }

    vector<positionInformation> * piForGenomicWindow = new vector<positionInformation> ();



    //string bamfilename = bamFileToOpen;

    //int bamdepth( const string &   bamfilename, const string & region, const string & bedfilename){
       









    // unsigned int totalBasesL=0;
    // unsigned int totalSitesL=0;
    

    string bamfilename = bamFileToOpen;


    int i, n, tid, reg_tid, beg, end, pos, *n_plp, baseQ = 0, mapQ = 0, min_len = MINLENGTHFRAGMENT;
    int prevPos=-1000;

    int all = 0; //status = EXIT_SUCCESS, 

    int max_depth = 20*MAXCOV;
    const bam_pileup1_t **plp;
    bool reg = !region.empty();

    void *bed = 0; // BED data structure
    if(!bedfilename.empty()){
	bed = bed_read(bedfilename.c_str()); // BED or position list file can be parsed now
    }

    bam_hdr_t *h = NULL; // BAM header of the 1st input
    aux_t **data;
    bam_mplp_t mplp;
    int last_pos = -1, last_tid = -1, ret;

    
     unsigned int totalBasesL=0;//=cv->getTotalBases();
     unsigned int totalSitesL=0; //=cv->getTotalSites();

    n =1;
    
    data = (aux_t **)calloc(1, sizeof(aux_t*)); // data[i] for the i-th input
    reg_tid = 0; beg = 0; end = INT_MAX;  // set the default region

    for (i = 0; i < n; ++i) {
	//cerr<<"i ="<<i<<endl;
	//for (i = 0; i < 1; ++i) {
        int rf;
        data[i] = (aux_t *)calloc(1, sizeof(aux_t));


	data[i]->fp = sam_open_format(bamfilename.c_str(), "r", NULL); // open BAM

        if (data[i]->fp == NULL) {
	    cerr<<"ERROR: Could not open BAM file "<<bamfilename<<""<<endl;
	    exit(1);

        }
        rf = SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR | SAM_SEQ;
        if (baseQ) rf |= SAM_QUAL;

        data[i]->min_mapQ = mapQ;                    // set the mapQ filter
        data[i]->min_len  = min_len;                 // set the qlen filter
        data[i]->hdr = sam_hdr_read(data[i]->fp);    // read the BAM header
        
	if (data[i]->hdr == NULL) {
	    cerr<<"ERROR: Could not read header for bamfile "<<bamfilename<<""<<endl;
	    exit(1);
        }
        if (reg) { // if a region is specified
	    hts_idx_t *idx = sam_index_load(data[i]->fp, bamfilename.c_str());  // load the index
            if (idx == NULL) {
                //print_error("depth", "can't load index for \"%s\"", argv[optind+i]);
		cerr<<"ERROR: Cannot load index for bamfile "<<bamfilename<<""<<endl;
		exit(1);
            }else{
		//cerr<<"index ok"<<endl;
	    }
            data[i]->iter = sam_itr_querys(idx, data[i]->hdr, region.c_str()); // set the iterator
            hts_idx_destroy(idx); // the index is not needed any more; free the memory
            if (data[i]->iter == NULL) {
		cerr<<"ERROR: Cannot parse region \""<<region<<"\""<<endl;
		exit(1);
            }else{
		//cerr<<"region ok"<<endl;
	    }
        }
    }

    h = data[0]->hdr; // easy access to the header of the 1st BAM
    //dataToWrite->refID = bam_get_seq(h,currentChunk->rangeGen.getChrName().c_str());
    dataToWrite->refID = sam_hdr_name2tid(h,currentChunk->rangeGen.getChrName().c_str());

    if (reg) {
        beg     = data[0]->iter->beg; // and to the parsed region coordinates
        end     = data[0]->iter->end;
        reg_tid = data[0]->iter->tid;
    }

    // the core multi-pileup loop
    mplp = bam_mplp_init(n, read_bamHET, (void**)data); // initialization
    if (0 < max_depth)
        bam_mplp_set_maxcnt(mplp,max_depth);  // set maximum coverage depth
    else if (!max_depth)
        bam_mplp_set_maxcnt(mplp,INT_MAX);
    n_plp = (int *)calloc(n, sizeof(int)); // n_plp[i] is the number of covering reads from the i-th BAM
    plp   = (const bam_pileup1_t **)calloc(n, sizeof(bam_pileup1_t*)); // plp[i] points to the array of covering reads (internal in mplp)
    //int seqptr=0;

#ifdef AROUNDINDELS
    bool prevIndel   =false;
    // int  prevLevels   [2*MAXCOV];
    // int  prevLevelsi =0;
    set<string> prevIndelSet;

    bool currIndel   =false;
    // int  currLevels   [2*MAXCOV];
    // int  currLevelsi =0;	    
    set<string> currIndelSet;

	    
#endif
    
    while ((ret=bam_mplp_auto(mplp, &tid, &pos, n_plp, plp)) > 0) { // come to the next covered position	
	//cerr<<endl<<(pos+1)<<" "<<" -----------------------"<<endl;
	
	if (pos < beg || pos >= end) continue; // out of range; skip
        if (tid >= h->n_targets) continue;     // diff number of @SQ lines per file?
	
	char refC = toupper(seq[ pos-currentChunk->rangeGen.getStartCoord()+1 ]);
	if(refC == 'N'){//skip unresolved
	    continue;
	}

	// #ifdef DEBUGHTS
	// 	cerr<<endl<<(pos+1)<<" "<<refc<<" -----------------------"<<endl;
	// #endif
	
        if (all) {
            while (tid > last_tid) {
                if (last_tid >= 0 && !reg) {
                    // Deal with remainder or entirety of last tid.
                    while (++last_pos < int(h->target_len[last_tid])) {
			//cerr<<(last_pos+1)<<endl;
			//exit(1);
			// // Horribly inefficient, but the bed API is an obfuscated black box.
                        if (bed && bed_overlap(bed, h->target_name[last_tid], last_pos, last_pos + 1) == 0)
                            continue;
			totalSitesL++;
                    }
                }
                last_tid++;
                last_pos = -1;
                if (all < 2)
                    break;
            }

            // Deal with missing portion of current tid
            while (++last_pos < pos) {
                if (last_pos < beg) continue; // out of range; skip
                if (bed && bed_overlap(bed, h->target_name[tid], last_pos, last_pos + 1) == 0)
                    continue;
		totalSitesL++;
            }

            last_tid = tid;
            last_pos = pos;
        }
	if (bed && bed_overlap(bed, h->target_name[tid], pos, pos + 1) == 0) continue;

	
#ifdef DEBUGHTS
	cout<<endl<<(pos+1)<<" "<<refC<<" -----------------------"<<endl;
#endif
	
#ifdef AROUNDINDELS
	//if indels are found, remove 
	if( (prevPos+1) == pos){
	    prevIndel     = currIndel;
	    prevIndelSet  = currIndelSet;
	    // prevLevelsi = currLevelsi;
	    // for(int ci=0;ci<currLevelsi;ci++){
	    //     prevLevels[ci] = currLevels[ci];
	    // }		
	}else{
	    prevIndel   =false;
	    //prevLevelsi =0;
	    prevIndelSet.clear();
	}
#endif

	//totalSitesL++;

	//fputs(h->target_name[tid], stdout); 
	//printf("\t%d\t", pos+1); // a customized printf() would be faster
	//cerr<<(pos+1)<<endl;
	//previous indels

        for (i = 0; i < n; ++i) { // base level filters have to go here
            int j;
	    int m = 0;//amount of invalid bases at site j that need to be removed from coverage calculations

	    
	    positionInformation piToAdd;

	    piToAdd.skipPosition                 = false;
	    piToAdd.posAlign                     = (pos+1);
	    piToAdd.refBase                      = refC;
	    double probMM=0;
	    int    basesRetained=0;
	    bool foundSites=false;

	    for(int n=0;n<4;n++){
		piToAdd.baseC[n]=0;
	    }
	    
#ifdef AROUNDINDELS
	    currIndel   =false;
	    // int  currLevels   [2*MAXCOV];
	    // int  currLevelsi =0;	    
	    currIndelSet.clear();
#endif
	    
            for (j = 0; j < n_plp[i]; ++j) {
                const bam_pileup1_t *p = plp[i] + j; // DON'T modfity plp[][] unless you really know
		// if( (pos+1) == 16050348){
		//     cerr<<bam1_qname(p->b)<<endl;	  
		// }
#ifdef DEBUGHTS
		cerr<< bam1_qname(p->b)<<" "<<j<<" "<<int((p->b)->core.flag)<<" "<<p->b->core.n_cigar<<" p="<<(p->cd.p)<<" i="<<int(p->cd.i)<<" f="<<float(p->cd.f)<<" "<<p->is_del<<" "<<p->is_refskip<<" "<<p->indel<<" "<<p->level<<" "<<p->aux<<" "<<endl;

#endif

#ifdef AROUNDINDELS	    
		//cerr<<"prevIndel1: "<<prevIndel<<" "<<prevIndelSet.size()<<endl;
#endif
		//prevIndel<<endl;
		
		//base level filters go here
		if (p->is_del || p->is_refskip){

#ifdef AROUNDINDELS	    
		    //currLevels[currLevelsi++] = p->level;
		    currIndel   = true;
		    //cerr<<"inserting "<<bam1_qname(p->b)<<endl;
		    currIndelSet.insert( bam1_qname(p->b) );
#endif
		    ++m; // having dels or refskips at tid:pos
		    continue;
		}

		//cerr<<"prevIndel2: "<<prevIndel<<endl;

		if ( p->indel != 0 ){// having dels or refskips at the next
		    ++m; 
		    continue;
		}

#ifdef AROUNDINDELS

		//cerr<<"prevIndel3: "<<prevIndel<<" "<<prevIndelSet.size()<<endl;
		
		if(prevIndel){
		    //cerr<<"prevIndel, testing  "<<bam1_qname(p->b)<<""<<endl;
		    
		    // for(set<string>::iterator it = prevIndelSet.begin(); it != prevIndelSet.end(); it++)  {
		    //  	cerr<<"record: "<<*it<<endl;
		    //  }
		    
		    if( prevIndelSet.find(bam1_qname(p->b)) != prevIndelSet.end() ){//skip reads previously flagged as indels. Since indels are few, this should be a quicker solution than to read ahead 
			//cerr<<"found  "<<bam1_qname(p->b)<<" in a previous indel"<<endl;			
			continue;
		    }
		    // bool skippos=false;
		    // for(int pi=0;pi<prevLevelsi;pi++){
		    // 	if(p->level == prevLevels[prevLevelsi]){
		    // 	    skippos=true;
		    // 	    break;
		    // 	}
		    // }
		    // if(skippos){
		    // 	continue;
		    // }
		}
#endif		



		// else{ 
		//     if(bam_get_qual(p->b)[p->qpos] < baseQ) 
		// 	++m; // low base quality
		// }
		//TODO reenable
		if( j>=MAXCOV){
		    break;
		}

		int bIndex = alphabetHTSLIB2idx[ bam_seqi(bam_get_seq(p->b),p->qpos) ];
		if(bIndex == -1){ continue; } //skip unresolved

		int   q    = MIN2( int(bam_get_qual(p->b)[p->qpos] ), MAXBASEQUAL);//offsetqual?
		int   m    = MIN2(  bam_mqual(p->b) , MAXMAPPINGQUAL );
		bool isRev = bam_is_rev(p->b);

		piToAdd.baseC[bIndex]++;
		probMM += likeMismatchProbMap[m]; 
		basesRetained++;


		foundSites=true;

		singleRead sr_;
		sr_.base    = uint8_t(bIndex);
		sr_.qual    = uint8_t(q);	    
		sr_.mapq    = uint8_t(m);
		sr_.lengthF = uint8_t( (p->b)->core.l_qseq );
		
		if(isRev){
		    sr_.pos5p = uint8_t( sr_.lengthF - p->qpos -1 );
		    sr_.base  = 3 - sr_.base;//complement
		}else{
		    sr_.pos5p = uint8_t( p->qpos ); 
		}
		sr_.isrv=isRev;

#ifdef DEBUGHTS
		//cerr<<isRev<<" "<<"ACGT"[int(sr_.base)]<<" "<<int(sr_.qual)<<" "<<int(sr_.mapq)<<" "<<int(m)<<" "<<int(sr_.lengthF)<<" "<<int(sr_.pos5p)<<" "<< bam1_qname(p->b)<<" "<<int((p->b)->core.flag)<<" "<<p->b->core.n_cigar<<" p="<<(p->cd.p)<<" i="<<int(p->cd.i)<<" f="<<float(p->cd.f)<<" "<<p->indel<<" "<<p->level<<endl;
		cerr<<"ADD "<<isRev<<" "<<int(sr_.base)<<" "<<int(sr_.qual)<<" "<<int(sr_.mapq)<<" "<<int(m)<<" "<<int(sr_.lengthF)<<" "<<int(sr_.pos5p)<<" "<< bam1_qname(p->b)<<" "<<int((p->b)->core.flag)<<endl;
#endif
		
		// if(isRev){
		//     sr_.pos5p = uint8_t(  pileupData.PileupAlignments[i].Alignment.Length-pileupData.PileupAlignments[i].PositionInAlignment-1 ); 
		//     sr_.base = 3 - sr_.base;//complement
		// }else{
		//     sr_.pos5p = uint8_t(  pileupData.PileupAlignments[i].PositionInAlignment ); 
		// }
		// sr_.isrv=isRev;
		totalBasesL++;
		piToAdd.readsVec.push_back(sr_);

//currLevelVec.push_back(p->level);
		//totalBasesL ++;

		
#ifdef DEBUGHTS
		if(0){
		// printf("%d,",p->b);
		// printf("%d,",bam_get_seq(p->b));
		// @discussion Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G,
		//  8 for T and 15 for N. Two bases are packed in one byte with the base
		//  at the higher 4 bits having smaller coordinate on the read. It is
		//  recommended to use bam_seqi() macro to get the base.
		//  */
		//coord
		//printf("%d,",p->qpos);
		cerr<<p->qpos<<"/"<<int(sr_.lengthF) <<",";
		//flag
		//printf("%d,",(p->b)->core.flag);
		cerr<<((p->b)->core.flag)<<",";
		//base
		//printf("%d,",bam_seqi(bam_get_seq(p->b),p->qpos));
		//printf("%c,",alphabetHTSLIB[ bam_seqi(bam_get_seq(p->b),p->qpos) ]);
		cerr<<alphabetHTSLIB[ bam_seqi(bam_get_seq(p->b),p->qpos) ]<<",";
		
		//int c = seq_nt16_int[bam_seqi(seq, p->qpos)];
		//qual
		//printf("%d,",bam_get_qual(p->b)[p->qpos]);
		cerr<<int(bam_get_qual(p->b)[p->qpos])<<",";
		//strand
		//printf("%d,",bam_is_rev(p->b));
		cerr<<bam_is_rev(p->b)<<",";
		
		//paired
		//printf("%d,",bam_is_paired(p->b));
		cerr<<bam_is_paired(p->b)<<",";
		cerr<<"D="<<p->is_del<<",RS="<<p->is_refskip<<",";
		//fail
		//printf("%d-",bam_is_failed(p->b));
		cerr<<bam_is_failed(p->b)<<",";
		cerr<<bam_mqual(p->b)<<"-";
		}
#endif
            }//end for each base at pos
            //printf("\tc=%d", n_plp[i] - m); // this the depth to output
	    //cout<<(n_plp[i] - m); // this the depth to output
	    // totalSitesL ++ ;

	    
	    if( foundSites ){
		piToAdd.avgMQ =  round(-10*log10(probMM/double(basesRetained)));
		totalSitesL++;
	    }

	    piForGenomicWindow->push_back(piToAdd);
	    //cerr<<"readsVec size= "<<piToAdd.readsVec.size()<<endl;
	    prevPos = pos;

	}//end for all pos


	// #ifdef DEBUGHTS
	// 	cerr<<endl<<" -----------------------"<<endl;
	// #endif
        
        //putchar('\n');
	//cout<<endl;
    }//end while mpileup
    
    if (ret < 0){ //status = EXIT_FAILURE;
	cerr<<"Problem parsing region:"<<region<<" in bamfile "<<bamfilename<<endl;
	exit(1);
    }
    free(n_plp); free(plp);
    bam_mplp_destroy(mplp);

    if (all) {
        // Handle terminating region
        if (last_tid < 0 && reg && all > 1) {
            last_tid = reg_tid;
            last_pos = beg-1;
        }
        while (last_tid >= 0 && last_tid < h->n_targets) {
            while (++last_pos < int(h->target_len[last_tid])) {
                if (last_pos >= end) break;
                if (bed && bed_overlap(bed, h->target_name[last_tid], last_pos, last_pos + 1) == 0)
                    continue;
		totalSitesL++;
            }
            last_tid++;
            last_pos = -1;
            if (all < 2 || reg)
                break;
        }
    }

    //depth_end:
    for (i = 0; i < n && data[i]; ++i) {
        bam_hdr_destroy(data[i]->hdr);
        if (data[i]->fp) sam_close(data[i]->fp);
        hts_itr_destroy(data[i]->iter);
        free(data[i]);
    }
    free(data); 
    free(seq);
    fai_destroy(fai);

    
    //free(reg);
    if (bed) bed_destroy(bed);

    
//     int i, n, tid, reg_tid, beg, end, pos, *n_plp, baseQ = 0, mapQ = 0, min_len = MINLENGTHFRAGMENT;
//     //int all = 0, status = EXIT_SUCCESS, nfiles, max_depth = -1;
//     int all = 0; //status = EXIT_SUCCESS, 
//     //int max_depth = -1;
//     int max_depth = MAXCOV;
//     const bam_pileup1_t **plp;
//     bool reg = !region.empty();
//     //bool bed = !bedfilename.empty();
//     //char *reg = 0; // specified region
//     void *bed = 0; // BED data structure
//     if(!bedfilename.empty()){
// 	bed = bed_read(bedfilename.c_str()); // BED or position list file can be parsed now
//     }

//     bam_hdr_t *h = NULL; // BAM header of the 1st input
//     aux_t **data;
//     bam_mplp_t mplp;
//     int last_pos = -1, last_tid = -1, ret;

//     n =1;

//     unsigned int totalBasesL=0;//=cv->getTotalBases();
//     unsigned int totalSitesL=0; //=cv->getTotalSites();

    
//     data = (aux_t **)calloc(1, sizeof(aux_t*)); // data[i] for the i-th input
//     reg_tid = 0; beg = 0; end = INT_MAX;  // set the default region
//     // bam_hdr_t *hin;

//     for (i = 0; i < n; ++i) {
// 	//cerr<<"i ="<<i<<endl;
// 	//for (i = 0; i < 1; ++i) {
//         int rf;
//         data[i] = (aux_t *)calloc(1, sizeof(aux_t));

// 	//cerr<<"bamfilename1 "<<bamfilename<<endl;

// 	data[i]->fp = sam_open_format(bamfilename.c_str(), "r", NULL); // open BAM
// 	//cerr<<"bamfilename2 "<<bamfilename<<endl;
        

//         if (data[i]->fp == NULL) {
// 	    cerr<<"Could not open \""<<bamfilename<<"\""<<endl;
// 	    exit(1);
//             //status = EXIT_FAILURE;
//             //goto depth_end;
//         }
// 	// hin = sam_hdr_read(fp[i]);
// 	// if (hin == NULL) {
// 	//     cerr<<"Could not read header from  \""<<bamfilename<<"\""<<endl;
// 	//     exit(1);
//         // }


//         rf = SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR | SAM_SEQ;
//         if (baseQ) rf |= SAM_QUAL;
//         // if (hts_set_opt(data[i]->fp, CRAM_OPT_REQUIRED_FIELDS, rf)) {
//         //     fprintf(stderr, "Failed to set CRAM_OPT_REQUIRED_FIELDS value\n");
//         //     return 1;
//         // }
//         // if (hts_set_opt(data[i]->fp, CRAM_OPT_DECODE_MD, 0)) {
//         //     fprintf(stderr, "Failed to set CRAM_OPT_DECODE_MD value\n");
//         //     return 1;
//         // }
//         data[i]->min_mapQ = mapQ;                    // set the mapQ filter
//         data[i]->min_len  = min_len;                 // set the qlen filter
//         data[i]->hdr = sam_hdr_read(data[i]->fp);    // read the BAM header
        
// 	if (data[i]->hdr == NULL) {
// 	    cerr<<"Could not read header for \""<<bamfilename<<"\""<<endl;
// 	    exit(1);
//             //status = EXIT_FAILURE;
//             //goto depth_end;
//         }
//         if (reg) { // if a region is specified
// 	//if (true) { // if a region is specified
//             // hts_idx_t *idx = sam_index_load(data[i]->fp, argv[optind+i]);  // load the index
// 	    hts_idx_t *idx = sam_index_load(data[i]->fp, bamfilename.c_str());  // load the index
// 	    // //int testingIDX = ti_get_tid( (const ti_index_t* )idx, currentChunk->rangeGen.getChrName().c_str()); // set the iterator
// 	    // //int testingIDX = get_tid( (const tbx_t* )idx, currentChunk->rangeGen.getChrName().c_str(),0); // set the iterator
// 	    // int32_t tempREFID = bam_get_tid(data[i]->hdr, currentChunk->rangeGen.getChrName().c_str() ); // set the iterator
// 	    // cout<<"tempREFID "<<tempREFID<<endl;
// 	    // dataToWrite->refID         =  tempREFID;
// 	    // exit(1);
//             if (idx == NULL) {
//                 //print_error("depth", "can't load index for \"%s\"", argv[optind+i]);
// 		cerr<<"cannot load index for \""<<bamfilename<<"\""<<endl;
// 		exit(1);
//                 //status = EXIT_FAILURE;
//                 //goto depth_end;
//             }else{
// 		//cerr<<"index ok"<<endl;
// 	    }
//             data[i]->iter = sam_itr_querys(idx, data[i]->hdr, region.c_str()); // set the iterator
//             hts_idx_destroy(idx); // the index is not needed any more; free the memory
//             if (data[i]->iter == NULL) {
//                 //print_error("depth", "can't parse region \"%s\"", reg);
// 		cerr<<"cannot parse region \""<<region<<"\""<<endl;
// 		exit(1);
//                 //status = EXIT_FAILURE;
//                 //goto depth_end;
//             }else{
// 		//cerr<<"region ok"<<endl;
// 	    }

// 	    if(verbose){
// 		rc = pthread_mutex_lock(&mutexCERR);
// 		checkResults("pthread_mutex_lock()\n", rc);
	
// 		cerr<<"Thread #"<<rankThread<<" setting region: "<<region<<"\t"<<"success!"<<endl;
	
// 		rc = pthread_mutex_unlock(&mutexCERR);
// 		checkResults("pthread_mutex_unlock()\n", rc);
// 	    }


//         }
//     }

//     h = data[0]->hdr; // easy access to the header of the 1st BAM
//     if (reg) {
//     //if(true){
//         beg     = data[0]->iter->beg; // and to the parsed region coordinates
//         end     = data[0]->iter->end;
//         reg_tid = data[0]->iter->tid;
// 	dataToWrite->refID = reg_tid ;
//     }

//     // the core multi-pileup loop
//     mplp = bam_mplp_init(n, read_bamHET, (void**)data); // initialization
//     if (0 < max_depth)
//         bam_mplp_set_maxcnt(mplp,max_depth);  // set maximum coverage depth
//     else if (!max_depth)
//         bam_mplp_set_maxcnt(mplp,INT_MAX);
//     n_plp = (int *)calloc(n, sizeof(int)); // n_plp[i] is the number of covering reads from the i-th BAM
//     plp   = (const bam_pileup1_t **)calloc(n, sizeof(bam_pileup1_t*)); // plp[i] points to the array of covering reads (internal in mplp)

//     while ((ret=bam_mplp_auto(mplp, &tid, &pos, n_plp, plp)) > 0) { // come to the next covered position
//         if (pos < beg || pos >= end) continue; // out of range; skip
//         if (tid >= h->n_targets) continue;     // diff number of @SQ lines per file?
// 	cout<<pos<<endl;
//         if (all) {
//             while (tid > last_tid) {
//                 if (last_tid >= 0 && !reg) {
//                     // Deal with remainder or entirety of last tid.
//                     while (++last_pos < int(h->target_len[last_tid])) {
//                         // // Horribly inefficient, but the bed API is an obfuscated black box.
//                         if (bed && bed_overlap(bed, h->target_name[last_tid], last_pos, last_pos + 1) == 0)
//                             continue;
//                         //fputs(h->target_name[last_tid], stdout); 
// 			//printf("\t%d", last_pos+1);
// 			//cout<<(last_pos+1)<<endl;
//                         // for (i = 0; i < n; i++)
//                         //     putchar('\t'), putchar('0');
//                         // putchar('\n');
//                     }
//                 }
//                 last_tid++;
//                 last_pos = -1;
//                 if (all < 2)
//                     break;
//             }

//             // Deal with missing portion of current tid
//             while (++last_pos < pos) {
//                 if (last_pos < beg) continue; // out of range; skip
//                 if (bed && bed_overlap(bed, h->target_name[tid], last_pos, last_pos + 1) == 0)
//                     continue;
//                 fputs(h->target_name[tid], stdout); 
// 		//printf("\t%d", last_pos+1);
//                 // for (i = 0; i < n; i++)
//                 //     putchar('\t'), putchar('0');
//                 // putchar('\n');
//             }

//             last_tid = tid;
//             last_pos = pos;
//         }
// 	if (bed && bed_overlap(bed, h->target_name[tid], pos, pos + 1) == 0) continue;
//         //fputs(h->target_name[tid], stdout); 
// 	//printf("\t%d\t", pos+1); // a customized printf() would be faster
//         for (i = 0; i < n; ++i) { // base level filters have to go here
//             int j;
// 	    int m = 0;//amount of invalid bases at site j that need to be removed from coverage calculations
	    
    
// 	    totalSitesL++; //=cv->getTotalSites();

//             for (j = 0; j < n_plp[i]; ++j) {
//                 const bam_pileup1_t *p = plp[i] + j; // DON'T modfity plp[][] unless you really know
// 		//base level filters go here
// 		if (p->is_del || p->is_refskip) 
// 		    ++m; // having dels or refskips at tid:pos
//                 else 
// 		    if(bam_get_qual(p->b)[p->qpos] < baseQ) 
// 			++m; // low base quality
	
// 		totalBasesL++;//=cv->getTotalBases();
// #ifndef NOTDEF
// 		// printf("%d,",p->b);
// 		// printf("%d,",bam_get_seq(p->b));
// 		// @discussion Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G,
// 		//  8 for T and 15 for N. Two bases are packed in one byte with the base
// 		//  at the higher 4 bits having smaller coordinate on the read. It is
// 		//  recommended to use bam_seqi() macro to get the base.
// 		//  */
// 		//coord
// 		//printf("%d,",p->qpos);
// 		cout<<p->qpos<<",";
// 		//flag
// 		//printf("%d,",(p->b)->core.flag);
// 		cout<<((p->b)->core.flag)<<",";
// 		//base
// 		//printf("%d,",bam_seqi(bam_get_seq(p->b),p->qpos));
// 		//printf("%c,",alphabet[ bam_seqi(bam_get_seq(p->b),p->qpos) ]);
// 		cout<<alphabet[ bam_seqi(bam_get_seq(p->b),p->qpos) ]<<",";
		
// 		//int c = seq_nt16_int[bam_seqi(seq, p->qpos)];
// 		//qual
// 		//printf("%d,",bam_get_qual(p->b)[p->qpos]);
// 		cout<<int(bam_get_qual(p->b)[p->qpos])<<",";
// 		//strand
// 		//printf("%d,",bam_is_rev(p->b));
// 		cout<<bam_is_rev(p->b)<<",";
		
// 		//paired
// 		//printf("%d,",bam_is_paired(p->b));
// 		cout<<bam_is_paired(p->b)<<",";
// 		cout<<"D="<<p->is_del<<",RS="<<p->is_refskip<<",";
// 		//fail
// 		//printf("%d-",bam_is_failed(p->b));
// 		cout<<bam_is_failed(p->b)<<"-";
// #endif
//             }
//             //printf("\tc=%d", n_plp[i] - m); // this the depth to output
// 	    //cout<<(n_plp[i] - m); // this the depth to output
// 	    totalBasesL += (n_plp[i] - m); // this the depth to output
// 	    totalSitesL ++ ;

//         }
//         //putchar('\n');
// 	//cout<<endl;
//     }
    
//     if (ret < 0){ //status = EXIT_FAILURE;
// 	cerr<<"Problem parsing region:"<<region<<" in bamfile "<<bamfilename<<endl;
// 	exit(1);
//     }
//     free(n_plp); free(plp);
//     bam_mplp_destroy(mplp);

//     if (all) {
//         // Handle terminating region
//         if (last_tid < 0 && reg && all > 1) {
//             last_tid = reg_tid;
//             last_pos = beg-1;
//         }
//         while (last_tid >= 0 && last_tid < h->n_targets) {
//             while (++last_pos < int(h->target_len[last_tid])) {
//                 if (last_pos >= end) break;
//                 if (bed && bed_overlap(bed, h->target_name[last_tid], last_pos, last_pos + 1) == 0)
//                     continue;
//                 fputs(h->target_name[last_tid], stdout); 
// 		printf("\t%d", last_pos+1);
//                 for (i = 0; i < n; i++)
//                     putchar('\t'), putchar('0');
//                 putchar('\n');
//             }
//             last_tid++;
//             last_pos = -1;
//             if (all < 2 || reg)
//                 break;
//         }
//     }

//     //depth_end:
//     for (i = 0; i < n && data[i]; ++i) {
//         bam_hdr_destroy(data[i]->hdr);
//         if (data[i]->fp) sam_close(data[i]->fp);
//         hts_itr_destroy(data[i]->iter);
//         free(data[i]);
//     }
//     free(data); 
    
//     //free(reg);
//     if (bed) bed_destroy(bed);
//     free(seq);
//     fai_destroy(fai);

    //dataToWrite->dataToWriteOut=new vector<PositionResult *>();
    // heteroComputerVisitor* hv = new heteroComputerVisitor(references,
    // 							  refID,
    // 							  currentChunk->rangeGen.getStartCoord(), 
    // 							  currentChunk->rangeGen.getEndCoord()  ,
    // 							  // dataToWrite->vecPositionResults,
    // 							  rankThread,
    // 							  piForGenomicWindow,
    // 							  &fastaReference);

    

    // PileupEngine pileup;
    // pileup.AddVisitor(hv);

    // BamAlignment al;
    // //unsigned int numReads=0;
    // //cerr<<MINLENGTHFRAGMENT<<" "<<MAXLENGTHFRAGMENT<<endl;
    // while ( reader.GetNextAlignment(al) ) {
    //     //cout<<"mainHeteroComputationThread al.Name="<<al.Name<<endl;

    // 	//if we have specified deamination rates	
    // 	if(specifiedDeam){
    // 	    if( al.IsPaired() ) continue; //skip paired-end reads because we cannot get proper deamination
    // 	    string rg;
    // 	    if(al.HasTag("RG")) {
    // 		if(!al.GetTag("RG",rg) ){ 
    // 		    cerr << "Heterozygosity computation: could not retrieve the RG tag from  "<<al.Name<<" for BAM file:" << bamFileToOpen <<" "<<currentChunk->rangeGen.getStartCoord()<<endl;
    // 		    exit(1);	
    // 		}
    // 	    }else{
    // 		rg="NA";
    // 	    }

    // 	    //The goal of the code below is to avoid having 
    // 	    //fragments for which we cannot quantify misincorporations
    // 	    //due to deamination at a given position
    // 	    if(rg2info.find(rg) == rg2info.end()){
    // 		cerr << "Heterozygosity computation: found an unknown RG tag from  "<<al.Name<<" for BAM file:" << bamFileToOpen <<" "<<currentChunk->rangeGen.getStartCoord()<<endl;
    // 		exit(1);	
    // 	    }else{
    // 		if(!rg2info[rg].isPe){ //if single end
    // 		    if(al.Length == rg2info[rg].maxReadLength){//probably reached the end of the read length, we will keep maxlength-1 and under		
    // 			//cerr<<"name "<<al.Name<<"\t"<<rg<<endl;
    // 			continue;
    // 		    }
    // 		}else{
    // 		    //since we skipped the paired reads, if we have reached here
    // 		    //we can add single-end fragments 		    
    // 		}
    // 	    }
    // 	}

    // 	if(al.Length>=int(MINLENGTHFRAGMENT) &&
    // 	   al.Length<=int(MAXLENGTHFRAGMENT) ){	   
    // 	    pileup.AddAlignment(al);
    // 	}// else{
    // 	//     cout<<"removed "<<al.Name<<endl;
    // 	// }
    // }

    // //clean up
    // pileup.Flush();
    // reader.Close();
    //fastaReference.Close();
    if(verbose){
	rc = pthread_mutex_lock(&mutexCERR);
	checkResults("pthread_mutex_lock()\n", rc);
	
	//cerr<<"Thread #"<<rankThread <<" done reading "<<thousandSeparator(hv->getTotalBases())<<" bases on "<<thousandSeparator(hv->getTotalSites())<<" sites average coverage: "<<(hv->getTotalSites()!=0 ? stringify(double(hv->getTotalBases())/double(hv->getTotalSites())) : "no sites found" )<<endl;
	cerr<<"Thread #"<<rankThread <<" done reading "<<thousandSeparator(totalBasesL)<<" bases on "<<thousandSeparator(totalSitesL)<<" sites average coverage: "<<(totalSitesL!=0 ? stringify(double(totalBasesL)/double(totalSitesL)) : "no sites found" )<<endl;
	
	rc = pthread_mutex_unlock(&mutexCERR);
	checkResults("pthread_mutex_unlock()\n", rc);
    }

//     if(totalSitesL == 0){
// 	delete piForGenomicWindow;

// #ifdef HETVERBOSE
// 	cerr<<"Thread #"<<rankThread <<" is re-starting"<<endl;
// #endif
// 	goto checkqueue;	   
//     }



    dataToWrite->hetEstResults = computeLL(piForGenomicWindow,
					   dataToWrite->vecPositionResults,
					   rankThread);   
    delete piForGenomicWindow;
	
    

    //call computeLL here

#ifdef HETVERBOSE
    cerr<<"Thread #"<<rankThread <<" is done with computations"<<endl;
#endif

    //////////////////////////////////////////////////////////////
    //                END   COMPUTATION                         //
    //////////////////////////////////////////////////////////////
	

    //COUNTERS
    rc = pthread_mutex_lock(&mutexQueueToWrite);
    checkResults("pthread_mutex_lock()\n", rc);
    
    delete currentChunk;    //delete old chunk

#ifdef DEBUGQUEUEDATATOWRITE
    //cout<<"range1 "<<currentChunk->rangeGen<<endl;
    cerr<<"dataToWrite range: "<<dataToWrite->rangeGen<<"rank: "<<dataToWrite->rank<<" #sites: "<<totalSitesL<<endl;
    //dataToWrite->refID         =  refID;

#endif

    queueDataTowrite.push(dataToWrite);
    if(usesuggestH){
	if(queueDataTowrite.size()<10){
	    sumhsuggest += dataToWrite->hetEstResults.h;
	    nhsuggest  ++;
	}
	    
	if(queueDataTowrite.size()==10){
	    sumhsuggest += dataToWrite->hetEstResults.h;
	    nhsuggest  ++;
	    usesuggestH=true;
	    suggestHval = (sumhsuggest / double(nhsuggest));
	}

    }
    rc = pthread_mutex_unlock(&mutexQueueToWrite);
    checkResults("pthread_mutex_unlock()\n", rc);

#ifdef HETVERBOSE
    cerr<<"Thread #"<<rankThread <<" is re-starting"<<endl;
#endif
    goto checkqueue;	   


    

#ifdef HETVERBOSE    
    cerr<<"Thread "<<rankThread<<" ended "<<endl;
#endif

    return NULL;

}// end mainHeteroComputationThread






queue< DataChunk * >  randomSubQueue(const queue< DataChunk * > queueDataToSubsample,unsigned int sizeToReturn,vector<GenomicRange> * rangesSelected){

    if( sizeToReturn > queueDataToSubsample.size()){
	cerr<<"Cannot subsample the queue to the size required"<<endl;
	exit(1);
    }

    queue< DataChunk *  >  toReturn=queueDataToSubsample;
    vector< DataChunk * >  myvectortemp;

    while(!toReturn.empty()){    
	DataChunk *dc = toReturn.front();
 	toReturn.pop();
	myvectortemp.push_back(dc);
    }

    random_shuffle ( myvectortemp.begin(), myvectortemp.end() );
    
    for(unsigned int i=0;i<sizeToReturn;i++){
	toReturn.push(myvectortemp[i]);
	rangesSelected->push_back(myvectortemp[i]->rangeGen);
    }

    return toReturn;
} // end randomSubQueue



queue< DataChunk * >  subFirstElemsQueue(const queue< DataChunk * > queueDataToSubsample,unsigned int sizeToReturn,vector<GenomicRange> * rangesSelected){

    if( sizeToReturn > queueDataToSubsample.size()){
	cerr<<"Cannot subsample the queue to the size required"<<endl;
	exit(1);
    }

    queue< DataChunk *  >  toReturn=queueDataToSubsample;
    vector< DataChunk * >  myvectortemp;

    while(!toReturn.empty()){    
	DataChunk *dc = toReturn.front();
 	toReturn.pop();
	myvectortemp.push_back(dc);
    }

    //random_shuffle ( myvectortemp.begin(), myvectortemp.end() );
    
    for(unsigned int i=0;i<sizeToReturn;i++){
	toReturn.push(myvectortemp[i]);
	rangesSelected->push_back(myvectortemp[i]->rangeGen);
    }

    return toReturn;
} // end randomSubQueue











//////////////////////////////////////////////////////////////////////
//                                                                  //
//                                                                  //
//                  BEGIN COVERAGE COMPUTATION                      //
//                                                                  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

// class coverageComputeVisitor : public PileupVisitor {
  
// public:
//     coverageComputeVisitor(const RefVector& references,unsigned int leftCoord, unsigned int rightCoord)
// 	: PileupVisitor()
// 	, m_references(references)
// 	, m_leftCoord(leftCoord)
// 	, m_rightCoord(rightCoord)
//     { 
// 	totalBases=0;
// 	totalSites=0;
	
//     }
//     ~coverageComputeVisitor(void) {}
  
//     // PileupVisitor interface implementation

    
//     void Visit(const PileupPosition& pileupData) {   
// 	//bool foundOneFragment=false;
// 	if(pileupData.Position < int(m_leftCoord)   || 
// 	   pileupData.Position > int(m_rightCoord) ){
// 	    return ;
// 	}
// 	//cout<<m_leftCoord<<"\t"<<m_rightCoord<<"\t"<<pileupData.Position<<endl;

// 	for(unsigned int i=0;i<pileupData.PileupAlignments.size();i++){
// 	    //skip reads that were QC failed
// 	    if(  pileupData.PileupAlignments[i].Alignment.IsFailedQC() ){
// 		continue;
// 	    }

// 	    if( pileupData.PileupAlignments[i].IsCurrentDeletion &&
// 	    	pileupData.PileupAlignments[i].IsNextInsertion ){
// 	    	continue;
// 	    }
// 	    //foundOneFragment=true;		  
// 	    totalBases++;
// 	}

// 	//if(foundOneFragment)
// 	totalSites++;
//     }
    
//     unsigned int getTotalBases() const{
//     	return totalBases;
//     }

//     unsigned int getTotalSites() const{
//     	return totalSites;
//     }

// private:
//     RefVector m_references;
//     //Fasta * m_fastaReference;
//     unsigned int m_leftCoord;
//     unsigned int m_rightCoord;
//     unsigned int totalBases;
//     unsigned int totalSites;

    
// };//end coverageComputeVisitor



static int read_bamCOV(void *data, bam1_t *b){ // read level filters better go here to avoid pileup
    //cerr<<"read_bamCOV"<<endl;
    aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
    int ret;
    while (1){
        ret = aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b);
        if ( ret<0 ) break;
	// int32_t qlen   = bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b));
	// cerr<<"qlen "<<qlen<<" "<<b->core.l_qseq<<endl;
	int32_t qlen   = b->core.l_qseq;

	if(qlen != b->core.l_qseq ){
	    cerr<<"name :"<<bam1_qname(b)<<endl;
	    exit(1);
	}
	int32_t isize  = b->core.n_cigar;
	
        if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;

	if( bam_is_paired( b) &&    //if a reads is paired and not properly paired, skip
	    !bam_is_properpair(b) ){ 
	    continue; 
	}

        if ( (int)b->core.qual < aux->min_mapQ ) continue;
        if ( aux->min_len && qlen < aux->min_len ) continue;
	uint8_t *rgptr = bam_aux_get(b, "RG");
	// cerr<<"rg1 "<<rgptr<<endl;
	// cout<<"isize "<<isize<<endl;
	string rg="UNKNOWN";
	if(rgptr){
	    rg = string( (const char*)(rgptr+1));
	}
	//cerr<<"rg2 "<<rg<<endl;	    

	if(rg2info.find(rg) == rg2info.end()){
	    rgInfo toadd;
	    toadd.isPe          = bam_is_paired( b); //al.IsPaired();
	    toadd.maxReadLength = MIN2( MAX2(qlen,isize), 255);
	    rg2info[rg]         = toadd;
	}else{
	    rg2info[rg].isPe          = rg2info[rg].isPe || bam_is_paired(b); //al.IsPaired();
	    rg2info[rg].maxReadLength = MIN2( MAX2( MAX2(qlen,isize), rg2info[rg].maxReadLength ), 255);
	}
	    
	// }else{
	//     //cerr<<"WARNING: Could not get RG tag from read "<<bam1_qname(b)<<endl;
	    
	// }
	if(specifiedDeam){
	    if( bam_is_paired( b) ) continue; //skip paired-end reads because we cannot get proper deamination
	}

        break;
    }
    return ret;
}


void *mainCoverageComputationThread(void * argc){

    int   rc;

#ifdef COVERAGETVERBOSE    
    int rankThread=0;
#endif

    rc = pthread_mutex_lock(&mutexRank);
    checkResults("pthread_mutex_lock()\n", rc);

    threadID2Rank[*(int *)pthread_self()]  = threadID2Rank.size()+1;
#ifdef COVERAGETVERBOSE    
    rankThread = threadID2Rank[*(int *)pthread_self()];
#endif
    
    rc = pthread_mutex_unlock(&mutexRank);
    checkResults("pthread_mutex_unlock()\n", rc);

 checkqueue:    
    // stackIndex=-1;
    //check stack

    
    rc = pthread_mutex_lock(&mutexQueueToRead);
    checkResults("pthread_mutex_lock()\n", rc);


    bool foundData=false;
#ifdef COVERAGETVERBOSE
    cerr<<"Thread coverage #"<<rankThread <<" started and is requesting data"<<endl;
#endif

    DataChunk * currentChunk;


    if(!queueDataForCoverage.empty()){    
 	foundData=true;
 	currentChunk = queueDataForCoverage.front();
 	queueDataForCoverage.pop();
#ifdef COVERAGETVERBOSE
 	cerr<<"Thread #"<<rankThread<<" is reading "<<currentChunk->rank<<endl;
#endif
	//cout<<"rank "<< &(currentChunk->dataToProcess) <<endl;
    }

    
  

    if(!foundData){
	rc = pthread_mutex_unlock(&mutexQueueToRead);
	checkResults("pthread_mutex_unlock()\n", rc);


	if(readDataDone){
#ifdef COVERAGETVERBOSE
	    cerr<<"Thread #"<<rankThread<<" is done"<<endl;
#endif
	    return NULL;	
	}else{

#ifdef COVERAGETVERBOSE
	    cerr<<"Thread #"<<rankThread<<" sleeping for "<<timeThreadSleep<<endl;
#endif
	    sleep(timeThreadSleep);
	    goto checkqueue;
	}

    }else{
	//release stack
	rc = pthread_mutex_unlock(&mutexQueueToRead);
	checkResults("pthread_mutex_unlock()\n", rc);
    }

    //////////////////////////////////////////////////////////////////
    //                BEGIN COV COMPUTATION                         //
    /////////////////////////////////////////////////////////////////

    //cout<<currentChunk->rangeGen<<endl;
#ifdef COVERAGETVERBOSE
    cerr<<"Thread #"<<rankThread<<" is reading "<<currentChunk->rangeGen<<endl;
#endif
    //sleep(10);


    















    unsigned int totalBasesL=0;
    unsigned int totalSitesL=0;
    

    string bamfilename = bamFileToOpen;
    string region      = currentChunk->rangeGen.getChrName()+":"+stringify(currentChunk->rangeGen.getStartCoord())+"-"+stringify(currentChunk->rangeGen.getEndCoord());

    faidx_t *fai = fai_load(fastaFile.c_str());
   
    if ( !fai ) {
         cerr << "ERROR: failed to open fasta file " <<fastaFile<<" and fasta index " << fastaIndex<<""<<endl;
         exit(1);
    }
    int seq_len;
    char *seq = fai_fetch(fai, region.c_str(), &seq_len);
    if ( seq_len < 0 ) {
	cerr << "ERROR: failed to retrieve sequence from fasta file " <<fastaFile<<" and fasta index " << fastaIndex<<" and region "<<region<<endl;
         exit(1);
    }


    
    if(seq_len != int(   currentChunk->rangeGen.getEndCoord() - currentChunk->rangeGen.getStartCoord() +1) ){
	cerr << "ERROR: failed to retrieve sequence from  file " <<fastaFile<<" and fasta index " << fastaIndex<<" asked for "<<int(   currentChunk->rangeGen.getEndCoord() - currentChunk->rangeGen.getStartCoord()+1 )<<" "<<seq_len<<endl;
	exit(1);	
    }

    int i, n, tid, reg_tid, beg, end, pos, *n_plp, baseQ = 0, mapQ = 0, min_len = MINLENGTHFRAGMENT;

    int all = 1; //status = EXIT_SUCCESS, 

    int max_depth = MAXCOV;
    const bam_pileup1_t **plp;
    bool reg = !region.empty();

    void *bed = 0; // BED data structure
    if(!bedfilename.empty()){
	bed = bed_read(bedfilename.c_str()); // BED or position list file can be parsed now
    }

    bam_hdr_t *h = NULL; // BAM header of the 1st input
    aux_t **data;
    bam_mplp_t mplp;
    int last_pos = -1, last_tid = -1, ret;

    n =1;
    
    data = (aux_t **)calloc(1, sizeof(aux_t*)); // data[i] for the i-th input
    reg_tid = 0; beg = 0; end = INT_MAX;  // set the default region

    for (i = 0; i < n; ++i) {
	//cerr<<"i ="<<i<<endl;
	//for (i = 0; i < 1; ++i) {
        int rf;
        data[i] = (aux_t *)calloc(1, sizeof(aux_t));


	data[i]->fp = sam_open_format(bamfilename.c_str(), "r", NULL); // open BAM

        if (data[i]->fp == NULL) {
	    cerr<<"ERROR: Could not open BAM file "<<bamfilename<<""<<endl;
	    exit(1);

        }
        rf = SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR | SAM_SEQ;
        if (baseQ) rf |= SAM_QUAL;

        data[i]->min_mapQ = mapQ;                    // set the mapQ filter
        data[i]->min_len  = min_len;                 // set the qlen filter
        data[i]->hdr = sam_hdr_read(data[i]->fp);    // read the BAM header
        
	if (data[i]->hdr == NULL) {
	    cerr<<"ERROR: Could not read header for bamfile "<<bamfilename<<""<<endl;
	    exit(1);
        }
        if (reg) { // if a region is specified
	    hts_idx_t *idx = sam_index_load(data[i]->fp, bamfilename.c_str());  // load the index
            if (idx == NULL) {
                //print_error("depth", "can't load index for \"%s\"", argv[optind+i]);
		cerr<<"ERROR: Cannot load index for bamfile "<<bamfilename<<""<<endl;
		exit(1);
            }else{
		//cerr<<"index ok"<<endl;
	    }
            data[i]->iter = sam_itr_querys(idx, data[i]->hdr, region.c_str()); // set the iterator
            hts_idx_destroy(idx); // the index is not needed any more; free the memory
            if (data[i]->iter == NULL) {
		cerr<<"ERROR: Cannot parse region \""<<region<<"\""<<endl;
		exit(1);
            }else{
		//cerr<<"region ok"<<endl;
	    }
        }
    }

    h = data[0]->hdr; // easy access to the header of the 1st BAM
    if (reg) {
        beg     = data[0]->iter->beg; // and to the parsed region coordinates
        end     = data[0]->iter->end;
        reg_tid = data[0]->iter->tid;
    }

    // the core multi-pileup loop
    mplp = bam_mplp_init(n, read_bamCOV, (void**)data); // initialization
    if (0 < max_depth)
        bam_mplp_set_maxcnt(mplp,max_depth);  // set maximum coverage depth
    else if (!max_depth)
        bam_mplp_set_maxcnt(mplp,INT_MAX);
    n_plp = (int *)calloc(n, sizeof(int)); // n_plp[i] is the number of covering reads from the i-th BAM
    plp   = (const bam_pileup1_t **)calloc(n, sizeof(bam_pileup1_t*)); // plp[i] points to the array of covering reads (internal in mplp)
    //int seqptr=0;
    
    while ((ret=bam_mplp_auto(mplp, &tid, &pos, n_plp, plp)) > 0) { // come to the next covered position	
	//cerr<<endl<<(pos+1)<<" "<<" -----------------------"<<endl;
	
	if (pos < beg || pos >= end) continue; // out of range; skip
        if (tid >= h->n_targets) continue;     // diff number of @SQ lines per file?
	
	char refc = toupper(seq[ pos-currentChunk->rangeGen.getStartCoord()+1 ]);
	if(refc == 'N'){//skip unresolved
	    continue;
	}

	
	//cerr<<endl<<(pos+1)<<" "<<refc<<" -----------------------"<<endl;

        if (all) {
            while (tid > last_tid) {
                if (last_tid >= 0 && !reg) {
                    // Deal with remainder or entirety of last tid.
                    while (++last_pos < int(h->target_len[last_tid])) {
			//cerr<<(last_pos+1)<<endl;
			//exit(1);
			// // Horribly inefficient, but the bed API is an obfuscated black box.
                        if (bed && bed_overlap(bed, h->target_name[last_tid], last_pos, last_pos + 1) == 0)
                            continue;
			totalSitesL++;
                    }
                }
                last_tid++;
                last_pos = -1;
                if (all < 2)
                    break;
            }

            // Deal with missing portion of current tid
            while (++last_pos < pos) {
                if (last_pos < beg) continue; // out of range; skip
                if (bed && bed_overlap(bed, h->target_name[tid], last_pos, last_pos + 1) == 0)
                    continue;
		totalSitesL++;
            }

            last_tid = tid;
            last_pos = pos;
        }
	if (bed && bed_overlap(bed, h->target_name[tid], pos, pos + 1) == 0) continue;

	//cout<<(pos+1)<<" -----------------------"<<endl;
	totalSitesL++;

	//fputs(h->target_name[tid], stdout); 
	//printf("\t%d\t", pos+1); // a customized printf() would be faster
	//cerr<<(pos+1)<<endl;
        for (i = 0; i < n; ++i) { // base level filters have to go here
            int j;
	    int m = 0;//amount of invalid bases at site j that need to be removed from coverage calculations

	    
	    //cout<<"-------------"<<endl;
	    
            for (j = 0; j < n_plp[i]; ++j) {
                const bam_pileup1_t *p = plp[i] + j; // DON'T modfity plp[][] unless you really know
		//base level filters go here
		if (p->is_del || p->is_refskip) 
		    ++m; // having dels or refskips at tid:pos
                else 
		    if(bam_get_qual(p->b)[p->qpos] < baseQ) 
			++m; // low base quality
		
#ifdef NOTDEF
		// printf("%d,",p->b);
		// printf("%d,",bam_get_seq(p->b));
		// @discussion Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G,
		//  8 for T and 15 for N. Two bases are packed in one byte with the base
		//  at the higher 4 bits having smaller coordinate on the read. It is
		//  recommended to use bam_seqi() macro to get the base.
		//  */
		//coord
		//printf("%d,",p->qpos);
		cerr<<p->qpos<<",";
		//flag
		//printf("%d,",(p->b)->core.flag);
		cerr<<((p->b)->core.flag)<<",";
		//base
		//printf("%d,",bam_seqi(bam_get_seq(p->b),p->qpos));
		//printf("%c,",alphabet[ bam_seqi(bam_get_seq(p->b),p->qpos) ]);
		cerr<<alphabet[ bam_seqi(bam_get_seq(p->b),p->qpos) ]<<",";
		
		//int c = seq_nt16_int[bam_seqi(seq, p->qpos)];
		//qual
		//printf("%d,",bam_get_qual(p->b)[p->qpos]);
		cerr<<int(bam_get_qual(p->b)[p->qpos])<<",";
		//strand
		//printf("%d,",bam_is_rev(p->b));
		cerr<<bam_is_rev(p->b)<<",";
		
		//paired
		//printf("%d,",bam_is_paired(p->b));
		cerr<<bam_is_paired(p->b)<<",";
		cerr<<"D="<<p->is_del<<",RS="<<p->is_refskip<<",";
		//fail
		//printf("%d-",bam_is_failed(p->b));
		cerr<<bam_is_failed(p->b)<<"-";
#endif
            }
            //printf("\tc=%d", n_plp[i] - m); // this the depth to output
	    //cout<<(n_plp[i] - m); // this the depth to output
	    totalBasesL += (n_plp[i] - m); // this the depth to output
	    // totalSitesL ++ ;

        }
        //putchar('\n');
	//cout<<endl;
    }
    
    if (ret < 0){ //status = EXIT_FAILURE;
	cerr<<"Problem parsing region:"<<region<<" in bamfile "<<bamfilename<<endl;
	exit(1);
    }
    free(n_plp); free(plp);
    bam_mplp_destroy(mplp);

    if (all) {
        // Handle terminating region
        if (last_tid < 0 && reg && all > 1) {
            last_tid = reg_tid;
            last_pos = beg-1;
        }
        while (last_tid >= 0 && last_tid < h->n_targets) {
            while (++last_pos < int(h->target_len[last_tid])) {
                if (last_pos >= end) break;
                if (bed && bed_overlap(bed, h->target_name[last_tid], last_pos, last_pos + 1) == 0)
                    continue;
		totalSitesL++;
            }
            last_tid++;
            last_pos = -1;
            if (all < 2 || reg)
                break;
        }
    }

    //depth_end:
    for (i = 0; i < n && data[i]; ++i) {
        bam_hdr_destroy(data[i]->hdr);
        if (data[i]->fp) sam_close(data[i]->fp);
        hts_itr_destroy(data[i]->iter);
        free(data[i]);
    }
    free(data); 
    free(seq);
    fai_destroy(fai);

    
    //free(reg);
    if (bed) bed_destroy(bed);

	
#ifdef COVERAGETVERBOSE    
    cerr<<"Thread #"<<rankThread <<" is done with computations"<<endl;
#endif
    //////////////////////////////////////////////////////////////
    //                END   COMPUTATION                         //
    //////////////////////////////////////////////////////////////
	

    //COUNTERS
    rc = pthread_mutex_lock(&mutexQueueToWrite);
    checkResults("pthread_mutex_lock()\n", rc);
    
//queueDataTowrite.push(currentChunk);
    totalBasesSum+=totalBasesL;
    totalSitesSum+=totalSitesL;

    rc = pthread_mutex_unlock(&mutexQueueToWrite);
    checkResults("pthread_mutex_unlock()\n", rc);

#ifdef COVERAGETVERBOSE    
    cerr<<"Thread #"<<rankThread <<" is re-starting"<<endl;
#endif

    goto checkqueue;	   


    

#ifdef COVERAGETVERBOSE        
    cerr<<"Thread "<<rankThread<<" ended "<<endl;
#endif

    return NULL;

} // mainCoverageComputationThread



















































































































// useminmidmax 0:min 1:mid 2:max

hmmRes runHMM(const string & outFilePrefix, const    vector<emissionUndef> & heteroEstResults, const int maxChains , const double fracChainsBurnin,const long double rohmu,const unsigned char useminmidmax, bool noROH){

    cerr<<"Creating HMM";
    // for(unsigned int i=0;i<heteroEstResults.size();i++){

    //   cout<<"h:"<<    heteroEstResults[i].h<<endl;
    //   cout<<"hlow:"<<    heteroEstResults[i].hlow<<endl;
    //   cout<<"hhigh:"<<    heteroEstResults[i].hhigh<<endl;
    //   cout<<"errb:"<<    heteroEstResults[i].errb<<endl;
    //   cout<<"weight:"<<        heteroEstResults[i].weight    <<endl;
    //   cout<<"undef:"<<    heteroEstResults[i].undef<<endl;
    //   cout<<"chrBreak:"<<    heteroEstResults[i].chrBreak<<endl;
    //   cout<<"rangeGen:"<<    heteroEstResults[i].rangeGen<<endl;
    // }


    
    if(useminmidmax == HMMCODEMIN) { cerr<<" for min. values"; }
    if(useminmidmax == HMMCODEMID) { cerr<<" for mid  values"; }
    if(useminmidmax == HMMCODEMAX) { cerr<<" for max. values"; }

    cerr<<"..";
    //write chains to output
    //bgzipWriterMCMC
    string headerHMMMCMC = "#llik\th\ts\tp\taccepted\tchains\tacptrate\n";
    //Internal::BgzfStream  bgzipWriterMCMC;
    BGZF *  bgzipWriterMCMC=NULL;


    string outFileSuffixMCMC;
    if(useminmidmax == HMMCODEMIN) { outFileSuffixMCMC = ".min"; }
    if(useminmidmax == HMMCODEMID) { outFileSuffixMCMC = ".mid"; }
    if(useminmidmax == HMMCODEMAX) { outFileSuffixMCMC = ".max"; }
    outFileSuffixMCMC += ".hmmmcmc.gz";

    bgzipWriterMCMC = bgzf_open(string(outFilePrefix+outFileSuffixMCMC).c_str(), "w");
    if (bgzipWriterMCMC == NULL) { // region invalid or reference name not found
	cerr<<"Cannot open file "<<(outFilePrefix+outFileSuffixMCMC)<<" in bgzip writer"<<endl;
	exit(1);	
    }

    // bgzipWriterMCMC.Open(outFilePrefix+outFileSuffixMCMC, IBamIODevice::WriteOnly);
    // if(!bgzipWriterMCMC.IsOpen()){
    // 	cerr<<"Cannot open file "<<(outFilePrefix+outFileSuffixMCMC)<<" in bgzip writer"<<endl;
    // 	exit(1);
    // }

    //computing the min/max segsites per chunk
    int minSegSitesPerChunk;
    int maxSegSitesPerChunk;
    if(sizeChunk == 1000000){
	minSegSitesPerChunk = minSegSitesPer1M;	
	maxSegSitesPerChunk = maxSegSitesPer1M;
    }else{
	minSegSitesPerChunk = int( (double(minSegSitesPer1M)/double(1000000))*double(sizeChunk) );
	maxSegSitesPerChunk = int( (double(maxSegSitesPer1M)/double(1000000))*double(sizeChunk) );
    }
    minSegSitesPerChunk=0;//to allow for ROHs
    //cerr<<"hmm "<<minSegSitesPerChunk<<" "<<maxSegSitesPerChunk<<endl;
    Hmm hmm (minSegSitesPerChunk,maxSegSitesPerChunk,sizeChunk,1000);
    
    //cerr<<"Begin running MCMC on HMM"<<endl;
    // cerr<<"generating a random set"<<endl;
    // vector<emission>       eTest = hmm.generateStates(250,sizeChunk);

    cerr<<".";
    //vector<emissionUndef>  eTestUndef;


    //BEGIN MCMC chain
    //init to random settings
    long double partition= (long double)(stepHMM);
    int accept=0;

    //likelihood
    long double x_i    ;
    long double x_i_1  ;
    
    //het rate
    //int sitesPer1M     = randomInt(200,1000);// between 1 and 1000 sites per million
    int sitesPer1M     = 1000;// 1000 heterozygous sites per million
    long double h_i    = double(       sitesPer1M )/double(1000000);
    long double h_i_1;

    long double hlower = double( minSegSitesPer1M )/double(1000000);
    long double hupper = double( maxSegSitesPer1M )/double(1000000);
    
    //number of non-recombining chunks per window
    long double s_i    = 100.0; //starting value, there are 100 non-recombining window in sizeChunk
    long double s_i_1;
    // cout<<"rhoWindow "<<rhoWindow<<endl;
    // for(int i=0;i<200;i++){
    // 	long double prho = s_i/(s_i+rho1M);
    // 	cout<<i<<"\t"<<prho<<"\t"<<gsl_ran_negative_binomial_pdf(i, prho ,s_i)<<endl;
	
    // }
    // exit(1);
    long double slower = 3.0;       // the whole sizeChunk is a non-recombining window
    long double supper = sizeChunk; // 1 base is a non-recombining window
    

    
    //transition rate
    long double pTlowerFirstHalf  =    0.5;//overestimate the transition probability 
    //long double pTlowerSecondHalf =    numeric_limits<double>::epsilon();
    long double pTlowerSecondHalf =    numeric_limits<double>::min();
	
    long double pTlower = pTlowerFirstHalf;	
    //long double pTupper = 1.0 - numeric_limits<double>::epsilon();
    long double pTupper = 1.0 - numeric_limits<double>::min();

    //generate a random initial probability between pTlower and pTupper
    //probability of transition
    //randomProb(false): srand() was already called once in main(), letting
    //libgab re-seed here would undo the seed set by --seed
    long double pT_i = randomProb(false)*(pTupper-pTlower) + pTlower;
    long double pT_i_1;
    //long double pTlower =       numeric_limits<double>::epsilon();
    if(	noROH ){
	pTlower  = 0.0;
	pT_i     = 0.0;
	pT_i_1   = 0.0;	
    }

    //the proposal distributions below draw from their own engine; seeding it
    //from rand() rather than a random_device ties it to the seed set in main()
    default_random_engine dre ( (default_random_engine::result_type)(rand()) );
    //int maxChains = 100000;
    //chain 0

    hmm.setHetRateForNonROH(h_i);
    if(noROH){
	hmm.setTransprob(0.0);
    }else{
	hmm.setTransprob(pT_i);
    }
    hmm.setNrwPerSizeChunk( (unsigned int)s_i );
    hmm.recomputeProbsNonROH();

    if(noROH){
	hmm.setHetRateForROH(0.0);	
    }else{
	hmm.setHetRateForROH(rohmu);
    }
    hmm.recomputeProbsROH();
    


    
    
    // cerr<<".";
    // for(unsigned int i=0;i<heteroEstResults.size();i++){
    //  	cerr<<"obs#"<<i<<" "<<heteroEstResults[i].chrBreak<<"\t"<<heteroEstResults[i].undef<<"\t"<<heteroEstResults[i].h<<"\t"<<heteroEstResults[i].hlow<<"\t"<<heteroEstResults[i].hhigh<<"\t"<<heteroEstResults[i].weight<<"\t"<<heteroEstResults[i].sites<<endl;
    // }
    
    // return 1;
    //x_i    =  forwardProb(&hmm, emittedH , sizeChunk);
    //cerr<<"test fwd"<<endl;
    fbreturnVal  tmpResFWD = forwardProbMissing(&hmm,heteroEstResults,sizeChunk,useminmidmax,verbose);
    x_i  = tmpResFWD.llik;
    // x_i    =  forwardProbUncertaintyMissing(&hmm, heteroEstResults , sizeChunk);
    // x_i    =  backwardProbUncertaintyMissing(&hmm, heteroEstResults , sizeChunk);
    // cerr<<"..done"<<endl;
    // cerr<<"x_i "<<x_i<<endl;
    // exit(1);
    //return 1;
    
    //cout<<setprecision(10)<<"\tinitial\t"<<h_i<<"\t"<<pT_i<<"\t"<<x_i<<"\t"<<endl;
    vector<long double> hvector;
    vector<long double> pvector;
    vector<long double> svector;
    
    cerr<<"..begin running MCMC on HMM using "<<thousandSeparator(maxChains)<<" chains.."<<endl;

    for(int chain=1;chain<=maxChains;chain++){

	//computing new state
	normal_distribution<long double> distribution_h(h_i,  (hupper-hlower)/partition  );
	h_i_1      = distribution_h(dre);
	if(h_i_1 <= hlower       ||  h_i_1 >= hupper     ){
	    h_i_1      = h_i;
	}


	normal_distribution<long double> distribution_pT(pT_i,(pTupper-pTlower)/partition  );
	if( !noROH ){
	    pT_i_1      = distribution_pT(dre);
	    if(pT_i_1 <= pTlower     ||  pT_i_1 >= pTupper     ){
		pT_i_1      = pT_i;
	    }
	}

	//normal_distribution<long double> distribution_s(s_i,  (supper-slower)/partition  );
	normal_distribution<long double> distribution_s(s_i,  2  );
	s_i_1      = distribution_s(dre);
	// long double prior_S;
	// long double prho = rho1M/(rho1M+rho1M);


	if(s_i_1 <= slower       ||  s_i_1 >= supper     ){
	    ///prior_S = gsl_ran_negative_binomial_pdf(s_i_1, prho , rho1M);
	    s_i_1      = s_i;
	}else{
	    //prior_S = gsl_ran_negative_binomial_pdf(s_i_1, prho , rho1M);
	}
	//prior_S = MAX2( numeric_limits<double>::min(), prior_S);
	//cout<<s_i<<"\t"<<s_i_1<<"\t"<<prior_S<<"\t"<<prho<<"\t"<<rho1M<<endl;

	//set new model parameters
	// cerr<<"old   h_i   "<<h_i<<" "<<" pT_i "<<" "<<pT_i<<" s_i  "<<s_i <<endl;
	// cerr<<"new   h_i_1 "<<h_i_1<<" "<<" pT_i_1 "<<" "<<pT_i_1<<" s_i_1  "<<s_i_1 <<endl;

	hmm.setHetRateForNonROH(h_i_1);
	hmm.setTransprob(pT_i_1);
	hmm.setNrwPerSizeChunk( (unsigned int)s_i_1 );
	hmm.recomputeProbsNonROH();
	hmm.recomputeProbsROH();

	
	//x_i_1=forwardProb(&hmm, emittedH , sizeChunk);
	//compute new likelihood
	//x_i_1=forwardProbUncertaintyMissing(&hmm, heteroEstResults , sizeChunk);
	// tmpResFWD = forwardProbUncertaintyMissing(&hmm,heteroEstResults, sizeChunk,useminmidmax);
	tmpResFWD = forwardProbMissing(&hmm,heteroEstResults, sizeChunk,useminmidmax);
 
	
	x_i_1     = tmpResFWD.llik; //+logl(prior_S);
	//cout<<prior_S<<" "<<logl(prior_S)<<endl;
	//	cerr<<"x_i_1 "<<x_i_1<<endl;
	if(chain>(maxChains/4)){
	    pTlower = pTlowerSecondHalf;
	}

	long double acceptance = MIN2( (long double)(1.0)  , expl(x_i_1-x_i) );
	if( (long double)(randomProb(false)) < acceptance){
	    h_i           =  h_i_1;
	    s_i           =  s_i_1;
	    pT_i          =  pT_i_1;
	    x_i           =  x_i_1;

	    if( (chain>=(maxChains*(1.0-fracChainsBurnin)))){
		hvector.push_back(h_i);
		pvector.push_back(pT_i);
		svector.push_back(s_i);

		string strToWrite = stringify( x_i )+"\t"+stringify( h_i )+"\t"+stringify( s_i )+"\t"+stringify( pT_i )+"\t"+stringify( accept )+"\t"+stringify( chain )+"\t"+stringify( double(accept)/double(chain) )+"\n";
		//bgzipWriterMCMC.Write(strToWrite.c_str(), strToWrite.size());
		if(bgzf_write(bgzipWriterMCMC,strToWrite.c_str(),strToWrite.size()) != int(strToWrite.size())){  cerr<<"Cannot write to bgzip stream"<<endl;   exit(1);  }
	    }
	    
	    accept++;
	    if(verbose)
		cout<<setprecision(10)<<"accepted jump from\t"<<h_i<<"\t"<<s_i<<"\t"<<pT_i<<"\t"<<x_i<<"\tto\t"<<h_i_1<<"\t"<<s_i_1<<"\t"<<pT_i_1<<"\t"<<x_i_1<<""<<"\t"<<acceptance<<" "<<accept<<" "<<chain<<" "<<double(accept)/double(chain)<<endl;
	    //cerr<<setprecision(10)<<"mcmc"<<mcmc<<"\taccepted\t"<<h_i<<"\t"<<pT_i<<"\t"<<x_i<<"\t"<<" "<<acceptance<<" "<<accept<<" "<<chain<<" "<<double(accept)/double(chain)<<endl;	    
	}else{
	  
	  if(verbose)
	      cout<<setprecision(10)<<"rejected jump from\t"<<h_i<<"\t"<<s_i<<"\t"<<pT_i<<"\t"<<x_i<<"\tto\t"<<h_i_1<<"\t"<<s_i_1<<"\t"<<pT_i_1<<"\t"<<x_i_1<<""<<"\t"<<acceptance<<" "<<accept<<" "<<chain<<" "<<double(accept)/double(chain)<<endl;	    
	}
	
	if(!verbose)
	  printprogressBarCerr( float(chain)/float(maxChains) );
	//chain++;
	//sleep(0.1);
    }
    cerr<<endl;
    cout<<setprecision(10)<<"mcmc"<<"\tfinal\t"<<h_i<<"\t"<<s_i<<"\t"<<pT_i<<"\t"<<x_i<<"\t"<<endl;

    if(bgzf_close(bgzipWriterMCMC) != 0 ){   cerr<<"Cannot close bgzip stream"<<endl;   exit(1);   }  
    //bgzipWriterMCMC.Close();    
	
	

    //hvector and pvector contain the values
    long double hSum=0.0;
    long double pSum=0.0;
    long double sSum=0.0;
    
    long double hAvg=0.0;
    long double pAvg=0.0;
    long double sAvg=0.0;
    
    long double hMin=1.0;    
    long double pMin=1.0;    
    long double sMin=supper+1;
	
    long double hMax=0.0;
    long double pMax=0.0;
    long double sMax=slower-1;

    

    
    for(unsigned int i=0;i<hvector.size();i++){
	hSum += hvector[i];
	pSum += pvector[i];
	sSum += svector[i];

	if(hvector[i] < hMin)
	    hMin = hvector[i] ;
	if(hvector[i] > hMax)
	    hMax = hvector[i] ;
	
	if(pvector[i] < pMin)
	    pMin = pvector[i] ;
	if(pvector[i] > pMax)
	    pMax = pvector[i] ;
	
	if(svector[i] < sMin)
	    sMin = svector[i] ;
	if(svector[i] > sMax)
	    sMax = svector[i] ;
	
    }

    hAvg = hSum/( (long double)hvector.size() );
    pAvg = pSum/( (long double)pvector.size() );
    sAvg = sSum/( (long double)svector.size() );


    //to remove
    // hMin = 0.00070;
    // hAvg = 0.0007467025205; 
    // hMax = 0.0008;

    // computing assignment prob
    //set average parameters
    //TODO remove
    // hAvg = 0.0007467025205;
    // pAvg = 0.0960255;
    
    //TODO remove
    //verbose=true;

    hmm.setHetRateForNonROH(hAvg);
    hmm.setTransprob(pAvg);
    hmm.setNrwPerSizeChunk( (unsigned int)sAvg );
    hmm.recomputeProbsNonROH(verbose);    
    hmm.recomputeProbsROH(verbose);

   
    
    fbreturnVal postprob = forwardBackwardProbMissing(&hmm, heteroEstResults , sizeChunk,useminmidmax,verbose);
    
    cerr<<"...HMM done"<<endl;
    //    exit(1);

    string outFileSuffixHMMpost;

    
    if(useminmidmax == HMMCODEMIN) { outFileSuffixHMMpost = ".min"; }
    if(useminmidmax == HMMCODEMID) { outFileSuffixHMMpost = ".mid"; }
    if(useminmidmax == HMMCODEMAX) { outFileSuffixHMMpost = ".max"; }
    string outFileSuffixHMMrohl = outFileSuffixHMMpost+".hmmrohl.gz";
    outFileSuffixHMMpost += ".hmmp.gz";


    

    string headerHMMpost =   "#CHROM\tBEGIN\tEND\tp[ROH]\tp[nonROH]\n";
    //Internal::BgzfStream  bgzipWriterHMMpost;
    BGZF * bgzipWriterHMMpost=NULL;
    BGZF * bgzipWriterHMMrohl=NULL;
    
    // bgzipWriterHMMpost.Open(outFilePrefix+outFileSuffixHMMpost, IBamIODevice::WriteOnly);
    // if(!bgzipWriterHMMpost.IsOpen()){
    // 	cerr<<"Cannot open file "<<(outFilePrefix+outFileSuffixHMMpost)<<" in bgzip writer"<<endl;
    // 	exit(1);
    // }
    if(	!noROH ){
	bgzipWriterHMMpost = bgzf_open(string(outFilePrefix+outFileSuffixHMMpost).c_str(), "w");
	if (bgzipWriterHMMpost == NULL) { // region invalid or reference name not found
	    cerr<<"Cannot open file "<<(outFilePrefix+outFileSuffixHMMpost)<<" in bgzip writer"<<endl;
	    exit(1);	
	}

	bgzipWriterHMMrohl = bgzf_open(string(outFilePrefix+outFileSuffixHMMrohl).c_str(), "w");
	if (bgzipWriterHMMrohl == NULL) { // region invalid or reference name not found
	    cerr<<"Cannot open file "<<(outFilePrefix+outFileSuffixHMMrohl)<<" in bgzip writer"<<endl;
	    exit(1);	
	}
	if(bgzf_write(bgzipWriterHMMpost,headerHMMpost.c_str(),headerHMMpost.size()) != int(headerHMMpost.size())){ cerr<<"Cannot write to bgzip stream"<<endl;   exit(1);  }
    }
    
    //bgzipWriterHMMpost.Write(headerHMMpost.c_str(), headerHMMpost.size());

    
    
    uint64_t rohSegments   =0;
    uint64_t nonrohSegments=0;
    uint64_t unsureSegments=0;

    unsigned int rohSegmentsContiguousCount=1;
    uint64_t rohSegmentsContiguousSum=0;
    uint64_t rohSegmentsContiguousSitesSum=0;
    string   rohSegmentsContiguousChr="";
    unsigned int rohSegmentsContiguousStart=0;
    unsigned int rohSegmentsContiguousEnd=0;
    vector<uint64_t> rohSegmentsContiguous;

    bool inROH=false;
    string strToWriterohl="#";
    
    for(unsigned int c=0;c<heteroEstResults.size();c++){
	string strToWrite=heteroEstResults[c].rangeGen.asBed()+"\t";//+"\t"+stringify(dataToWrite->hetEstResults.sites)+"\t";

	//The bp tally below has to use the actual span of this window rather than
	//sizeChunk: the last window of every chr/scaffold is shorter than sizeChunk and,
	//when --bed is used, the windows can be of any length while sizeChunk stays at
	//its default.
	const uint64_t lengthWindow = uint64_t( heteroEstResults[c].rangeGen.getEndCoord() -
						heteroEstResults[c].rangeGen.getStartCoord() + 1 );

	if(heteroEstResults[c].chrBreak){
	    if(inROH && c!=0){//was already in ROH
	      rohSegmentsContiguous.push_back(rohSegmentsContiguousSum);
	      if(strToWriterohl=="#"){   strToWriterohl = stringify( rohSegmentsContiguousCount )+"\t"+stringify( rohSegmentsContiguousChr )+"\t"+stringify( rohSegmentsContiguousStart )+"\t"+stringify( rohSegmentsContiguousEnd )+"\t"+stringify( rohSegmentsContiguousSum )+"\t"+stringify( rohSegmentsContiguousSitesSum ); }else{  strToWriterohl += "\n"+ stringify( rohSegmentsContiguousCount )+"\t"+stringify( rohSegmentsContiguousChr )+"\t"+stringify( rohSegmentsContiguousStart )+"\t"+stringify( rohSegmentsContiguousEnd )+"\t"+stringify( rohSegmentsContiguousSum )+"\t"+stringify( rohSegmentsContiguousSitesSum );  }
	      rohSegmentsContiguousSum=0;//reset
	      rohSegmentsContiguousSitesSum=0;
	      rohSegmentsContiguousChr="";
	      rohSegmentsContiguousStart=0;
	      rohSegmentsContiguousEnd=0;
	      rohSegmentsContiguousCount++;
	    }
	    inROH=false;
	}

	if(heteroEstResults[c].undef){ 
	    strToWrite+="NA\tNA\n";
	    
	    if(inROH && c!=0){//was already in ROH	      
		rohSegmentsContiguous.push_back(rohSegmentsContiguousSum);
		if(strToWriterohl=="#"){  strToWriterohl =  stringify( rohSegmentsContiguousCount )+"\t"+stringify( rohSegmentsContiguousChr )+"\t"+stringify( rohSegmentsContiguousStart )+"\t"+stringify( rohSegmentsContiguousEnd )+"\t"+stringify( rohSegmentsContiguousSum )+"\t"+stringify( rohSegmentsContiguousSitesSum ); }else{ strToWriterohl += "\n"+ stringify( rohSegmentsContiguousCount )+"\t"+stringify( rohSegmentsContiguousChr )+"\t"+stringify( rohSegmentsContiguousStart )+"\t"+stringify( rohSegmentsContiguousEnd )+"\t"+stringify( rohSegmentsContiguousSum )+"\t"+stringify( rohSegmentsContiguousSitesSum );  }
		rohSegmentsContiguousSum=0;//reset
		rohSegmentsContiguousSitesSum=0;
		rohSegmentsContiguousChr="";
		rohSegmentsContiguousStart=0;
		rohSegmentsContiguousEnd=0;
		rohSegmentsContiguousCount++;
	    }
	    inROH=false;

	}else{
	    strToWrite+=stringify( exp(postprob.m[0][c]) )+"\t"+stringify( exp(postprob.m[1][c]) )+"\n";
	}

	if(     exp(postprob.m[0][c]) > 0.9){//confident in ROH
	    rohSegments        += lengthWindow;


	    if(inROH){//was already in ROH
		rohSegmentsContiguousSum += lengthWindow;
		rohSegmentsContiguousSitesSum += heteroEstResults[c].sites;
		if(rohSegmentsContiguousChr      != heteroEstResults[c].rangeGen.getChrName()){
		    cerr<<"The current chr "<<heteroEstResults[c].rangeGen.getChrName()<<" is different than the one in the same ROH:"<< rohSegmentsContiguousChr<<" please email the programmer with your .hEst file."<<endl;
		}
		//rohSegmentsContiguousStart    = heteroEstResults[c].rangeGen.getStartCoord(); start stays as is
		rohSegmentsContiguousEnd      = heteroEstResults[c].rangeGen.getEndCoord();
	    }else{//new ROH
		rohSegmentsContiguousSum      = lengthWindow;
		rohSegmentsContiguousSitesSum = heteroEstResults[c].sites;
		rohSegmentsContiguousChr      = heteroEstResults[c].rangeGen.getChrName();
		rohSegmentsContiguousStart    = heteroEstResults[c].rangeGen.getStartCoord();
		rohSegmentsContiguousEnd      = heteroEstResults[c].rangeGen.getEndCoord();
	    }
	    inROH=true;

	}else{
	    if( exp(postprob.m[1][c]) > 0.9){//confident not in ROH
		nonrohSegments += lengthWindow;
		
		if(inROH){//was already in ROH

		  rohSegmentsContiguous.push_back(rohSegmentsContiguousSum);
		  if(strToWriterohl=="#"){  strToWriterohl =  stringify( rohSegmentsContiguousCount )+"\t"+stringify( rohSegmentsContiguousChr )+"\t"+stringify( rohSegmentsContiguousStart )+"\t"+stringify( rohSegmentsContiguousEnd )+"\t"+stringify( rohSegmentsContiguousSum )+"\t"+stringify( rohSegmentsContiguousSitesSum ); }else{ strToWriterohl += "\n"+ stringify( rohSegmentsContiguousCount )+"\t"+stringify( rohSegmentsContiguousChr )+"\t"+stringify( rohSegmentsContiguousStart )+"\t"+stringify( rohSegmentsContiguousEnd )+"\t"+stringify( rohSegmentsContiguousSum )+"\t"+stringify( rohSegmentsContiguousSitesSum );  }
		  rohSegmentsContiguousSum=0;//reset
		  rohSegmentsContiguousSitesSum=0;
		  rohSegmentsContiguousChr="";
		  rohSegmentsContiguousStart=0;
		  rohSegmentsContiguousEnd=0;
		  rohSegmentsContiguousCount++;
		}else{//do nothing
		    rohSegmentsContiguousSum  = 0;//superfluous
		    rohSegmentsContiguousSitesSum=0;
		    rohSegmentsContiguousChr="";
		    rohSegmentsContiguousStart=0;
		    rohSegmentsContiguousEnd=0;		    
		}
		inROH=false;

	    }else{
		unsureSegments += lengthWindow;
	    }
	}
	//bgzipWriterHMMpost.Write(strToWrite.c_str(), strToWrite.size());
	if( !noROH ){
	    if( bgzf_write(bgzipWriterHMMpost,strToWrite.c_str(),strToWrite.size()) != int(strToWrite.size())){	cerr<<"Cannot write to bgzip stream"<<endl;   exit(1);  }
	}
    }

    if(inROH ){//was already in ROH
	rohSegmentsContiguous.push_back(rohSegmentsContiguousSum);
	if(strToWriterohl=="#"){  strToWriterohl =  stringify( rohSegmentsContiguousCount )+"\t"+stringify( rohSegmentsContiguousChr )+"\t"+stringify( rohSegmentsContiguousStart )+"\t"+stringify( rohSegmentsContiguousEnd )+"\t"+stringify( rohSegmentsContiguousSum )+"\t"+stringify( rohSegmentsContiguousSitesSum ); }else{ strToWriterohl += "\n"+ stringify( rohSegmentsContiguousCount )+"\t"+stringify( rohSegmentsContiguousChr )+"\t"+stringify( rohSegmentsContiguousStart )+"\t"+stringify( rohSegmentsContiguousEnd )+"\t"+stringify( rohSegmentsContiguousSum )+"\t"+stringify( rohSegmentsContiguousSitesSum );  }
    }

    if(strToWriterohl=="#"){  strToWriterohl = "none found\n"; }else{ strToWriterohl = "#ROH_ID\tCHROM\tBEGIN\tEND\tROH_LENGTH\tVALIDATED_SITES\n"+strToWriterohl+"\n"; }
    
    if(	!noROH ){
	if( bgzf_write(bgzipWriterHMMrohl,strToWriterohl.c_str(),strToWriterohl.size()) != int(strToWriterohl.size())){	cerr<<"Cannot write to bgzip stream"<<endl;   exit(1);  }

	if(bgzf_close(bgzipWriterHMMpost) != 0 ){   cerr<<"Cannot close bgzip stream"<<endl;   exit(1);   }  
	if(bgzf_close(bgzipWriterHMMrohl) != 0 ){   cerr<<"Cannot close bgzip stream"<<endl;   exit(1);   }  

	cerr<<"Written trace to "<<(outFilePrefix+outFileSuffixHMMpost)<<endl;
	cerr<<"Written ROH lengths to "<<(outFilePrefix+outFileSuffixHMMrohl)<<endl;
    }
    hmmRes toreturn;

    toreturn.hAvg = hAvg;
    toreturn.pAvg = pAvg;
    toreturn.sAvg = sAvg;
    
    toreturn.hMin = hMin;
    toreturn.pMin = pMin;
    toreturn.sMin = sMin;
    
    toreturn.hMax = hMax;
    toreturn.pMax = pMax;
    toreturn.sMax = sMax;


    toreturn.rohSegments    = rohSegments   ;
    toreturn.nonrohSegments = nonrohSegments;
    toreturn.unsureSegments = unsureSegments;
    
    double rohSegmentsContiguousS=0;
    double rohSegmentsContiguousN=0;

    for(unsigned int i=0;i<rohSegmentsContiguous.size();i++){
      rohSegmentsContiguousS += double(  rohSegmentsContiguous[i] );
      rohSegmentsContiguousN += 1.0;
    }

    if(rohSegmentsContiguousN!=0)
      toreturn.avgLengthROHSegments =     (rohSegmentsContiguousS / rohSegmentsContiguousN);
    else
      toreturn.avgLengthROHSegments =     0.0;

    toreturn.postprob       = postprob;

    
    return toreturn;
}






//! Main method
/*!
  The main:
    calls initScores(), 
    captures the arguments
    reads the deamination profiles
*/

int main (int argc, char *argv[]) {
    setlocale(LC_ALL, "POSIX");
    
    

    // 	    //cout<<"\tobs="<<eTest[i].idx<<"\t"<<eTest[i].p <<"\tpredicted="<< hp.seq[i] << endl;
    // }
    // cout<<returnGitHubVersion(string(argv[0]),"..")<<endl;
    // return 1;
    
    ////////////////////////////////////
    //        Variables               //
    ////////////////////////////////////


    string cwdProg=getCWD(argv[0]);
    string deam5pfreqE  = getFullPath(cwdProg+"../deaminationProfile/none.prof");
    string deam3pfreqE  = getFullPath(cwdProg+"../deaminationProfile/none.prof");
    //string illuminafreq = getFullPath(cwdProg+"../illuminaProf/null.prof");
    string illuminafreq = getFullPath(cwdProg+"../illuminaProf/error.prof");
    
    string dnafreqFile  = getFullPath(cwdProg+"../DNAprof/default_30_20_20_30.freq");

    // cout<<deam5pfreqE<<endl;
    // cout<<deam3pfreqE<<endl;

    // return 1;

    //no contaminant deamination for now
    // string deam5pfreqC = getCWD(argv[0])+"deaminationProfile/none.prof";
    // string deam3pfreqC = getCWD(argv[0])+"deaminationProfile/none.prof";

    // vector<substitutionRates>    deam5PsubE;
    // vector<substitutionRates>    deam3PsubE;
    // vector<substitutionRates>    deam5PsubC;
    // vector<substitutionRates>    deam3PsubC;


    int    numberOfThreads   = 1;
    string outFilePrefix;
    bool   outFilePrefixFlag=false;

    string sampleName        = "sample";
    // bool   useVCFoutput      = false;
    int lastOpt=1;

    double lambdaCov=0;
    bool   lambdaCovSpecified=false;
    TStoTVratio=2.1;
    // string genoFileAsInput    ="";
    // bool   genoFileAsInputFlag=false;
    bool wroteEverything=false;
    int lastWrittenChunk=-1;   
    string headerHest="#CHROM\tBEGIN\tEND\tVALIDSITES\th\terr\thLow\thHigh\n";    
    //Internal::BgzfStream  bgzipWriterInfo;
    BGZF *  bgzipWriterInfo=NULL;
    
    string stringinfo ;
    vector<emissionUndef> heteroEstResults;
    //Internal::BgzfStream  bgzipWriterHest;
    BGZF * bgzipWriterHest=NULL;


    
    unsigned int      rank  = 0;
    int          lastRank   =-1;
    vector<GenomicRange> v ;
    vector<GenomicRange> allGenomic ;

    GenomicWindows     rw ;
    int    bpToExtract;
    int                   rc;
    // RefVector  references;
    // BamReader reader;
    aux_t *data   =NULL;//bam reader
    hts_idx_t *idx=NULL; //bam index

    
    string previousChrWritten="###";
    bool skipToHMM =false;
    bool skipTheHMM=false;
    bool noROH=false;


    ifstream myFileFAI;
    //string filenameFAI;
    string headerVCFFile;
    //    Internal::BgzfStream  bgzipWriterGL;
    BGZF *  bgzipWriterGL=NULL;
    int maxChains   =  50000;
    double fracChainsBurnin   =  0.1;

    //seed for every source of randomness in rohan (the MCMC mostly); taken
    //from the wall clock unless --seed is given, and reported in the summary
    //so that a given run can be reproduced
    unsigned int randomSeed;
    {
	timeval timeForSeed;
	gettimeofday(&timeForSeed, NULL);
	randomSeed = (unsigned int)( (timeForSeed.tv_sec * 1000) + (timeForSeed.tv_usec / 1000) );
    }

    string bedFile;

    string autosomeFile;
    bool ignoreExistingRGINFO=false;
    bool shuffleWindCoverage=true;


    string argvProg;
    size_t posSlash;

    ////////////////////////////////////
    // BEGIN Parsing arguments        //
    ////////////////////////////////////
    
    long double rohmu=0.00001;
    
    const string usage=string("\nThis program co-estimates heterozygosity rates (theta) and large runs of homozygosity\n")+
	"for modern and ancient samples\n\n"+
	string(argv[0])+                        
	" [options] [fasta file] [bam file]  "+"\n\n"+
	"\twhere:\n"+
	"\t\t\t\t\t[fasta file]\t\tThe fasta file used for alignement\n"
	"\t\t\t\t\t[bam file]\t\tThe aligned and indexed BAM file\n"+
	"\n\n"
	
	"\n\tI/O options:\n"+
	"\t\t"+"-o"+","+"--out"      + "\t\t"     + "[out prefix]" +"\t\t"+"Output prefix  (default: none)"+"\n"+
	"\t\t"+""  +"" +"--name"     + "\t\t\t"   + "[name]"       +"\t\t\t"+"Sample name (default: "+sampleName+")"+"\n"+
	//"\t\t"+""  +""+"--vcf"     + "\t\t\t" +    ""          +"\t\t\t"+"Use VCF as output format (default: "+booleanAsString(useVCFoutput)+")"+"\n"+
	//"\t\t"+""+"\t"+"--ingeno"  + "\t\t"   +    "[infile]" +"\t\t"+"Read likelihoods in BGZIP and start comp. from there (default: none)"+"\n"+
	"\t\t"+"-v"+","+"--verbose"  +"\t\t"      + ""             +"\t\t\t"+"Print extensive info  (default: "+booleanAsString(verbose)+")"+"\n"+  
	//"\t\t"+"-f"+","+""           +"\t\t"      + ""             +"\t\t\t"+"Overwrite any .rginfo.gz (default: "+booleanAsString(ignoreExistingRGINFO)+")"+"\n"+  
	//"\t\t"+""+""+"--nogl"       +"\t\t\t"      + ""             +"\t\t\t"+"Do not output genotype likelihoods  (default: "+booleanAsString(!outputgenol)+")"+"\n"+  
	"\t\t"+""+""+"--gl"       +"\t\t\t"      + ""             +"\t\t\t"+"Produce genotype likelihoods  (default: "+booleanAsString(outputgenol)+")"+"\n"+  
			      
	"\n\tComputation options:\n"+	
	"\t\t"+"-t"+"" +""           +"\t\t\t"    + "[threads]" +"\t\t"+"Number of threads to use (default: "+stringify(numberOfThreads)+")"+"\n"+
	"\t\t"+""  +"" +"--phred64"  +"\t\t\t"    + ""          +"\t\t"+"Use PHRED 64 as the offset for QC scores (default : PHRED33)"+"\n"+
	"\t\t"+""  +"" +"--size"     +"\t\t\t"    + "[bp]"      +"\t\t\t"+"Size of windows in bp  (default: "+thousandSeparator(sizeChunk)+")"+"\n"+	      
	"\t\t"+""  +"" +"--bed"      +"\t\t\t"    + "[bed file]"+"\t\t"+"Do not automatically generate genomic windows, use the regions in this bed file  (default: none)"+"\n"+
	"\t\t"+""  +"" +"     "      +"\t\t\t"    + "          "+"\t\t"+"by default, the windows are automatically generated using the --size parameter"+"\n"+
	"\t\t"+""  +"" +"     "      +"\t\t\t"    + "          "+"\t\t"+"can be combined with --size, if --size is not given, the median length of the"+"\n"+
	"\t\t"+""  +"" +"     "      +"\t\t\t"    + "          "+"\t\t"+"windows in the bed file is used as the window size"+"\n"+

	"\t\t"+""  +"" +"--map"      +"\t\t\t"    + "[bed file]"+"\t\t"+"Use a mappability filter to filter on a per site basis  (default: none)"+"\n"+
   ///"\t\t"+""  +"" +"--first"      +"\t\t\t"    + ""+"\t\t"+"Do not shuffle the windows for coverage computations (default: "+booleanAsString(!shuffleWindCoverage)+")"+"\n"+	      
	"\t\t"+""  +"" +"--tstv"     +"\t\t\t"    + "[tstv]"  +"\t\t\t"+"Ratio of transitions to transversions  (default: "+stringify(TStoTVratio)+")"+"\n"+
	"\t\t"+""  +"" +"--tvonly"     +"\t\t\t"    + ""  +"\t\t"+"Only consider transversions  (default: "+booleanAsString(tvonly)+")"+"\n"+
	"\t\t"+""  +"" +""             +"\t\t\t"    + ""  +"\t\t\t"+"recommended for highly damaged samples where damage cannot be accurately quantified "+"\n"+
	//"\t\t"+""  +"" +""             +"\t\t\t"    + ""  +"\t\t\t"+"heterozygosity estimates will be multiplied by ([tstv]+1)"+"\n"+
	//"\t\t"+""  +"" +"--noroh"    +"\t\t\t"    + ""        +"\t\t\t"+"Do not allow any region to be flagged as ROH (default: "+booleanAsString(noROH)+")"+"\n"+	

	"\t\t"+""  +"" +"--auto"     +"\t\t\t"    + "[file]"  +"\t\t\t"+"Use only the chromosome/scaffolds in this file   (default: use every chromosome)"+"\n"+
	"\t\t"+""  +"" +""           +"\t\t\t"    + ""        +"\t\t\t"+"this is done to avoid including sex chromosomes in the calculation"+"\n"+
	"\t\t"+""  +"" +"--rohmu"    +"\t\t\t"    + "[rate]"  +"\t\t\t"+"Use this value as the expected theta in ROHs   (default: "+stringify(rohmu)+")"+"\n"+
	"\t\t"+""  +"" +""           +"\t\t\t"    + ""        +"\t\t\t"+"be careful when using this option as it can inflate the background estimate for theta"+"\n"+
	"\t\t"+""  +"" +"--cov"      +"\t\t\t"    + "mincov,maxcov"+"\t\t"+"Ignore the prior probabilities on coverage, simply use these cutoffs  (default: prior probabilities are used)"+"\n"+
	"\t\t"+""  +"" +"--seed"     +"\t\t\t"    + "[seed]"       +"\t\t\t"+"Seed for the random number generator  (default: from the wall clock)"+"\n"+
	"\t\t"+""  +"" +"     "      +"\t\t\t"    + "      "       +"\t\t\t"+"give the seed reported in [out prefix].summary.txt to reproduce a run"+"\n"+

	"\n\tPlease modify these parameters carefully and read what they mean:\n"+	      
	"\t\t"+""  +"" +"--minss1M"    +"\t\t"    + "[rate]"  +"\t\t\t"+"Minimum number of segregating sites per 1M  (default: "+stringify(minSegSitesPer1M)+")"+"\n"+
	"\t\t"+""  +"" +"--maxss1M"    +"\t\t"    + "[rate]"  +"\t\t\t"+"Maximum number of segregating sites per 1M  (default: "+stringify(maxSegSitesPer1M)+")"+"\n"+
	//			      "\t\t"+""  +""+"--lambda"     +"\t\t"    + "[lambda]" +"\t\t"+"Skip coverage computation, specify lambda manually  (default: "+booleanAsString(lambdaCovSpecified)+")"+"\n"+	      
	"\n\t\tHMM:\n"+

	"\t\t"+""  +"" +"--step"     +"\t\t\t"    + "[step]"  +"\t\t\t"+"Steps used for MCMC sampling (default: "+thousandSeparator(stepHMM)+")"+"\n"+ 
	"\t\t"+""  +"" +"--chains"   +"\t\t"      + "[chains]"+"\t\t"+"Number of chains for MCMC  (default: "+thousandSeparator(maxChains)+")"+"\n"+  
	"\t\t"+""  +"" +"--hmm"      +"\t\t\t"    + ""        +"\t\t\t"+"Skip the computation of local het. rates,              (default: "+stringify(skipToHMM)+")"+"\n"+  
	"\t\t"+""  +"" +""           +"\t\t\t"    + ""        +"\t\t\t"+"read the previous het. rates [out prefix].hEst.gz and skip to HMM"+"\n"+
	"\t\t"+""  +"" +"--nohmm"    +"\t\t\t"    + ""        +"\t\t\t"+"Skip the HMM to estimate global het rates (default: "+booleanAsString(skipTheHMM)+")"+"\n"+	
	"\t\t"+""  +"" +"--burnin"   +"\t\t"      + "[frac ]" +"\t\t\t"+"Fraction of the number of chains for MCMC"+"\n"+  
	"\t\t"+""  +"" +""   +"\t\t"      + "" +"\t\t\t\t"+"to use as burning  (default: "+stringify(fracChainsBurnin)+")"+"\n"+  

			      
	// "\n\tSample options:\n"+
	// "\t\t"+""  +""+"--cont"  +"\t\t\t"    +  "[cont rate:0-1]" +"\t\t"+"Present-day human contamination rate (default: "+stringify(contrate)+")"+"\n"+
	// // "\t\t"+"--phred64" +"\t\t\t\t"+"Use PHRED 64 as the offset for QC scores (default : PHRED33)"+"\n"+
	       		      
	"\n\tDeamination and error options:\n"+                                   
	"\t\t"+""  +"" +"--deam5p"   +"\t\t"      +"[.prof]"  +"\t\t\t"+"5p deamination frequency for the endogenous\n\t\t\t\t\t\t\t\t(default: "+deam5pfreqE+")"+"\n"+
	"\t\t"+""  +"" +"--deam3p"   +"\t\t"      +"[.prof]"  +"\t\t\t"+"3p deamination frequency for the endogenous\n\t\t\t\t\t\t\t\t(default: "+deam3pfreqE+")"+"\n"+
	// "\t\t"+"-deam5pc [.prof file]" +"\t\t"+"5p deamination frequency for the contaminant (default: "+deam5pfreqC+")"+"\n"+
	// "\t\t"+"-deam3pc [.prof file]" +"\t\t"+"3p deamination frequency for the contaminant (default: "+deam3pfreqC+")"+"\n"+			      
	"\t\t"+""  +"" +"--err"      +"\t\t\t"    +"[.prof]"  +"\t\t\t"+"Illumina error profile (default: "+illuminafreq+")"+"\n"+
	"\t\t"+""  +"" +"--base"     +"\t\t\t"    +"[.freq]"  +"\t\t\t"+"Frequency of DNA bases in the genome (default: "+dnafreqFile+")"+"\n"+
	
	"";

    if( (argc== 1) ||
        (argc== 2 && string(argv[1]) == "-h") ||
        (argc== 2 && string(argv[1]) == "-help") ||
        (argc== 2 && string(argv[1]) == "--help") ){
        cerr<<usage<<endl;
        return 1;
    }

    bool specifiedSIZE=false;

    for(int i=1;i<(argc);i++){ 

	//cout<<i<<"\t"<<string(argv[i])<<endl;

        if(string(argv[i])[0] != '-' ){
            lastOpt=i;
            break;
        }


	// if( string(argv[i]) == "--ingeno" ){
	//     genoFileAsInput     = string(argv[i+1]);
	//     genoFileAsInputFlag = true;
        //     i++;
        //     continue;
	// }
	
        if( string(argv[i]) == "--step"  ){
	    stepHMM=destringify<long double>(argv[i+1]);
            i++;
            continue;
        }

        if( string(argv[i]) == "--chains"  ){
	    maxChains=destringify<int>(argv[i+1]);
            i++;
            continue;
        }

        if( string(argv[i]) == "--burnin"  ){
	    fracChainsBurnin=destringify<double>(argv[i+1]);
            i++;
            continue;
        }

	if( (string(argv[i]) == "-v") || 
	    (string(argv[i]) == "--verbose") ){
	    verbose=true;
            continue;
	}
	
	if( (string(argv[i]) == "-f") ){
	    ignoreExistingRGINFO=true;
            continue;
	}

        if( string(argv[i]) == "--hmm"  ){
	    skipToHMM=true;
            continue;
        }

        if( string(argv[i]) == "--nohmm"  ){
	    skipTheHMM=true;
            continue;
        }

        if( string(argv[i]) == "--noroh"  ){
	    noROH=true;
            continue;
        }

	
        if( string(argv[i]) == "--cov"  ){
	    covcutoffs=true;
	    string sc = string(argv[i+1]);
	    vector<string> scf=allTokens(sc,',');
	    if(scf.size() != 2 ){
		cerr<<"Error: for --cov, the format is mincov,maxcov found "<<scf.size()<<" fields "<<endl;
		return 1;
	    }

	    covcutoffsmin = destringify<unsigned int>( scf[0] );
	    covcutoffsmax = destringify<unsigned int>( scf[1] );
	    i++;
            continue;
        }


        if( string(argv[i]) == "--size"  ){
	    sizeChunk=destringify<unsigned int>(argv[i+1]);	    
	    specifiedSIZE=true;
            i++;
            continue;
        }

        if( string(argv[i]) == "--rohmu"  ){
	    rohmu=destringify<long double>(argv[i+1]);
            i++;
            continue;
        }

        if( string(argv[i]) == "--seed"  ){
	  randomSeed=destringify<unsigned int>(argv[i+1]);
	  i++;
	  continue;
        }

        if( string(argv[i]) == "--minss1M"  ){
	  minSegSitesPer1M=destringify<int>(argv[i+1]);
	  i++;
	  continue;
        }

        if( string(argv[i]) == "--maxss1M"  ){
	  maxSegSitesPer1M=destringify<int>(argv[i+1]);
	  i++;
	  continue;
        }



        if( string(argv[i]) == "--auto"  ){
	    autosomeFile=string(argv[i+1]);
            i++;
            continue;
        }

        if( string(argv[i]) == "--bed"  ){
	    bedFile=string(argv[i+1]);
            i++;
            continue;
        }

        if( string(argv[i]) == "--map"  ){
	    bedfilename=string(argv[i+1]);
            i++;
            continue;
        }


        // if( string(argv[i]) == "--vcf"  ){
	//     useVCFoutput=true;
        //     continue;
        // }

        if( string(argv[i]) == "--name"  ){
            sampleName=string(argv[i+1]);
            i++;
            continue;
        }
       
        if( string(argv[i]) == "-t"  ){
            numberOfThreads=destringify<int>(argv[i+1]);
            i++;
            continue;
        }

        if( string(argv[i]) == "--cont"  ){
            contrate=destringify<long double>(argv[i+1]);
            i++;
            continue;
        }

        if( string(argv[i]) == "--alpha"  ){
            alpha=destringify<long double>(argv[i+1]);
            i++;
            continue;
        }

        if( string(argv[i]) == "--beta"  ){
            beta=destringify<long double>(argv[i+1]);
            i++;
            continue;
        }


        if( string(argv[i]) == "-o"    ||
	    string(argv[i]) == "--out" ){
            outFilePrefix=string(argv[i+1]);
	    outFilePrefixFlag=true;
            i++;
            continue;
        }


        // if( string(argv[i]) == "--nogl" ){
	//     outputgenol=false;
        //     continue;
        // }

        if( string(argv[i]) == "--gl" ){
	    outputgenol=true;
            continue;
        }
	
        if( string(argv[i]) == "--first" ){
	    shuffleWindCoverage=false;
            continue;
        }

       	    
        if(string(argv[i]) == "--phred64"  ){
            offsetQual=64;
            continue;
        }

        if(string(argv[i]) == "--lambda"  ){
	    lambdaCovSpecified=true;
            lambdaCov         =destringify<double>(argv[i+1]);
	    i++;
            continue;
        }

        if(string(argv[i]) == "--tstv"  ){
            TStoTVratio         =destringify<long double>(argv[i+1]);
	    i++;
            continue;
        }

        if(string(argv[i]) == "--tvonly"  ){
	    tvonly  = true;
            continue;
        }


        if(string(argv[i]) == "--deam5p"  ){
            deam5pfreqE=string(argv[i+1]);
            i++;
	    specifiedDeam=true;
            continue;
        }

        if(string(argv[i]) == "--deam3p"  ){
            deam3pfreqE=string(argv[i+1]);
            i++;
	    specifiedDeam=true;
            continue;
        }

	if(string(argv[i]) == "--err"  ){
	    illuminafreq=string(argv[i+1]);
	    i++;
	    continue;
	}
	
	if(string(argv[i]) == "--base"  ){
	    dnafreqFile=string(argv[i+1]);
	    i++;
	    continue;
	}


	cerr<<"Error: unknown option "<<string(argv[i])<<endl;
	return 1;
    }

    //--bed and --size can be combined: --bed gives the windows to analyze whereas --size
    //tells the HMM how large a window is. If --size is omitted, it is inferred from the
    //median length of the windows in the BED file (see below).

    pthread_t             threadCov[numberOfThreads];//coverage threads
    pthread_t             threadHet[numberOfThreads];//het      threads
    
    if(!outFilePrefixFlag){
	cerr<<"The output file has to be defined"<<endl;
	return 1;	
    }// else{

    // }

    fastaFile         = string(argv[lastOpt]);
    bamFileToOpen     = string(argv[lastOpt+1]);
    fastaIndex        = fastaFile+".fai";
    
    set<string> listAutosomes;

    if( !autosomeFile.empty()){
       
	ifstream myFileAUTO;
	myFileAUTO.open(autosomeFile.c_str(), ios::in);
	string line;
	if (myFileAUTO.is_open()){
	    while ( getline (myFileAUTO,line)){     
		listAutosomes.insert(line);
	    }
	    myFileAUTO.close();	   
	}else{
	    cerr << "Error in reading list of autosomes: Unable to open file "<<autosomeFile<<endl;
	    exit(1);
	}
	
    }


    

    if(skipToHMM){//skip het computations
	goto beginhmm;
    }
    
    cerr<<"Parsing arguments ...";

    
    if( !isFile(fastaFile) ){
	cerr<<"The fasta file "<<fastaFile<<" does not exists"<<endl;
	return 1;	
    }

    if( !isFile(fastaIndex) ){
	cerr<<"The fasta file "<<fastaFile<<"  does not have an index: "<<fastaIndex<<endl;
	return 1;	
    }

    if( !bedFile.empty() ){
	if( !isFile(bedFile) ){
	    cerr<<"The supplied BED file "<<bedFile<<"  does not exist"<<endl;
	    return 1;	
	}
    }

    if(	   (0>fracChainsBurnin) ||
	   (fracChainsBurnin>1) ){
	cerr<<"The supplied fraction of burnins  "<<fracChainsBurnin<<"  should be between 0 and 1"<<endl;
	return 1;	
    }
    //    if(outFilePrefixFlag)
    
    // if( !strEndsWith(outFilePrefix,".gz")){
    //     cerr<<"The output file "<<outFilePrefix<<" must end with .gz"<<endl;
    //     return 1;	
    // }



    //Testing BAM file

    //bam_hdr_t *h = NULL; // BAM header of the 1st input

    data = (aux_t *)calloc(1, sizeof(aux_t));


    data->fp = sam_open_format(bamFileToOpen.c_str(), "r", NULL); // open BAM

    if(data->fp == NULL) {
	cerr<<"ERROR: Could not open input BAM file "<<bamFileToOpen<<""<<endl;
	exit(1);
    }

    idx = sam_index_load(data->fp, bamFileToOpen.c_str());  // load the index
    if (idx == NULL) {
	cerr<<"ERROR: Cannot load index for bamfile "<<bamFileToOpen<<""<<endl;
	exit(1);
    }
    data->hdr = sam_hdr_read(data->fp);    // read the BAM header

    // for(int rid=0;rid<100;rid++){
    // 	cout<<rid<<endl;
    // 	cout<<data->hdr->target_name[rid]<<endl;
    // }

    
    //bam_get_tid(data->hdr, currentChunk->rangeGen.getChrName().c_str() );
    
    // if ( !reader.Open(bamFileToOpen) ) {
    // 	cerr << "Could not open input BAM file:" << bamFileToOpen <<endl;
    // 	exit(1);
    // }

    // reader.LocateIndex();

    // if(!reader.HasIndex()){
    // 	cerr << "The BAM file: " << bamFileToOpen <<" does not have an index"<<endl;
    // 	exit(1);
    // }

    // // retrieve reference data
    // references = reader.GetReferenceData();


    // reader.Close();


    
    ////////////////////////////////////
    //   END Parsing arguments        //
    ////////////////////////////////////

    //the single seeding for the whole program: every call to randomProb() uses
    //randomProb(false) so that libgab never re-seeds behind our back, and the
    //MCMC's default_random_engine is seeded from rand()
    srand( randomSeed );
    cerr<<"Random seed: "<<randomSeed<<endl;

















    ////////////////////////////////////
    // BEGIN Initializing scores      //
    ////////////////////////////////////
    initScores();
    ////////////////////////////////////
    //    END Initializing scores     //
    ////////////////////////////////////




       


    ///////////////////////////////////////
    // BEGIN read Illumina error profile //
    ///////////////////////////////////////    


    readIlluminaError(illuminafreq,illuminaErrorsProb);
    for(int b1=0;b1<4;b1++){
	for(int b2=0;b2<4;b2++){
	    int dinucIndex = b1 *4+ b2;
	    illuminaErrorsProbMatrix.p[b1][b2] = illuminaErrorsProb.s[dinucIndex];
	}
    }

#ifdef DEBUGILLUMINAFREQ
    for(unsigned int i=0;i<16;i++){
    	cerr<<i<<"\t"<<illuminaErrorsProb.s[i]<<endl;	
    }

    for(int b1=0;b1<4;b1++)
	cerr<<"ACGT"[b1]<<"\t";
    cerr<<endl;
    for(int b1=0;b1<4;b1++){
	cerr<<"ACGT"[b1]<<"\t";
	for(int b2=0;b2<4;b2++){
	    cerr<<illuminaErrorsProbMatrix.p[b1][b2]<<"\t";
	}
	cerr<<endl;
    }
    cerr<<endl;
    
    return 1;
#endif

    cerr<<"...done"<<endl;
    ///////////////////////////////////////
    // END read Illumina error profile //
    ///////////////////////////////////////    

    //    return 1;













    /////////////////////////////
    // BEGIN  Compute coverage //
    /////////////////////////////
    
    cerr<<"Computing average coverage.."<<endl;

    rc=0;


    //If the user gave an explicit list of chr/scaffolds via --auto, that list is
    //authoritative and the automatic sex chromosome exclusion is turned off. Otherwise
    //sequences named like a sex chromosome are skipped.
    rw = GenomicWindows  (fastaIndex, !autosomeFile.empty() );



    if( !bedFile.empty() ){
	v = readBEDfile(bedFile);

	//sizeChunk is the scale the HMM works at: it converts rates of heterozygosity into
	//expected counts per window and normalizes the weight given to each window. When the
	//windows come from a BED file, it therefore has to reflect the size of those windows
	//and not the default of 1Mbp, otherwise every window is down-weighted by the ratio of
	//the two and the HMM ends up driven by its prior rather than by the data. If the user
	//did not say what size to use, take the median length of the windows in the BED.
	if( !specifiedSIZE && !v.empty() ){
	    vector<unsigned int> lengthWindows;
	    for(unsigned int i=0;i<v.size();i++){
		lengthWindows.push_back( v[i].getEndCoord() - v[i].getStartCoord() + 1 );
	    }
	    sort( lengthWindows.begin(), lengthWindows.end() );
	    sizeChunk = lengthWindows[ lengthWindows.size()/2 ];
	    if(sizeChunk == 0){ sizeChunk = 1; }
	    cerr<<"Windows were read from "<<bedFile<<", using the median length of these windows ("<<thousandSeparator(sizeChunk)<<"bp) as window size, use --size to override"<<endl;
	}
    }

    bpToExtract       = sizeChunk;

    if( bedFile.empty() ){
	v = rw.getGenomicWindows(bpToExtract,0);
    }
    allGenomic = rw.getGenomicWindows(bpToExtract,0);
    
    if( !autosomeFile.empty()){

	vector<GenomicRange> newV;
	vector<GenomicRange> newAllgenomic;
	
	
	for(unsigned int i=0;i<v.size();i++){
	    if( listAutosomes.find( v[i].getChrName() ) == listAutosomes.end() ){//not found
		
	    }else{//found
		newV.push_back(          v[i] );
		newAllgenomic.push_back( v[i] );		
	    }
	}

	v=newV;
	allGenomic = newAllgenomic;
    }

    if( v.size() == 0 ){    
	cerr<<"No range found using these chr/loci settings"<<endl; 
	return 1;
    }
    
    //unsigned int sizeGenome = 0;

    for(unsigned int i=0;i<v.size();i++){
#ifdef DEBUGFIRSTWINDOWS
	if(i==DEBUGFIRSTWINDOWS)
	    break;
#endif

	//cout<<"genomic region #"<<i<<" "<<v[i]<<endl;
	DataChunk * currentChunk = new DataChunk();
	
	currentChunk->rangeGen   = v[i];
	currentChunk->rank       = rank;
	lastRank                 = rank;
	//sizeGenome             += v[i].getLength();

	
	queueDataToprocess.push(currentChunk);
	rank++;
    }



    rank  = 0;
    for(unsigned int i=0;i<allGenomic.size();i++){
#ifdef DEBUGFIRSTWINDOWS
	if(i==DEBUGFIRSTWINDOWS)
	    break;
#endif

	//cout<<"genomic region #"<<i<<" "<<v[i]<<endl;
	DataChunk * currentChunk = new DataChunk();
	
	currentChunk->rangeGen   = allGenomic[i];
	currentChunk->rank       = rank;
	//lastRank                 = rank;
	//sizeGenome             += v[i].getLength();

	
	queueDataAllGenomic.push(currentChunk);
	rank++;
    }




    readDataDone=true;
    //    return 1;

    // vector<GenoResults *> vectorGenoResults;


    if( ignoreExistingRGINFO || (!isFile(outFilePrefix+".rginfo.gz")) ){

	// if(!genoFileAsInputFlag){
	    int bpToComputeCoverage = 10000000;
	    int genomicRegionsToUse = bpToComputeCoverage/bpToExtract;
	    if( genomicRegionsToUse > int(queueDataToprocess.size())){
		genomicRegionsToUse = int(queueDataToprocess.size());
	    }
    
	    vector<GenomicRange> * rangesSelected=new vector<GenomicRange>();
	    
	    //renabled
	    if(shuffleWindCoverage){//we shuffle
		queueDataForCoverage = randomSubQueue(     queueDataToprocess, genomicRegionsToUse,rangesSelected);
	    }else{//we just take the first elements
		//queueDataForCoverage = subFirstElemsQueue( queueDataToprocess,genomicRegionsToUse);
		queueDataForCoverage = subFirstElemsQueue( queueDataAllGenomic,genomicRegionsToUse,rangesSelected);
	    }
	    //queueDataForCoverage = subFirstElemsQueue( queueDataToprocess,genomicRegionsToUse);
	    unsigned int queueDataForCoverageOrigsize = queueDataForCoverage.size();
	    //cerr<<queueDataForCoverageOrigsize<<endl;

	    pthread_mutex_init(&mutexQueueToRead,   NULL);
	    pthread_mutex_init(&mutexQueueToWrite,  NULL);
	    pthread_mutex_init(&mutexRank ,         NULL);
	    pthread_mutex_init(&mutexCERR ,         NULL);
	    pthread_mutex_init(&mutexCOUNTCOV ,     NULL);


	    for(int i=0;i<numberOfThreads;i++){
		rc = pthread_create(&threadCov[i], NULL, mainCoverageComputationThread, NULL);
		checkResults("pthread_create()\n", rc);
	    }

	    if(verbose){	    
		cerr<<"Creating threads for coverage calculation, need to process="<<queueDataForCoverage.size()<<" out of a total of "<<queueDataToprocess.size()<<endl;
	    }


	    //TODO to code?
	    // rc = pthread_mutex_lock(&mutexCERR);
	    // checkResults("pthread_mutex_lock()\n", rc);

	    // if(verbose){
	    // 	cerr<<"Thread#"<<threadID<<" starting pre-computations size of window: "<<thousandSeparator(piForGenomicWindow->size())<<endl;
	    // }
    
	    // rc = pthread_mutex_unlock(&mutexCERR);
	    // checkResults("pthread_mutex_unlock()\n", rc);

	    // pthread_mutex_t  mutexCOUNTCOV              = PTHREAD_MUTEX_INITIALIZER;
	    // unsigned int queueDataForCoverageOrigsizeSum;

	    while(queueDataForCoverage.size()!=0){
		if(verbose){	    
		    cerr<<getDateString()<<" "<<getTimeString()<<" # of slices left to process: "<<queueDataForCoverage.size()<<"/"<<queueDataForCoverageOrigsize<<endl;
		}else{
		    printprogressBarCerr( float(queueDataForCoverageOrigsize- queueDataForCoverage.size() ) / float(queueDataForCoverageOrigsize) );
		}
		sleep(timeThreadSleep);
	    }
	    printprogressBarCerr( float(queueDataForCoverageOrigsize- queueDataForCoverage.size() ) / float(queueDataForCoverageOrigsize) );//print 100%
	
	    //waiting for threads to finish
	    for (int i=0; i <numberOfThreads; ++i) {	
		rc = pthread_join(threadCov[i], NULL);
		checkResults("pthread_join()\n", rc);
	    }

	    if(!verbose)
		cerr<<endl;//flush the progress bar
	
	    cerr<<"..done"<<endl;

	    pthread_mutex_destroy(&mutexRank);
	    pthread_mutex_destroy(&mutexQueueToRead);
	    pthread_mutex_destroy(&mutexQueueToWrite);
	    pthread_mutex_destroy(&mutexCERR);
	    pthread_mutex_destroy(&mutexCOUNTCOV);



	    //    cout<<"Final" <<" "<<totalBasesSum<<"\t"<<totalSitesSum<<"\t"<<double(totalBasesSum)/double(totalSitesSum)<<endl;
	
	    //The read groups are recorded before reads are filtered out, so rg2info tells us
	    //whether the data is paired-end even when every single read has been discarded.
	    bool foundPEreadgroup=false;
	    bool foundSEreadgroup=false;
	    for(map<string,rgInfo>::iterator it=rg2info.begin(); it!=rg2info.end(); ++it){
		if(it->second.isPe){ foundPEreadgroup=true; }else{ foundSEreadgroup=true; }
	    }

	    rateForPoissonCov    = ((long double)totalBasesSum)/((long double)totalSitesSum);
	    if(totalSitesSum ==0 ){
		//By far the most common reason for ending up with no data at all is asking for
		//deamination rates on paired-end data: those reads are all skipped (see below)
		//and nothing is left. Say so instead of blaming the BAM file.
		if(specifiedDeam && foundPEreadgroup && !foundSEreadgroup){
		    cerr<<"No data was retained: every read in this BAM file is paired-end and deamination profiles were given via --deam5p/--deam3p."<<endl;
		    cerr<<"The deamination rates are indexed on the distance to the ends of the original molecule, which cannot be recovered from a read"<<endl;
		    cerr<<"that is only one half of a pair, ROHan therefore skips paired-end reads when deamination profiles are used."<<endl;
		    cerr<<"Please merge/collapse your overlapping read pairs (e.g. leeHom, AdapterRemoval --collapse, fastp -m or SeqPrep) and align the"<<endl;
		    cerr<<"collapsed reads as single-end. If you do not need the damage correction, run without --deam5p/--deam3p instead."<<endl;
		    return 1;
		}
		//cerr<<"No data was found for the entire BAM file for the region "<<endl;
		cerr<<"No data was found for the entire BAM file for the regions you defined, please verify the BAM file. "<<endl;
		cerr<<"The following regions had been selected. "<<endl;
		for(unsigned int i=0;i<rangesSelected->size();i++){
		    cerr<<(i+1)<<"\t"<<rangesSelected->at(i)<<endl;
		}
		cerr<<"If you are running tests, you can use the hidden option --first it will select the first genomic windows. "<<endl;
		return 1;
	    }

	    //Mixed data: the single-end reads got us here but the paired-end ones were all
	    //dropped, warn as this silently removes part of the data.
	    if(specifiedDeam && foundPEreadgroup){
		cerr<<"WARNING: deamination profiles were given via --deam5p/--deam3p and some read groups are paired-end."<<endl;
		cerr<<"WARNING: paired-end reads are skipped as their deamination rates cannot be placed, only the single-end reads are used."<<endl;
		cerr<<"WARNING: merge/collapse your read pairs and align them as single-end to use all of the data."<<endl;
	    }

	    cerr<<"Final computation:" <<" bases="<<totalBasesSum<<"\tsites="<<totalSitesSum<<"\tlambda coverage="<<rateForPoissonCov<<endl;
	    delete(rangesSelected);
	    //pthread_exit(NULL);	
	    //	cerr<<"Lambda coverage: " <<rateForPoissonCov<<endl;	

	    stringinfo= stringify(rateForPoissonCov)+"\n";       
	
	bgzipWriterInfo = bgzf_open( string(outFilePrefix+".rginfo.gz").c_str(), "w");
	if (bgzipWriterInfo == NULL) { // region invalid or reference name not found
	    cerr<<"Cannot open file "<<(outFilePrefix+".rginfo.gz")<<" in bgzip writer"<<endl;
	    exit(1);	
	}

	// bgzipWriterInfo.Open(outFilePrefix+".rginfo.gz", IBamIODevice::WriteOnly);
	// if(!bgzipWriterInfo.IsOpen()){
	//     cerr<<"Cannot open file "<<(outFilePrefix+".rginfo.gz")<<" in bgzip writer"<<endl;
	//     return 1;
	// }
	
	
	for(map<string,rgInfo>::iterator it=rg2info.begin(); it!=rg2info.end(); ++it){
	    string s = stringify(it->first)+"\t"+stringify(it->second.isPe)+"\t"+stringify(it->second.maxReadLength)+"\n";
	    
	    if(it->second.isPe)
		MAXLENGTHFRAGMENT   = MAX2( MAXLENGTHFRAGMENT , (unsigned int)it->second.maxReadLength );
	    else
		MAXLENGTHFRAGMENT   = MAX2( MAXLENGTHFRAGMENT , (unsigned int)(it->second.maxReadLength) );

	    MAXLENGTHFRAGMENT = MIN2( 255, MAXLENGTHFRAGMENT );//as the length is stored on 8 bits (for now)
	    
	    //cerr<<"RG:\t" <<s<<endl;
	    cerr<<"RG:\t"<<stringify(it->first)<<"\t"<<(it->second.isPe?"PE":"SE")<<"\t"<<stringify(it->second.maxReadLength)<<endl;
	    stringinfo+=s;
	}

	if(bgzf_write(bgzipWriterInfo,stringinfo.c_str(),stringinfo.size()) != int(stringinfo.size())){ cerr<<"Cannot write to bgzip stream"<<endl;   exit(1);  }
	//bgzipWriterInfo.Write(stringinfo.c_str(), stringinfo.size());
	

	//bgzipWriterInfo.Close();
	if(bgzf_close(bgzipWriterInfo) != 0 ){   cerr<<"Cannot close bgzip stream"<<endl;   exit(1);   }  

	if(!lambdaCovSpecified){

	}else{
	    rateForPoissonCov=lambdaCov;
	    cerr<<"Using user-specified lambda coverage="<<rateForPoissonCov<<endl;
	}

	//exit(1);
    }else{ //not lambdaCovSpecified
	//use the rate specified via the command line

	//rateForPoissonCov = lambdaCov;
	igzstream infoFilePreComputed;
	string filen = outFilePrefix+".rginfo.gz";
	cerr<<"Found read group file: "<<filen<<endl;

	infoFilePreComputed.open(filen.c_str(), ios::in);

	if (infoFilePreComputed.good()){
	    string l;
	    //retrieving coverage
	    getline (infoFilePreComputed,l);
	    rateForPoissonCov = destringify<long double>(l);
	    if(lambdaCovSpecified){
		rateForPoissonCov = lambdaCov;
	    }
	    cerr<<"lambda coverage="<<rateForPoissonCov<<endl;
	    while ( getline (infoFilePreComputed,l)){		
		vector<string> tempv = allTokens(l,'\t');		
		if(tempv.size() != 3){
		    cerr<<"ERROR line "<<l<<" in file "<<filen<<" does not have 3 fields, has "<<tempv.size()<<" fields instead"<<endl;
		    return 1;
		}
		
		rgInfo toadd;
		toadd.isPe          =                 (tempv[1]=="1");
		toadd.maxReadLength = destringify<int>(tempv[2]);
		if(toadd.isPe)
		    MAXLENGTHFRAGMENT   = MAX2( MAXLENGTHFRAGMENT , (unsigned int)toadd.maxReadLength );
		else
		    MAXLENGTHFRAGMENT   = MAX2( MAXLENGTHFRAGMENT , (unsigned int)(toadd.maxReadLength) );

		MAXLENGTHFRAGMENT = MIN2( 255, MAXLENGTHFRAGMENT );//as the length is stored on 8 bits (for now)
		
		rg2info[ tempv[0] ] = toadd;

		cerr<<"RG:\t"<<tempv[0]<<"\t"<<(toadd.isPe?"PE":"SE")<<"\t"<<stringify(toadd.maxReadLength)<<endl;

	    }
	    infoFilePreComputed.close();

	}else{
	    cerr << "Unable to open info file: "<<filen<<endl;
	    return 1;
	}

    }
    //    return 1;
    
    //     long double rateForPoissonCovFloor = floorl(rateForPoissonCov);
    //     long double rateForPoissonCovCeil  =  ceill(rateForPoissonCov);

    
    // #ifdef DEBUGCOV
    //     cerr<<"rateForPoissonCovFloor "<<rateForPoissonCovFloor<<endl;
    //     cerr<<"rateForPoissonCovCeil "<<rateForPoissonCovCeil<<endl;
    // #endif
    
    cov2ProbSite = new vector<long double> (MAXCOV+1,0);
    //#ifdef CORRECTCOV

    argvProg=string(argv[0]);
    posSlash=argvProg.find_last_of("/");
    if(posSlash == string::npos){
	//do nothing, leave argvProg as is
    }else{
        argvProg=argvProg.substr(posSlash+1);
    }

    populatedCoverateVector( (cwdProg+argvProg) , cov2ProbSite , rateForPoissonCov, MAXCOV );
    //#endif


    for(int i=0;i<=MAXCOV;i++){
	// long double florPPMF=poisson_pmfl( (long double)(i), rateForPoissonCovFloor)/poisson_pmfl( rateForPoissonCovFloor, rateForPoissonCovFloor);
	// long double ceilPPMF=poisson_pmfl( (long double)(i), rateForPoissonCovCeil) /poisson_pmfl( rateForPoissonCovCeil,  rateForPoissonCovCeil);

#ifdef DEBUGCOV
	// cerr<<"florPPMF "<<i<<" = "<<poisson_pmfl( (long double)(i), rateForPoissonCovFloor)<<"/"<<poisson_pmfl( rateForPoissonCovFloor, rateForPoissonCovFloor)<<" "<<florPPMF<<" w= "<<(1.0-(rateForPoissonCov-rateForPoissonCovFloor))<<endl;
	// cerr<<"ceilPPMF "<<i<<" = "<<poisson_pmfl( (long double)(i), rateForPoissonCovCeil) <<"/"<<poisson_pmfl( rateForPoissonCovCeil,  rateForPoissonCovCeil)<<" "<<ceilPPMF<<" w= "<<(1.0-(rateForPoissonCovCeil-rateForPoissonCov))<<endl;
	cerr<<"cov2ProbSite["<<i<<"] = "<<cov2ProbSite->at(i)<<endl;
#endif

    }

    //return 1;
    

    cerr<<"Results\tbp="<<totalBasesSum<<"\tsites="<<totalSitesSum<<"\tlambda="<<rateForPoissonCov<<endl;
    // cerr<<MAXLENGTHFRAGMENT<<endl;
    //return 1;

    // for(int i=0;i<20;i++){
    // 	cout<<i<<"\t"<<pdfPoisson( (long double)i, rateForPoissonCov)/pdfPoisson( rateForPoissonCov, rateForPoissonCov)<<endl;
    // }


    // for(int i=0;i<100;i++){
    // 	cout<<i<<"\t"<<pdfPoisson( (long double)i, 20)/pdfPoisson( 20, 20)<<endl;
    // }

    cerr<<"..done"<<endl;

    //return 1;

    // pthread_exit(NULL);

    ////////////////////////////
    // END   Compute coverage //
    ////////////////////////////

















    cerr<<"Begin pre-computations "<<endl;
    ////////////////////////////////////////////////////////////////////////
    //
    // BEGIN DEAMINATION PROFILE
    //
    ////////////////////////////////////////////////////////////////////////
    
    cerr<<"Begin computing deamination substitution rates."<<endl;
    
    initDeamProbabilities(deam5pfreqE,deam3pfreqE);

    cerr<<"done"<<endl;

    ////////////////////////////////////////////////////////////////////////
    //
    // END  DEAMINATION PROFILE
    //
    ////////////////////////////////////////////////////////////////////////







    ////////////////////////////////////////////////////////////////////////
    //
    // BEGIN DNA BASE FREQUENCY
    //
    ////////////////////////////////////////////////////////////////////////
    cerr<<"Begin computing default base frequencies post-deamination."<<endl;
    //initDefaultBaseFreq_(dnafreqFile);
    initDefaultBaseFreq(dnafreqFile);
    cerr<<"done"<<endl;

    ////////////////////////////////////////////////////////////////////////
    //
    // END DNA BASE FREQUENCY
    //
    ////////////////////////////////////////////////////////////////////////






    ////////////////////////////////////////////////////////////////////////
    //
    // BEGIN COMPUTE P[OBSERVED BASE|THEORITICAL BASE]
    //
    ////////////////////////////////////////////////////////////////////////

    cerr<<"Begin computing substitution probabilities."<<endl;
    initLikelihoodScores();
    cerr<<"done"<<endl;

    // exit(1);
    // initLikelihoodScores();
    // cerr<<"..done"<<endl;
    ////////////////////////////////////////////////////////////////////////
    //
    // END COMPUTE P[OBSERVED BASE|THEORITICAL BASE]
    //
    ////////////////////////////////////////////////////////////////////////

    cerr<<"pre-computations are done"<<endl;









    // doneReading=true;    

    ///////////////////////
    //  Compute hetero   //
    ///////////////////////
    
    threadID2Rank=map<unsigned int, int> ();


    cerr<<"Creating "<<numberOfThreads<<" thread(s) for heterozygosity calculation"<<endl;
    pthread_mutex_init(&mutexQueueToRead,   NULL);
    pthread_mutex_init(&mutexQueueToWrite,  NULL);
    pthread_mutex_init(&mutexRank,          NULL);
    pthread_mutex_init(&mutexCERR,          NULL);
    

    for(int i=0;i<numberOfThreads;i++){
	rc = pthread_create(&threadHet[i], NULL, mainHeteroComputationThread, NULL);
	checkResults("pthread_create()\n", rc);
    }
    //    return 1;
    // 	//threads are running here

    // unsigned int originalSize = queueDataToprocess.size();
     // while(queueDataToprocess.size()!=0){
     // 	 cout<<"# of slices left to process: "<<queueDataToprocess.size()<<"/"<<originalSize<<endl;
     // 	 sleep(timeThreadSleep);
     // }
    cerr<<"done creating threads "<<endl;
    //cerr<<"outFilePrefixFlag ="<<outFilePrefixFlag<<endl;
    ///////////////////
    //Writing data out/
    ///////////////////
    // sleep(1000000000000);
    // return 1;
    //writing h estimates

    bgzipWriterHest = bgzf_open(string(outFilePrefix+".hEst.gz").c_str(), "w");

     if (bgzipWriterHest == NULL) { // region invalid or reference name not found
	cerr<<"Cannot open file "<<(outFilePrefix+".hEst.gz")<<" in bgzip writer"<<endl;
	exit(1);	
     }
    // bgzipWriterHest.Open(outFilePrefix+".hEst.gz", IBamIODevice::WriteOnly);
    // if(!bgzipWriterHest.IsOpen()){
    // 	cerr<<"Cannot open file "<<(outFilePrefix+".hEst.gz")<<" in bgzip writer"<<endl;
    // 	return 1;
    // }

     if(bgzf_write(bgzipWriterHest,headerHest.c_str(),headerHest.size()) != int(headerHest.size())) { cerr<<"Cannot write to bgzip stream"<<endl;   exit(1);  }

     //bgzipWriterHest.Write(headerHest.c_str(), headerHest.size());

    

#ifndef DEBUGFIRSTWINDOWS
    
    
    if(outputgenol){

	bgzipWriterGL = bgzf_open(string(outFilePrefix+".vcf.gz").c_str(), "w");

	if (bgzipWriterGL == NULL) { // region invalid or reference name not found
	    cerr<<"Cannot open file "<<(outFilePrefix+".vcf.gz")<<" in bgzip writer"<<endl;
	    exit(1);	
	}

	// bgzipWriterGL.Open(outFilePrefix+".vcf.gz", IBamIODevice::WriteOnly);
	// if(!bgzipWriterGL.IsOpen()){
	//     cerr<<"Cannot open file "<<(outFilePrefix+".vcf.gz")<<" in bgzip writer"<<endl;
	//     return 1;
	// }
    
    
	// if(outFileSiteLLFlag){
	// 	bgzipWriter.Open(outFileSiteLL, IBamIODevice::WriteOnly);
	// 	if(!bgzipWriter.IsOpen()){
	// 	    cerr<<"Cannot open file "<<outFileSiteLL<<" in bgzip writer"<<endl;
	// 	    return 1;
	// 	}

    
	headerVCFFile=string("")+
	    "##fileformat=VCFv4.2\n";
	headerVCFFile+="##ROHanVersion="+returnGitHubVersion(string(argv[0]),"..")+"\n";
	headerVCFFile+="##reference="+fastaFile+"\n";
    }


    // 	if(useVCFoutput){
    // 	    headerOutFile="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+sampleName+"\n";	
    // 	}else{
    // 	    headerOutFile="#CHROM\tPOS\tA\tC\tG\tT\tGENO\tGENOS\tQualL\tCovL\tAA\tAC\tAG\tAT\tCC\tCG\tCT\tGG\tGT\tTT\n";	
    // 	}

    //filenameFAI = fastaIndex;

    myFileFAI.open(fastaIndex.c_str(), ios::in);

   if (myFileFAI.is_open()){
       string l;
       while ( getline (myFileFAI,l)){
	   vector<string> tem = allTokens(l,'\t');
	   if(tem.size()!=5){
	       cerr<<"ERROR line "<<l<<" does not have 5 fields, has "<<tem.size()<<" fields instead"<<endl;
	       return 1;
	   }
	   headerVCFFile+="##contig=<ID="+tem[0]+",length="+tem[1]+">\n";
       }
       myFileFAI.close();
   }else{
       cerr << "Unable to open fasta fai file "<<fastaIndex<<endl;
       return 1;
    }


    headerVCFFile+="##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele at this location\">\n";
    headerVCFFile+="##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Root-mean-square mapping quality of covering reads\">\n";
    headerVCFFile+="##INFO=<ID=AC4,Number=4,Type=Integer,Description=\"Count of each 4 bases\">\n";

    headerVCFFile+="##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    headerVCFFile+="##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth\">\n";
    headerVCFFile+="##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n";    
    headerVCFFile+="##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled likelihoods for genotypes\">\n";
    headerVCFFile+="##FORMAT=<ID=PL10,Number=G,Type=Integer,Description=\"Phred-scaled likelihoods for all 10 genotypes\">\n";

    headerVCFFile+="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+sampleName+"\n";	

    if(outputgenol){
	//bgzipWriterGL.Write(headerVCFFile.c_str(),headerVCFFile.size());
	if(bgzf_write(bgzipWriterGL,headerVCFFile.c_str(),headerVCFFile.size()) != int(headerVCFFile.size())) { cerr<<"Cannot write to bgzip stream"<<endl;   exit(1);  }

    }
    // }
#endif

   
    //#ifdef LATER
    //cerr<<"test wroteEverything="<<wroteEverything<<endl;
   cerr<<""<<numberOfThreads<<" thread(s) are running in the background, results will be written at it becomes available"<<endl;
   
   if(outputgenol){	
	cerr<<"Writing genotype data to:        "<<outFilePrefix+".vcf.gz" << endl;
   }
   cerr<<"Writing local het. estimates to: "<<outFilePrefix+".hEst.gz"<< endl;
   if(!verbose)
       printprogressBarCerr( 0 );

   //outFilePrefix+".hEst.gz", IBamIODevice::WriteOnly);
   
    while(!wroteEverything){

	//threads are running here
	rc = pthread_mutex_lock(&mutexQueueToWrite);
	checkResults("pthread_mutex_lock()\n", rc);
	
	bool wroteData=false;
	if(!queueDataTowrite.empty()){
	
	    DataToWrite *  dataToWrite= queueDataTowrite.top();

	    if( lastWrittenChunk == (dataToWrite->rank-1) ){ 	    //correct order
		queueDataTowrite.pop();



		
		rc = pthread_mutex_lock(&mutexCERR);
		checkResults("pthread_mutex_lock()\n", rc);

		// if(verbose){
		if(verbose)
		    cerr<<getDateString()<<" "<<getTimeString()<<" writing chunk#"<<dataToWrite->rank<<" with "<<thousandSeparator(dataToWrite->vecPositionResults->size())<<" records, "<<(queueDataTowrite.size())<<" chunk(s) left in queue"<<endl;
		
		// }
		
		rc = pthread_mutex_unlock(&mutexCERR);
		checkResults("pthread_mutex_unlock()\n", rc);

		rc = pthread_mutex_unlock(&mutexQueueToWrite);
		checkResults("pthread_mutex_unlock()\n", rc);


		
		string strToWrite=dataToWrite->rangeGen.asBed()+"\t"+stringify(dataToWrite->hetEstResults.sites)+"\t";
		emissionUndef hetResToAdd;
		
		hetResToAdd.rangeGen=dataToWrite->rangeGen;
		
		if(     previousChrWritten != dataToWrite->rangeGen.getChrName()){
		    if(previousChrWritten == "###"){//first chr
			hetResToAdd.chrBreak   = false;//not a true chr break
		    }else{//true new chr
			hetResToAdd.chrBreak   = true;			
		    }
		    previousChrWritten = dataToWrite->rangeGen.getChrName();
		}else{//same chr
		    hetResToAdd.chrBreak   = false;		
		}
		
		//hetResToAdd
		//heteroEstResults.push_back( dataToWrite->hetEstResults );
		
		if(dataToWrite->hetEstResults.hasConverged){
		    hetResToAdd.undef  = false;		
		    long double h      = dataToWrite->hetEstResults.h    ;
		    long double hLow   = dataToWrite->hetEstResults.hLow ;
		    long double hHigh  = dataToWrite->hetEstResults.hHigh;
		    long double errb   = dataToWrite->hetEstResults.errb ;
		    
		    hetResToAdd.h      = h;
		    hetResToAdd.hlow   = hLow;
		    hetResToAdd.hhigh  = hHigh;
		    
		    if(h<0)      h     = 0;
		    if(hLow<0)   hLow  = 0;
		    if(hHigh<0)  hHigh = 0;
		    
		    strToWrite+=stringify( h )+"\t"+stringify( errb )+"\t"+stringify( hLow )+"\t"+stringify( hHigh )+"\n";
		}else{
		    hetResToAdd.undef  = true;		
		    strToWrite+="NA\tNA\tNA\tNA\n";
		}
		
		hetResToAdd.weight = ( (long double)(dataToWrite->hetEstResults.sites) ) / ( (long double)(sizeChunk) );
		hetResToAdd.sites  = dataToWrite->hetEstResults.sites;

		hetResToAdd.undef  = (hetResToAdd.undef || isnan(hetResToAdd.h) || isnan(hetResToAdd.hlow) || isnan(hetResToAdd.hhigh) );
		
		heteroEstResults.push_back(hetResToAdd);
				
		// cerr<<"writing "<<strToWrite<<endl;
		//bgzipWriterHest.Write(strToWrite.c_str(), strToWrite.size());
		if(bgzf_write(bgzipWriterHest,strToWrite.c_str(),strToWrite.size()) != int(strToWrite.size())) { cerr<<"Cannot write to bgzip stream"<<endl;   exit(1);  }

		//#ifdef LATER
		//sizeGenome+=dataToWrite->vecPositionResults->size();

#ifndef DEBUGFIRSTWINDOWS
		if(outputgenol){
		    strToWrite="";
		    // cerr<<"SIZE "<<dataToWrite->vecPositionResults->size()<<endl;
		    for(unsigned int i=0;i<dataToWrite->vecPositionResults->size();i++){
			//cerr<<i<<endl;
			// cerr<<"\t"<<dataToWrite->vecPositionResults->at(i)->toString(&references,dataToWrite->refID)<<endl;

			strToWrite += dataToWrite->vecPositionResults->at(i)->toString( data->hdr,dataToWrite->refID);
			//strToWrite += dataToWrite->vecPositionResults->at(i)->toString(&references,dataToWrite->refID);
			//cout<<i<<"\t"<<dataToWrite->vecPositionResults->at(i)->toString(references);
			if( (i%500) == 499){
			    //if(outFileSiteLLFlag){
			    if(bgzf_write(bgzipWriterGL,strToWrite.c_str(),strToWrite.size())  != int(strToWrite.size())){ cerr<<"Cannot write to bgzip stream"<<endl;   exit(1);  }

			    //bgzipWriterGL.Write(strToWrite.c_str(), strToWrite.size());
			    //}
			    strToWrite="";
			}
		    
			// GenoResults * toadd =  new GenoResults( dataToWrite->vecPositionResults->at(i) );
			// vectorGenoResults.push_back(toadd);
		    }

		    if(!strToWrite.empty()){
			// if(outFileSiteLLFlag){
			// 	if(outFileSiteLLFlag){ 
			if(bgzf_write(bgzipWriterGL,strToWrite.c_str(),strToWrite.size()) != int(strToWrite.size())){ cerr<<"Cannot write to bgzip stream"<<endl;   exit(1);  }

			//bgzipWriterGL.Write(strToWrite.c_str(), strToWrite.size());
			
			// 	}
			// }
		    }
		}
		//#endif		    
#endif
		
		wroteData=true;		
		lastWrittenChunk=dataToWrite->rank;


		if(!verbose)
		    printprogressBarCerr( float(lastWrittenChunk)/float(lastRank) );
		
		if(dataToWrite->rank == lastRank)
		    wroteEverything=true;	

		delete dataToWrite;
	    }else{
		//do nothing, we have to wait for the chunk with the right rank
		rc = pthread_mutex_unlock(&mutexQueueToWrite);
		checkResults("pthread_mutex_unlock()\n", rc);	
	    }

	}else{//end if queue not empty
	    rc = pthread_mutex_unlock(&mutexQueueToWrite);
	    checkResults("pthread_mutex_unlock()\n", rc);
	}

	if(!wroteData)
	    sleep(timeSleepWrite);
    }
    //#endif

    pthread_mutex_destroy(&mutexRank);
    pthread_mutex_destroy(&mutexQueueToRead);
    pthread_mutex_destroy(&mutexQueueToWrite);
    pthread_mutex_destroy(&mutexCERR);

    ///////////////////////
    //end Writing data out/
    ///////////////////////

    if(bgzf_close(bgzipWriterHest) != 0 ){   cerr<<"Cannot close bgzip stream"<<endl;   exit(1);   }

    //bgzipWriterHest.Close();

#ifndef DEBUGFIRSTWINDOWS
    if(outputgenol){
	//bgzipWriterGL.Close();
	if(bgzf_close(bgzipWriterGL) != 0 ){   cerr<<"Cannot close bgzip stream"<<endl;   exit(1);   }  
    }
#endif

    cerr<<"done writing data"<<endl;




    
    //////////////////////////////////
    //                              //
    // BEGIN HMM                    //
    //                              //
    //////////////////////////////////

    if(skipTheHMM){
	//goto endhmm;
	delete cov2ProbSite;
	
	cerr<<"ROHan finished succesfully"<<endl;
	
	return 0;
    }



    
    
 beginhmm:
    //heteroEstResults contains the results
    if(skipToHMM){//skip het computations read from file
	string hEstFileToRead = outFilePrefix+".hEst.gz";
	cerr<<"Reading previous est. of h from "<<hEstFileToRead<<".";
	igzstream hEstFileSt;

	hEstFileSt.open(hEstFileToRead.c_str(), ios::in);
	string  previousChrWritten = "###";
	if (hEstFileSt.good()){
	    vector<string> fields;
	    string line;
	    
	    getline (hEstFileSt,line) ; 		//header
	    
	    //cerr<<line<<endl;
	    while (getline (hEstFileSt,line) ){
		if(verbose)	
		    cerr<<line<<endl;	    
		fields = allTokens(line,'\t');

		if(fields.size() != 8){
		    cerr << "ERROR: line from previous h est. does not have 9 fields ("<<(fields.size()+1)<<") "<<line<<endl;
		    return 1;
		}

	    	emissionUndef hetResToAdd;		
		hetResToAdd.rangeGen.setChrName(                                fields[0]    );
		hetResToAdd.rangeGen.setStartCoord(  destringify<unsigned int>( fields[1])+1 );
		hetResToAdd.rangeGen.setEndCoord(    destringify<unsigned int>( fields[2])   );
		unsigned int sitesDefinedLine      = destringify<unsigned int>( fields[3]    );
		
		if(     previousChrWritten != fields[0]){
		    if(previousChrWritten == "###"){//first chr
			hetResToAdd.chrBreak   = false;//not a true chr break
		    }else{//true new chr
			hetResToAdd.chrBreak   = true;			
		    }
		    previousChrWritten = fields[0];
		}else{//same chr
		    hetResToAdd.chrBreak   = false;		
		}
	    
		if(fields[4] != "NA"){ //if data is not missing

		    if(fields[5] == "-nan"){//numerical errors			
			hetResToAdd.undef  = true;  //undefined 

			string htouse = fields[4];
			if(htouse == "-nan"){
			    htouse = "0";
			}
			
			long double h      = destringify<long double>( htouse );
			long double errb   = destringify<long double>( "0" );
			long double hLow   = destringify<long double>( htouse );			
			long double hHigh  = destringify<long double>( htouse );

			hetResToAdd.h      = h;
			hetResToAdd.errb   = errb;
			
			if(hLow == 0){  //went under 0
			    hLow = h-errb;
			}
			
			hetResToAdd.hlow   = hLow;
			hetResToAdd.hhigh  = hHigh;
			
		    }else{//good to go

			hetResToAdd.undef  = false;		
			long double h      = destringify<long double>( fields[4] );
			long double errb   = destringify<long double>( fields[5] );
			long double hLow   = destringify<long double>( fields[6] );
			long double hHigh  = destringify<long double>( fields[7] );
		    
			hetResToAdd.h      = h;
			hetResToAdd.errb   = errb;
		    
			if(hLow == 0){  //went under 0
			    hLow = h-errb;
			}
		    
			hetResToAdd.hlow   = hLow;
			hetResToAdd.hhigh  = hHigh;
		    }	    
		    // if(h<0)      h     = 0;
		    // if(hLow<0)   hLow  = 0;
		    // if(hHigh<0)  hHigh = 0;		    
		}else{
		    hetResToAdd.undef  = true; //window is undefined
		
		}
		hetResToAdd.undef  = (hetResToAdd.undef || isnan(hetResToAdd.h) || isnan(hetResToAdd.hlow) || isnan(hetResToAdd.hhigh) );
		hetResToAdd.weight = ( (long double)(sitesDefinedLine) ) / ( (long double)(sizeChunk) );
		hetResToAdd.sites  = sitesDefinedLine;
		//cerr<<hetResToAdd.sites<<endl;

		if(verbose)
		    cerr<<hetResToAdd.chrBreak<<"\t"<<hetResToAdd.undef<<"\t"<<hetResToAdd.h<<"\t"<<hetResToAdd.hlow<<"\t"<<hetResToAdd.hhigh<<"\t"<<hetResToAdd.weight<<endl;
		
		heteroEstResults.push_back(hetResToAdd);
	    }           
	    hEstFileSt.close();
	}else{//not able to open
	    cerr << "Unable to open file "<<hEstFileToRead<<endl;
	    exit(1);
	}
	cerr<<"..done"<<endl;
	
    } //end if skipToHMM


    
    
    //exit(1);
    cerr<<"Running HMM again to obtain global theta outside of ROHs"<<endl;
    //lower
    noROH=false;
    hmmRes hmmResmin=runHMM(outFilePrefix,heteroEstResults,maxChains,fracChainsBurnin,rohmu,HMMCODEMIN,noROH);
    cerr<<"min h est. "<<hmmResmin.hAvg<<" hMin "<<hmmResmin.hMin<<" hMax "<<hmmResmin.hMax<<" s "<<hmmResmin.sAvg<<" sMin "<<hmmResmin.sMin<<" sMax "<<hmmResmin.sMax<<" p avg. "<<hmmResmin.pAvg<<" pMin "<<hmmResmin.pMin<<" pMax "<<hmmResmin.pMax<<" rohS "<<hmmResmin.rohSegments<<" nonrohS "<<hmmResmin.nonrohSegments<<" unsure "<<hmmResmin.unsureSegments<<endl;

    //mid
    hmmRes hmmResmid=runHMM(outFilePrefix,heteroEstResults,maxChains,fracChainsBurnin,rohmu,HMMCODEMID,noROH);
    cerr<<"mid h est. "<<hmmResmid.hAvg<<" hMin "<<hmmResmid.hMin<<" hMax "<<hmmResmid.hMax<<" s "<<hmmResmid.sAvg<<" sMin "<<hmmResmid.sMin<<" sMax "<<hmmResmid.sMax<<" p avg. "<<hmmResmid.pAvg<<" pMin "<<hmmResmid.pMin<<" pMax "<<hmmResmid.pMax<<" rohS "<<hmmResmid.rohSegments<<" nonrohS "<<hmmResmid.nonrohSegments<<" unsure "<<hmmResmid.unsureSegments<<endl;


    //upper
    hmmRes hmmResmax=runHMM(outFilePrefix,heteroEstResults,maxChains,fracChainsBurnin,rohmu,HMMCODEMAX,noROH);
    cerr<<"max h est. "<<hmmResmax.hAvg<<" hMin "<<hmmResmax.hMin<<" hMax "<<hmmResmax.hMax<<" s "<<hmmResmax.sAvg<<" sMin "<<hmmResmax.sMin<<" sMax "<<hmmResmax.sMax<<" p avg. "<<hmmResmax.pAvg<<" pMin "<<hmmResmax.pMin<<" pMax "<<hmmResmax.pMax<<" rohS "<<hmmResmax.rohSegments<<" nonrohS "<<hmmResmax.nonrohSegments<<" unsure "<<hmmResmax.unsureSegments<<endl;

    long double hAvg = MID3(  hmmResmin.hMin , hmmResmid.hMin  ,  hmmResmax.hMin );
    long double hMin = MIN3(  hmmResmin.hMin , hmmResmid.hMin  ,  hmmResmax.hMin );
    long double hMax = MAX3(  hmmResmin.hMax , hmmResmid.hMax  ,  hmmResmax.hMax );

    cerr<<"Running HMM again to obtain global theta both within and outside ROHs"<<endl;
    noROH=true;
    hmmRes hmmResminNOROH=runHMM(outFilePrefix,heteroEstResults,maxChains,fracChainsBurnin,rohmu,HMMCODEMIN,noROH);
    cerr<<"min h est. "<<hmmResminNOROH.hAvg<<" hMin "<<hmmResminNOROH.hMin<<" hMax "<<hmmResminNOROH.hMax<<" s "<<hmmResminNOROH.sAvg<<" sMin "<<hmmResminNOROH.sMin<<" sMax "<<hmmResminNOROH.sMax<<" p avg. "<<hmmResminNOROH.pAvg<<" pMin "<<hmmResminNOROH.pMin<<" pMax "<<hmmResminNOROH.pMax<<" rohS "<<hmmResminNOROH.rohSegments<<" nonrohS "<<hmmResminNOROH.nonrohSegments<<" unsure "<<hmmResminNOROH.unsureSegments<<endl;

    //mid
    hmmRes hmmResmidNOROH=runHMM(outFilePrefix,heteroEstResults,maxChains,fracChainsBurnin,rohmu,HMMCODEMID,noROH);
    cerr<<"mid h est. "<<hmmResmidNOROH.hAvg<<" hMin "<<hmmResmidNOROH.hMin<<" hMax "<<hmmResmidNOROH.hMax<<" s "<<hmmResmidNOROH.sAvg<<" sMin "<<hmmResmidNOROH.sMin<<" sMax "<<hmmResmidNOROH.sMax<<" p avg. "<<hmmResmidNOROH.pAvg<<" pMin "<<hmmResmidNOROH.pMin<<" pMax "<<hmmResmidNOROH.pMax<<" rohS "<<hmmResmidNOROH.rohSegments<<" nonrohS "<<hmmResmidNOROH.nonrohSegments<<" unsure "<<hmmResmidNOROH.unsureSegments<<endl;


    //upper
    hmmRes hmmResmaxNOROH=runHMM(outFilePrefix,heteroEstResults,maxChains,fracChainsBurnin,rohmu,HMMCODEMAX,noROH);
    cerr<<"max h est. "<<hmmResmaxNOROH.hAvg<<" hMin "<<hmmResmaxNOROH.hMin<<" hMax "<<hmmResmaxNOROH.hMax<<" s "<<hmmResmaxNOROH.sAvg<<" sMin "<<hmmResmaxNOROH.sMin<<" sMax "<<hmmResmaxNOROH.sMax<<" p avg. "<<hmmResmaxNOROH.pAvg<<" pMin "<<hmmResmaxNOROH.pMin<<" pMax "<<hmmResmaxNOROH.pMax<<" rohS "<<hmmResmaxNOROH.rohSegments<<" nonrohS "<<hmmResmaxNOROH.nonrohSegments<<" unsure "<<hmmResmaxNOROH.unsureSegments<<endl;


    //return 1;
    //////////////////////////////////
    //                              //
    //  END HMM                     //
    //                              //
    //////////////////////////////////

    long double hAvgNOROH = MID3(  hmmResminNOROH.hMin , hmmResmidNOROH.hMin  ,  hmmResmaxNOROH.hMin );
    long double hMinNOROH = MIN3(  hmmResminNOROH.hMin , hmmResmidNOROH.hMin  ,  hmmResmaxNOROH.hMin );
    long double hMaxNOROH = MAX3(  hmmResminNOROH.hMax , hmmResmidNOROH.hMax  ,  hmmResmaxNOROH.hMax );



    






    
    
    //////////////////////////////////
    //                              //
    // BEGIN h global               //
    //                              //
    //////////////////////////////////
    //cerr<<"writing plot"<<endl;

    double heightChr=45;
    vector<PdfWriter *> vecOfPdfWritersH;
    PdfWriter * pdfwriterH =new PdfWriter(outFilePrefix+".het.pdf",heightChr);

    if( pdfwriterH->drawFrame(fastaIndex,sizeChunk,0, &listAutosomes, !autosomeFile.empty()) == 1 ){
	cerr<<"ERROR writing frame to pdf file:"<<(outFilePrefix+".het.pdf")<<endl;
	return 1;
    }
    
    vecOfPdfWritersH.push_back(pdfwriterH);
    int numberPagesPDFhet =1;
    int numberOfChrPerPage = int( floor( double(pdfwriterH->getPageHeight()) / double( pdfwriterH->getTotalHeightChr()  )) );

    if(pdfwriterH->getTotalNumChrToDraw() > numberOfChrPerPage){//multipage
	numberPagesPDFhet  = int( ceil( double(pdfwriterH->getTotalNumChrToDraw())/double(numberOfChrPerPage)  ) );
	//cerr<<"numberPagesPDFhet "<<numberPagesPDFhet<<endl;
	
	vecOfPdfWritersH[0]->setFname(    outFilePrefix+".het_"+  stringify(1)+"_"+stringify(numberPagesPDFhet)+".pdf");
	//cerr<<"0 "<<(outFilePrefix+".het_"+  stringify(1)+"_"+stringify(numberPagesPDFhet)+".pdf")<<endl;
	for(int i=1;i<numberPagesPDFhet;i++){
	    PdfWriter * pdfwriterH_ =new PdfWriter(outFilePrefix+".het_"+stringify(i+1)+"_"+stringify(numberPagesPDFhet)+".pdf",heightChr);
	    //cerr<<"0 "<<(outFilePrefix+".het_"+  stringify(i+1)+"_"+stringify(numberPagesPDFhet)+".pdf")<<endl;
	    vecOfPdfWritersH.push_back(pdfwriterH_);

	    if( pdfwriterH_->drawFrame(fastaIndex,sizeChunk,i*numberOfChrPerPage, &listAutosomes, !autosomeFile.empty()) == 1 ){
		cerr<<"ERROR writing frame to pdf file:"<<(outFilePrefix+".het_"+stringify(i+1)+"_"+stringify(numberPagesPDFhet)+".pdf")<<endl;
		return 1;
	    }

	    
	}
    }
    //cerr<<numberOfChrPerPage<<endl;
    //return 1;
		 
    
    // cerr<<pdfwriterH->getPageHeight()<<endl;
    // cerr<<pdfwriterH->getPageWidth()<<endl;
    // cerr<<pdfwriterH->getHeightLabel()<<endl;
    // cerr<<pdfwriterH->getHeightChr()<<endl;
    // cerr<<pdfwriterH->getTotalNumChrToDraw()<<endl;

    // delete(pdfwriterH);
    // cerr<<numberOfChrPerPage<<endl;
    
    //exit(1);
    
    //pdfwriter.drawHorizontalLine(100,100,102);
    long double minHFoundPlotting = 0.0;
    //long double maxHFoundPlotting = numeric_limits<double>::epsilon();
    long double maxHFoundPlotting = numeric_limits<double>::min();
    
    for(unsigned int c=0;c<heteroEstResults.size();c++){
	if( heteroEstResults[c].hhigh > maxHFoundPlotting){
	    maxHFoundPlotting = heteroEstResults[c].hhigh;
	}
    }
    //cerr<<"maxHFoundPlotting "<<maxHFoundPlotting<<endl;
    //maxHFoundPlotting=0.0015;
    maxHFoundPlotting = double( maxSegSitesPer1M )/double(1000000);
    //vector<PdfWriter *> vecOfPdfWriters;
    for(unsigned int i=0;i<vecOfPdfWritersH.size();i++){
	if( vecOfPdfWritersH[i]->drawYLabels(minHFoundPlotting,maxHFoundPlotting,true) == 1 ){
	    cerr<<"ERROR writing y labels to pdf file:"<<vecOfPdfWritersH[i]->getFname()<<endl;
	    return 1;
	}
    }

    for(unsigned int c=0;c<heteroEstResults.size();c++){
	if(heteroEstResults[c].undef)
	    continue;

	for(unsigned int i=0;i<vecOfPdfWritersH.size();i++){
	    //cerr<<heteroEstResults[c].rangeGen.getChrName()<<endl;


	    if( vecOfPdfWritersH[i]->chrIspresent(heteroEstResults[c].rangeGen.getChrName())){//very inefficient..

		//cerr<<"present "<<heteroEstResults[c].rangeGen.getChrName()<<" i="<<i<<":"<<vecOfPdfWritersH[i]->chrIspresent(heteroEstResults[c].rangeGen.getChrName())<<endl;
		if(    vecOfPdfWritersH[i]->drawHEst(heteroEstResults[c].rangeGen,
						     (heteroEstResults[c].h ),
						     // ( (heteroEstResults[c].h )-6.0188e-05),
						     // ( (heteroEstResults[c].h )+6.0188e-05),				   				   
						     MAX2(heteroEstResults[c].hlow,0),
						     heteroEstResults[c].hhigh,
						     minHFoundPlotting,//0.0,//double( minSegSitesPer1M )/double(1000000),
						     maxHFoundPlotting, // 0.00500    = 4*2e-8*62500
						     sizeChunk
		)  != 0 ){
		    cerr<<"ERROR writing data point#"<<c<<" "<<heteroEstResults[c].rangeGen<<" to pdf file:"<<(vecOfPdfWritersH[i]->getFname())<<endl;
		    return 1;
		}
		//}
	    }
	}

    }
    //produce plot libharoutFilePrefix+".vcf.gz"u?


    // endhmm:
    //write out h estimate
    for(unsigned int i=0;i<vecOfPdfWritersH.size();i++){
	vecOfPdfWritersH[i]->drawGlobalHEst(hAvg,
					   hMin,
					   hMax,			     
					   double( minSegSitesPer1M )/double(1000000),
					   maxHFoundPlotting); // 0.00500    = 4*2e-8*62500

    }
    for(unsigned int i=0;i<vecOfPdfWritersH.size();i++){
	delete(vecOfPdfWritersH[i]);
    }

    //#ifdef NOTDEF    



    

    //exit(1);

    vector<PdfWriter *> vecOfPdfWritersHMM;
    PdfWriter *  pdfwriterHMM = new PdfWriter (outFilePrefix+".hmm.pdf",heightChr);
    
    if( pdfwriterHMM->drawFrame(fastaIndex,sizeChunk, 0,&listAutosomes, !autosomeFile.empty()) == 1 ){
	cerr<<"ERROR writing frame to pdf file:"<<(outFilePrefix+".hmm.pdf")<<endl;
	return 1;
    }
    vecOfPdfWritersHMM.push_back(pdfwriterHMM);
    int numberPagesPDFhetHMM =1;
    int numberOfChrPerPageHMM = int( floor( double(pdfwriterHMM->getPageHeight()) / double( pdfwriterHMM->getTotalHeightChr()  )) );

    if(pdfwriterHMM->getTotalNumChrToDraw() > numberOfChrPerPageHMM){//multipage
	numberPagesPDFhetHMM  = int( ceil( double(pdfwriterHMM->getTotalNumChrToDraw())/double(numberOfChrPerPageHMM)  ) );

	
	vecOfPdfWritersHMM[0]->setFname(    outFilePrefix+".hmm_"+  stringify(1)+"_"+stringify(numberPagesPDFhetHMM)+".pdf");

	for(int i=1;i<numberPagesPDFhetHMM;i++){
	    PdfWriter * pdfwriterHMM_ =new PdfWriter(outFilePrefix+".hmm_"+stringify(i+1)+"_"+stringify(numberPagesPDFhet)+".pdf",heightChr);

	    vecOfPdfWritersHMM.push_back(pdfwriterHMM_);

	    if( pdfwriterHMM_->drawFrame(fastaIndex,sizeChunk,i*numberOfChrPerPageHMM, &listAutosomes, !autosomeFile.empty()) == 1 ){
		cerr<<"ERROR writing frame to pdf file:"<<(outFilePrefix+".hmm_"+stringify(i+1)+"_"+stringify(numberPagesPDFhet)+".pdf")<<endl;
		return 1;
	    }
	    
	}
    }
    
    //pdfwriter.drawHorizontalLine(100,100,102);
    long double minPPlotting = 0.0;
    long double maxPPlotting = 1.0;   

    for(unsigned int i=0;i<vecOfPdfWritersHMM.size();i++){
  
	if( vecOfPdfWritersHMM[i]->drawYLabels(minPPlotting,maxPPlotting,false) == 1 ){
	    cerr<<"ERROR writing y labels to pdf file:"<<(vecOfPdfWritersHMM[i]->getFname())<<endl;
	    return 1;
	}
    }
    
    for(unsigned int c=0;c<heteroEstResults.size();c++){
    //for(unsigned int c=0;c<150;c++){
	//cerr<<c<<" plotting 0="<<exp(postprob.m[0][c])<<" 1="<<exp(postprob.m[1][c])<<endl;
	//long double pROH = expl(postprob.m[0][c]);

	// if(heteroEstResults[c].undef)
	//     continue;
	for(unsigned int i=0;i<vecOfPdfWritersHMM.size();i++){

	    if( vecOfPdfWritersHMM[i]->chrIspresent(heteroEstResults[c].rangeGen.getChrName())){//very inefficient..
		//cerr<<"present "<<heteroEstResults[c].rangeGen<<" "<<i<<endl;
		
		if(    vecOfPdfWritersHMM[i]->drawHMM(heteroEstResults[c].rangeGen,
						      exp(hmmResmin.postprob.m[1][c]),
						      exp(hmmResmin.postprob.m[1][c]),
						      exp(hmmResmin.postprob.m[1][c]),
						      minPPlotting,//0.0,//double( minSegSitesPer1M )/double(1000000),
						      maxPPlotting, // 0.00500    = 4*2e-8*62500
						      sizeChunk,
						      !heteroEstResults[c].undef,
						      heteroEstResults[c].chrBreak,						      
						      HMMCODEMIN)  != 0 ){
		    cerr<<"ERROR writing data point#"<<c<<" "<<heteroEstResults[c].rangeGen<<" to pdf file:"<<(vecOfPdfWritersHMM[i]->getFname())<<endl;
		    return 1;
		}
		
		if(  vecOfPdfWritersHMM[i]->drawHMM(heteroEstResults[c].rangeGen,
						    exp(hmmResmid.postprob.m[1][c]),
						    exp(hmmResmid.postprob.m[1][c]),
						    exp(hmmResmid.postprob.m[1][c]),
						    minPPlotting,//0.0,//double( minSegSitesPer1M )/double(1000000),
						    maxPPlotting, // 0.00500    = 4*2e-8*62500
						    sizeChunk,
						    !heteroEstResults[c].undef,
						    heteroEstResults[c].chrBreak,						      
						    HMMCODEMID)  != 0 ){
		    cerr<<"ERROR writing data point#"<<c<<" "<<heteroEstResults[c].rangeGen<<" to pdf file:"<<(vecOfPdfWritersHMM[i]->getFname())<<endl;
		    return 1;
		}


		
		if(  vecOfPdfWritersHMM[i]->drawHMM(heteroEstResults[c].rangeGen,
						    exp(hmmResmax.postprob.m[1][c]),
						    exp(hmmResmax.postprob.m[1][c]),
						    exp(hmmResmax.postprob.m[1][c]),
						    minPPlotting,//0.0,//double( minSegSitesPer1M )/double(1000000),
						    maxPPlotting, // 0.00500    = 4*2e-8*62500
						    sizeChunk,
						    !heteroEstResults[c].undef,
						    heteroEstResults[c].chrBreak,						      
						    HMMCODEMAX)  != 0 ){
		    cerr<<"ERROR writing data point#"<<c<<" "<<heteroEstResults[c].rangeGen<<" to pdf file:"<<(vecOfPdfWritersHMM[i]->getFname())<<endl;
		    return 1;
		}

	    }
	}	
    }

    for(unsigned int i=0;i<vecOfPdfWritersHMM.size();i++){
	delete(vecOfPdfWritersHMM[i]);
    }
    
    // toreturn.rohSegments    = rohSegments   ;
    // toreturn.nonrohSegments = nonrohSegments;
    // toreturn.unsureSegments = unsureSegments;
    
    // toreturn.postprob       = postprob;





// #ifdef notdef
    ofstream fileSummary;
    string filenameSummary = outFilePrefix+".summary.txt";
    fileSummary.open(filenameSummary.c_str());
    
    if (fileSummary.is_open()){

	fileSummary << "Command line:";
	for(int i=0;i<argc;i++){
	    fileSummary<<" "<<argv[i]<<" ";
	}
	fileSummary <<endl;
	fileSummary<<"Github version: "<< returnGitHubVersion(string(argv[0]),"..") <<" "<<endl;
	fileSummary<<"Random seed: "<< randomSeed <<" (re-run with --seed "<<randomSeed<<" to reproduce)"<<endl;
	
	
	fileSummary << "Genome-wide theta outside ROH:\t"<<hAvg     <<"\t"<<hMin     <<"\t"<<hMax     <<endl;
	fileSummary << "Genome-wide theta inc. ROH:\t"   <<hAvgNOROH<<"\t"<<hMinNOROH<<"\t"<<hMaxNOROH<<endl;
	
	//fileSummary << "\t" <<"total\tfraction in %"<<endl;
	fileSummary << "Segments unclassified    :\t"<<MID3( hmmResmin.unsureSegments, hmmResmid.unsureSegments, hmmResmax.unsureSegments) <<" ("<<MIN3( hmmResmin.unsureSegments, hmmResmid.unsureSegments, hmmResmax.unsureSegments)<<","<<MAX3( hmmResmin.unsureSegments, hmmResmid.unsureSegments, hmmResmax.unsureSegments)<<")"<<endl;
	double ratio1=(double(hmmResmid.unsureSegments)/double(hmmResmid.rohSegments+hmmResmid.nonrohSegments+hmmResmid.unsureSegments));
	double ratio2=MIN2(1.0,
			   double(MIN3(hmmResmin.unsureSegments,hmmResmid.unsureSegments,hmmResmax.unsureSegments))/
			   double(MAX3(hmmResmin.rohSegments+hmmResmin.nonrohSegments+hmmResmin.unsureSegments,
				       hmmResmid.rohSegments+hmmResmid.nonrohSegments+hmmResmid.unsureSegments,
				       hmmResmax.rohSegments+hmmResmax.nonrohSegments+hmmResmax.unsureSegments)));
	double ratio3=MIN2(1.0,double(MAX3(hmmResmin.unsureSegments,hmmResmid.unsureSegments,hmmResmax.unsureSegments))/
			   double(MIN3(hmmResmin.rohSegments+hmmResmin.nonrohSegments+hmmResmin.unsureSegments,
				       hmmResmid.rohSegments+hmmResmid.nonrohSegments+hmmResmid.unsureSegments,
				       hmmResmax.rohSegments+hmmResmax.nonrohSegments+hmmResmax.unsureSegments)));
	fileSummary << "Segments unclassified (%):\t"<<100*MID3(ratio1,ratio2,ratio3)<<" ("<<100*MIN3(ratio1,ratio2,ratio3)<<","<<100*MAX3(ratio1,ratio2,ratio3)<<")" <<endl;

	// fileSummary << "Segments unclassified (%):\t"<<100* (double(hmmResmid.unsureSegments)/double(hmmResmid.rohSegments+hmmResmid.nonrohSegments+hmmResmid.unsureSegments))<<" ("<<
	//   100* MIN2(1.0,
	// 	   double(MIN3(hmmResmin.unsureSegments,hmmResmid.unsureSegments,hmmResmax.unsureSegments))/
	// 	   double(MAX3(hmmResmin.rohSegments+hmmResmin.nonrohSegments+hmmResmin.unsureSegments,
	// 		      hmmResmid.rohSegments+hmmResmid.nonrohSegments+hmmResmid.unsureSegments,
	// 		       hmmResmax.rohSegments+hmmResmax.nonrohSegments+hmmResmax.unsureSegments)))<<","<<
	//   100* MIN2(1.0,double(MAX3(hmmResmin.unsureSegments,hmmResmid.unsureSegments,hmmResmax.unsureSegments))/
	// 	  double(MIN3(hmmResmin.rohSegments+hmmResmin.nonrohSegments+hmmResmin.unsureSegments,
	// 		      hmmResmid.rohSegments+hmmResmid.nonrohSegments+hmmResmid.unsureSegments,
	// 		      hmmResmax.rohSegments+hmmResmax.nonrohSegments+hmmResmax.unsureSegments)))<<")"
	// 	    <<endl;
	
	// fileSummary << "\t"<<"total\tfraction(defined)\tfraction(total)"<<endl;

	fileSummary << "Segments in ROH          :\t"      << MID3( hmmResmin.rohSegments, hmmResmid.rohSegments, hmmResmax.rohSegments)    <<" ("<<MIN3( hmmResmin.rohSegments, hmmResmid.rohSegments, hmmResmax.rohSegments)<<","<<MAX3( hmmResmin.rohSegments, hmmResmid.rohSegments, hmmResmax.rohSegments)<<")"<<endl;
	ratio1= double(hmmResmid.rohSegments)/double(hmmResmid.rohSegments+hmmResmid.nonrohSegments);
	ratio2= MIN2(1.0,double(MIN3(hmmResmin.rohSegments,hmmResmid.rohSegments,hmmResmax.rohSegments))/
		     double( MAX3( (hmmResmin.rohSegments+hmmResmin.nonrohSegments) , (hmmResmid.rohSegments+hmmResmin.nonrohSegments) , (hmmResmax.rohSegments+hmmResmin.nonrohSegments))));
	ratio3= MIN2(1.0,double(MAX3(hmmResmin.rohSegments,hmmResmid.rohSegments,hmmResmax.rohSegments))/
		     double( MIN3( (hmmResmin.rohSegments+hmmResmin.nonrohSegments),(hmmResmid.rohSegments+hmmResmin.nonrohSegments),(hmmResmax.rohSegments+hmmResmin.nonrohSegments))));

	fileSummary << "Segments in ROH(%)       :\t"      <<100*MID3(ratio1,ratio2,ratio3)<<" ("<<100*MIN3(ratio1,ratio2,ratio3)<<","<<100*MAX3(ratio1,ratio2,ratio3)<<")" <<endl;
	// fileSummary << "Segments in ROH(%)       :\t"      <<100*double(hmmResmid.rohSegments)/double(hmmResmid.rohSegments+hmmResmid.nonrohSegments)<<" ("<<
	//   100*MIN2(1.0,double(MIN3(hmmResmin.rohSegments,hmmResmid.rohSegments,hmmResmax.rohSegments))/
	// 	  double( MAX3( (hmmResmin.rohSegments+hmmResmin.nonrohSegments) , (hmmResmid.rohSegments+hmmResmin.nonrohSegments) , (hmmResmax.rohSegments+hmmResmin.nonrohSegments))))
	// 	    <<","<<
	//   100*MIN2(1.0,double(MAX3(hmmResmin.rohSegments,hmmResmid.rohSegments,hmmResmax.rohSegments))/
	// 	  double( MIN3( (hmmResmin.rohSegments+hmmResmin.nonrohSegments),(hmmResmid.rohSegments+hmmResmin.nonrohSegments),(hmmResmax.rohSegments+hmmResmin.nonrohSegments))))
	//     <<")"<<endl;

	fileSummary << "Segments in non-ROH      :\t"  << MID3( hmmResmin.nonrohSegments, hmmResmid.nonrohSegments, hmmResmax.nonrohSegments) <<" ("<<MIN3( hmmResmin.nonrohSegments, hmmResmid.nonrohSegments, hmmResmax.nonrohSegments)<<","<<MAX3( hmmResmin.nonrohSegments, hmmResmid.nonrohSegments, hmmResmax.nonrohSegments)<<")"<<endl;


	ratio1= double(hmmResmid.nonrohSegments)/double(hmmResmid.rohSegments+hmmResmid.nonrohSegments);
	ratio2= MIN2(1.0,double(MIN3(hmmResmin.nonrohSegments,hmmResmid.nonrohSegments,hmmResmax.nonrohSegments))/
		     double( MAX3( (hmmResmin.rohSegments+hmmResmin.nonrohSegments),(hmmResmid.rohSegments+hmmResmin.nonrohSegments),(hmmResmax.rohSegments+hmmResmin.nonrohSegments))));
	ratio3= MIN2(1.0,double(MAX3(hmmResmin.nonrohSegments,hmmResmid.nonrohSegments,hmmResmax.nonrohSegments))/
		     double( MIN3( (hmmResmin.rohSegments+hmmResmin.nonrohSegments),(hmmResmid.rohSegments+hmmResmin.nonrohSegments),(hmmResmax.rohSegments+hmmResmin.nonrohSegments))));

	fileSummary << "Segments in non-ROH (%)  :\t"  <<100*MID3(ratio1,ratio2,ratio3)<<" ("<<100*MIN3(ratio1,ratio2,ratio3)<<","<<100*MAX3(ratio1,ratio2,ratio3)<<")" <<endl;
	// fileSummary << "Segments in non-ROH (%)  :\t"  <<100*double(hmmResmid.nonrohSegments)/double(hmmResmid.rohSegments+hmmResmid.nonrohSegments)<<" ("<<	    
	//   100*MIN2(1.0,double(MIN3(hmmResmin.nonrohSegments,hmmResmid.nonrohSegments,hmmResmax.nonrohSegments))/
	// 	  double( MAX3( (hmmResmin.rohSegments+hmmResmin.nonrohSegments),(hmmResmid.rohSegments+hmmResmin.nonrohSegments),(hmmResmax.rohSegments+hmmResmin.nonrohSegments))))
	// 	    <<","<<
	//   100*MIN2(1.0,double(MAX3(hmmResmin.nonrohSegments,hmmResmid.nonrohSegments,hmmResmax.nonrohSegments))/
	// 	  double( MIN3( (hmmResmin.rohSegments+hmmResmin.nonrohSegments),(hmmResmid.rohSegments+hmmResmin.nonrohSegments),(hmmResmax.rohSegments+hmmResmin.nonrohSegments))))
	// 	    <<")"<<endl;


	fileSummary << "Avg. length of ROH       :\t"      <<MID3( hmmResmin.avgLengthROHSegments, hmmResmid.avgLengthROHSegments, hmmResmax.avgLengthROHSegments)   <<" ("<<MIN3( hmmResmin.avgLengthROHSegments, hmmResmid.avgLengthROHSegments, hmmResmax.avgLengthROHSegments)<<","<<MAX3( hmmResmin.avgLengthROHSegments, hmmResmid.avgLengthROHSegments, hmmResmax.avgLengthROHSegments)<<")"<<endl;


	
    }else{
	cerr << "Unable to print to file "<<filenameSummary<<endl;
    }
    fileSummary.close();
    cerr<<"final summary written to "<<filenameSummary<<endl;
    
// #endif


    
    //////////////////////////////////
    //                              //
    //  END h global                //
    //                              //
    //////////////////////////////////




    //BEGIN CLEANUP and SHUTDOWN
    
    delete cov2ProbSite;
    if(!skipToHMM){
	deleteLikelihoodScores();
    }

    cerr<<"ROHan finished succesfully"<<endl;
    
    return 0;
}
    







    

//     //#endif
// #ifdef TESTHMMSIMS
//     //calculating min/max number of sites per sizeChunk
//     int minSegSitesPerChunk;
//     int maxSegSitesPerChunk;
//     if(sizeChunk == 1000000){
// 	//those parameters are defined in the header, they correspond to the lowest possible h in a non-ROH region
// 	//and the highest
// 	minSegSitesPerChunk = minSegSitesPer1M;	
// 	maxSegSitesPerChunk = maxSegSitesPer1M;
//     }else{
// 	minSegSitesPerChunk = int( (double(minSegSitesPer1M)/double(1000000))*double(sizeChunk) );
// 	maxSegSitesPerChunk = int( (double(maxSegSitesPer1M)/double(1000000))*double(sizeChunk) );
//     }
    
//     Hmm hmm (minSegSitesPerChunk,maxSegSitesPerChunk,sizeChunk);
    
//     cerr<<"Begin running HMM"<<endl;
//     cerr<<"generating a random set"<<endl;
//     vector<emission>       eTest = hmm.generateStates(250,sizeChunk);

    
//     vector<emissionUndef>  eTestUndef;
//     cerr<<"done"<<endl;
//     //TODO add multiple bootstraps
//     //- generate a number between hLower and hUpper
//     // or modified probEmmission to reflect the uncertainty
//     vector<long double> emittedH;
    
//     //TODO add a way to avoid missing data,
//     //in the log-likelihood computation for missing data and chr breaks
//     //for(unsigned int hWindow=0;hWindow<heteroEstResults.size();hWindow++){
    
//     // for(unsigned int hWindow=0;hWindow<250;hWindow++){
//     // 	// if(dataToWrite->hetEstResults.hasConverged){
//     // 	//     emittedH.push_back( heteroEstResults[i].h );
//     // 	// }
//     // }
//     long double uncertainty = 0.000000000001;
//     //                        0.00072
//     //                        0.00005
//     //                        0.00001
//     //long double uncertainty = 0.00001000;
//     for(unsigned int i=0;i<eTest.size();i++){

// 	emittedH.push_back( eTest[i].p );
// 	emissionUndef e;
// 	e.hlow     = eTest[i].p - uncertainty;
// 	if(e.hlow < 0)
// 	    e.hlow=0.0;
// 	e.hhigh    = eTest[i].p + uncertainty;
	
// 	e.undef    = false;
// 	e.chrBreak = false;
// 	cout<<i<<"\t"<<eTest[i].idx<<"\t"<<eTest[i].p<<"\t"<<e.hlow<<" "<<e.hhigh<<" "<<e.undef<<" "<<e.chrBreak <<endl;	
// 	eTestUndef.push_back( e );
//     }

    
//     cerr<<"ok"<<endl;
//     //BEGIN MCMC chain
//     //init to random settings
//     long double partition= (long double)(stepHMM);
//     int accept=0;
//     long double x_i    ;
//     long double x_i_1  ;
    
//     //het rate
//     //int sitesPer1M     = randomInt(200,1000);// between 1 and 1000 sites per million
//     int sitesPer1M     = 1000;// 1000 sites per million
//     long double h_i    = double(       sitesPer1M )/double(1000000);
//     long double h_i_1;

//     long double hlower = double( minSegSitesPer1M )/double(1000000);
//     long double hupper = double( maxSegSitesPer1M )/double(1000000);
    
//     //transition rate
//     // long double pTlowerFirstHalf  =    numeric_limits<double>::epsilon();//TODO put back 0.5 
//     // long double pTlowerSecondHalf =    numeric_limits<double>::epsilon();
//     long double pTlowerFirstHalf  =    numeric_limits<double>::min();//TODO put back 0.5 
//     long double pTlowerSecondHalf =    numeric_limits<double>::min();
	
//     long double pTlower = pTlowerFirstHalf;	
//     //long double pTupper = 1.0 - numeric_limits<double>::epsilon();
//     long double pTupper = 1.0 - numeric_limits<double>::min();

//     long double pT_i = randomProb()*(pTupper-pTlower) + pTlower;
//     long double pT_i_1;
//     //long double pTlower =       numeric_limits<double>::epsilon();
        
//     random_device rd;
//     default_random_engine dre (rd());
//     //int maxChains = 100000;
//     //int maxChains   =  50000;
//     //chain 0

//     hmm.setHetRateForNonROH(h_i);
//     hmm.setTransprob(pT_i);
//     //x_i    =  forwardProb(&hmm, emittedH , sizeChunk);
//     x_i    =  forwardProbUncertaintyMissing(&hmm, eTestUndef , sizeChunk);
//     cout<<setprecision(10)<<"\tinitial\t"<<h_i<<"\t"<<pT_i<<"\t"<<x_i<<"\t"<<endl;    
//     for(int chain=1;chain<=maxChains;chain++){

// 	//computing new state
// 	normal_distribution<long double> distribution_h(h_i,(hupper-hlower)/partition  );
// 	h_i_1      = distribution_h(dre);

// 	if(h_i_1 <= hlower     ||  h_i_1 >= hupper     ){
// 	    h_i_1      = h_i;
// 	}


// 	normal_distribution<long double> distribution_pT(pT_i,(pTupper-pTlower)/partition  );
// 	pT_i_1      = distribution_pT(dre);

// 	if(pT_i_1 <= pTlower     ||  pT_i_1 >= pTupper     ){
// 	    pT_i_1      = pT_i;
// 	}
       
// 	hmm.setHetRateForNonROH(h_i_1);
// 	hmm.setTransprob(pT_i_1);

// 	//x_i_1=forwardProb(&hmm, emittedH , sizeChunk);
// 	x_i_1=forwardProbUncertaintyMissing(&hmm, eTestUndef , sizeChunk);

// 	if(chain>(maxChains/4)){
// 	    pTlower = pTlowerSecondHalf;
// 	}

// 	long double acceptance = min( (long double)(1.0)  , expl(x_i_1-x_i) );
// 	if( (long double)(randomProb()) < acceptance){
// 	    h_i           =  h_i_1;
// 	    pT_i          =  pT_i_1;
// 	    x_i           =  x_i_1;
// 	    accept++;
// 	    //cout<<setprecision(10)<<"accepted jump from\t"<<h_i<<"\t"<<pT_i<<"\t"<<x_i<<"\tto\t"<<h_i_1<<"\t"<<pT_i_1<<"\t"<<x_i_1<<""<<"\t"<<acceptance<<" "<<accept<<" "<<chain<<" "<<double(accept)/double(chain)<<endl;
// 	    //cerr<<setprecision(10)<<"mcmc"<<mcmc<<"\taccepted\t"<<h_i<<"\t"<<pT_i<<"\t"<<x_i<<"\t"<<" "<<acceptance<<" "<<accept<<" "<<chain<<" "<<double(accept)/double(chain)<<endl;	    
// 	}else{
// 	    //cout<<setprecision(10)<<"rejected jump from\t"<<h_i<<"\t"<<pT_i<<"\t"<<x_i<<"\tto\t"<<h_i_1<<"\t"<<pT_i_1<<"\t"<<x_i_1<<""<<"\t"<<acceptance<<" "<<accept<<" "<<chain<<" "<<double(accept)/double(chain)<<endl;
// 	}
// 	//chain++;
// 	//sleep(0.1);
//     }
// cout<<setprecision(10)<<"mcmc"<<"\tfinal\t"<<h_i<<"\t"<<pT_i<<"\t"<<x_i<<"\t"<<endl;

	
//     //cerr<<"Baum Welch"<<endl;
//     //baum_welch(&hmm,&emittedH);
	
//     cerr<<"HMM done"<<endl;
// #endif






// //! A method to initialize the likelihood of observing data to avoid recomputation
// /*!
//   This method is called by the main after calling initDefaultBaseFreq
// */
// void initLikelihoodScores(){







// //  // 1D: pos away from 5'/3' end
// //  // 2D: mapping quality 
// //  // 3D: base qual
// // vector< vector< vector<diNucleotideProb> > >  pos2mpq2BaseQual2SubMatrix5p;
// // vector< vector< vector<diNucleotideProb> > >  pos2mpq2BaseQual2SubMatrix3p;
    
//     for(int l=0;l<( int(MAXLENGTHFRAGMENT/2) +1);l++){//for each position

// 	mpq2bsq2submatrix  mpq2BaseQualSubMatrix5p;
// 	mpq2bsq2submatrix  mpq2BaseQualSubMatrix3p;
	
// 	for(int mq=0;mq<=MAXMAPPINGQUAL;mq++){ //for each mapping quality 




// 	    //we need 4 theoretical bases X 4 observed bases X MAXMAPPINGQUAL mapping quality X MAXLENGTHFRAGMENT positions on the fragment X base quality
// 	    vector<diNucleotideProb>  baseQual2SubMatrix5p;//FOR EACH QUAL SCORE
// 	    vector<diNucleotideProb>  baseQual2SubMatrix3p;//FOR EACH QUAL SCORE



// 	    for(int q=0;q<=MAXBASEQUAL;q++){ //for each base quality
		
// 		//vector<diNucleotideProb> pos2SubMatrix;//FOR EACH QUAL SCORE
// 		//vector< vector< vector<diNucleotideProb> > > mpq2Length2BaseQual2Pos2SubMatrix;//FOR EACH FRAGMENT LENGTH
		
		
// 		//we know the probability of P(b1->b2) after deamination is vector<diNucleotideProb> sub5pDiNuc; vector<diNucleotideProb> sub3pDiNuc;
		
// 		diNucleotideProb toAddForBaseQual5p;
// 		diNucleotideProb toAddForBaseQual3p;

// 		diNucleotideProb toAddForBaseQual5p_;
// 		diNucleotideProb toAddForBaseQual3p_;
		
		
// 		for(int bTheo=0;bTheo<4;bTheo++){
// 		    for(int bObs=0;bObs<4;bObs++){
// 			toAddForBaseQual5p_.p[bTheo][bObs]=0.0;
// 			toAddForBaseQual3p_.p[bTheo][bObs]=0.0;
// 		    }
// 		}
		    
// 		//q is base quality 
// 		for(int bTheo=0;bTheo<4;bTheo++){               // each possible theoritical base
// 		    for(int bpstDeam=0;bpstDeam<4;bpstDeam++){	// each possible deaminated base		
// 			for(int bObs=0;bObs<4;bObs++){	        // each possible observed base
			    
// 			    //cerr<<"tripleloop1\tl="<<l<<"\tmq="<<mq<<"\tq="<<q<<"\t"<<sub5pDiNuc[    l].p[bTheo][bpstDeam] <<"\t"<<likeMatchProb[q] <<"\t"<<   defaultSubMatchMatrix.p[bpstDeam][bObs]<<"\t"<<likeMismatchProb[q] <<"\t"<<illuminaErrorsProbMatrix.p[bpstDeam][bObs]<<"\t"<<toAddForBaseQual5p_.p[bTheo][bObs]<<"\t"<<toAddForBaseQual3p_.p[bTheo][bObs]<<"\t"<<sub5pDiNuc[    l].p[bTheo][bpstDeam]<<"*"<<(		    likeMatchProb[q]    *   defaultSubMatchMatrix.p[bpstDeam][bObs]				    +				    likeMismatchProb[q] * illuminaErrorsProbMatrix.p[bpstDeam][bObs]				)<<"\t"<<				sub3pDiNuc[    l].p[bTheo][bpstDeam]<<" * "<<(				    likeMatchProb[q]    *   defaultSubMatchMatrix.p[bpstDeam][bObs]				    +				    likeMismatchProb[q] * illuminaErrorsProbMatrix.p[bpstDeam][bObs]				)<<endl;

// 			    toAddForBaseQual5p_.p[bTheo][bObs] +=        
// 				sub5pDiNuc[    l].p[bTheo][bpstDeam] * (
// 									likeMatchProb[q]    *   defaultSubMatchMatrix.p[bpstDeam][bObs]
// 									+
// 									likeMismatchProb[q] * illuminaErrorsProbMatrix.p[bpstDeam][bObs]
// 									);
			    
// 			    toAddForBaseQual3p_.p[bTheo][bObs] += 
// 				sub3pDiNuc[    l].p[bTheo][bpstDeam] * (
// 									likeMatchProb[q]    *   defaultSubMatchMatrix.p[bpstDeam][bObs]
// 									+
// 									likeMismatchProb[q] * illuminaErrorsProbMatrix.p[bpstDeam][bObs]
// 									);

// 			    // cerr<<"tripleloop2\tl="<<l<<"\tmq="<<mq<<"\tq="<<q<<"\t"<<sub5pDiNuc[    l].p[bTheo][bpstDeam] <<"\t"<<likeMatchProb[q] <<"\t"<<   defaultSubMatchMatrix.p[bpstDeam][bObs]<<"\t"<<likeMismatchProb[q] <<"\t"<<illuminaErrorsProbMatrix.p[bpstDeam][bObs]<<"\t"<<toAddForBaseQual5p_.p[bTheo][bObs]<<"\t"<<toAddForBaseQual3p_.p[bTheo][bObs]<<endl;

// 			}//end each bObs
// 		    }//end each bpstDeam
// 		}//for each bTheo


		
		
// 		for(int bTheo=0;bTheo<4;bTheo++){               // each possible theoretical base
// 		    for(int bObs=0;bObs<4;bObs++){	        // each possible observed base
			
// 			//0.5 since this is the prior prob of having sampled a given chromosome
// #ifdef PRECOMPUTELOG
// 			toAddForBaseQual5p.p[bTheo][bObs] = logl(0.5 * (likeMatchProbMap[mq]*toAddForBaseQual5p_.p[bTheo][bObs] + likeMismatchProbMap[mq] * defaultDNA5p[l].f[bObs]) );			
// 			toAddForBaseQual3p.p[bTheo][bObs] = logl(0.5 * (likeMatchProbMap[mq]*toAddForBaseQual3p_.p[bTheo][bObs] + likeMismatchProbMap[mq] * defaultDNA3p[l].f[bObs]) );

// // 			if(mq==37 && q==38 && bObs==3 && bTheo == 0 && l==5 ){
			    
// // 			    cerr<<setprecision(20)<<"doubleloop\tl="<<l<<"\tmq="<<mq<<"\tq="<<q<<"\t"<<toAddForBaseQual5p.p[bTheo][bObs]<<"\t=\t"<<likeMatchProbMap[mq]<<"*"<<toAddForBaseQual5p_.p[bTheo][bObs] <<"+"<< likeMismatchProbMap[mq]<<"*"<< defaultDNA5p[l].f[bObs]

// // 				<<"="<<(likeMatchProbMap[mq]*toAddForBaseQual5p_.p[bTheo][bObs] + likeMismatchProbMap[mq] * defaultDNA5p[l].f[bObs] )
// // 				<<"\t"<<logl(likeMatchProbMap[mq]*toAddForBaseQual5p_.p[bTheo][bObs] + likeMismatchProbMap[mq] * defaultDNA5p[l].f[bObs] )			
// // 				<<"\t"<<logl(likeMatchProbMap[mq]*toAddForBaseQual5p_.p[bTheo][bObs] + likeMismatchProbMap[mq] * defaultDNA5p[l].f[bObs] )			
// // 				<<"\t"<<toAddForBaseQual5p.p[bTheo][bObs]<<"\t"<<expl(toAddForBaseQual5p.p[bTheo][bObs])
// // <<endl;
// // 			}

// #else
// 			// without logl
// 			toAddForBaseQual5p.p[bTheo][bObs] =  likeMatchProbMap[mq]*toAddForBaseQual5p_.p[bTheo][bObs] + likeMismatchProbMap[mq] * defaultDNA5p[l].f[bObs];			
// 			toAddForBaseQual3p.p[bTheo][bObs] =  likeMatchProbMap[mq]*toAddForBaseQual3p_.p[bTheo][bObs] + likeMismatchProbMap[mq] * defaultDNA3p[l].f[bObs];
			
			
// 			// if(mq==37 && q==38 && bObs==3 && bTheo == 0 && l==5 ){
			    
// 			//     cerr<<setprecision(20)<<"doubleloop\tl="<<l<<"\tmq="<<mq<<"\tq="<<q<<"\t"<<toAddForBaseQual5p.p[bTheo][bObs]<<"\t=\t"<<likeMatchProbMap[mq]<<"*"<<toAddForBaseQual5p_.p[bTheo][bObs] <<"+"<< likeMismatchProbMap[mq]<<"*"<< defaultDNA5p[l].f[bObs]<<"\t"<<toAddForBaseQual5p.p[bTheo][bObs]<<endl;
// 			// }
			
// #endif
// 			// cerr<<"doubleloop\tl="<<l<<"\tmq="<<mq<<"\tq="<<q<<"\t"<<toAddForBaseQual5p.p[bTheo][bObs]<<"\t=\t"<<likeMatchProbMap[mq]<<"*"<<toAddForBaseQual5p_.p[bTheo][bObs] <<"+"<< likeMismatchProbMap[mq]<<"*"<< defaultDNA5p[l].f[bObs]<<endl;
// 			// cerr<<"doubleloop\tl="<<l<<"\tmq="<<mq<<"\tq="<<q<<"\t"<<toAddForBaseQual3p.p[bTheo][bObs]<<"\t=\t"<<likeMatchProbMap[mq]<<"*"<<toAddForBaseQual3p_.p[bTheo][bObs] <<"+"<< likeMismatchProbMap[mq]<<"*"<< defaultDNA3p[l].f[bObs]<<endl;
			
// 		    }//for each bObs
// 		}//for each bTheo, bpstDeam and bObs
		
// 		baseQual2SubMatrix5p.push_back(toAddForBaseQual5p);
// 		baseQual2SubMatrix3p.push_back(toAddForBaseQual3p);
		
		
// #ifdef DEBUGINITLIKELIHOODSCORES
// 		cerr<<"ppos = "<<l<<" mq = "<<mq<<" ("<<likeMatchProbMap[mq]<<" "<<likeMismatchProbMap[mq]<<" )  q = "<<q<<endl;
// 		cerr<<"5'----------"<<endl;
// 		for(int nuc1=0;nuc1<4;nuc1++){
// 		    cerr<<"ACGT"[nuc1]<<"\t";		    
// 		    long double s=0;
// 		    for(int nuc2=0;nuc2<4;nuc2++){
// 			cerr<<toAddForBaseQual5p_.p[nuc1][nuc2]<<"\t";			
// 			s+=toAddForBaseQual5p_.p[nuc1][nuc2];
// 		    }
// 		    cerr<<"\t"<<s<<"\t";
// 		    s=0;
		    
// 		    for(int nuc2=0;nuc2<4;nuc2++){
// 			cerr<<defaultDNA5p[l].f[nuc2]<<"\t";			
// 			s+=defaultDNA5p[l].f[nuc2];
// 		    }

// 		    cerr<<"\t"<<s<<"\t";
// 		    s=0;
		    
// 		    for(int nuc2=0;nuc2<4;nuc2++){
// 			cerr<<toAddForBaseQual5p.p[nuc1][nuc2]<<"\t";
// 			s+=toAddForBaseQual5p.p[nuc1][nuc2];
// 		    }

// 		    cerr<<"\t"<<s<<endl;
// 		}
// 		cerr<<endl;
// 		cerr<<"3'----------"<<endl;
// 		for(int nuc1=0;nuc1<4;nuc1++){
// 		    cerr<<"ACGT"[nuc1]<<"\t";
// 		    long double s=0;

// 		    for(int nuc2=0;nuc2<4;nuc2++){
// 			cerr<<toAddForBaseQual3p_.p[nuc1][nuc2]<<"\t";
// 			s+=toAddForBaseQual3p_.p[nuc1][nuc2];
// 		    }
// 		    cerr<<"\t"<<s<<"\t";
// 		    s=0;
		    
// 		    for(int nuc2=0;nuc2<4;nuc2++){
// 			cerr<<defaultDNA3p[l].f[nuc2]<<"\t";
// 			s+=defaultDNA3p[l].f[nuc2];			
// 		    }
// 		    cerr<<"\t"<<s<<"\t";
// 		    s=0;
		    
// 		    for(int nuc2=0;nuc2<4;nuc2++){
// 			cerr<<toAddForBaseQual3p.p[nuc1][nuc2]<<"\t";
// 			s+=toAddForBaseQual3p.p[nuc1][nuc2];
// 		    }
// 		    cerr<<"\t"<<s<<endl;
// 		}
// 		cerr<<endl;
// #endif
		    
// 	    }//for each base qual
	    


// 	    mpq2BaseQualSubMatrix5p.push_back(baseQual2SubMatrix5p);
// 	    mpq2BaseQualSubMatrix3p.push_back(baseQual2SubMatrix3p);
	    
	    
// 	}// for each mapping quality

// 	pos2mpq2BaseQual2SubMatrix5p.push_back(mpq2BaseQualSubMatrix5p);
// 	pos2mpq2BaseQual2SubMatrix3p.push_back(mpq2BaseQualSubMatrix3p);

	
//     }//for each position
    


//     //if(mq==37 && q==38 && bObs==0 && bTheo == 0 && l==5 ){
//     //    cerr<<setprecision(20)<<"TEST1 "<<pos2mpq2BaseQual2SubMatrix5p[  12  ][37][38].p[0][3]<<"\t"<<pos2mpq2BaseQual2SubMatrix5p[  12  ][37][18].p[0][3]<<endl;
// 	//cerr<<"doubleloop\tl="<<l<<"\tmq="<<mq<<"\tq="<<q<<"\t"<<toAddForBaseQual5p.p[bTheo][bObs]<<"\t=\t"<<likeMatchProbMap[mq]<<"*"<<toAddForBaseQual5p_.p[bTheo][bObs] <<"+"<< likeMismatchProbMap[mq]<<"*"<< defaultDNA5p[l].f[bObs]<<"\t"<<toAddForBaseQual5p.p[bTheo][bObs]<<endl;
//     //}
    

//     //vector< vector< mpq2bsq2submatrix * > > length2pos2mpq2bsq2submatrix;
//     //dummy values
//     for(int L=0;L<int(MINLENGTHFRAGMENT);L++){//for each fragment length
// 	vector< mpq2bsq2submatrix * > vectorToAdd;	
// 	length2pos2mpq2bsq2submatrix.push_back(vectorToAdd);
//     }

//     for(int L=int(MINLENGTHFRAGMENT);L<=int(MAXLENGTHFRAGMENT);L++){//for each fragment length
// 	vector< mpq2bsq2submatrix * > vectorToAdd;
// 	for(int l=0;l<L;l++){//for each pos
// 	    //cerr<<l<<" "<<(L-l-1)<<" "<<L<<endl;	   
	    
// 	    if( l<(L/2) ){//use 5' substitutions
// 		vectorToAdd.push_back( &pos2mpq2BaseQual2SubMatrix5p[  l  ] );
// 		//cerr<< "5' size "<<pos2mpq2BaseQual2SubMatrix5p[  l  ].size() <<" sizemq0 "<<pos2mpq2BaseQual2SubMatrix5p[  l  ][0].size()<<endl;
// 	    }else{        //use 3' substitutions
// 		//cerr<< "3' size "<<pos2mpq2BaseQual2SubMatrix3p[  l  ].size() <<" sizemq0 "<<pos2mpq2BaseQual2SubMatrix3p[  l  ][0].size()<<endl;
// 		vectorToAdd.push_back( &pos2mpq2BaseQual2SubMatrix3p[L-l-1] );
// 	    }

// 	}
// 	length2pos2mpq2bsq2submatrix.push_back(vectorToAdd);
//     }//for each fragment length
    
//     // cerr<<setprecision(20)<<"TEST2 AT 38 "<<pos2mpq2BaseQual2SubMatrix5p[  12  ][37][38].p[0][3]<<"\tL2P="<<length2pos2mpq2bsq2submatrix[113][12]->at(37)[38].p[0][3]<<endl ;
//     // cerr<<setprecision(20)<<"TEST2 AT 18 "<<pos2mpq2BaseQual2SubMatrix5p[  12  ][37][18].p[0][3]<<"\tL2P="<<length2pos2mpq2bsq2submatrix[113][12]->at(37)[18].p[0][3]<<endl ;

// #ifdef DEBUGINITLIKELIHOODSCORES2

//     for(int L=MINLENGTHFRAGMENT;L<=MAXLENGTHFRAGMENT;L++){//for each fragment length

	
// 	for(int l=0;l<L;l++){//for each pos
// 	    // cerr<<"pos "<<l<<"/"<<L<<endl;//<<"\t"<<length2pos2mpq2bsq2submatrix[L][l]->size()<<endl;
	    
// 	    if( l == DEBUGINITLIKELIHOODSCORES2){

// 		for(int mq=0;mq<=MAXMAPPINGQUAL;mq++){ //for each mapping quality 
// 		    //cerr<<"mq="<<mq<<"\t"<<endl;//length2pos2mpq2bsq2submatrix[L][l]->at(mq).size()<<endl;

// 		    for(int q=0;q<=MAXBASEQUAL;q++){ //for each base quality
// 			cerr<<"pos "<<l<<"/"<<L<<" mq="<<mq<<" bq="<<q<<endl;

// 			//cerr<<" mq = "<<mq<<" ("<<likeMatchProbMap[mq]<<" "<<likeMismatchProbMap[mq]<<" )  q = "<<q<<endl;
// 			for(int nuc1=0;nuc1<4;nuc1++){
// 			    cerr<<"ACGT"[nuc1]<<"\t";		    
// 			    long double s=0;//sum of probs

// 			    for(int nuc2=0;nuc2<4;nuc2++){
// 				cerr<<2*exp(length2pos2mpq2bsq2submatrix[L][l]->at(mq)[q].p[nuc1][nuc2])<<"\t";
// 				//length2pos2mpq2bsq2submatrix[L][l]->at(mq)
// 				//toAddForBaseQual5p.p[nuc1][nuc2]<<"\t";
// 				s+=length2pos2mpq2bsq2submatrix[L][l]->at(mq)[q].p[nuc1][nuc2];
// 			    }

// 			    cerr<<"\tsum="<<s<<endl;
// 			}
// 			cerr<<endl;
		    
// 		    }//for each base qual	   		    

// 		}// for each mapping quality	   
// 	    }
// 	}//for each pos
//     }//for each fragment length

// #endif    

//     //cerr<<"done "<<endl;

//     //exit(1);

// }// END initLikelihoodScores()






//OLD HMM code

// Hmm hmm;
// //cout<<"Hmm:"<<hmm<<endl;
// vector<double> emittedH;
// //cerr<<"generating states"<<endl;
    
// vector<emission>  eTest = hmm.generateStates(250,sizeChunk);

// //cerr<<"test data generated:"<<endl;
    
// for(unsigned int i=0;i<eTest.size();i++){
//     cerr<<i<<"\t"<<eTest[i].idx<<"\t"<<eTest[i].p<<endl;
//     emittedH.push_back( eTest[i].p );
//  }

// // cerr<<"forward "<<endl;
// // vector< vector<double> > f=forward(&hmm, emittedH , sizeChunk);//, const int n) {
// // cerr<<"backward "<<endl;
// // vector< vector<double> > b=backward(&hmm, emittedH, sizeChunk);//, const int n) {

// // cerr<<"forward-backward"<<endl;
// // double sumProbData=forwardBackward(&hmm, emittedH, sizeChunk);//, const int n) {
// // cerr<<"forward-backward "<<sumProbData<<endl;

// // cerr<<"forward Prob()"<<endl;
// // //for(double d=0.000
// // // for (long double hHMM = 0.0001;hHMM < 0.001; hHMM += 0.00005){	
// // // 	//cout<<"test "<<hHMM<<endl;
// // // 	hmm.setHetRateForNonROH(hHMM);
// // // 	double fwdProb=forwardProb(&hmm, emittedH , sizeChunk);
// // // 	cout<<hHMM<<"\t"<<fwdProb<<endl;
// // // }

// // // for (long double pTrans = 0.01;pTrans < 0.99; pTrans += 0.01){	
// // // 	//cout<<"test "<<hHMM<<endl;
// // // 	hmm.setHetRateForNonROH(0.0007);
// // // 	hmm.setTransprob(pTrans);
// // // 	double fwdProb=forwardProb(&hmm, emittedH , sizeChunk);
// // // 	cout<<pTrans<<"\t"<<fwdProb<<endl;
// // // }
// cout<<"MCMC"<<endl;

// double step = 1000;
// //init to random settings
// long double partition= (long double)(step);
// int accept=0;
// long double x_i    ;
// long double x_i_1  ;

// //het rate
// //int sitesPer1M     = randomInt(200,1000);// between 1 and 1000 sites per million
// int sitesPer1M     = 1000;// 1000 sites per million
// long double h_i    = double(sitesPer1M)/double(1000000);
// long double h_i_1;
// long double hlower = double( 200)/double(1000000);
// long double hupper = double(1000)/double(1000000);
    
// //transition rate
// long double pTlowerFirstHalf  =    0.5;
// long double pTlowerSecondHalf =    numeric_limits<double>::epsilon();
	
// long double pTlower = pTlowerFirstHalf;	
// long double pTupper = 1.0 - numeric_limits<double>::epsilon();

// long double pT_i = randomProb()*(pTupper-pTlower) + pTlower;
// long double pT_i_1;
// //long double pTlower =       numeric_limits<double>::epsilon();
        
// random_device rd;
// default_random_engine dre (rd());
// int maxChains = 100000;
// //chain 0

// hmm.setHetRateForNonROH(h_i);
// hmm.setTransprob(pT_i);
// x_i    =  forwardProb(&hmm, emittedH , sizeChunk);
// cerr<<setprecision(10)<<"mcmc"<<mcmc<<"\tinitial\t"<<h_i<<"\t"<<pT_i<<"\t"<<x_i<<"\t"<<endl;    
// for(int chain=1;chain<=maxChains;chain++){

//     //computing new state
//     normal_distribution<long double> distribution_h(h_i,(hupper-hlower)/partition  );
//     h_i_1      = distribution_h(dre);

//     if(h_i_1 <= hlower     ||  h_i_1 >= hupper     ){
// 	h_i_1      = h_i;
//     }


//     normal_distribution<long double> distribution_pT(pT_i,(pTupper-pTlower)/partition  );
//     pT_i_1      = distribution_pT(dre);

//     if(pT_i_1 <= pTlower     ||  pT_i_1 >= pTupper     ){
// 	pT_i_1      = pT_i;
//     }
       
//     hmm.setHetRateForNonROH(h_i_1);
//     hmm.setTransprob(pT_i_1);

//     x_i_1=forwardProb(&hmm, emittedH , sizeChunk);

//     if(chain>(maxChains/2)){
// 	pTlower = pTlowerSecondHalf;
//     }

//     long double acceptance = min( (long double)(1.0)  , expl(x_i_1-x_i) );
//     if( (long double)(randomProb()) < acceptance){
// 	h_i           =  h_i_1;
// 	pT_i          =  pT_i_1;
// 	x_i           =  x_i_1;
// 	accept++;
// 	cerr<<setprecision(10)<<"accepted jump from\t"<<h_i<<"\t"<<pT_i<<"\t"<<x_i<<"\tto\t"<<h_i_1<<"\t"<<pT_i_1<<"\t"<<x_i_1<<""<<"\t"<<acceptance<<" "<<accept<<" "<<chain<<" "<<double(accept)/double(chain)<<endl;
// 	//cerr<<setprecision(10)<<"mcmc"<<mcmc<<"\taccepted\t"<<h_i<<"\t"<<pT_i<<"\t"<<x_i<<"\t"<<" "<<acceptance<<" "<<accept<<" "<<chain<<" "<<double(accept)/double(chain)<<endl;
	    
//     }else{
// 	cerr<<setprecision(10)<<"rejected jump from\t"<<h_i<<"\t"<<pT_i<<"\t"<<x_i<<"\tto\t"<<h_i_1<<"\t"<<pT_i_1<<"\t"<<x_i_1<<""<<"\t"<<acceptance<<" "<<accept<<" "<<chain<<" "<<double(accept)/double(chain)<<endl;
//     }
//     chain++;
//     //sleep(0.1);
//  }
// cerr<<setprecision(10)<<"mcmc"<<mcmc<<"\tfinal\t"<<h_i<<"\t"<<pT_i<<"\t"<<x_i<<"\t"<<endl;

	
// //cerr<<"Baum Welch"<<endl;
// //baum_welch(&hmm,&emittedH);
	
// cerr<<"testing viterbi algorithm"<<endl;
	
// hmmpath hp=viterbi(&hmm, emittedH, sizeChunk);//, const int n) {
// // cout<<"done"<<endl;
















// void runHMM(const string & outFilePrefix, const    vector<emissionUndef> & heteroEstResults, const int maxChains , const double fracChainsBurnin,const unsigned char usemin){

//     cerr<<"Creating HMM..";
//     //write chains to output
//     //bgzipWriterMCMC
//     string headerHMMMCMC = "#llik\th\ts\tp\taccepted\tchains\tacptrate\n";
//     Internal::BgzfStream  bgzipWriterMCMC;
//     bgzipWriterMCMC.Open(outFilePrefix+".hmmmcmc.gz", IBamIODevice::WriteOnly);
//     if(!bgzipWriterMCMC.IsOpen()){
// 	cerr<<"Cannot open file "<<(outFilePrefix+".hmmmcmc.gz")<<" in bgzip writer"<<endl;
// 	exit(1);
//     }

//     //computing the min/max segsites per chunk
//     int minSegSitesPerChunk;
//     int maxSegSitesPerChunk;
//     if(sizeChunk == 1000000){
// 	minSegSitesPerChunk = minSegSitesPer1M;	
// 	maxSegSitesPerChunk = maxSegSitesPer1M;
//     }else{
// 	minSegSitesPerChunk = int( (double(minSegSitesPer1M)/double(1000000))*double(sizeChunk) );
// 	maxSegSitesPerChunk = int( (double(maxSegSitesPer1M)/double(1000000))*double(sizeChunk) );
//     }
    
//     Hmm hmm (minSegSitesPerChunk,maxSegSitesPerChunk,sizeChunk,1000);
    
//     //cerr<<"Begin running MCMC on HMM"<<endl;
//     // cerr<<"generating a random set"<<endl;
//     // vector<emission>       eTest = hmm.generateStates(250,sizeChunk);

//     cerr<<".";
//     //vector<emissionUndef>  eTestUndef;


//     //BEGIN MCMC chain
//     //init to random settings
//     long double partition= (long double)(stepHMM);
//     int accept=0;

//     //likelihood
//     long double x_i    ;
//     long double x_i_1  ;
    
//     //het rate
//     //int sitesPer1M     = randomInt(200,1000);// between 1 and 1000 sites per million
//     int sitesPer1M     = 1000;// 1000 heterozygous sites per million
//     long double h_i    = double(       sitesPer1M )/double(1000000);
//     long double h_i_1;

//     long double hlower = double( minSegSitesPer1M )/double(1000000);
//     long double hupper = double( maxSegSitesPer1M )/double(1000000);

//     //number of non-recombining chunks per window
//     long double s_i    = 100.0; //starting value, there are 100 non-recombining window in sizeChunk
//     long double s_i_1;

//     long double slower = 1.0;       // the whole sizeChunk is a non-recombining window
//     long double supper = sizeChunk; // 1 base is a non-recombining window
    

    
//     //transition rate
//     long double pTlowerFirstHalf  =    0.5;//overestimate the transition probability 
//     //long double pTlowerSecondHalf =    numeric_limits<double>::epsilon();
//     long double pTlowerSecondHalf =    numeric_limits<double>::min();
	
//     long double pTlower = pTlowerFirstHalf;	
//     //long double pTupper = 1.0 - numeric_limits<double>::epsilon();
//     long double pTupper = 1.0 - numeric_limits<double>::min();

//     //generate a random initial probability between pTlower and pTupper
//     //probability of transition
//     long double pT_i = randomProb()*(pTupper-pTlower) + pTlower;
//     long double pT_i_1;
//     //long double pTlower =       numeric_limits<double>::epsilon();
        
//     random_device rd;
//     default_random_engine dre (rd());
//     //int maxChains = 100000;
//     //chain 0

//     hmm.setHetRateForNonROH(h_i);
//     hmm.setTransprob(pT_i);
//     hmm.setNrwPerSizeChunk( (unsigned int)s_i );
//     hmm.recomputeProbs();
    
//     cerr<<".";
//     // for(unsigned int i=0;i<heteroEstResults.size();i++){
//     // 	cerr<<"obs#"<<i<<" "<<heteroEstResults[i].chrBreak<<"\t"<<heteroEstResults[i].undef<<"\t"<<heteroEstResults[i].h<<"\t"<<heteroEstResults[i].hlow<<"\t"<<heteroEstResults[i].hhigh<<"\t"<<heteroEstResults[i].weight<<endl;
//     // }
//     // return 1;
//     //x_i    =  forwardProb(&hmm, emittedH , sizeChunk);
//     fbreturnVal  tmpResFWD = forwardProbUncertaintyMissing(&hmm, heteroEstResults , sizeChunk,true);
//     x_i  = tmpResFWD.llik;
//     // x_i    =  forwardProbUncertaintyMissing(&hmm, heteroEstResults , sizeChunk);
//     // x_i    =  backwardProbUncertaintyMissing(&hmm, heteroEstResults , sizeChunk);
//     cerr<<"..done"<<endl;
//     //cerr<<"x_i "<<x_i<<endl;
//     //return 1;
//     //cout<<setprecision(10)<<"\tinitial\t"<<h_i<<"\t"<<pT_i<<"\t"<<x_i<<"\t"<<endl;
//     vector<long double> hvector;
//     vector<long double> pvector;
//     vector<long double> svector;
    
//     cerr<<"Begin running MCMC on HMM using "<<thousandSeparator(maxChains)<<" chains"<<endl;

//     for(int chain=1;chain<=maxChains;chain++){

// 	//computing new state
// 	normal_distribution<long double> distribution_h(h_i,  (hupper-hlower)/partition  );
// 	h_i_1      = distribution_h(dre);
// 	if(h_i_1 <= hlower       ||  h_i_1 >= hupper     ){
// 	    h_i_1      = h_i;
// 	}


// 	normal_distribution<long double> distribution_pT(pT_i,(pTupper-pTlower)/partition  );
// 	pT_i_1      = distribution_pT(dre);
// 	if(pT_i_1 <= pTlower     ||  pT_i_1 >= pTupper     ){
// 	    pT_i_1      = pT_i;
// 	}

// 	//normal_distribution<long double> distribution_s(s_i,  (supper-slower)/partition  );
// 	normal_distribution<long double> distribution_s(s_i,  2  );
// 	s_i_1      = distribution_s(dre);
// 	if(s_i_1 <= slower       ||  s_i_1 >= supper     ){
// 	    s_i_1      = s_i;
// 	}


// 	//set new model parameters
// 	// cerr<<"old   h_i   "<<h_i<<" "<<" pT_i "<<" "<<pT_i<<" s_i  "<<s_i <<endl;
// 	// cerr<<"new   h_i_1 "<<h_i_1<<" "<<" pT_i_1 "<<" "<<pT_i_1<<" s_i_1  "<<s_i_1 <<endl;

// 	hmm.setHetRateForNonROH(h_i_1);
// 	hmm.setTransprob(pT_i_1);
// 	hmm.setNrwPerSizeChunk( (unsigned int)s_i_1 );
// 	hmm.recomputeProbs();

	
// 	//x_i_1=forwardProb(&hmm, emittedH , sizeChunk);
// 	//compute new likelihood
// 	//x_i_1=forwardProbUncertaintyMissing(&hmm, heteroEstResults , sizeChunk);
// 	tmpResFWD = forwardProbUncertaintyMissing(&hmm, heteroEstResults , sizeChunk,true);

// 	x_i_1     = tmpResFWD.llik;
// 	//	cerr<<"x_i_1 "<<x_i_1<<endl;
// 	if(chain>(maxChains/4)){
// 	    pTlower = pTlowerSecondHalf;
// 	}

// 	long double acceptance = min( (long double)(1.0)  , expl(x_i_1-x_i) );
// 	if( (long double)(randomProb()) < acceptance){
// 	    h_i           =  h_i_1;
// 	    s_i           =  s_i_1;
// 	    pT_i          =  pT_i_1;
// 	    x_i           =  x_i_1;

// 	    if( (chain>=(maxChains*(1.0-fracChainsBurnin)))){
// 		hvector.push_back(h_i);
// 		pvector.push_back(pT_i);
// 		svector.push_back(s_i);

// 		string strToWrite = stringify( x_i )+"\t"+stringify( h_i )+"\t"+stringify( s_i )+"\t"+stringify( pT_i )+"\t"+stringify( accept )+"\t"+stringify( chain )+"\t"+stringify( double(accept)/double(chain) )+"\n";
// 		bgzipWriterMCMC.Write(strToWrite.c_str(), strToWrite.size());
// 	    }
	    
// 	    accept++;
// 	    if(verbose)
// 	      cout<<setprecision(10)<<"accepted jump from\t"<<h_i<<"\t"<<s_i<<"\t"<<pT_i<<"\t"<<x_i<<"\tto\t"<<h_i_1<<"\t"<<s_i_1<<"\t"<<pT_i_1<<"\t"<<x_i_1<<""<<"\t"<<acceptance<<" "<<accept<<" "<<chain<<" "<<double(accept)/double(chain)<<endl;
// 	    //cerr<<setprecision(10)<<"mcmc"<<mcmc<<"\taccepted\t"<<h_i<<"\t"<<pT_i<<"\t"<<x_i<<"\t"<<" "<<acceptance<<" "<<accept<<" "<<chain<<" "<<double(accept)/double(chain)<<endl;	    
// 	}else{
	  
// 	  if(verbose)
// 	    cout<<setprecision(10)<<"rejected jump from\t"<<h_i<<"\t"<<s_i<<"\t"<<pT_i<<"\t"<<x_i<<"\tto\t"<<h_i_1<<"\t"<<s_i_1<<"\t"<<pT_i_1<<"\t"<<x_i_1<<""<<"\t"<<acceptance<<" "<<accept<<" "<<chain<<" "<<double(accept)/double(chain)<<endl;	    
// 	}
	
// 	if(!verbose)
// 	  printprogressBarCerr( float(chain)/float(maxChains) );
// 	//chain++;
// 	//sleep(0.1);
//     }
//     cerr<<endl;
//     cout<<setprecision(10)<<"mcmc"<<"\tfinal\t"<<h_i<<"\t"<<s_i<<"\t"<<pT_i<<"\t"<<x_i<<"\t"<<endl;

//     bgzipWriterMCMC.Close();    
	
//     //cerr<<"Baum Welch"<<endl;
//     //baum_welch(&hmm,&emittedH);
	

//     //hvector and pvector contain the values
//     long double hSum=0.0;
//     long double pSum=0.0;
//     long double sSum=0.0;
    
//     long double hAvg=0.0;
//     long double pAvg=0.0;
//     long double sAvg=0.0;
    
//     long double hMin=1.0;    
//     long double pMin=1.0;    
//     long double sMin=supper+1;
	
//     long double hMax=0.0;
//     long double pMax=0.0;
//     long double sMax=slower-1;

    

    
//     for(unsigned int i=0;i<hvector.size();i++){
// 	hSum += hvector[i];
// 	pSum += pvector[i];
// 	sSum += svector[i];

// 	if(hvector[i] < hMin)
// 	    hMin = hvector[i] ;
// 	if(hvector[i] > hMax)
// 	    hMax = hvector[i] ;
	
// 	if(pvector[i] < pMin)
// 	    pMin = pvector[i] ;
// 	if(pvector[i] > pMax)
// 	    pMax = pvector[i] ;
	
// 	if(svector[i] < sMin)
// 	    sMin = svector[i] ;
// 	if(svector[i] > sMax)
// 	    sMax = svector[i] ;
	
//     }

//     hAvg = hSum/( (long double)hvector.size() );
//     pAvg = pSum/( (long double)pvector.size() );
//     sAvg = sSum/( (long double)svector.size() );


//     //to remove
//     // hMin = 0.00070;
//     // hAvg = 0.0007467025205; 
//     // hMax = 0.0008;

//     // computing assignment prob
//     //set average parameters
//     //TODO remove
//     // hAvg = 0.0007467025205;
//     // pAvg = 0.0960255;
    
//     hmm.setHetRateForNonROH(hAvg);
//     hmm.setTransprob(pAvg);
//     hmm.setNrwPerSizeChunk( (unsigned int)sAvg );
//     hmm.recomputeProbs();

    
//     fbreturnVal postprob = forwardBackwardProbUncertaintyMissing(&hmm, heteroEstResults , sizeChunk,true);

//     cerr<<"HMM done"<<endl;


// }
