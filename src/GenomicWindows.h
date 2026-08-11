/*
 * GenomicWindows
 * Date: Aug-17-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#ifndef GenomicWindows_h
#define GenomicWindows_h

#include <stdlib.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sys/time.h>

#include "libgab.h"
#include "GenomicRange.h"
/* #include "RandomGenomicCoord.h" */


using namespace std;

typedef struct{
    string name;
    uint64_t startIndexChr;
    uint64_t endIndexChr;
    uint64_t length;
} chrinfo;


//Returns true if the sequence name is that of a sex chromosome. The match is on
//the whole name, not on a substring: scaffold/contig names routinely contain an
//"X" or a "Y" (e.g. the GenBank WGS accessions PVHY01000001.1) and those must not
//be mistaken for sex chromosomes.
inline bool isSexChrName(const string & chrName){
    return (chrName == "X"     || chrName == "Y"     ||
	    chrName == "x"     || chrName == "y"     ||
	    chrName == "chrX"  || chrName == "chrY"  ||
	    chrName == "chrx"  || chrName == "chry"  ||
	    chrName == "chr_X" || chrName == "chr_Y" ||
	    chrName == "chr_x" || chrName == "chr_y" );
}

class GenomicWindows{
 private:
    unsigned int genomeLength;
    vector<chrinfo> chrFound;
    bool allowSexChr;

public:
    GenomicWindows( );
    GenomicWindows(string fastaIndex,bool allowSexChr=false );
    ~GenomicWindows();

    vector<GenomicRange> getGenomicWindows(int windowSize,int overlap=0);
    vector<GenomicRange> getGenomicWindowsChr(string chrName, int windowSize,int overlap=0);


    vector<GenomicRange> getGenomeWide();
    vector<GenomicRange> getChr(string chrName);




};
#endif
