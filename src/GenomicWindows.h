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
