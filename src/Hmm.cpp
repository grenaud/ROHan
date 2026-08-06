/*
 * Hmm
 * Date: Dec-05-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "Hmm.h"


Hmm::Hmm(int minSegSitesPerChunk_,int maxSegSitesPerChunk_,int sizeChunk,unsigned int nrwPerSizeChunk){
    //cerr<<"Hmm constr"<<endl;
    minSegSitesPerChunk =  minSegSitesPerChunk_;
    maxSegSitesPerChunk =  maxSegSitesPerChunk_;

    for(int n=0;n<NBRSTATES;n++){
	if(n==0)
	    hmmstates[n] = new  HmmState (n,0.000000012,minSegSitesPerChunk,maxSegSitesPerChunk,sizeChunk,nrwPerSizeChunk);
	if(n==1)
	    hmmstates[n] = new  HmmState (n,0.00072    ,minSegSitesPerChunk,maxSegSitesPerChunk,sizeChunk,nrwPerSizeChunk);
    }

    for(int n=0;n<NBRSTATES;n++){
	if(n==0)
	    hmmstates[n]->setSecond(hmmstates[n+1]);
	if(n==1)
	    hmmstates[n]->setSecond(hmmstates[n-1]);
    }
    // for(int m=0;m<1000;m++)
    // 	cout<<m<<"\tp="<<hmmstates[1]->probEmission(m,1000000)<<endl;
    long double pTrans = 0.1;
    trans = new long double*[NBRSTATES];
    
    for(int i = 0; i < NBRSTATES; ++i)
	trans[i] = new long double[NBRSTATES];
    
    for(int n=0;n<NBRSTATES;n++){
	probTrans[n] =   pTrans;
	probStay[n]  = 1-pTrans;
    }
    trans[0][0]  = 1-pTrans;	
    trans[0][1]  =   pTrans;
    trans[1][0]  =   pTrans;
    trans[1][1]  = 1-pTrans;	

    
    for(int n=0;n<NBRSTATES;n++){
	startingState[n] = 1/(  (long double)NBRSTATES );
    }
    // HmmState roh    (0.0007);
    // HmmState normal (0.000000012);
    // roh.setSecond(&normal);
    // normal.setSecond(&roh);

    // hmmstates[0]=roh;
    // hmmstates[1]=normal;

    // double pTrans = 0.1;
    // probTrans[0]=pTrans;
    // probTrans[1]=pTrans;

    // probStay[0]=1-pTrans;
    // probStay[1]=1-pTrans;
    
    //    cout<<"test"<<endl;

    // for(int m=0;m<1500;m++){
    // 	cout<<m<<"\t"<<hmmstates[1]->probEmission(m,1000000)<<endl;       
    // }
    //exit(1);
}

Hmm::Hmm(const Hmm & other){
    if(NBRSTATES != other.getNumberStates()){
	cerr<<"Discrepency between number of states"<<endl;
	exit(1);
    }

    for(int n=0;n<NBRSTATES;n++){
	startingState[n] = other.startingState[n];
	probTrans[n]     = other.probTrans[n];
	probStay[n]      = other.probStay[n];
	
	hmmstates[n]     = new  HmmState ( *other.hmmstates[n] );
    }
    
     for(int n=0;n<NBRSTATES;n++){
	if(n==0)
	    hmmstates[n]->setSecond(hmmstates[n+1]);
	if(n==1)
	    hmmstates[n]->setSecond(hmmstates[n-1]);
    }
}

Hmm::~Hmm(){

    for(int n=0;n<NBRSTATES;n++){
	delete hmmstates[n];
    }

    for(int i = 0; i < NBRSTATES; ++i)
      delete 	trans[i] ;
    delete trans;
}


int  Hmm::getMinSegSitesPerChunk(){
    return minSegSitesPerChunk;
}
int  Hmm::getMaxSegSitesPerChunk(){
    return maxSegSitesPerChunk;
}
 

vector<emission> Hmm::generateStates(unsigned int N,unsigned int total) const{
    vector<emission> vecToReturn;
    //int initState = randomInt(0,NBRSTATES-1);
    int initState = 1;
    HmmState * currrentState = hmmstates[ initState ];
    
    for(unsigned int i=0;i<N;i++){
	//cerr<<"generateStates "<<i<<endl;
	unsigned int d_ = currrentState->randomEmission(total);
	emission eToAdd;
	eToAdd.p     =  ( (long double)d_ )/( (long double)total );
	eToAdd.idx   = currrentState->getIdx();
	
	vecToReturn.push_back( eToAdd );
	//randomProb(false): rohan seeds rand() once in main(), a re-seed here
	//would undo the seed set by --seed
	long double pt = randomProb(false);
	//cout<<"state#"<<currrentState->getIdx()<<"\t"<<d<<"\t"<<pt<<"\t"<<probTrans[currrentState->getIdx()]<<endl;	

	//if(0)
	if(pt < probTrans[currrentState->getIdx()]){//transit
	    if(currrentState->getIdx() == 0)
		currrentState = hmmstates[1];
	    else
		if(currrentState->getIdx()     == 1)
		    currrentState  = hmmstates[0];
	}else{
	    //do nothing
	}
    }
    return vecToReturn;
}

int  Hmm::getNumberStates() const{
    return NBRSTATES;
}


long double Hmm::getTrans(int i,int j) const{
    return trans[i][j];
}

void Hmm::setHetRateForNonROH(long double newH){
    hmmstates[1]->setH(newH);
}

void Hmm::setHetRateForROH(long double newH){
    hmmstates[0]->setH(newH);
}

void Hmm::setNrwPerSizeChunk(unsigned int nrwPerSizeChunk){
    hmmstates[0]->setNrwPerSizeChunk(nrwPerSizeChunk);
    hmmstates[1]->setNrwPerSizeChunk(nrwPerSizeChunk);
}

void Hmm::recomputeProbsROH(bool verbose){
    hmmstates[0]->recomputeProbs(verbose);
}

void Hmm::recomputeProbsNonROH(bool verbose){
    hmmstates[1]->recomputeProbs(verbose);
}

void Hmm::setTransprob(long double newTrans){
    long double pTrans = newTrans;    
    for(int n=0;n<NBRSTATES;n++){
	probTrans[n] =   pTrans;
	probStay[n]  = 1-pTrans;
    }
    trans[0][0]  = 1-pTrans;	
    trans[0][1]  =   pTrans;
    trans[1][0]  =   pTrans;
    trans[1][1]  = 1-pTrans;	

}

std::ostream& operator<<(ostream& os, const Hmm& hmm)  {
    os<<"# of states"<<hmm.getNumberStates()<<endl;
    os<<"--------"<<endl;

    for(int s=0;s<hmm.getNumberStates();s++){
    	os<<"p start      at state#"<<s<<" = "<<hmm.startingState[s]<<endl;
    }
    os<<"--------"<<endl;
    
    for(int s=0;s<hmm.getNumberStates();s++){
    	os<<"p transition at state#"<<s<<" = "<<hmm.probTrans[s]<<endl;
    }
    os<<"--------"<<endl;
    for(int s=0;s<hmm.getNumberStates();s++){
    	os<<"p staying    at state#"<<s<<" = "<<hmm.probStay[s]<<endl;
    }
    os<<"--------"<<endl;    
    for(int s=0;s<hmm.getNumberStates();s++){
    	os<<"h            at state#"<<s<<" = "<<hmm.hmmstates[s]->getH()<<endl;
    }
    os<<"--------"<<endl;
    for(int s=0;s<hmm.getNumberStates();s++){
    	os<<"theta        at state#"<<s<<" = "<<hmm.hmmstates[s]->getTheta()<<endl;
    }
    os<<"--------"<<endl;
    for(int s=0;s<hmm.getNumberStates();s++){
    	os<<"rate geom.   at state#"<<s<<" = "<<hmm.hmmstates[s]->getRateGeom()<<endl;
    }
    

    return os;
}




