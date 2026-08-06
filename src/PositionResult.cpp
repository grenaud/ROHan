/*
 * PositionResult
 * Date: Mar-21-2016 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "PositionResult.h"

PositionResult::PositionResult(){
    //cerr<<"Constructor addr: "<<this<<endl;
}

PositionResult::~PositionResult(){
    //cerr<<"Destructor  addr: "<<this<<endl;
}

pair<char,char> PositionResult::hetIndex2Bases() const{
    if(mostLikelyGenoHetIdx==1)
	return make_pair<char,char>('A','C');
    if(mostLikelyGenoHetIdx==2)
	return make_pair<char,char>('A','G');
    if(mostLikelyGenoHetIdx==3)
	return make_pair<char,char>('A','T');
    if(mostLikelyGenoHetIdx==5)
	return make_pair<char,char>('C','G');
    if(mostLikelyGenoHetIdx==6)
	return make_pair<char,char>('C','T');
    if(mostLikelyGenoHetIdx==8)
	return make_pair<char,char>('G','T');
    cerr<<"PositionResult: wrong state hetIndex2Bases() "<<mostLikelyGenoHetIdx<<endl;
    exit(1);
    
}

int PositionResult::bases2hetIndex(char c1, char c2) const{
    if(c1>c2){
	char c=c1;
	c1=c2;
	c2=c;
    }

    if(c1=='A' && c2=='C')
	return 1;

    if(c1=='A' && c2=='G')
	return 2;

    if(c1=='A' && c2=='T')
	return 3;

    if(c1=='C' && c2=='G')
	return 5;

    if(c1=='C' && c2=='T')
	return 6;

    if(c1=='G' && c2=='T')
	return 8;

    cerr<<"PositionResult: wrong state bases2hetIndex() c1="<<c1<<" c2="<<c2<<endl;
    exit(1);
    
}


char PositionResult::homoIndex2Base() const{
    if(mostLikelyGenoIdx==0)
	return 'A';
    if(mostLikelyGenoIdx==4)
	return 'C';
    if(mostLikelyGenoIdx==7)
	return 'G';
    if(mostLikelyGenoIdx==9)
	return 'T';
    cerr<<"PositionResult: wrong state homoIndex2Base()"<<mostLikelyGenoIdx<<endl;
    exit(1);
    
}


int PositionResult::base2HomoIndex(const char b) const{
    if(	b == 'A')
	return 0;
    if(	b == 'C')
	return 4;
    if( b == 'G')
	return 7;
    if(	b == 'T')
	return 9;

    cerr<<"PositionResult: wrong state base2HomoIndex() b="<<b<<"#"<<endl;
    exit(1);    
}

//string PositionResult::toString(const RefVector  * references, const int & refID) const{
string PositionResult::toString(const bam_hdr_t * references,const int & refID) const{

    //cerr<<"toString: "<<endl;
    stringstream s;
    s.precision(13);
    // cerr<<refID<<endl;
    //string toReturn="";

    // cerr<<references->at(refID).RefName<<"\t";
    // cerr<<pos<<"\t";

    //s<<references->at(refID).RefName<<"\t";
    s<<references->target_name[refID]<<"\t";
    s<<pos<<"\t";
    //ID
    s<<"."<<"\t";
    
    // for(int n=0;n<4;n++)
    // 	s<<baseC[n]<<"\t";
    //REF
    s<<refB<<"\t";
    //cerr<<refB<<"\t";

    string altB;
    string gt="./.";
    stringstream pl;
    //cerr<<"mostLikelyGenoIdx="<<mostLikelyGenoIdx<<"\tmostLikelyGenoHetIdx="<<mostLikelyGenoHetIdx<<"\t";
    
    if(mostLikelyGenoIdx ==  mostLikelyGenoHetIdx){//heterozygous
	bool hasRefAsAlt=false;
	pair<char,char> hetBase = hetIndex2Bases();
	hasRefAsAlt = (refB == hetBase.first) || (refB == hetBase.second) ;
	
	if(hasRefAsAlt){
	    gt="0/1";
	    if( refB == hetBase.first )
		altB = stringify(hetBase.second);
	    else
		altB = stringify(hetBase.first);	    

	    //pl is 00, 01, 11
	    pl<<setprecision(0) << fixed << -1.0*ll[base2HomoIndex(refB)]<<","
		<<setprecision(0) << fixed << -1.0*ll[mostLikelyGenoIdx]<<","
		<<setprecision(0) << fixed << -1.0*ll[base2HomoIndex(altB[0])];

	}else{
	    gt="1/2";
	    altB = stringify(hetBase.first)+","+stringify(hetBase.second);

	    pl<<setprecision(0) << fixed << -1.0*ll[base2HomoIndex(altB[0])]<<","
	      <<setprecision(0) << fixed << -1.0*ll[mostLikelyGenoIdx]<<","
	      <<setprecision(0) << fixed << -1.0*ll[base2HomoIndex(altB[2])];
		
	}
    }else{//homozygous
	char bhomo=homoIndex2Base();
	bool hasRefAsAlt=false;
	//cerr<<"bhomo "<<bhomo<<" "<<"\t";

	if(bhomo != refB){ //is homo not reference
	    altB = stringify(bhomo);
	    gt="1/1";

	    pl<<setprecision(0) << fixed << -1.0*ll[base2HomoIndex(refB)]<<","
	      <<setprecision(0) << fixed << -1.0*ll[bases2hetIndex(refB,altB[0])]<<","
	      <<setprecision(0) << fixed << -1.0*ll[mostLikelyGenoIdx];

	}else{//is homo reference
	    gt="0/0";
	    pair<char,char> hetBase = hetIndex2Bases();//we are homozygous but finding the most likely heterozygous state
	    
	    hasRefAsAlt = (refB == hetBase.first) || (refB == hetBase.second) ;

	    if(hasRefAsAlt){//if the most likely het state has the ref
		if( refB == hetBase.first )
		    altB = stringify(hetBase.second);
		else
		    altB = stringify(hetBase.first);
		pl<<setprecision(0) << fixed << -1.0*ll[mostLikelyGenoIdx]<<","
		  <<setprecision(0) << fixed << -1.0*ll[bases2hetIndex(refB,altB[0])]<<","
		  <<setprecision(0) << fixed << -1.0*ll[base2HomoIndex(altB[0])];
		
	    }else{ //extreme corner case if homo ref but the most likely het. state does not have the ref. likely due to a statistical tie
		//neither the hetBase.first nor the hetBase.secondselecting the a random 
		//randomProb(false): rohan seeds rand() once in main(), a re-seed
		//here would undo the seed set by --seed
		if( randomProb(false) )
		    altB = stringify(hetBase.second);
		else
		    altB = stringify(hetBase.first);
		
		pl<<setprecision(0) << fixed << -1.0*ll[mostLikelyGenoIdx]<<","
		  <<setprecision(0) << fixed << -1.0*ll[bases2hetIndex(refB,altB[0])]<<","
		  <<setprecision(0) << fixed << -1.0*ll[base2HomoIndex(altB[0])];
		
		// cerr<<"ERROR at site "<<references->at(refID).RefName<<":"<<pos<<"\t"<<mostLikelyGenoIdx<<"\t"<<refB<<"\t"<<altB<<endl;
		// pl<<"NA"<<","
		//   <<"NA"<<","
		//   <<"NA"<<"1="<<hetBase.first<<"#2="<<hetBase.second<<"#";
	    }

	}
	//cerr<<endl;
    }

    //find if has ref
    s<<altB<<"\t";



    s<< setprecision(1) << fixed << -1.0*gq;
    s<<"\t.\t";//the . is the filter
    //INFO
    s<<"MQ="<<avgMQ<<";";
    s<<"AC4="<<baseC[0]<<","<<baseC[1]<<","<<baseC[2]<<","<<baseC[3];
    s<<"\t";//INFO field

    s<<"GT:DP:GQ:PL:PL10\t";//format field
    
    s<<gt<<":";
    s<<dp<<":";
    s<< setprecision(1) << fixed << -1.0*gq <<":";
    s<<pl.str()<<":";

    for(int g=0;g<9;g++)
	s<<setprecision(0) << fixed << -1.0*ll[g]<<",";
    s<<setprecision(0) << fixed << -1.0*ll[9]<<"";

    s<<"\n";

    //cout<<"toString "<<toReturn<<endl;
    return s.str();

    // s<<refC<<"\t";
    // s<<altC<<"\t";

    // if(geno==0){
    // 	s<<"0/0\t";
    // }else{
    // 	if(geno==1){
    // 	    s<<"0/1\t";
    // 	}else{
    // 	    if(geno==2){
    // 		s<<"1/1\t";
    // 	    }else{
    // 		cerr<<"Internal error for genotype ="<<geno<<endl;
    // 		exit(1);
    // 	    }	    
    // 	}
    // }

    // s<<genoS[0]<<genoS[1]<<"\t";

    // s<<lqual<<"\t";
    // s<<llCov<<"\t";

    //for(int g=0;g<10;g++)
    // 	s<<ll[g]<<"\t";	

    // // s<<rrll<<"\t";
    // // s<<rall<<"\t";
    // // s<<aall<<"\t";


}

//TODO code to VCF
// string PositionResult::toString(const RefVector  references) const{
//     //cerr<<"Constructor addr: "<<this<<endl;
//     string toReturn="";
//     toReturn += ""+references[refID].RefName+"\t";
//     toReturn += ""+stringify(pos)+"\t";
//     toReturn += ""+stringify(refB)+"\t";
//     toReturn += ""+stringify(altB)+"\t";
//     toReturn += ".\t"; //ID
//     toReturn += "0\t"; //QUAL   
//     toReturn += "\t";

//     toReturn += "\n";

//     //cout<<"toString "<<toReturn<<endl;
//     return toReturn;
// }

// ostream & operator << (ostream & os, const PositionResult & pr){


//     os<<pr.toString();


//     return os;
// }



    // if(refB=='A'){
    // 	hasRefAsAlt = hasRefAsAlt || (mostLikelyGenoHetIdx==1) || (mostLikelyGenoHetIdx==2) || (mostLikelyGenoHetIdx==3) ;
    // }else{
    // if(refB=='C'){
    // 	hasRefAsAlt = hasRefAsAlt || (mostLikelyGenoHetIdx==1) || (mostLikelyGenoHetIdx==5) || (mostLikelyGenoHetIdx==6) ;
    // }else{
    // if(refB=='G'){
    // 	hasRefAsAlt = hasRefAsAlt || (mostLikelyGenoHetIdx==2) || (mostLikelyGenoHetIdx==5) || (mostLikelyGenoHetIdx==8) ;
    // }else{
    // if(refB=='T'){
    // 	hasRefAsAlt = hasRefAsAlt || (mostLikelyGenoHetIdx==3) || (mostLikelyGenoHetIdx==6) || (mostLikelyGenoHetIdx==8) ;
    // }else{
    // }}}}
