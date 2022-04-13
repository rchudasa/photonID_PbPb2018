// This code is based on Aravid Sugunan's code 
// located on https://github.com/ats2008/BsMMGAnalysis/tree/photnIDdev/PhotonID/test
// Updated by Ruchi Chudasama on 13 April 2022
//
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "sstream"
#include "fstream"
#include "iostream"
using namespace std;


#include "TreeMaker.h" 

int  main(int argc,char *argv[])
{

    if(argc<2)
    {
        std::cout<<" Usage : \n"
                 <<"         ./main.exe <configfile.cfg> <input> <output> <DO_GEN>\n\n";
        exit(1);
    }
    int val(0);
    
    if(argc > 3)
    {
        val=atoi(argv[4]);        
    }

    std::cout<<"\n VAL = "<<val<<"\n";
    string cfgFile(argv[1]);
    string inFile(argv[2]);
    string outFile(argv[3]);
    
    TreeMaker aTreeeMaker;
    aTreeeMaker.Init(cfgFile);
    aTreeeMaker.SetupAnalysis(true, inFile, outFile);
    
    if(val == 0) aTreeeMaker.genParticleSCMaker();
    if(val == 1) aTreeeMaker.Pi0ParticleSCMaker();
    if(val == 2) aTreeeMaker.DataSCMaker();
    aTreeeMaker.SaveFile();

}
