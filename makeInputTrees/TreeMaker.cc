void TreeMaker::DataSCMaker()
{


    AddSCTree("dataAllSCTree");
    
    
    std::cout<<"\nBegining Analysis Script !";
    if (maxEvents >0 ) maxEvents = nentries > maxEvents ? maxEvents : nentries;
    cout<<"\nProcessing total "<<maxEvents<<" events \n\n";
   
    Long64_t EventCount=0;
    Long64_t EventCountWithCand=0;
    Long64_t nCands=0;
    
    Long64_t nb = 0,nbytes=0 ;

    auto t_start = std::chrono::high_resolution_clock::now();
    auto t_end = std::chrono::high_resolution_clock::now();

    Double_t dr=-1.0;
    Double_t drMin;

    Bool_t hasACand=false;
    
    for (Long64_t jentry=0; jentry<maxEvents; jentry++)
    {  
       Long64_t ientry_evt = ntupleRawTree.LoadTree(jentry);
       if (ientry_evt < 0) break;
       nb = ntupleRawTree.fChain->GetEntry(jentry);   nbytes += nb;
       
       if(jentry%reportEvery == 0 )
       {
             t_end = std::chrono::high_resolution_clock::now();
             std::cout<<"Processing Entry in event loop : "<<jentry<<" / "<<maxEvents<<"  [ "<<100.0*jentry/maxEvents<<"  % ]  "
                      << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
                      <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( maxEvents - jentry)/(1e-9 + jentry)* 0.001
                      <<std::endl;
       
       }
       
       EventCount++;
       hasACand=false;
       for( Int_t j =0 ;j< ntupleRawTree.nPho ;j++)
       {
            if(abs(ntupleRawTree.phoEta->at(j)) < scAbsEtaMin ) continue;
            if(abs(ntupleRawTree.phoEta->at(j)) > scAbsEtaMax ) continue;
        
            drMin=-1.0;
            hasACand=true;
            fillSCVariablesToOutTree(j,"dataAllSCTree");
            nCands++;
       }

       if(hasACand) EventCountWithCand++;

    }

    std::cout<<" Number of Evnets processed        : "<<EventCount<<"\n";
    std::cout<<" Number of Evnets with candidates  : "<<EventCountWithCand<<"\n";
    std::cout<<" Number of candidates              : "<<nCands<<"\n";

}
void TreeMaker::Pi0ParticleSCMaker()
{
    AddSCTree("mergedPi0_SCTree");
    AddSCTree("leadPi0_SCTree");
    AddSCTree("subLeadPi0Gamma_SCTree");
    Double_t dr;
    
    std::cout<<"\nBegining Pi0 Analysis Script !";
    if (maxEvents >0 ) maxEvents = nentries > maxEvents ? maxEvents : nentries;
    cout<<"\nProcessing total "<<maxEvents<<" events , reporting every "<<reportEvery<<" events\n\n";
   
    Long64_t EventCount=0;
    Long64_t EventCountWithCand=0;
    Long64_t nCands(0),nMergedCands(0),nLeadGammaCands(0),nSubLeadGammaCands(0);
    Long64_t nCandsFromInclusion(0), nCandsFormSameSC(0);
    Long64_t nb = 0,nbytes=0 ;

    auto t_start = std::chrono::high_resolution_clock::now();
    auto t_end = std::chrono::high_resolution_clock::now();
    bool goodRunLumi = false;

    Double_t drMin;
    Int_t scMatchIdx,scG1MatchIdx,scG2MatchIdx,tempI;
    Int_t g1MCIdx,g2MCIdx;
    Bool_t foundmatch=false;
    TLorentzVector g1,g2,pi0;
    Double_t g1g2DR,drMinG1,drMinG2,tempD;
    Bool_t isMerged;

    th1fStore["gen_p0daugterGammaGammaDR"] = new TH1F("gen_p0daugterGammaGammaDR","gen_p0daugterGammaGammaDR",100,0.0,2.0);
    th1fStore["gen_Pi0Gamma1DR"] = new TH1F("gen_Pi0Gamma1DR","gen_Pi0Gamma1DR",100,0.0,2.0);
    th1fStore["gen_Pi0Gamma2DR"] = new TH1F("gen_Pi0Gamma2DR","gen_Pi0Gamma2DR",100,0.0,2.0);
    

    //for (Long64_t jentry=0; jentry<2; jentry++)
    for (Long64_t jentry=0; jentry<maxEvents; jentry++)
    {  
       eventGenMultiplicity=0;
       Long64_t ientry_evt = ntupleRawTree.LoadTree(jentry);
       if (ientry_evt < 0) break;
       nb = ntupleRawTree.fChain->GetEntry(jentry);   nbytes += nb;
       
       if(jentry%reportEvery == 0 )
       {
             t_end = std::chrono::high_resolution_clock::now();
             std::cout<<"Processing Entry in event loop : "<<jentry<<" / "<<maxEvents<<"  [ "<<100.0*jentry/maxEvents<<"  % ]  "
                      << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
                      <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( maxEvents - jentry)/(1e-9 + jentry)* 0.001
                      <<std::endl;
       
       }
       
       EventCount++;
       foundmatch=false;
        
       if(isMC)
       for(int i=0;i < ntupleRawTree.nMC ; i++)
       {
           //std::cout << "\n"; 
           //std::cout<<"Entry:"<<jentry << "\t" << i<<"/"<<ntupleRawTree.nMC<<"pdg id : "<<(ntupleRawTree.mcPID)->at(i)<<" pt : "<<(ntupleRawTree.mcPt)->at(i)<<" eta : "<<(ntupleRawTree.mcEta)->at(i)<<" phi : "<<ntupleRawTree.mcPhi->at(i)<<"\n";
            if(ntupleRawTree.mcPID->at(i) != genParticlePDGID ) continue;
            if(ntupleRawTree.mcPt->at(i)  < genParticlePtMin ) continue;
            if(genParticleIsStable) if( ntupleRawTree.mcStatus->at(i) != 1 ) continue;
                
            pi0.SetPtEtaPhiM(ntupleRawTree.mcPt->at(i) , ntupleRawTree.mcEta->at(i) , ntupleRawTree.mcPhi->at(i) , ntupleRawTree.mcMass->at(i) );
            drMin=1e9;
            for(int idx=0;idx < ntupleRawTree.nMC ; idx++)
            {
                if(ntupleRawTree.mcPID->at(idx) != 22 ) continue;
                g1.SetPtEtaPhiM(ntupleRawTree.mcPt->at(idx) , ntupleRawTree.mcEta->at(idx) , ntupleRawTree.mcPhi->at(idx) , ntupleRawTree.mcMass->at(idx) );
                 //std::cout<<"\t\tpdg id : "<<(ntupleRawTree.mcPID)->at(idx)<<" pt : "<<(ntupleRawTree.mcPt)->at(idx)<<" eta : "<<(ntupleRawTree.mcEta)->at(idx)<<" phi : "<<ntupleRawTree.mcPhi->at(idx)<<"\n";
                for(int jdx=idx+1;jdx < ntupleRawTree.nMC ; jdx++)
                {
                    if(ntupleRawTree.mcPID->at(jdx) != 22 ) continue;
                    g2.SetPtEtaPhiM(ntupleRawTree.mcPt->at(jdx) , ntupleRawTree.mcEta->at(jdx) , ntupleRawTree.mcPhi->at(jdx) , ntupleRawTree.mcMass->at(jdx) );
                    dr=pi0.DeltaR(g1+g2);
                    //std::cout<<"\t\t\t\tpdg id : "<<(ntupleRawTree.mcPID)->at(jdx)<<" pt : "<<(ntupleRawTree.mcPt)->at(jdx)<<" eta : "<<(ntupleRawTree.mcEta)->at(jdx)<<" phi : "<<ntupleRawTree.mcPhi->at(jdx)<<" dr : "<<dr<<"\n";
                    if(dr<drMin)
                    {
                       drMin=dr;
                       g1MCIdx = idx; 
                       g2MCIdx = jdx; 
                    }
                    
                }
            }


            if( drMin < 0.001 ) {
           //std::cout<<"Entry:"<<jentry << "\t g1 ID:" << g1MCIdx<<"/"<<" pt : "<<(ntupleRawTree.mcPt)->at(g1MCIdx)<<" eta : "<<(ntupleRawTree.mcEta)->at(g1MCIdx)<<" phi : "<<ntupleRawTree.mcPhi->at(g1MCIdx)<<"\n";
           //std::cout<<"Entry:"<<jentry << "\t g2 ID:" << g2MCIdx<<"/"<<" pt : "<<(ntupleRawTree.mcPt)->at(g2MCIdx)<<" eta : "<<(ntupleRawTree.mcEta)->at(g2MCIdx)<<" phi : "<<ntupleRawTree.mcPhi->at(g2MCIdx)<<"\n";

              g1.SetPtEtaPhiM(ntupleRawTree.mcPt->at(g1MCIdx) , ntupleRawTree.mcEta->at(g1MCIdx) , ntupleRawTree.mcPhi->at(g1MCIdx) , ntupleRawTree.mcMass->at(g1MCIdx) );
              g2.SetPtEtaPhiM(ntupleRawTree.mcPt->at(g2MCIdx) , ntupleRawTree.mcEta->at(g2MCIdx) , ntupleRawTree.mcPhi->at(g2MCIdx) , ntupleRawTree.mcMass->at(g2MCIdx) );
              g1g2DR=g1.DeltaR(g2);
              th1fStore["gen_p0daugterGammaGammaDR"]->Fill(g1g2DR);
              th1fStore["gen_Pi0Gamma1DR"]->Fill(pi0.DeltaR(g1));
              th1fStore["gen_Pi0Gamma2DR"]->Fill(pi0.DeltaR(g2));
            }
            else
            {
                continue;
            }

 //           continue;
            
            drMinG1=drGenMatchMin;
            drMinG2=drGenMatchMin;

            double detaMin1 = 0.1;
	    double detaMin2 = 0.1;
            scG1MatchIdx=-1;
            scG2MatchIdx=-1;
            
            isMerged=false;

          // std::cout << jentry << " Number of superclusters :" << ntupleRawTree.nPho << "\n";
            for( Int_t j =0 ;j< ntupleRawTree.nPho ;j++)
            {
                if(abs(ntupleRawTree.phoEta->at(j)) < scAbsEtaMin ) continue;
                if(abs(ntupleRawTree.phoEta->at(j)) > scAbsEtaMax ) continue;

                    // match each Supercluster to both MC photons
                    // matching with first Gen photon
		    Double_t drSC1 =getDR(ntupleRawTree.mcEta->at(g1MCIdx),ntupleRawTree.mcPhi->at(g1MCIdx), ntupleRawTree.phoEta->at(j),ntupleRawTree.phoPhi->at(j));
 		    Double_t deta1 = getDETA(ntupleRawTree.mcEta->at(g1MCIdx),ntupleRawTree.phoEta->at(j));
 		    //std::cout<<"\t\t: "<< j <<"/" << ntupleRawTree.nPho << " SC et: " << (ntupleRawTree.phoEt)->at(j) << " SC eta : "<<(ntupleRawTree.phoEta)->at(j)<<" SC phi : "<<ntupleRawTree.phoPhi->at(j)<<" dr1 : "<<drSC1 << " deta1:  " << deta1 << "\n";
                    if(drSC1<drMinG1)
                     {
 			if(deta1 < detaMin1) {
		        detaMin1 = deta1;
                        drMinG1=drSC1;
                        scG1MatchIdx=j;
			}
                    }

                    // matching with second Gen photon
                    Double_t drSC2=getDR(ntupleRawTree.mcEta->at(g2MCIdx),ntupleRawTree.mcPhi->at(g2MCIdx), ntupleRawTree.phoEta->at(j),ntupleRawTree.phoPhi->at(j));
 		    Double_t deta2 = getDETA(ntupleRawTree.mcEta->at(g2MCIdx),ntupleRawTree.phoEta->at(j));
 		   // std::cout<<"\t\t: "<< j <<"/" << ntupleRawTree.nPho << " SC et: " << (ntupleRawTree.phoEt)->at(j) << " SC eta : "<<(ntupleRawTree.phoEta)->at(j)<<" SC phi : "<<ntupleRawTree.phoPhi->at(j)<<" dr2 : "<<drSC2<< " deta2:  " << deta2 << "\n";
                    if(drSC2<drMinG2)
                    {
		        if(deta2 < detaMin2){
			detaMin2 = deta2;
                        drMinG2=drSC2;
                        scG2MatchIdx=j;
			}
                    }
            } //supercluster loop
            
	    // check the Index of the Supercluster matched to gen photon , if the same index is matched with 2 gen photon then the SC is merged.
            if( (scG1MatchIdx==scG2MatchIdx) and (scG2MatchIdx !=-1) )
            {   
 		//std::cout<<"Entry: "<< jentry << " " << scG1MatchIdx << "  " << scG2MatchIdx << "Both SCs are merged" <<  "\n";
                isMerged = true;
                nCandsFormSameSC++;
            }
            else
            {
                if( scG1MatchIdx > -1 )
                {
                    if(scG2MatchIdx > -1)
                    {
                        if( ntupleRawTree.phoEt->at(scG1MatchIdx) < ntupleRawTree.phoEt->at(scG2MatchIdx))
                        {
                            tempI = scG1MatchIdx;
                            scG1MatchIdx=scG2MatchIdx;
                            scG2MatchIdx=tempI;

                            tempD = drMinG1 ;
                            drMinG1 = drMinG2;
                            drMinG2= tempD;

                            tempI=g1MCIdx;
                            g1MCIdx=g2MCIdx;
                            g2MCIdx=tempI;
                        }
                        isMerged=false;
                    }
                }
                else
                {
                    if(scG2MatchIdx > -1)
                    {
                            tempI = scG1MatchIdx;
                            scG1MatchIdx=scG2MatchIdx;
                            scG2MatchIdx=tempI;

                            tempD = drMinG1 ;
                            drMinG1 = drMinG2;
                            drMinG2= tempD;

                            tempI=g1MCIdx;
                            g1MCIdx=g2MCIdx;
                            g2MCIdx=tempI;
                    }
                }

                if(scG1MatchIdx > -1 and scG2MatchIdx < 0 )
                {
                    if( 0.5*( sqrt( 
                                ntupleRawTree.phoSCEtaWidth->at(scG1MatchIdx)*ntupleRawTree.phoSCEtaWidth->at(scG1MatchIdx) 
                                + ntupleRawTree.phoSCPhiWidth->at(scG1MatchIdx)*ntupleRawTree.phoSCPhiWidth->at(scG1MatchIdx) 
                                ))  
                          > 
                           getDR(ntupleRawTree.phoEta->at(scG1MatchIdx) , ntupleRawTree.phoPhi->at(scG1MatchIdx) , ntupleRawTree.mcEta->at(g2MCIdx) , ntupleRawTree.mcPhi->at(g2MCIdx) ) 
                      ) 
                    {
 			   std::cout<<"Entry: "<< jentry << " " << scG1MatchIdx << "  " << scG2MatchIdx << "Inclusion  merged" <<  "\n";
                            nCandsFromInclusion++;
                           isMerged=true; 
                    }
                    else
                    {
                           isMerged=false;
                    }
                }
                else
                {
                      isMerged=false;
                }
                
            }

            if(isMerged) 
            { 
                
                fillSCVariablesToOutTree(scG1MatchIdx,"mergedPi0_SCTree");
                nMergedCands++;
                foundmatch=true;
            }
            else
            {
                if(scG1MatchIdx > -1 ) 
                {
                    
                    fillSCVariablesToOutTree(scG1MatchIdx,"leadPi0_SCTree");
                    nLeadGammaCands++;
                    foundmatch=true;
                }
                if(scG2MatchIdx > -1 ) 
                {
                    
                    fillSCVariablesToOutTree(scG2MatchIdx,"subLeadPi0Gamma_SCTree");
                    nSubLeadGammaCands++;
                    foundmatch=true;
                }
                if(scG1MatchIdx > -1 and scG2MatchIdx > -1 ) 
                {
                    nCands++;
                }
            }
        
       }
       if(foundmatch==true)     EventCountWithCand++;
       
       for( Int_t j =0 ;j< ntupleRawTree.nPho ;j++)
       {
            if(abs(ntupleRawTree.phoEta->at(j)) < scAbsEtaMin ) continue;
            if(abs(ntupleRawTree.phoEta->at(j)) > scAbsEtaMax ) continue;

            drMin=2.99;
            if(isMC)
            for(int i=0;i < ntupleRawTree.nMC ; i++)
            {
                 if( ntupleRawTree.mcPID->at(i) != genParticlePDGID ) continue;
                 if(genParticleIsStable) if( ntupleRawTree.mcStatus->at(i) != 1 ) continue;
                 if(ntupleRawTree.mcPt->at(i)  < genParticlePtMin ) continue;
                 dr=getDR(ntupleRawTree.mcEta->at(i),ntupleRawTree.mcPhi->at(i), ntupleRawTree.phoEta->at(j),ntupleRawTree.phoPhi->at(j));
                 if(dr<drMin)
                 {
                        drMin=dr;
                 }
            }
            
       }
      
    }

    std::cout<<" Number of Evnets processed                      : "<<EventCount<<"\n";
    std::cout<<" Number of Evnets with candidates                : "<<EventCountWithCand<<"\n";
    std::cout<<" Number of Events with merged candidates         : "<<nMergedCands<<"\n";
    std::cout<<"            single SC Matched candidates         : "<<nCandsFormSameSC<<"\n";
    std::cout<<"          Incusive SC Matched candidates         : "<<nCandsFromInclusion<<"\n";
    std::cout<<" Number of Events with lead candidates     [NM]  : "<<nLeadGammaCands<<"\n";
    std::cout<<" Number of Events with sub lead candidates [NM]  : "<<nSubLeadGammaCands<<"\n";
    std::cout<<" Number of Events with both candidates     [NM]  : "<<nCands<<"\n";
}

void TreeMaker::genParticleSCMaker()
{

    AddSCTree("genMatchedSCTree");
    
    Double_t dr;
    
    std::cout<<"\nBegining Analysis Script !";
    if (maxEvents >0 ) maxEvents = nentries > maxEvents ? maxEvents : nentries;
    cout<<"\nProcessing total "<<maxEvents<<" events \n\n";
   
    Long64_t EventCount=0;
    Long64_t EventCountWithCand=0;
    Long64_t nCands=0;
    
    Long64_t nb = 0,nbytes=0 ;

    auto t_start = std::chrono::high_resolution_clock::now();
    auto t_end = std::chrono::high_resolution_clock::now();
    bool goodRunLumi = false;

    Double_t drMin;
    Int_t scMatchIdx;
    Bool_t foundmatch=false;
    for (Long64_t jentry=0; jentry<maxEvents; jentry++)
    {  
       eventGenMultiplicity=0;
       Long64_t ientry_evt = ntupleRawTree.LoadTree(jentry);
       if (ientry_evt < 0) break;
       nb = ntupleRawTree.fChain->GetEntry(jentry);   nbytes += nb;
       
       std::cout << "Entry:" << jentry << "**********************Load tree "<< ntupleRawTree.nMC << std::endl;

       if(jentry%reportEvery == 0 )
       {
             t_end = std::chrono::high_resolution_clock::now();
             std::cout<<"Processing Entry in event loop : "<<jentry<<" / "<<maxEvents<<"  [ "<<100.0*jentry/maxEvents<<"  % ]  "
                      << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
                      <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( maxEvents - jentry)/(1e-9 + jentry)* 0.001
                      <<std::endl;
       
       }
       
       EventCount++;
       foundmatch=false;

       if(isMC)
       for(int i=0;i < ntupleRawTree.nMC ; i++)
       {

	   if(abs(ntupleRawTree.mcEta->at(i)) < scAbsEtaMin ) continue;
           if(abs(ntupleRawTree.mcEta->at(i)) > scAbsEtaMax ) continue;
            std::cout<<"\tpdg id : "<<(ntupleRawTree.mcPID)->at(i)<<" pt : "  \
                                    <<(ntupleRawTree.mcPt)->at(i)<<" eta : "<<(ntupleRawTree.mcEta)->at(i)<<" phi : "<<ntupleRawTree.mcPhi->at(i) \
                                    <<" status : "<<ntupleRawTree.mcStatus->at(i) \
                                    <<"\n";
            if(ntupleRawTree.mcPID->at(i) != genParticlePDGID )      continue;
            if(genParticleIsStable)   if( ntupleRawTree.mcStatus->at(i) != 1 )          continue;
            
            drMin=drGenMatchMin;
            scMatchIdx=-1;
            for( Int_t j =0 ;j< ntupleRawTree.nPho ;j++)
            {
                if(abs(ntupleRawTree.phoEta->at(j)) < scAbsEtaMin ) continue;
                if(abs(ntupleRawTree.phoEta->at(j)) > scAbsEtaMax ) continue;

                    dr=getDR(ntupleRawTree.mcEta->at(i),ntupleRawTree.mcPhi->at(i), ntupleRawTree.phoEta->at(j),ntupleRawTree.phoPhi->at(j));
                    if(dr<drMin)
                    {
                        drMin=dr;
                        scMatchIdx=j;
                    }
            }
            

            if( scMatchIdx > -1)
            {
                cout << jentry << "  matched " << endl;
                fillSCVariablesToOutTree(scMatchIdx,"genMatchedSCTree");
 	        cout<< "filled sc variables" << endl;
                nCands++;
                foundmatch=true;
            }
       }
       if(foundmatch==true)     EventCountWithCand++;
       
 	std::cout << jentry << "  just before photon loop " << std::endl;
       for( Int_t j =0 ;j< ntupleRawTree.nPho ;j++)
       {
		               
            if(abs(ntupleRawTree.phoEta->at(j)) < scAbsEtaMin ) continue;
            if(abs(ntupleRawTree.phoEta->at(j)) > scAbsEtaMax ) continue;

            drMin=2.99;
            if(isMC)
            for(int i=0;i < ntupleRawTree.nMC ; i++)
            {
                 if( ntupleRawTree.mcPID->at(i) != genParticlePDGID ) continue;
                 if(genParticleIsStable) if(not ntupleRawTree.mcStatus->at(i) != 1 ) continue;
                 dr=getDR(ntupleRawTree.mcEta->at(i),ntupleRawTree.mcPhi->at(i), ntupleRawTree.phoEta->at(j),ntupleRawTree.phoPhi->at(j));
                 if(dr<drMin)
                 {
                        drMin=dr;
                 }
            }
           
            

       }

    }

    std::cout<<" Number of Evnets processed        : "<<EventCount<<"\n";
    std::cout<<" Number of Evnets with candidates  : "<<EventCountWithCand<<"\n";
    std::cout<<" Number of candidates              : "<<nCands<<"\n";
}


void TreeMaker::setupOutputSCTree()
{
}


void TreeMaker::AddSCTree(TString SCTreeName)
{
    auto outSC_Tree = new TTree(SCTreeName,"variables for the PhotonID from SC");

    Int_t idx(0), offset(candidateMapDouble["SCTreeStorage"]);
    outSC_Tree->Branch("phoE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoEt"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoEta"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoPhi"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoSCE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoSCEt"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoSCRawE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoSCEta"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoSCPhi"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoSCEtaWidth"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoSCPhiWidth"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoSCBrem"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoSCnBC"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoESEn"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoPSCE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoPSCRawE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoPSCEta"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoPSCPhi"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoPSCEtaWidth"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoPSCPhiWidth"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoPSCBrem"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoPSCnHits"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoPSCflags"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoPSCinClean"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoPSCinUnClean"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoPSCnBC"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoPESEn"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoIsPFPhoton"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoIsStandardPhoton"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoHasPixelSeed"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoHasConversionTracks"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoHadTowerOverEm"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoHoverE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoHoverEValid"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoSigmaIEtaIEta"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoR9"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoE1x5"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoE2x5"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoE3x3"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoE5x5"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoMaxEnergyXtal"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoSigmaEtaEta"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoSigmaIEtaIEta_2012"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoR9_2012"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoE1x5_2012"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoE2x5_2012"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoE3x3_2012"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoE5x5_2012"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoMaxEnergyXtal_2012"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoSigmaEtaEta_2012"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoHadTowerOverEm1"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoHadTowerOverEm2"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoHoverE1"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoHoverE2"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoSigmaIEtaIPhi"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoSigmaIPhiIPhi"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoR1x5"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoR2x5"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoE2nd"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoETop"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoEBottom"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoELeft"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoERight"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoE1x3"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoE2x2"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoE2x5Max"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoE2x5Top"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoE2x5Bottom"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoE2x5Left"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoE2x5Right"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoSigmaIEtaIPhi_2012"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoSigmaIPhiIPhi_2012"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoR1x5_2012"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoR2x5_2012"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoE2nd_2012"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoETop_2012"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoEBottom_2012"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoELeft_2012"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoERight_2012"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoE1x3_2012"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoE2x2_2012"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoE2x5Max_2012"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoE2x5Top_2012"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoE2x5Bottom_2012"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoE2x5Left_2012"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoE2x5Right_2012"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoBC1E"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoBC1Ecorr"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoBC1Eta"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoBC1Phi"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoBC1size"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoBC1flags"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoBC1inClean"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoBC1inUnClean"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("phoBC1rawID"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pho_seedTime"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfcIso1"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfcIso2"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfcIso3"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfcIso4"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfcIso5"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfpIso1"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfpIso2"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfpIso3"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfpIso4"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfpIso5"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfnIso1"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfnIso2"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfnIso3"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfnIso4"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfnIso5"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfpIso1subSC"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfpIso2subSC"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfpIso3subSC"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfpIso4subSC"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfpIso5subSC"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfcIso1subUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfcIso2subUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfcIso3subUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfcIso4subUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfcIso5subUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfpIso1subUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfpIso2subUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfpIso3subUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfpIso4subUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfpIso5subUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfnIso1subUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfnIso2subUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfnIso3subUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfnIso4subUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfnIso5subUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfpIso1subSCsubUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfpIso2subSCsubUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfpIso3subSCsubUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfpIso4subSCsubUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfpIso5subSCsubUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfcIso1pTgt1p0subUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfcIso2pTgt1p0subUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfcIso3pTgt1p0subUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfcIso4pTgt1p0subUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfcIso5pTgt1p0subUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfcIso1pTgt2p0subUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfcIso2pTgt2p0subUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfcIso3pTgt2p0subUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfcIso4pTgt2p0subUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfcIso5pTgt2p0subUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfcIso1pTgt3p0subUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfcIso2pTgt3p0subUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfcIso3pTgt3p0subUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfcIso4pTgt3p0subUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    outSC_Tree->Branch("pfcIso5pTgt3p0subUE"	,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
    treeStore[SCTreeName]=outSC_Tree;
}

void TreeMaker::fillSCVariablesToOutTree(Int_t scIDX,TString SCTreeName)
{
  Int_t idx(0) , offset(candidateMapDouble["SCTreeStorage"]);
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoEt->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoEta->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoPhi->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoSCE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoSCEt->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoSCRawE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoSCEta->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoSCPhi->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoSCEtaWidth->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoSCPhiWidth->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoSCBrem->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoSCnBC->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoESEn->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoPSCE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoPSCRawE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoPSCEta->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoPSCPhi->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoPSCEtaWidth->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoPSCPhiWidth->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoPSCBrem->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoPSCnHits->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoPSCflags->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoPSCinClean->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoPSCinUnClean->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoPSCnBC->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoPESEn->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoIsPFPhoton->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoIsStandardPhoton->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoHasPixelSeed->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoHasConversionTracks->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoHadTowerOverEm->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoHoverE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoHoverEValid->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoSigmaIEtaIEta->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoR9->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoE1x5->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoE2x5->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoE3x3->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoE5x5->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoMaxEnergyXtal->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoSigmaEtaEta->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoSigmaIEtaIEta_2012->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoR9_2012->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoE1x5_2012->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoE2x5_2012->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoE3x3_2012->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoE5x5_2012->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoMaxEnergyXtal_2012->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoSigmaEtaEta_2012->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoHadTowerOverEm1->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoHadTowerOverEm2->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoHoverE1->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoHoverE2->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoSigmaIEtaIPhi->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoSigmaIPhiIPhi->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoR1x5->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoR2x5->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoE2nd->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoETop->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoEBottom->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoELeft->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoERight->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoE1x3->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoE2x2->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoE2x5Max->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoE2x5Top->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoE2x5Bottom->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoE2x5Left->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoE2x5Right->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoSigmaIEtaIPhi_2012->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoSigmaIPhiIPhi_2012->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoR1x5_2012->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoR2x5_2012->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoE2nd_2012->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoETop_2012->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoEBottom_2012->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoELeft_2012->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoERight_2012->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoE1x3_2012->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoE2x2_2012->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoE2x5Max_2012->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoE2x5Top_2012->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoE2x5Bottom_2012->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoE2x5Left_2012->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoE2x5Right_2012->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoBC1E->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoBC1Ecorr->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoBC1Eta->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoBC1Phi->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoBC1size->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoBC1flags->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoBC1inClean->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoBC1inUnClean->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.phoBC1rawID->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pho_seedTime->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfcIso1->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfcIso2->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfcIso3->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfcIso4->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfcIso5->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfpIso1->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfpIso2->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfpIso3->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfpIso4->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfpIso5->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfnIso1->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfnIso2->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfnIso3->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfnIso4->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfnIso5->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfpIso1subSC->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfpIso2subSC->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfpIso3subSC->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfpIso4subSC->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfpIso5subSC->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfcIso1subUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfcIso2subUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfcIso3subUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfcIso4subUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfcIso5subUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfpIso1subUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfpIso2subUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfpIso3subUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfpIso4subUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfpIso5subUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfnIso1subUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfnIso2subUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfnIso3subUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfnIso4subUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfnIso5subUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfpIso1subSCsubUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfpIso2subSCsubUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfpIso3subSCsubUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfpIso4subSCsubUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfpIso5subSCsubUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfcIso1pTgt1p0subUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfcIso2pTgt1p0subUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfcIso3pTgt1p0subUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfcIso4pTgt1p0subUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfcIso5pTgt1p0subUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfcIso1pTgt2p0subUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfcIso2pTgt2p0subUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfcIso3pTgt2p0subUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfcIso4pTgt2p0subUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfcIso5pTgt2p0subUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfcIso1pTgt3p0subUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfcIso2pTgt3p0subUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfcIso3pTgt3p0subUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfcIso4pTgt3p0subUE->at(scIDX); idx+=1 ;
  storageArrayDouble[idx + offset ] =ntupleRawTree.pfcIso5pTgt3p0subUE->at(scIDX); idx+=1 ;
  //AssignPFVariables(scIDX);
  treeStore[SCTreeName]->Fill();
}

/*
  void TreeMaker::AssignPFVariables( Int_t scIDX  )
  {
  storageFloat["scPF_nGammaInDr0p1"]=0;
  storageFloat["scPF_nEleGammaInDr0p1"]=0;
  storageFloat["scPF_nCHadronInDr0p1"]=0;
  storageFloat["scPF_nNHadronInDr0p1"]=0;
  storageFloat["scPF_nOtherInDr0p1"]=0;
  storageFloat["scPF_sumPtGammaInDr0p1"]=0;
  storageFloat["scPF_sumPtElectronInDr0p1"]=0;
  storageFloat["scPF_sumPtCHadronInDr0p1"]=0;
  storageFloat["scPF_sumPtNHadronInDr0p1"]=0;
  storageFloat["scPF_sumPtOtherInDr0p1"]=0;

  storageFloat["scPF_nGammaInDr0p3"]=0;
  storageFloat["scPF_nEleGammaInDr0p3"]=0;
    storageFloat["scPF_nCHadronInDr0p3"]=0;
    storageFloat["scPF_nNHadronInDr0p3"]=0;
    storageFloat["scPF_nOtherInDr0p3"]=0;
    storageFloat["scPF_sumPtGammaInDr0p3"]=0;
    storageFloat["scPF_sumPtElectronInDr0p3"]=0;
    storageFloat["scPF_sumPtCHadronInDr0p3"]=0;
    storageFloat["scPF_sumPtNHadronInDr0p3"]=0;
    storageFloat["scPF_sumPtOtherInDr0p3"]=0;

    storageFloat["scPF_nGammaInDr0p5"]=0;
    storageFloat["scPF_nEleGammaInDr0p5"]=0;
    storageFloat["scPF_nCHadronInDr0p5"]=0;
    storageFloat["scPF_nNHadronInDr0p5"]=0;
    storageFloat["scPF_nOtherInDr0p5"]=0;
    storageFloat["scPF_sumPtGammaInDr0p5"]=0;
    storageFloat["scPF_sumPtElectronInDr0p5"]=0;
    storageFloat["scPF_sumPtCHadronInDr0p5"]=0;
    storageFloat["scPF_sumPtNHadronInDr0p5"]=0;
    storageFloat["scPF_sumPtOtherInDr0p5"]=0;


    
    Float_t eta0(ntupleRawTree.phoEta->at(scIDX));
    Float_t phi0(ntupleRawTree.phoPhi->at(scIDX));
    Float_t dr;
    for(Int_t i =0;i< ntupleRawTree.nPFCandidates;i++)
    {   
               dr=getDR(eta0,phi0,ntupleRawTree.pf_eta[i],ntupleRawTree.pf_phi[i]);
               if(dr<0.1)
               {
                    if(ntupleRawTree.pf_id[i] == 4.0)         storageFloat["scPF_nGammaInDr0p1"]+=1;
                    else if(ntupleRawTree.pf_id[i] == 2.0)    storageFloat["scPF_nEleGammaInDr0p1"]+=1;
                    else if(ntupleRawTree.pf_id[i] == 1.0)    storageFloat["scPF_nCHadronInDr0p1"]+=1;
                    else if(ntupleRawTree.pf_id[i] == 5.0)    storageFloat["scPF_nNHadronInDr0p1"]+=1;
                    else                                 storageFloat["scPF_nOtherInDr0p1"]+=1;

                    if(ntupleRawTree.pf_id[i] == 4.0)         storageFloat["scPF_sumPtGammaInDr0p1"]   +=ntupleRawTree.pf_pt[i];
                    else if(ntupleRawTree.pf_id[i] == 2.0)    storageFloat["scPF_sumPtElectronInDr0p1"]+=ntupleRawTree.pf_pt[i];
                    else if(ntupleRawTree.pf_id[i] == 1.0)    storageFloat["scPF_sumPtCHadronInDr0p1"] +=ntupleRawTree.pf_pt[i];
                    else if(ntupleRawTree.pf_id[i] == 5.0)    storageFloat["scPF_sumPtNHadronInDr0p1"] +=ntupleRawTree.pf_pt[i];
                    else                                      storageFloat["scPF_sumPtOtherInDr0p1"]   +=ntupleRawTree.pf_pt[i];
               }

               if(dr<0.3)
               {
                    if(ntupleRawTree.pf_id[i] == 4.0)         storageFloat["scPF_nGammaInDr0p3"]+=1;
                    else if(ntupleRawTree.pf_id[i] == 2.0)    storageFloat["scPF_nEleGammaInDr0p3"]+=1;
                    else if(ntupleRawTree.pf_id[i] == 1.0)    storageFloat["scPF_nCHadronInDr0p3"]+=1;
                    else if(ntupleRawTree.pf_id[i] == 5.0)    storageFloat["scPF_nNHadronInDr0p3"]+=1;
                    else                                 storageFloat["scPF_nOtherInDr0p3"]+=1;

                    if(ntupleRawTree.pf_id[i] == 4.0)         storageFloat["scPF_sumPtGammaInDr0p3"]   +=ntupleRawTree.pf_pt[i];
                    else if(ntupleRawTree.pf_id[i] == 2.0)    storageFloat["scPF_sumPtElectronInDr0p3"]+=ntupleRawTree.pf_pt[i];
                    else if(ntupleRawTree.pf_id[i] == 1.0)    storageFloat["scPF_sumPtCHadronInDr0p3"] +=ntupleRawTree.pf_pt[i];
                    else if(ntupleRawTree.pf_id[i] == 5.0)    storageFloat["scPF_sumPtNHadronInDr0p3"] +=ntupleRawTree.pf_pt[i];
                    else                                      storageFloat["scPF_sumPtOtherInDr0p3"]   +=ntupleRawTree.pf_pt[i];


               }
               
               if(dr<0.5)
               {
                    if(ntupleRawTree.pf_id[i] == 4.0)         storageFloat["scPF_nGammaInDr0p5"]+=1;
                    else if(ntupleRawTree.pf_id[i] == 2.0)    storageFloat["scPF_nEleGammaInDr0p5"]+=1;
                    else if(ntupleRawTree.pf_id[i] == 1.0)    storageFloat["scPF_nCHadronInDr0p5"]+=1;
                    else if(ntupleRawTree.pf_id[i] == 5.0)    storageFloat["scPF_nNHadronInDr0p5"]+=1;
                    else                                 storageFloat["scPF_nOtherInDr0p5"]+=1;

                    if(ntupleRawTree.pf_id[i] == 4.0)         storageFloat["scPF_sumPtGammaInDr0p5"]   +=ntupleRawTree.pf_pt[i];
                    else if(ntupleRawTree.pf_id[i] == 2.0)    storageFloat["scPF_sumPtElectronInDr0p5"]+=ntupleRawTree.pf_pt[i];
                    else if(ntupleRawTree.pf_id[i] == 1.0)    storageFloat["scPF_sumPtCHadronInDr0p5"] +=ntupleRawTree.pf_pt[i];
                    else if(ntupleRawTree.pf_id[i] == 5.0)    storageFloat["scPF_sumPtNHadronInDr0p5"] +=ntupleRawTree.pf_pt[i];
                    else                                      storageFloat["scPF_sumPtOtherInDr0p5"]   +=ntupleRawTree.pf_pt[i];

               }
    }    

}*/

