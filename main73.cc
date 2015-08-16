// main12.cc is a part of the PYTHIA event generator.
// Copyright (C) 2010 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. 
// It illustrates how Les Houches Event File input can be used in Pythia8.
// It uses the ttsample.lhe input file, the latter only with 100 events.
// Other samples could be generated as illustrated by main53.f.


#include <fstream>
#include <cstdlib>

#include "Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

using namespace Pythia8; 
int main() {

  std::ofstream output("data.in");
  std::ofstream output2("data2.in");
  // Number of events to print.
  int nPrint = 1;             
  double pi = acos(-1.0);
  int asym[3];
  for(int i=0;i<3;i++)
    asym[i]=0;

  // Generator           
  Pythia pythia;                            

  // Stick with default values, so do not bother with a separate file
  // for changes. However, do one change, to show readString in action.
  pythia.readString("PartonLevel:ISR = on"); 
  pythia.readString("PartonLevel:FSR = on"); 
  pythia.readString("PartonLevel:MI = off"); 
  pythia.readString("HadronLevel:Hadronize = on"); 

  // Initialize Les Houches Event File run. List initialization information.
  pythia.init("higgsll.lhef");      

  // Allow for possibility of a few faulty events.
  int nAbort = 10;
  int iAbort = 0;

  int nEvent = 0;
  int iEvent;
  // Begin event loop; generate until none left in input file.     
  for (iEvent = 0; ; ++iEvent) {

    int njets = 0;
    int nlep = 0;
    int leptons[10];

    

    // Generate events, and check whether generation failed.
    if (!pythia.next()) {

      // If failure because reached end of file then exit event loop.
      if (pythia.info.atEndOfFile()) {
	cout<<"End of file reached"<<endl;
	break; 
      }

      // First few failures write off as "acceptable" errors, then quit.
      if (++iAbort < nAbort) continue;
      break;
    }
  
    // List first few events: Les Houches, hard process and complete.
    if (iEvent < nPrint) {     
      pythia.LHAeventList();               
      pythia.info.list();          
      pythia.process.list();          
      pythia.event.list();           
    }                           
    
    // Fastjet analysis - select algorithm and parameters
    double Rparam = 0.4;
    fastjet::Strategy               strategy = fastjet::Best;
    fastjet::RecombinationScheme    recombScheme = fastjet::E_scheme;
    fastjet::JetDefinition         *jetDef = NULL;
    jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, Rparam,
					recombScheme, strategy);
    // Fastjet input
    std::vector <fastjet::PseudoJet> fjInputs;
    
    // Keep track of missing ET
    Vec4 missingETvec;

    // Loop over event record to decide what to pass to FastJet
    for (int i = 0; i < pythia.event.size(); i++) {
      // Final state only
      if (!pythia.event[i].isFinal()) continue;

      // No neutrinos
      if (pythia.event[i].idAbs() == 12 || pythia.event[i].idAbs() == 14 ||
          pythia.event[i].idAbs() == 16)     continue;

      //No isolated leptons
      if ((pythia.event[i].idAbs() == 11 || pythia.event[i].idAbs() == 13)
	  && pythia.event[i].pT() > 10.0 && abs(pythia.event[i].eta()) < 2.5){
	double ptsum = 0.0;
	for(int k = 0; k < pythia.event.size(); k++){
	  double sep,deta,dphi;
	  if(k != i && pythia.event[k].isFinal() && pythia.event[k].pT() > 0.5){
	    deta = pythia.event[k].eta() - pythia.event[i].eta();
	    dphi = pythia.event[k].phi() - pythia.event[i].phi();
	    if(dphi > pi) dphi = dphi-2*pi;
	    if(dphi < -pi) dphi = 2*pi+dphi;
	    sep = sqrt(pow2(deta) + pow2(dphi));
	    if(sep < 0.2) ptsum += pythia.event[k].pT();
	  }
	}
	if(ptsum < 10.0) {
	  leptons[nlep]=i;
	  nlep++;
	  // Still contributes to missing ET
	  missingETvec += pythia.event[i].p();
	  continue;
	}
      }

      // Only |eta| < 2.5
      if (fabs(pythia.event[i].eta()) > 2.5) continue;

      // Missing ET
      missingETvec += pythia.event[i].p();

      // Store as input to Fastjet
      fjInputs.push_back( fastjet::PseudoJet (pythia.event[i].px(),
					      pythia.event[i].py(), pythia.event[i].pz(),
					      pythia.event[i].e() ) );
    }

    if (fjInputs.size() == 0) {
      cout << "Error: event "<<iEvent<<" with no final state particles" << endl;
      continue;
    }



    // Sort the leptons
    if(nlep > 1){
      if( pythia.event[leptons[0]].pT()<pythia.event[leptons[1]].pT()){
	int itmp   = leptons[1];
	leptons[1] = leptons[0];
	leptons[0] = itmp;
      }
    }

    // Start Cuts
    // Cut 2: Need two isolated leptons with correct pT
    if(nlep < 2 || pythia.event[leptons[0]].pT() < 40.0 
       || pythia.event[leptons[1]].pT() < 30.0) continue;

    // Cut3: Missing ET cut
    double missingET = missingETvec.pT();
    if(missingET < 30.0) continue;

    nEvent++;

    // Run Fastjet algorithm

    vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
    fastjet::ClusterSequence clustSeq(fjInputs, *jetDef);

    // Extract inclusive jets sorted by pT (note minimum pT of 20.0 GeV)
    inclusiveJets = clustSeq.inclusive_jets(20.0);
    sortedJets    = sorted_by_pt(inclusiveJets);  
    
    njets = sortedJets.size();
    if(njets > 0){
      for (int i = 0; i < sortedJets.size()-1; i++) {
	if(abs(sortedJets[i].eta()) > 2.5) njets=i;
	if(sortedJets[i].perp() < 20.0) njets=i;
	break;
      }
    }
  
    //Reconstruct the W
    double mj1j2[4], minv,minvmin,j1min,j2min;
    minvmin = 0.0;
    if(njets>1){
      for (int i = 0; i < sortedJets.size()-1; i++) {
	for (int j = i+1; j < sortedJets.size(); j++) {
	  mj1j2[0]=sortedJets[i].e()+sortedJets[j].e();
	  mj1j2[1]=sortedJets[i].px()+sortedJets[j].px();
	  mj1j2[2]=sortedJets[i].py()+sortedJets[j].py();
	  mj1j2[3]=sortedJets[i].pz()+sortedJets[j].pz();
	  
	  minv = sqrt(pow2(mj1j2[0])-pow2(mj1j2[1])-pow2(mj1j2[2])-pow2(mj1j2[3]));
	  if(abs(minv-80.4)<abs(minvmin-80.4)) {
	    j1min = i;
	    j2min = j;
	    minvmin = minv;
	  }
	}
      }
    }

    // Start output

    output<<iEvent<<"  "<<njets<<"  ";
    output<<(pythia.event[leptons[0]].pT() - pythia.event[leptons[1]].pT())<<"  ";

    double sep,deta,dphi,sepmin = 10.0;
    int lep1=0 ,lep2=1;

    // if(njets>1){
    //   // Find the one closest to the W
    //   double phiW = atan(abs(mj1j2[2])/abs(mj1j2[1]));
    //   if(mj1j2[1]<0 && mj1j2[2]>0) phiW = pi - phiW;
    //   if(mj1j2[1]<0 && mj1j2[2]<0) phiW = pi + phiW;
    //   if(mj1j2[1]>0 && mj1j2[2]<0) phiW = 2*pi - phiW;
    //   deta = pythia.event[leptons[0]].eta()-etaW;
    //   dphi = pythia.event[leptons[0]].phi()-phiW;
    //   if(dphi > pi) dphi = dphi-2*pi;
    //   if(dphi < -pi) dphi = 2*pi+dphi;
    //   double sep1 = sqrt(pow2(deta)+pow2(dphi));

    //   deta = pythia.event[leptons[1]].eta()-etaW;
    //   dphi = pythia.event[leptons[1]].phi()-phiW;
    //   if(dphi > pi) dphi = dphi-2*pi;
    //   if(dphi < -pi) dphi = 2*pi+dphi;

    //   if(abs(dphi)>pi) dphi=2*pi-abs(dphi);
    //   double sep2 = sqrt(pow2(deta)+pow2(dphi));

    //   if(sep1<sep2){
    // 	lep1=1;
    // 	lep2=0;
    //   }
    // }
    if(njets > 0 ){
      for (int i = 0 ; i < njets ; i++){
	deta = pythia.event[leptons[0]].eta()-sortedJets[i].eta();
	dphi = pythia.event[leptons[0]].phi()-sortedJets[i].phi();
	if(dphi > pi) dphi = dphi-2*pi;
	if(dphi < -pi) dphi = 2*pi+dphi;
	sep = sqrt(pow2(deta)+pow2(dphi));
	if(sep<sepmin) {
    	lep1=1;
    	lep2=0;
    	sepmin=sep;
	}
      }
      for (int i = 0 ; i < njets ; i++){
	deta = pythia.event[leptons[1]].eta()-sortedJets[i].eta();
	dphi = pythia.event[leptons[1]].phi()-sortedJets[i].phi();
	if(dphi > pi) dphi = dphi-2*pi;
	if(dphi < -pi) dphi = 2*pi+dphi;
	sep = sqrt(pow2(deta)+pow2(dphi));
	if(sep<sepmin) {
	  lep1=0;
	  lep2=1;
	  sepmin=sep;
	}
      }
    }

    dphi = pythia.event[leptons[lep1]].phi()-pythia.event[leptons[lep2]].phi();
    if(dphi > pi) dphi = dphi-2*pi;
    if(dphi < -pi) dphi = 2*pi+dphi;

    if(dphi<0) asym[0]--;
    else asym[0]++;
    
    if(abs(dphi)<pi/2) asym[1]++;
    else asym[1]--;

    output<<dphi<<endl;

    deta = -99.0;
    double etaW;
    if(njets>1) {
      mj1j2[0]=sortedJets[j1min].e()+sortedJets[j2min].e();
      mj1j2[1]=sortedJets[j1min].px()+sortedJets[j2min].px();
      mj1j2[2]=sortedJets[j1min].py()+sortedJets[j2min].py();
      mj1j2[3]=sortedJets[j1min].pz()+sortedJets[j2min].pz();
      etaW = 0.5 * log((mj1j2[0]+mj1j2[3])/(mj1j2[0]-mj1j2[3]));
      deta = abs(pythia.event[leptons[0]].eta()-etaW)-abs(pythia.event[leptons[1]].eta()-etaW);

      if(deta<0) asym[2]--;
      else asym[2]++;

      output2<<deta<<"  "<<minvmin<<endl;
    }


  }

  cout<<"Efficiency of cuts: "<<(nEvent*1.0/iEvent)<<endl;
  for(int i=0;i<3;i++)
    cout<<"Asymmetry["<<i<<"]="<<asym[i]*1.0/nEvent<<endl;
  //pythia.statistics();
  // Done.                           
  return 0;
  

}
