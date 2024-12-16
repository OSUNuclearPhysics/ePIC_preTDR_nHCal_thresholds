#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

//R__LOAD_LIBRARY(libfmt.so)
#include "fmt/core.h"

#include "ROOT/RDataFrame.hxx"
#include <TH1D.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TMath.h>
#include <TVector3.h>

#include "TROOT.h"
#include "TRandom.h"
#include "TH3.h"


#include "DD4hep/Detector.h"
#include "DDRec/CellIDPositionConverter.h"

#include "podio/ROOTReader.h"
#include "podio/EventStore.h"
#include "podio/CollectionIDTable.h"
#include "podio/ObjectID.h"

#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/MCParticleCollectionData.h"
#include "edm4hep/MCParticle.h"
#include "edm4hep/MCParticleData.h"

#include "edm4hep/SimCalorimeterHitCollectionData.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/SimCalorimeterHitData.h"
#include "edm4hep/SimCalorimeterHit.h"

#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollectionData.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/CalorimeterHitData.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitObj.h"

//dd4hep::sim::Geant4Calorimeter::Hit

#include "edm4eic/ClusterCollection.h"
#include "edm4eic/Cluster.h"
#include "edm4eic/ClusterData.h"

#include "edm4eic/CalorimeterHit.h"
#include "edm4eic/CalorimeterHitCollectionData.h"
#include "edm4eic/CalorimeterHitCollection.h"
#include "edm4eic/CalorimeterHitData.h"
#include "edm4eic/CalorimeterHit.h"
#include "edm4eic/CalorimeterHitObj.h"


#include <edm4eic/vector_utils_legacy.h>
#include <edm4hep/Vector3f.h>

///#include "eicd/Vector3f.h"

///#include "eicd/VectorXYZ.h"
//#include "eicd/Cluster.h"
//#include "eicd/ClusterData.h"
//#include "edm4hep/Vector3f.h"
//#include "Vector3D.h"
//include "eic/Vector3D.h"
//#include <eic/vector_utils.h>
//#include "dd4pod/CalorimeterHitData.h"

#include "FileList.h"
#include "HistogramsSim.h"

#pragma link C++ class vector<edm4hep::MCParticleData>+;
#pragma link C++ class vector<eicd::ClusterData>+;
#pragma link C++ class vector<podio::ObjectID>+;
#pragma link C++ class vector<edm4hep::SimCalorimeterHitData>+;
#pragma link C++ class vector<edm4eic::CalorimeterHitData>+;


using namespace std;
using namespace ROOT;
using namespace TMath;
//using namespace eicd;
//using namespace edm4eic;
using namespace edm4hep;

int MakeEvent(podio::EventStore *store, unsigned ev);


bool printEvNum = true;
bool debug = true;

int main(int argc, char **argv)
{
	TString list = "data/filesSim_hcal_only_test_good.list";
	TString ofname = "output/output_test.root";
	long nevents = -1;


	//TChain *chain = new TChain("events"); // "events" - edm4hep

	//if(!openFileList(chain, list)) return 0;

	auto store = new podio::EventStore();
	podio::ROOTReader *reader = new podio::ROOTReader();

    cout<<"ROOTReader created"<<endl;

	std::vector<std::string> filenames = openList(list);
	if(filenames.size() != 0) reader->openFiles(filenames);
	else {
		cout<<"Can't open file list! Exiting."<<endl;
		return 0;
	}

	TFile *output = new TFile(ofname, "recreate");
	output->cd();

	CreateHistogamsSim();

	//TTree *tree = (TTree*)chain;
	//tree->Print("toponly");

	//auto store = new podio::EventStore();
	store->setReader(reader);

	reader->getCollectionIDTable()->print();

	unsigned nEvents = reader->getEntries();
	cout<<"Number of events = "<<nEvents<<endl;

	if(nevents>0) nEvents = nevents;
	

	for(unsigned ev=0; ev<nEvents; ++ev) {

		//tree->GetEntry(ev);
		//reader->goToEvent(ev);
		//reader->readEvent(); // use it! Only breaks when trying to access HcalEndcapNHits for some files
		if(printEvNum) std::cout<<"reading event "<<ev<<std::endl;

		MakeEvent(store, ev); // execute for each event

		if(debug) std::cout<<"Calling clear event"<<std::endl;

		store->clear();
		reader->endOfEvent();

	} // event loop

	output->Write();

	return 1;

}


int MakeEvent(podio::EventStore *store, unsigned ev)
{

	// Get collections

	//edm4hep::SimCalorimeterHitCollection nHCal_hits_store;
	auto& nHCal_hits_store  = store->get<edm4hep::SimCalorimeterHitCollection>("HcalEndcapNHits");
	auto& MCParticles_store  = store->get<edm4hep::MCParticleCollection>("MCParticles");


	// MC particle loop

	if(!MCParticles_store.isValid())
		cout<<"SimCalorimeterHitCollection does not exist!"<<endl;

	if(debug) cout<<"MCParticles size = "<<MCParticles_store.size()<<endl;

	for (unsigned mc_iter = 0; mc_iter < MCParticles_store.size(); ++mc_iter) {

		MCParticle mcpart =  MCParticles_store.at(mc_iter);


		TVector3 mcMom(mcpart.getMomentum().x, mcpart.getMomentum().y, mcpart.getMomentum().z);
		TVector3 mcMomEnd(mcpart.getMomentumAtEndpoint ().x, mcpart.getMomentumAtEndpoint ().y, mcpart.getMomentumAtEndpoint ().z);
		TVector3 mcStart(mcpart.getVertex().x, mcpart.getVertex().y, mcpart.getVertex().z);
		TVector3 mcEnd(mcpart.getEndpoint().x, mcpart.getEndpoint().y, mcpart.getEndpoint().z);

		if(debug)
		{
			cout<<"MCParticle px = "<<mcMom.x()<<endl;
			cout<<"MCParticle py = "<<mcMom.y()<<endl;
			cout<<"MCParticle pz = "<<mcMom.z()<<endl;
		}

		h_MCpart_mass->Fill(mcpart.getMass());
		h_MCpart_charge->Fill(mcpart.getCharge());
		h_MCpart_E->Fill(mcpart.getEnergy());
		h_MCpart_p->Fill(mcMom.Mag());
		h_MCpart_pT->Fill(mcMom.Pt());

		h_MCpart_mom_x->Fill(mcMom.x());
		h_MCpart_mom_y->Fill(mcMom.y());
		h_MCpart_mom_z->Fill(mcMom.z());

		h_MCpart_eta->Fill(mcMom.Eta());
		h_MCpart_etaphi->Fill(mcMom.Eta(), mcMom.Phi());

		h_MCpart_xy->Fill(mcStart.x(), mcStart.y());
		h_MCpart_zr->Fill(mcStart.z(), mcStart.Pt());

		h_MCpart_end_p->Fill(mcMomEnd.Mag());
		h_MCpart_end_pT->Fill(mcMomEnd.Pt());

		h_MCpart_posEnd_xy->Fill(mcEnd.x(), mcEnd.y());
		h_MCpart_posEnd_zr->Fill(mcEnd.z(), mcEnd.Pt());


		// Generated MC particles
		if(!mcpart.isCreatedInSimulation())
		{
			h_MCpart_gen_mass->Fill(mcpart.getMass());
			h_MCpart_gen_charge->Fill(mcpart.getCharge());
			h_MCpart_gen_E->Fill(mcpart.getEnergy());
			h_MCpart_gen_p->Fill(mcMom.Mag());
			h_MCpart_gen_pT->Fill(mcMom.Pt());

			h_MCpart_gen_eta->Fill(mcMom.Eta());
			h_MCpart_gen_etaphi->Fill(mcMom.Eta(), mcMom.Phi());

			h_MCpart_gen_xy->Fill(mcStart.x(), mcStart.y());
			h_MCpart_gen_zr->Fill(mcStart.z(), mcStart.Pt());

			h_MCpart_gen_end_p->Fill(mcMomEnd.Mag());
			h_MCpart_gen_end_pT->Fill(mcMomEnd.Pt());

			h_MCpart_gen_posEnd_xy->Fill(mcEnd.x(), mcEnd.y());
			h_MCpart_gen_posEnd_zr->Fill(mcEnd.z(), mcEnd.Pt());
		}


		// Secondary MC particles
		if(mcpart.isCreatedInSimulation())
		{
			h_MCpart_sec_mass->Fill(mcpart.getMass());
			h_MCpart_sec_charge->Fill(mcpart.getCharge());
			h_MCpart_sec_E->Fill(mcpart.getEnergy());
			h_MCpart_sec_p->Fill(mcMom.Mag());
			h_MCpart_sec_pT->Fill(mcMom.Pt());

			h_MCpart_sec_eta->Fill(mcMom.Eta());
			h_MCpart_sec_etaphi->Fill(mcMom.Eta(), mcMom.Phi());

			h_MCpart_sec_xy->Fill(mcStart.x(), mcStart.y());
			h_MCpart_sec_zr->Fill(mcStart.z(), mcStart.Pt());

			h_MCpart_sec_end_p->Fill(mcMomEnd.Mag());
			h_MCpart_sec_end_pT->Fill(mcMomEnd.Pt());

			h_MCpart_sec_posEnd_xy->Fill(mcEnd.x(), mcEnd.y());
			h_MCpart_sec_posEnd_zr->Fill(mcEnd.z(), mcEnd.Pt());

			bool parentIsPrimary = false;

			podio::RelationRange<edm4hep::MCParticle> MCparents = mcpart.getParents();

			for (int parentIter = 0; parentIter < MCparents.size(); ++parentIter) {

				if(!MCparents.at(parentIter).isCreatedInSimulation())	parentIsPrimary = true;
			}

			// MC particles - 1st generation daughters
			if(parentIsPrimary)
			{
				h_MCpart_1stgen_daughter_mass->Fill(mcpart.getMass());
				h_MCpart_1stgen_daughter_charge->Fill(mcpart.getCharge());
				h_MCpart_1stgen_daughter_E->Fill(mcpart.getEnergy());
				h_MCpart_1stgen_daughter_p->Fill(mcMom.Mag());
				h_MCpart_1stgen_daughter_pT->Fill(mcMom.Pt());

				h_MCpart_1stgen_daughter_eta->Fill(mcMom.Eta());
				h_MCpart_1stgen_daughter_etaphi->Fill(mcMom.Eta(), mcMom.Phi());

				h_MCpart_1stgen_daughter_xy->Fill(mcStart.x(), mcStart.y());
				h_MCpart_1stgen_daughter_zr->Fill(mcStart.z(), mcStart.Pt());

				h_MCpart_1stgen_daughter_end_p->Fill(mcMomEnd.Mag());
				h_MCpart_1stgen_daughter_end_pT->Fill(mcMomEnd.Pt());

				h_MCpart_1stgen_daughter_posEnd_xy->Fill(mcEnd.x(), mcEnd.y());
				h_MCpart_1stgen_daughter_posEnd_zr->Fill(mcEnd.z(), mcEnd.Pt());
			}

		}



	} // MCParticles loop

	if(!nHCal_hits_store.isValid())
		cout<<"SimCalorimeterHitCollection does not exist!"<<endl;

	if(debug) cout<<"SimCalorimeterHitCollection size = "<<nHCal_hits_store.size()<<endl;

		for (unsigned hit_iter = 0; hit_iter < nHCal_hits_store.size(); ++hit_iter) {

			SimCalorimeterHit hit_s =  nHCal_hits_store[hit_iter];

			if(!hit_s.isAvailable())
				cout<<"SimCalorimeterHit does not exist! index = "<<hit_s<<endl;


			auto contrib = hit_s.getContributions();

			//if(contrib==NULL)
			//	cout<<"Contributions vector does not exist!"<<endl;

			//if(debug) cout<<"contributions size = "<<contrib.size()<<endl;

			for (unsigned c = 0; c < contrib.size(); ++c) {

				//if(contrib[c]==NULL)
				if(!contrib.at(c).isAvailable())
					cout<<"Contribution does not exist! index = "<<c<<endl;

				//if(debug) cout<<"hit time = "<<contrib.at(c).getTime()<<endl;

			}


		} // HcalEndcapNHits loop


	return 1;
}
