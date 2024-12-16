#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

R__LOAD_LIBRARY(libfmt.so)
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


int readEvents(TString list = "data/filesSim_hcal_only_test.list", TString ofname = "output/output_test.root", long nevents = -1)
{

	bool printEvNum = true;
	bool debug = true;

	TChain *chain = new TChain("events"); // "events" - edm4hep

	if(!openFileList(chain, list)) return 0;

	auto store = new podio::EventStore();
	podio::ROOTReader *reader = new podio::ROOTReader();

    cout<<"ROOTReader created"<<endl;

	std::vector<std::string> filenames = openList(list);
	if(filenames.size() != 0) reader->openFiles(filenames);
	else {
		cout<<"Can't open file list! Exiting."<<endl;
		return 0;
	}

	//TFile *input = new TFile(filenames[0].data(), "open");
	//TTree *tree = (TTree*)input->Get("events");
	//tree->Print();
	TTree *tree = (TTree*)chain;
	tree->Print("toponly");

	//auto store = new podio::EventStore();
	store->setReader(reader);

	reader->getCollectionIDTable()->print();

	unsigned nEvents = reader->getEntries();
	cout<<"Number of events = "<<nEvents<<endl;

	if(nevents>0) nEvents = nevents;
	

	vector<edm4hep::SimCalorimeterHitData> *nHCal_hitscoll = 0;

    //tree->SetBranchAddress("HcalEndcapNHits", &nHCal_hitscoll);


	for(unsigned ev=0; ev<nEvents; ++ev) {

		tree->GetEntry(ev);
		//reader->goToEvent(ev);
		reader->readEvent(); // use it! Only breaks when trying to access HcalEndcapNHits for some files
		if(printEvNum) std::cout<<"reading event "<<ev<<std::endl;

		//store->endOfEvent();

		if(nHCal_hitscoll!=NULL)
		{
			if(debug) cout<<"SimCalorimeterHit from tree size = "<<nHCal_hitscoll->size()<<endl;

			for (unsigned bhchit = 0; bhchit < nHCal_hitscoll->size(); ++bhchit) {

				SimCalorimeterHitData nHCal_hits = nHCal_hitscoll->at(bhchit);

				//if(debug) cout<<"nHCal_hitscoll energy = "<<nHCal_hits.energy<<endl;

			}

		}
		edm4hep::SimCalorimeterHitCollection nHCal_hits_store;
		//auto& nHCal_hits_store  = store->get<edm4hep::SimCalorimeterHitCollection>("HcalEndcapNHits");
		auto& MCParticles_store  = store->get<edm4hep::MCParticleCollection>("MCParticles");

		for (unsigned mc_iter = 0; mc_iter < MCParticles_store.size(); ++mc_iter) {

			MCParticle mcpart =  MCParticles_store.at(mc_iter);

			if(debug) cout<<"MCParticle px = "<<mcpart.getMomentum().x<<endl;
			if(debug) cout<<"MCParticle py = "<<mcpart.getMomentum().y<<endl;
			if(debug) cout<<"MCParticle pz = "<<mcpart.getMomentum().z<<endl;
		}

		if(!nHCal_hits_store.isValid())
			cout<<"SimCalorimeterHitCollection does not exist!"<<endl;

		if(debug) cout<<"SimCalorimeterHit size = "<<nHCal_hits_store.size()<<endl;

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

			if(debug) std::cout<<"Calling clear event"<<std::endl;

		store->clear();
		reader->endOfEvent();

	} // event loop



	return 1;

}
