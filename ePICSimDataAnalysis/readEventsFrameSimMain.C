#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

//R__LOAD_LIBRARY(libfmt.so)
//R__LOAD_LIBRARY(fmt)
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
/*
#include "podio/ROOTReader.h"
#include "podio/EventStore.h"
#include "podio/CollectionIDTable.h"
#include "podio/ObjectID.h"
*/
//R__LOAD_LIBRARY(podio)
//R__LOAD_LIBRARY(podioRootIO)
#include "podio/ROOTFrameReader.h"
#include "podio/Frame.h"
/*
R__LOAD_LIBRARY(edm4hep)
R__LOAD_LIBRARY(edm4eic)
R__LOAD_LIBRARY(edm4hepDict)
R__LOAD_LIBRARY(edm4eicDict)
*/
#include "edm4hep/utils/kinematics.h"

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

#include "edm4eic/ReconstructedParticleCollection.h"
#include "edm4eic/MCRecoParticleAssociationCollection.h"

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

int readEventsFrameSim(TString list = "data/filesSim_neutrons.list", TString ofname = "output/output_test.root", long nevents = -1);
int MakeEvent(podio::Frame *frame, unsigned ev);


bool printEvNum = true;
bool debug = false;


int main(int argc, char **argv)
{
	readEventsFrameSim();

	return 1;
}

int readEventsFrameSim(TString list, TString ofname, long nevents)
{

	//TChain *chain = new TChain("events"); // "events" - edm4hep

	//if(!openFileList(chain, list)) return 0;

	//auto store = new podio::EventStore();
	//podio::ROOTReader *reader = new podio::ROOTReader();

	auto readerAll = new podio::ROOTFrameReader();
    cout<<"ROOTReader created"<<endl;

	std::vector<std::string> filenames = openList(list);
	if(filenames.size() != 0) readerAll->openFiles(filenames);
	else {
		cout<<"Can't open file list! Exiting."<<endl;
		return 0;
	}

	cout<<"Number of all events = "<<readerAll->getEntries("events")<<endl;
	delete readerAll;

	//reader.openFile("/gpfs/mnt/gpfs02/eic/palsp/simdir/100k/neutron/epic/170_45/neutron_gun_5GeV_170_45_sci_tiles_version2_new_geo_341062_0.edm4hep.root");
	//reader.openFile("/gpfs/mnt/gpfs02/eic/palsp/simdir/100k/neutron/epic/170_45/neutron_gun_5GeV_170_45_sci_tiles_version2_new_geo_341062_0.edm4hep.root");
	//reader.openFile("/gpfs/mnt/gpfs02/eic/palsp/simdir/100k/neutron/epic/170_45/neutron_gun_5GeV_170_45_sci_tiles_version2_new_geo_341062_1.edm4hep.root");
	//reader.openFile("/gpfs/mnt/gpfs02/eic/palsp/simdir/100k/neutron/epic/170_45/neutron_gun_5GeV_170_45_sci_tiles_version2_new_geo_341062_2.edm4hep.root");


	TFile *output = new TFile(ofname, "recreate");
	output->cd();

	CreateHistogamsSim();

	unsigned ev_all = 0;

	for (int rd = 0; rd < filenames.size(); ++rd) {

	auto reader = new podio::ROOTFrameReader();
	reader->openFile(filenames.at(rd));

	//TTree *tree = (TTree*)chain;
	//tree->Print("toponly");

	//auto store = new podio::EventStore();
	//store->setReader(reader);

	//reader->getCollectionIDTable()->print();
	//auto cat_vec = reader->getAvailableCategories();

	cout<<"Running podio version = "<<podio::version::build_version<<endl;
	cout<<"File created with podio version = "<<reader->currentFileVersion()<<endl;


    fmt::print("\n{:=^50}\n", " COLLECTIONS ");
    for(auto cat_vec : reader->getAvailableCategories())
      fmt::print("  {}\n", cat_vec);

	//cout<<"Data model definition = "<<reader.getDatamodelDefinition("events")<<endl;

    fmt::print("\n{:=^50}\n", " DATA MODELS  ");
    for(auto mod_vec : reader->getAvailableDatamodels())
      fmt::print("  {}\n", mod_vec);

	unsigned nEvents = reader->getEntries("events");
	//cout<<"Number of events = "<<nEvents<<endl;
	cout<<"Number of events in a file = "<<nEvents<<endl;

	if(nevents>0) nEvents = nevents;
	

	for(unsigned ev=0; ev<nEvents; ++ev) {

		//reader.readEventMetaData(reader.getCategoryInfo("events"));

		// next event
		auto frame = new podio::Frame(reader->readNextEntry("events"));
		//auto frame = reader.readNextEntry("events");
		//auto frame = reader.readEntry("events", ev);

		//if(printEvNum) std::cout<<"reading event "<<ev<<std::endl;
		if(printEvNum) std::cout<<"reading event "<<ev_all<<std::endl;

		MakeEvent(frame, ev_all); // execute for each event

		ev_all++;
		delete frame;

		if(debug) std::cout<<"End of event"<<std::endl;

	} // event loop

	delete reader;
	} // readers loop

	output->Write();

	return 1;

}


int MakeEvent(podio::Frame *frame, unsigned ev)
{

    if(ev==0) {
      fmt::print("\n{:=^50}\n", " COLLECTIONS ");
      for(auto coll : frame->getAvailableCollections())
        fmt::print("  {}\n", coll);
    }

	// Get collections

	//edm4hep::SimCalorimeterHitCollection nHCal_hits_frame;
	auto& nHCal_hits_frame  = frame->get<edm4hep::SimCalorimeterHitCollection>("HcalEndcapNHits");
	auto& MCParticles_frame  = frame->get<edm4hep::MCParticleCollection>("MCParticles");

	int nPion_p = 0;
	int nPion_n = 0;
	int nKaon_p = 0;
	int nKaon_n = 0;
	int nProton_p = 0;
	int nProton_n = 0;
	int nElectron_p = 0;
	int nElectron_n = 0;

	int nNeutron = 0;

	int nMCpart_gen = 0;
	int nMCpart_sec = 0;


	// MC particle loop

	if(!MCParticles_frame.isValid())
		cout<<"SimCalorimeterHitCollection does not exist!"<<endl;

	if(debug) cout<<"MCParticles size = "<<MCParticles_frame.size()<<endl;

	for (unsigned mc_iter = 0; mc_iter < MCParticles_frame.size(); ++mc_iter) {

		MCParticle mcpart =  MCParticles_frame.at(mc_iter);


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

		if(mcpart.getPDG() == 211) nPion_p++;
		if(mcpart.getPDG() == -211) nPion_n++;
		if(mcpart.getPDG() == 321) nKaon_p++;
		if(mcpart.getPDG() == -321) nKaon_n++;
		if(mcpart.getPDG() == 2212) nProton_p++;
		if(mcpart.getPDG() == -2212) nProton_n++;
		if(mcpart.getPDG() == -11) nElectron_p++;
		if(mcpart.getPDG() == 11) nElectron_n++;

		if(mcpart.getPDG() == 2112) nNeutron++;


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

			nMCpart_gen++;
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

			nMCpart_sec++;

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


	h_MCpart->Fill(MCParticles_frame.size());

	h_MCpart_nPion_p->Fill(nPion_p);
	h_MCpart_nPion_n->Fill(nPion_n);
	h_MCpart_nKaon_p->Fill(nKaon_p);
	h_MCpart_nKaon_n->Fill(nKaon_n);
	h_MCpart_nProton_p->Fill(nProton_p);
	h_MCpart_nProton_n->Fill(nProton_n);
	h_MCpart_nElectron_p->Fill(nElectron_p);
	h_MCpart_nElectron_n->Fill(nElectron_n);

	h_MCpart_nNeutron->Fill(nNeutron);

	h_MCpart_nGen->Fill(nMCpart_gen);
	h_MCpart_nSec->Fill(nMCpart_sec);




	if(!nHCal_hits_frame.isValid())
		cout<<"SimCalorimeterHitCollection does not exist!"<<endl;

	if(debug) cout<<"SimCalorimeterHitCollection size = "<<nHCal_hits_frame.size()<<endl;

	h_nHCal_nhits->Fill(nHCal_hits_frame.size());


		for (unsigned hit_iter = 0; hit_iter < nHCal_hits_frame.size(); ++hit_iter) {

			SimCalorimeterHit hit_nHCal =  nHCal_hits_frame[hit_iter];

			if(!hit_nHCal.isAvailable())
				cout<<"SimCalorimeterHit does not exist! index = "<<hit_nHCal<<endl;


			h_nHCal_hit_E->Fill(hit_nHCal.getEnergy());

			h_nHCal_hit_pos_x->Fill(hit_nHCal.getPosition().x);
			h_nHCal_hit_pos_y->Fill(hit_nHCal.getPosition().y);
			h_nHCal_hit_pos_z->Fill(hit_nHCal.getPosition().z);
			h_nHCal_hit_pos_xy->Fill(hit_nHCal.getPosition().x, hit_nHCal.getPosition().y);








			auto contrib = hit_nHCal.getContributions();

			//if(contrib==NULL)
			if(contrib.size()==0)
				cout<<"Contributions vector does not exist!"<<endl;

			if(debug) cout<<"contributions size = "<<contrib.size()<<endl;

			for (unsigned c = 0; c < contrib.size(); ++c) {

				if(contrib[c]==NULL)
				if(!contrib.at(c).isAvailable())
					cout<<"Contribution does not exist! index = "<<c<<endl;

				if(debug) cout<<"hit time = "<<contrib.at(c).getTime()<<endl;

			}


		} // HcalEndcapNHits loop


	return 1;
}
