/*
 * EICutil.h
 *
 *  Created on: 11 maj 2023
 *      Author: Khaless
 */

#ifndef EICUTIL_H_
#define EICUTIL_H_

#include <vector>
#include <map>

#include "edm4hep/MCParticle.h"
#include "edm4hep/MCParticleData.h"

#include "edm4eic/CalorimeterHitData.h"
#include "edm4eic/CalorimeterHit.h"
#include "edm4hep/RawCalorimeterHitData.h"
#include "edm4hep/RawCalorimeterHit.h"
#include "edm4hep/SimCalorimeterHitData.h"
#include "edm4hep/SimCalorimeterHit.h"

using namespace std;

edm4hep::MCParticle GetMCParticle(edm4hep::MCParticleData mcpart_data);
int GetMCParentData(edm4hep::MCParticleData mcpart_data, vector<edm4hep::MCParticleData> *MCParticles_data, vector<podio::ObjectID> *MCparents_data, vector<edm4hep::MCParticleData> *parents);
int GetMCDaughtersData(edm4hep::MCParticleData mcpart_data, vector<edm4hep::MCParticleData> *MCParticles_data, vector<podio::ObjectID> *MCdaughters_data, vector<edm4hep::MCParticleData> *daughters);

edm4eic::CalorimeterHit GetCaloHit(edm4eic::CalorimeterHitData hitData);
edm4hep::RawCalorimeterHit GetCaloRawHit(edm4hep::RawCalorimeterHitData hitData);
edm4hep::SimCalorimeterHit GetCaloSimHit(edm4hep::SimCalorimeterHit hitData);

int GetCaloHitContributionsData(edm4hep::SimCalorimeterHitData hit_data, vector<edm4hep::CaloHitContributionData> *contribs_data, vector<podio::ObjectID> *relations_data, vector<edm4hep::CaloHitContributionData> *contributions);
edm4hep::CaloHitContribution GetCaloHitContribution(edm4hep::CaloHitContributionData calo_data);

int GetMCParticleDataFromCaloHitContributions(vector<edm4hep::MCParticleData> *MCParticles_data, vector<podio::ObjectID> *MCparents_data, vector<edm4hep::MCParticleData> *parents);
int GetMCParticleIdFromCaloHitContributions(vector<podio::ObjectID> *MCparents_data, vector<unsigned> *parents);

int CreateMCParticleToCaloHitContributionMap(vector<edm4hep::MCParticleData> *MCParticles_data, vector<podio::ObjectID> *MCparents_data, map<int, vector<int>> *map);



edm4hep::MCParticle GetMCParticle(edm4hep::MCParticleData mcpart_data)
{
	edm4hep::MCParticle mcpart(mcpart_data.PDG, mcpart_data.generatorStatus, mcpart_data.simulatorStatus, mcpart_data.charge, mcpart_data.time, mcpart_data.mass, mcpart_data.vertex, mcpart_data.endpoint, mcpart_data.momentum, mcpart_data.momentumAtEndpoint, mcpart_data.spin, mcpart_data.colorFlow);

	return mcpart;
}

int GetMCParentsData(edm4hep::MCParticleData mcpart_data, vector<edm4hep::MCParticleData> *MCParticles_data, vector<podio::ObjectID> *MCparents_data, vector<edm4hep::MCParticleData> *parents)
{
	for (int i = mcpart_data.parents_begin; i < mcpart_data.parents_end; ++i) {

    	unsigned parentID = MCparents_data->at(i).index;

		edm4hep::MCParticleData mcpart_data_parent = MCParticles_data->at(parentID);

		parents->push_back(mcpart_data_parent);
	}

	return 1.0;
}
int GetMCDaughtersData(edm4hep::MCParticleData mcpart_data, vector<edm4hep::MCParticleData> *MCParticles_data, vector<podio::ObjectID> *MCdaughters_data, vector<edm4hep::MCParticleData> *daughters)
{
	for (int i = mcpart_data.daughters_begin; i < mcpart_data.daughters_end; ++i) {

    	unsigned parentID = MCdaughters_data->at(i).index;

		edm4hep::MCParticleData mcpart_data_daughter = MCParticles_data->at(parentID);

		daughters->push_back(mcpart_data_daughter);
	}

	return 1.0;
}

edm4eic::CalorimeterHit GetCaloHit(edm4eic::CalorimeterHitData hitData)
{
	edm4eic::CalorimeterHit hit(hitData.cellID, hitData.energy, hitData.energyError, hitData.time, hitData.timeError, hitData.position, hitData.dimension, hitData.sector, hitData.layer, hitData.local);

	return hit;
}

edm4hep::RawCalorimeterHit GetCaloRecHit(edm4hep::RawCalorimeterHitData hitData)
{
	edm4hep::RawCalorimeterHit hit(hitData.cellID, hitData.amplitude, hitData.timeStamp);

	return hit;
}


edm4hep::SimCalorimeterHit GetCaloSimHit(edm4hep::SimCalorimeterHitData hitData)
{
	edm4hep::SimCalorimeterHit hit(hitData.cellID, hitData.energy, hitData.position);

	return hit;
}


int GetCaloHitContributionsData(edm4hep::SimCalorimeterHitData hit_data, vector<edm4hep::CaloHitContributionData> *contribs_data, vector<podio::ObjectID> *relations_data, vector<edm4hep::CaloHitContributionData> *contributions)
{
	for (int i = hit_data.contributions_begin; i < hit_data.contributions_end; ++i) {

    	unsigned parentID = relations_data->at(i).index;

		edm4hep::CaloHitContributionData contrib_data = contribs_data->at(parentID);

		contributions->push_back(contrib_data);
	}

	return 1.0;
}


edm4hep::CaloHitContribution GetCaloHitContribution(edm4hep::CaloHitContributionData calo_data)
{
	edm4hep::CaloHitContribution caloHitCont(calo_data.PDG, calo_data.energy, calo_data.time, calo_data.stepPosition);

	return caloHitCont;
}


int GetMCParticleDataFromCaloHitContributions(vector<edm4hep::MCParticleData> *MCParticles_data, vector<podio::ObjectID> *MCparents_data, vector<edm4hep::MCParticleData> *parents)
{
	//for (int i = calo_data.parents_begin; i < calo_data.parents_end; ++i) {
	for (int i = 0; i < MCparents_data->size(); ++i) {

	unsigned parentID = MCparents_data->at(i).index;

	edm4hep::MCParticleData mcpart_data_parent = MCParticles_data->at(parentID);

	parents->push_back(mcpart_data_parent);
	}

	return 1.0;
}


int GetMCParticleIdFromCaloHitContributions(vector<podio::ObjectID> *MCparents_data, vector<unsigned> *parents)
{
	//for (int i = calo_data.parents_begin; i < calo_data.parents_end; ++i) {
	for (int i = 0; i < MCparents_data->size(); ++i) {

	unsigned parentID = MCparents_data->at(i).index;

	//edm4hep::MCParticleData mcpart_data_parent = MCParticles_data->at(parentID);

	parents->push_back(parentID);
	}

	return 1.0;
}


int CreateMCParticleToCaloHitContributionMap(vector<edm4hep::MCParticleData> *MCParticles_data, vector<podio::ObjectID> *MCparents_data, map<int, vector<int>> *map_vec)
{

	for (int i = 0; i < MCparents_data->size(); ++i) {

	unsigned parentID = MCparents_data->at(i).index;

	edm4hep::MCParticleData mcpart_data_parent = MCParticles_data->at(parentID);

	map<int, vector<int>>::iterator it;
	it = map_vec->find(parentID);

	if (it != map_vec->end()) it->second.push_back(i);
	else{

		vector<int> vec = {i};

		map_vec->insert(it, pair<int, vector<int>>(parentID, vec));
		//cout<<"EICutil inserted = "<<parentID<<","<<i<<endl;
		}
/*
	it = map_vec->find(parentID);

	cout<<"EICutil size[id="<<it->first<<"] ="<<it->second.size()<<endl;
	for (int j = 0; j < it->second.size(); ++j) {

		cout<<it->second[j]<<", ";
	}
	cout<<endl;*/

	}

	return 1.0;
}

#endif /* EICUTIL_H_ */
