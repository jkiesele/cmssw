/*
 * TiclTrackstersDNNTree.cc
 *
 *  Created on: 28 Nov 2020
 *      Author: jkiesele
 */


#include <Validation/HGCalValidation/interface/TiclTrackstersDNNTree.h>
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "SimGeneral/CaloAnalysis/interface/SimClusterMerger.h"

void TiclTrackstersDNNTree::assignTruth(const std::vector<SimCluster> &scs,
        const std::vector<SimTrack>& simtracks,
        const std::vector<SimVertex>& simvertices,
        const std::vector<const HGCRecHit *>& rechits){


    if(!layerClusters_)
        throw edm::Exception(edm::errors::NullPointerError,
                        "First assign layer clusters before assigning the truth");

    std::vector<std::map<DetId, float> >  trackstersHitsAndE(tracksters_.size());
    for(const auto& t: tracksters_){

        //build combined hits and energies
        auto trksClusterIdcs = t->vertices();
        std::map<DetId, float> thisHitsAndE;
        for(const auto &clidx: trksClusterIdcs){
            auto haf = layerClusters_->at(clidx).hitsAndFractions();
            float clenergy = layerClusters_->at(clidx).energy();
            for(const auto hf: haf){
                auto hit=hf.first;
                auto fidx = thisHitsAndE.find(hit);
                if(fidx == thisHitsAndE.end())
                    thisHitsAndE[hit]=hf.second*clenergy;
                else
                    fidx->second += hf.second*clenergy;
            }
        }
        trackstersHitsAndE.push_back(thisHitsAndE);
    }
    //create rechit detid map
    std::map<DetId, const HGCRecHit*> rhmap;
    for(const auto& rh: rechits)
        rhmap[rh->id()]=rh;

    //now for simclusters, only detID needed
    std::vector<std::vector<DetId> > scDetIDs(scs.size());
    std::vector<float> scDepEnergies;
    for(size_t i=0;i<scs.size();i++){
        const auto& sc = scs.at(i);

        float depEnergy=0;
        auto haf = sc.hits_and_fractions();

        std::vector<DetId> thisScDetids;
        for(const auto& hf: haf){
            float rhenergy = rhmap[hf.first]->energy();
            depEnergy += rhenergy*hf.second;
            thisScDetids.push_back(hf.first);
        }
        scDetIDs.at(i)=thisScDetids;
        scDepEnergies.push_back(depEnergy);
    }

    //create SimCluster history here (for now!)
    std::vector<std::vector<const SimTrack*> > scHistory(scs.size());
    for(size_t i=0;i<scs.size();i++){
        const auto& sc = scs.at(i);
        auto& hist = scHistory.at(i);
        if(sc.eventId().bunchCrossing() || sc.eventId().event())
            continue; //this is PU?, history will not be covered for now

        const SimTrack * parent = &sc.g4Tracks().at(0);
        while(true){
            parent=parentSimTrack(parent,simtracks,simvertices);
            if(parent){
                hist.push_back(parent);
            }
            else{
                break;
            }
        }
    }

    //match tracksters and simclusters by intersection energies
    std::vector<int> maskedsc;
    std::vector<std::vector<size_t> > tobemergedSc(tracksters_.size());
    for (size_t i_t = 0; i_t < tracksters_.size(); i_t++) {
        const auto& thafmap = trackstersHitsAndE.at(i_t);

        float overlapenergy=0;
        for (size_t j_sc; j_sc < scs.size(); j_sc++) {
            if(std::find(maskedsc.begin(),maskedsc.end(),j_sc) != maskedsc.end())
                continue;

            for(const auto& scHit: scDetIDs.at(j_sc)){
                auto fidx = thafmap.find(scHit);
                if(fidx != thafmap.end()){
                    overlapenergy+= fidx->second;
                }
            }
            if(overlapenergy / scDepEnergies.at(j_sc) > fractionThresh_){
                if(fractionThresh_>0.5)//masking only works then
                    maskedsc.push_back(j_sc);//used now
                tobemergedSc.at(i_t).push_back(j_sc);
            }
        }
    }

    trueEnergies_.resize(tracksters_.size());
    trueIds_.resize(tracksters_.size());
    trueIdValid_.resize(tracksters_.size(),0);

    SimClusterMerger scmerger(simtracks);
    for (size_t i_t = 0; i_t < tracksters_.size(); i_t++) {
        const auto& matchedSc = tobemergedSc.at(i_t);
        for(const auto& i_sc: matchedSc){
            scmerger.addSimCluster(scs.at(i_sc),scHistory.at(i_sc));
        }
        scmerger.merge();
        std::array<float, TICLTRACKSTER_N_PARTICLETYPES> pids{};//all zeroed
        if(scmerger.pdgIdValid()){
            trueIdValid_.at(i_t)=1;
            int pididx = (int)ticl::Trackster::pdgIdToType(scmerger.mergedPdgID());
            pids[pididx]= 1.;
            trueIds_.at(i_t)=pids; //copy
            pids[pididx]= 0.;
        }
        else{
            trueIdValid_.at(i_t)=0;
            pids[(int)ticl::Trackster::ParticleType::unknown]=1.;//just to have one default
            trueIds_.at(i_t)=pids;
            pids[(int)ticl::Trackster::ParticleType::unknown]=0.;
        }
        trueEnergies_.at(i_t) = scmerger.mergedMomentum().energy();

        scmerger.clear();
    }




}

void TiclTrackstersDNNTree::fill(){
    /*
     * change from event-wise to trackster-wise tree filling
     *
     */

    for(size_t i_t=0;i_t<tracksters_.size();i_t++){

        trTrueEnergy_ = trueEnergies_.at(i_t);
        trTrueId_ = trueIds_.at(i_t);
        trTrueIdValid_ = trueIdValid_.at(i_t);

        //row split logic
        int rs_start = clusterRowSplits_.at(i_t);
        int rs_stop = clusterRowSplits_.at(i_t+1);

        //use to flatten cluster features (nested vectors in TTrees are slow)
        trClusterFeatures_.clear();
        auto tmp = std::vector<std::vector<float> > (clusterFeatures_.begin()+rs_start,
                clusterFeatures_.begin()+rs_stop);
        for(const auto& t:tmp)
            trClusterFeatures_.insert(trClusterFeatures_.begin(), t.begin(),t.end());

        trGlobalFeatures_=globalFeatures_.at(i_t+1);

        //pure validation branches do *not* use for training!
        val_trEta_ = tracksters_.at(i_t)->barycenter().eta();
        val_trPhi_ = tracksters_.at(i_t)->barycenter().phi();

        tree_->Fill();

    }


}

void TiclTrackstersDNNTree::createBranches(){

    tree_->Branch("true_energy",&trTrueEnergy_);
    tree_->Branch("true_id",&trTrueId_);
    tree_->Branch("true_id_valid",&trTrueIdValid_);
    tree_->Branch("n_layerclusters",&trNlayerClusters_);
    tree_->Branch("layercluster_feat",&trClusterFeatures_);
    tree_->Branch("global_feat",&trGlobalFeatures_);


    //pure validation branches do *not* use for training!
    tree_->Branch("val_trEta",&val_trEta_);
    tree_->Branch("val_trPhi",&val_trPhi_);


}


const SimTrack *  TiclTrackstersDNNTree::parentSimTrack(
        const SimTrack * t,
        const std::vector<SimTrack>& simtracks,
        const std::vector<SimVertex>& simvertices)const{

    auto origvert = t->vertIndex(); //origin
    if(simvertices.at(origvert).noParent()){
        return NULL;
    }
    return &simtracks.at(simvertices.at(origvert).parentIndex());

}
