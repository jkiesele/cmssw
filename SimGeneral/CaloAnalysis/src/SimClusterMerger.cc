/*
 * SimClusterMerger.cc
 *
 *  Created on: 26 Nov 2020
 *      Author: jkiesele
 */

#include "../interface/SimClusterMerger.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include <algorithm>

SimClusterMerger::SimClusterMerger():eventSimtracks_(0),
pidRelEnTresh_(0.9), pidValidThresh_(0.9){
    clear(); //set defaults
}

SimClusterMerger::SimClusterMerger(const std::vector<SimTrack>& v):
        pidRelEnTresh_(0.9), pidValidThresh_(0.9){
    clear();
    setEventSimTracks(v);
}

void SimClusterMerger::addSimCluster(const SimCluster& sc, std::vector<const SimTrack*> scHistory){

    simclusters_.push_back(&sc);
    std::vector<const SimTrack*> thishist;
    if(!eventSimtracks_){
        histories_.push_back(thishist);
        return;
    }
    histories_.push_back(scHistory);

}


std::vector<const SimCluster* > SimClusterMerger::getAssoClusters(const SimTrack * t,
       const std::vector<const SimCluster* >& masked)const{

    std::vector<const SimCluster* >  out;
    for(size_t i=0;i<simclusters_.size();i++){
        const auto& sc=simclusters_.at(i);
        const auto hist = histories_.at(i);
        if(std::find(masked.begin(),masked.end(),sc)!=masked.end())
            continue;
        if(std::find(hist.begin(),hist.end(),t) != hist.end())
            out.push_back(sc);
    }
    return out;
}

void SimClusterMerger::merge(){
    //sets
    /*
     *
    mergedMomentum_.SetXYZT(0,0,0,0);
    mergedPosition_.SetXYZT(0,0,0,0);
    idvalid_=false;
    mergedId_=0;
     */

    //determine total momentum and position
    mergedMomentum_.SetXYZT(0,0,0,0);
    for(const auto& sc: simclusters_){
        mergedMomentum_+=sc->caloSurfaceMomentum();
        mergedPosition_+=sc->caloSurfacePosition() * sc->energy();
    }
    mergedPosition_ /= mergedMomentum_.energy();
    idvalid_=false;
    mergedId_=0;

    if(!eventSimtracks_)
        return;


    //ID is a bit more complicated
    //search for common ancestors in the history, determine groups

    //look at all simtracks and see how many and which clusters have them in common.

    // start from the back (so from PV), and see how many SC have that track in common.
    // only if at least two have that track in common (otherwise one could have gone somewhere else)
    // still check the next track. If the number does not decrease, take that next track
    // (one must have gone somewhere else)

    // mask all of them for the next track
    // do that for each not yet masked SC

    std::vector<const SimCluster* >  masked;
    std::vector<std::pair<const SimTrack*, std::vector<const SimCluster* > > > groups;

    for(size_t i=0;i<simclusters_.size();i++){
        const auto& sc = simclusters_.at(i);

        if(std::find(masked.begin(),masked.end(),sc)!=masked.end())
            continue;

        const auto& hist = histories_.at(i);
        std::vector<const SimCluster* > assoClusters;

        for(size_t j=0;j<hist.size();j++){
            const auto& st = hist.at(hist.size()-j-1);
            auto asso = getAssoClusters(st,masked);

            if(asso.size()<assoClusters.size()){//further up the chain, less are associated
                //save and break
                masked.insert(masked.end(),assoClusters.begin(),assoClusters.end());
                groups.emplace_back(st, assoClusters);
                break;//and take assoClusters
            }
            else{
                assoClusters=asso;
            }
        }
    }

    //check if the ID can be valid:
    float tracked_energy=0;
    for(const auto& group: groups)
        tracked_energy += group.first->momentum().E();

    float totalenergy=mergedMomentum_.energy();
    if(tracked_energy/totalenergy > pidValidThresh_)
        idvalid_=true;
    else{
        idvalid_=false;
        mergedId_=0;
    }
    if(idvalid_)
        for(const auto& group: groups){
            if(group.first->momentum().E()/totalenergy > pidRelEnTresh_){
                mergedId_ = group.first->type();
                break;
            }
        }

    //FIXME debug: print everything
    std::cout << "merged "<< simclusters_.size()  <<" simclusters to simcluster with energy "
            << mergedMomentum_.energy() << " at " << mergedPosition_ << " and pdgid " <<  mergedId_
            << " which is valid: "<<idvalid_ << std::endl;

}

SimCluster SimClusterMerger::constructMergedCluster()const{
    if(!simclusters_.size())
        return SimCluster();

    SimCluster sc(*simclusters_.at(0));
    sc.setCaloSurfaceMomentum(mergedMomentum_);
    sc.setCaloSurfacePosition(mergedPosition_);
    sc.setPdgId(mergedId_);

    for(const auto & s: simclusters_){
        const auto hafs = s->hits_and_fractions();
        for(const auto& haf: hafs)
            sc.addDuplicateRecHitAndFraction(haf.first,haf.second);
    }
    return sc;
}

void SimClusterMerger::clear(){

    simclusters_.clear();
    histories_.clear();

    mergedMomentum_.SetXYZT(0,0,0,0);
    mergedPosition_.SetXYZT(0,0,0,0);
    idvalid_=false;
    mergedId_=0;

}

