/*
 * SimClusterMerger.h
 *
 *  Created on: 26 Nov 2020
 *      Author: jkiesele
 */

#ifndef SIMGENERAL_CALOANALYSIS_INTERFACE_SIMCLUSTERMERGER_H_
#define SIMGENERAL_CALOANALYSIS_INTERFACE_SIMCLUSTERMERGER_H_


#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/SimClusterFwd.h"


#include <vector>
#include <algorithm>
//for now here, needs to be moved

class SimClusterMerger{

    /*
     * change a bit to
     * addSC(sc, vecto<int> its history)
     * and have simtracks stored. no need for keeping full history
     */

public:
    SimClusterMerger();
    SimClusterMerger(const std::vector<SimTrack>& v);

    void setEventSimTracks(const std::vector<SimTrack>& v){
        eventSimtracks_=&v;
    }

    /**
     * Add a simcluster to the ones to be merged.
     * If a history is provided, and the simcluster can be identified
     * by an index in that history, also add the history index.
     */
    void addSimCluster(const SimCluster&, std::vector<const SimTrack*> scHistory);

    void merge();



    const math::XYZTLorentzVectorF& mergedMomentum()const{return mergedMomentum_;}
    const math::XYZTLorentzVectorF& mergedPosition()const{return mergedPosition_;}
    /**
     * if the contribution of clusters without a valid history is above 1-pidRelEnTresh_,
     * e.g. because pileup information has not been stored,
     * the pdgId cannot be determined unambiguously.
     * In this case, this function returns false.
     *
     * This should not be confused with an ambiguous ID (here defined as pdgID=0),
     * e.g. when two showers without a common ancestor fully overlap, but have tracable history.
     */
    bool pdgIdValid()const{return idvalid_;}
    int mergedPdgID()const{return mergedId_;}

    /**
     * Unless this function is called, only momentum, position, and pdgId are determined.
     * Only by calling this function, also the hits are merged.
     */
    SimCluster constructMergedCluster()const;

    void clear();

private:
    const std::vector<SimTrack>* eventSimtracks_;


    float pidRelEnTresh_, pidValidThresh_;

    std::vector<const SimCluster*> simclusters_;
    std::vector<std::vector<const SimTrack *> > histories_;

    math::XYZTLorentzVectorF mergedMomentum_;
    math::XYZTLorentzVectorF mergedPosition_;
    int mergedId_;
    bool idvalid_;

    //helper
    std::vector<const SimCluster* > getAssoClusters(const SimTrack * t,
                    const std::vector<const SimCluster* >& masked=std::vector<const SimCluster* >())const;

};






#endif /* SIMGENERAL_CALOANALYSIS_INTERFACE_SIMCLUSTERMERGER_H_ */
