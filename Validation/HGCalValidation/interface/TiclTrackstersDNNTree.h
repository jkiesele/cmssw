/*
 * TiclTrackstersDNNTree.h
 *
 *  Created on: 28 Nov 2020
 *      Author: jkiesele
 */

#ifndef VALIDATION_HGCALVALIDATION_INTERFACE_TICLTRACKSTERSDNNTREE_H_
#define VALIDATION_HGCALVALIDATION_INTERFACE_TICLTRACKSTERSDNNTREE_H_

#include "RecoHGCal/TICL/interface/TrackstersFeatureExtractor.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimGeneral/CaloAnalysis/interface/SimClusterMerger.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"

#include "TTree.h"

// Makes sure inputs to tree are consistent with inputs to DNN
// TTree bindings
// does not own tree (comes from TFileService)
// single thread
/*
 *
 * The class works on event level.
 * The tree on trackster level
 *
 */
class TiclTrackstersDNNTree: public ticl::TrackstersFeatureExtractor{
public:
    TiclTrackstersDNNTree( TTree *t) :
        TrackstersFeatureExtractor(), fractionThresh_(0.6), tree_(t),
        trTrueEnergy_(0), trTrueIdValid_(0), trNlayerClusters_(0){
        createBranches();
    }

    void setTruthFractionThreshold(float thresh){
        fractionThresh_=thresh;
    }

    void assignTruth(const std::vector<SimCluster>& scs,
            const std::vector<SimTrack>& simtracks,
            const std::vector<SimVertex>& simvertices,
            const std::vector<const HGCRecHit *>& rechits);

    void fill(); //fill and clear

private:
    //cannot live without at Tree attached to it
    TiclTrackstersDNNTree() :
        TrackstersFeatureExtractor(), fractionThresh_(0.9), tree_(0),
        trTrueEnergy_(0), trTrueIdValid_(0), trNlayerClusters_(0),val_trEta_(0),val_trPhi_(0){
    }

    void createBranches();

    const SimTrack * parentSimTrack(
            const SimTrack* t,
            const std::vector<SimTrack>& simtracks,
            const std::vector<SimVertex>& simvertices)const;

    // truth
    float fractionThresh_;

    std::vector<float> trueEnergies_;
    std::vector<std::array<float, TICLTRACKSTER_N_PARTICLETYPES> > trueIds_;
    std::vector<int> trueIdValid_;

    TTree* tree_;

    //trackster-wise properties
    float trTrueEnergy_;
    std::array<float, TICLTRACKSTER_N_PARTICLETYPES> trTrueId_;
    int trTrueIdValid_;
    int trNlayerClusters_; //to flatten
    std::vector<float> trClusterFeatures_;
    std::vector<float> trGlobalFeatures_;


    /////pure validation quantities. These should *not* be used for the training, there can be overlap.
    float val_trEta_;
    float val_trPhi_;
    //etc...


};






#endif /* VALIDATION_HGCALVALIDATION_INTERFACE_TICLTRACKSTERSDNNTREE_H_ */
