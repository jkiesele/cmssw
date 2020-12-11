/*
 * TrackstersFeatureExtractor.h
 *
 *  Created on: 10 Nov 2020
 *      Author: jkiesele
 */

#ifndef RECOHGCAL_TICL_PLUGINS_TRACKSTERSFEATUREEXTRACTOR_H_
#define RECOHGCAL_TICL_PLUGINS_TRACKSTERSFEATUREEXTRACTOR_H_


#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include <vector>
#include <string>


namespace ticl {

/*
 * Class that fully comprises the input to the DNN model
 * as well as how to interpret the output.
 *
 * workflow:
 * - add tracksters
 * - trigger calculate
 * - use arrays (for training tree writing or inference)
 *
 *
 * also includes how to interpret the output of the model
 */
class TrackstersFeatureExtractor {
public:

    TrackstersFeatureExtractor(
            const std::vector<reco::CaloCluster> &layerClusters) :
            layerClusters_(&layerClusters) {
    }

    TrackstersFeatureExtractor();

    void setLayerClusters(const std::vector<reco::CaloCluster> &layerClusters){
        layerClusters_=&layerClusters;
    }

    //maybe not all should be added?
    void addTrackster(const Trackster& t){
        tracksters_.push_back(&t);
    }
    //implements all features extraction

    void calculateFeatures();

    void clear();

    tensorflow::NamedTensorList createTFInputs(const std::string& clusterName,
            const std::string& rowSplitName,const std::string& globalsName) const;


    const float& getRegressedEnergy(const std::vector<tensorflow::Tensor>& dnnout, const size_t& tracksteridx) const;
    std::array<float, TICLTRACKSTER_N_PARTICLETYPES>  getPIDProbs(const std::vector<tensorflow::Tensor>& dnnout, const size_t& tracksteridx) const;


private:

    tensorflow::Tensor makeTFClusterInput() const;
    tensorflow::Tensor makeTFRowSplitInput() const;
    tensorflow::Tensor makeTFGlobalInput() const;



    std::vector<float> calculateClusterFeatures(const reco::CaloCluster&, const Trackster*) const;
    std::vector<float> calculateTracksterFeatures(const Trackster*) const;


protected:
    //also used for the tree writer
    std::vector<const Trackster*> tracksters_;
    const std::vector<reco::CaloCluster> *layerClusters_;

    std::vector<long long int> makeClusterArrShape()const;
    std::vector<long long int> makeRowSplitArrShape()const;
    std::vector<long long int> makeGlobalArrShape()const;

    std::vector<std::vector<float> >clusterFeatures_;
    std::vector<int> clusterRowSplits_;
    std::vector<std::vector<float> > globalFeatures_;

};


}//ticl






#endif /* RECOHGCAL_TICL_PLUGINS_TRACKSTERSFEATUREEXTRACTOR_H_ */
