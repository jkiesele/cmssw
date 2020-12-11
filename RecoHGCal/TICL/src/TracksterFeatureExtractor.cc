/*
 * TracksterFeatureExtractor.h
 *
 *  Created on: 10 Nov 2020
 *      Author: jkiesele
 */

#include "RecoHGCal/TICL/interface/TrackstersFeatureExtractor.h"
#include "FWCore/Utilities/interface/EDMException.h"

namespace ticl{

TrackstersFeatureExtractor::TrackstersFeatureExtractor(): layerClusters_(0){
}

void TrackstersFeatureExtractor::calculateFeatures(){
    /*
     * Do not calculate anything here, just fill
     */
    clusterFeatures_.clear();
    clusterRowSplits_.clear();
    globalFeatures_.clear();

    if(!layerClusters_)
        throw edm::Exception(edm::errors::NullPointerError,
                        "First assign layer clusters before calculating features");

    size_t clustercount=0;//for row splits
    clusterRowSplits_.push_back(clustercount);//TF format, begins with 0
    for(const auto& trackster: tracksters_){

        for(const auto& clidx: trackster->vertices()){
            auto cluster = layerClusters_->at(clidx);
            auto clusterfeats = calculateClusterFeatures(cluster,trackster);
            clusterFeatures_.push_back(clusterfeats);
        }
        clustercount+=trackster->vertices().size();
        clusterRowSplits_.push_back(clustercount);//TF format, ends with total size

        globalFeatures_.push_back(calculateTracksterFeatures(trackster));//move

    }

}

std::vector<float>  TrackstersFeatureExtractor::calculateClusterFeatures(
        const reco::CaloCluster &cl, const Trackster* tr) const {
    /*
     * This function should contain the calculation and extraction of
     * all layercluster features used.
     *
     * This should be normalised information, in extent and energy
     *
     */

    std::vector<float>  out;

    float tracksterenergy=tr->raw_energy(); //just to scale similarly for all tracksters
    float tracksterlength=tr->eigenvalues().at(0); //just to scale similarly for all tracksters

    out.push_back(cl.energy()/tracksterenergy);
    out.push_back((cl.position().x()-tr->barycenter().x())/tracksterlength);
    out.push_back((cl.position().y()-tr->barycenter().y())/tracksterlength);
    out.push_back((cl.position().z()-tr->barycenter().z())/tracksterlength);
    //...

    return out;


}

std::vector<float>  TrackstersFeatureExtractor::calculateTracksterFeatures(
        const Trackster *tr) const {
    /*
     * This function should contain the calculation and extraction of
     * all trackster features used
     */
    std::vector<float> out;
    out.push_back(tr->raw_energy());
    out.push_back(tr->eigenvalues().at(0));
    out.push_back(tr->eigenvalues().at(1));
    out.push_back(tr->eigenvalues().at(2));
    out.push_back(tr->eigenvectors(0).x());
    out.push_back(tr->eigenvectors(0).y());
    out.push_back(tr->eigenvectors(0).z());
    out.push_back(tr->eigenvectors(1).x());
    out.push_back(tr->eigenvectors(1).y());
    out.push_back(tr->eigenvectors(1).z());
    out.push_back(tr->eigenvectors(2).x());
    out.push_back(tr->eigenvectors(2).y());
    out.push_back(tr->eigenvectors(2).z());

    return out;
}


void TrackstersFeatureExtractor::clear(){
    tracksters_.clear();
    clusterFeatures_.clear();
    clusterRowSplits_.clear();
    globalFeatures_.clear();
}




////////// TF interface

tensorflow::NamedTensorList TrackstersFeatureExtractor::createTFInputs(const std::string& clusterName,
        const std::string& rowSplitName,const std::string& globalsName) const{

    if(!globalFeatures_.size())
        throw edm::Exception(edm::errors::LogicError,
                        "First assign calculate features before creating TF inputs");

    auto clusterInput=makeTFClusterInput();
    auto rsInput = makeTFRowSplitInput();
    auto globalInput = makeTFGlobalInput();

    return tensorflow::NamedTensorList( { { clusterName, clusterInput }, {
            rowSplitName, rsInput }, { globalsName, globalInput } });

}


std::vector<long long int> TrackstersFeatureExtractor::makeClusterArrShape()const{
    long long int nclf = 0;
    if(clusterFeatures_.size())
        nclf=clusterFeatures_.at(0).size();
    return {(long long int)clusterFeatures_.size(), nclf};
}

std::vector<long long int> TrackstersFeatureExtractor::makeRowSplitArrShape()const{
    return {(long long int)clusterRowSplits_.size()};
}

std::vector<long long int> TrackstersFeatureExtractor::makeGlobalArrShape()const{
    long long int nglf = 0;
    if(globalFeatures_.size())
        nglf=globalFeatures_.at(0).size();
    return {(long long int)globalFeatures_.size(),nglf};
}


tensorflow::Tensor TrackstersFeatureExtractor::makeTFClusterInput() const{

    tensorflow::TensorShape shape(makeClusterArrShape());
    tensorflow::Tensor t(tensorflow::DT_FLOAT, shape);

    for (size_t i = 0; i < clusterFeatures_.size(); i++) {
        for (size_t j = 0; j < clusterFeatures_.at(i).size(); j++)
            t.tensor<float, 2>()(i, j) = clusterFeatures_[i][j];
    }
    return t;
}
tensorflow::Tensor TrackstersFeatureExtractor::makeTFRowSplitInput() const{

    tensorflow::TensorShape shape(makeRowSplitArrShape());
    tensorflow::Tensor t(tensorflow::DT_INT64, shape);

    for (size_t i = 0; i < clusterRowSplits_.size(); i++) {
        t.tensor<int, 1>()(i) = clusterRowSplits_[i];
    }

    return t;
}
tensorflow::Tensor TrackstersFeatureExtractor::makeTFGlobalInput() const{

    tensorflow::TensorShape shape(makeGlobalArrShape());

    tensorflow::Tensor t(tensorflow::DT_FLOAT, shape);

    for (size_t i = 0; i < globalFeatures_.size(); i++) {
        for (size_t j = 0; j < globalFeatures_.at(i).size(); j++)
            t.tensor<float, 2>()(i, j) = globalFeatures_[i][j];
    }

    return t;
}



}//ticl
