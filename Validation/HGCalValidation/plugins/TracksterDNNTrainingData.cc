// -*- C++ -*-

//
// Original Author:  Jan Kieseler
//         Created:  Wed, 22 Nov 2020 11:57:20 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

#include "../interface/TiclTrackstersDNNTree.h"



class TracksterDNNTrainingData : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit TracksterDNNTrainingData(const edm::ParameterSet&);
  ~TracksterDNNTrainingData();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;


  edm::EDGetTokenT<std::vector<ticl::Trackster>> tracksterToken_;
  edm::EDGetTokenT<reco::CaloClusterCollection> layerClusterToken_;
  edm::EDGetTokenT<std::vector<SimCluster>> simClusterToken_;
  edm::EDGetTokenT<std::vector<SimVertex>> simVertexToken_;
  edm::EDGetTokenT<std::vector<SimTrack>> simTrackToken_;
  std::vector<edm::EDGetTokenT<HGCRecHitCollection>> recHitTokens_;


  edm::Service<TFileService> fs_;
  TTree* outTree_;
  TiclTrackstersDNNTree * tracksterTree_;

  float truthFractionThresh_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TracksterDNNTrainingData::TracksterDNNTrainingData(const edm::ParameterSet& iConfig)
    : tracksterToken_(consumes<std::vector<ticl::Trackster>>(iConfig.getUntrackedParameter<edm::InputTag>("tracksters"))),
      layerClusterToken_(consumes<reco::CaloClusterCollection>(iConfig.getUntrackedParameter<edm::InputTag>("layerclusters"))),
      simClusterToken_(consumes<std::vector<SimCluster>>(iConfig.getUntrackedParameter<edm::InputTag>("simclusters"))),
      simVertexToken_(consumes<std::vector<SimVertex>>(iConfig.getUntrackedParameter<edm::InputTag>("simvertices"))),
      simTrackToken_(consumes<std::vector<SimTrack>>(iConfig.getUntrackedParameter<edm::InputTag>("simtracks"))),
      outTree_(0),
      tracksterTree_(0),
      truthFractionThresh_(0.6){

    for (edm::InputTag& recHitCollection : iConfig.getParameter<
            std::vector<edm::InputTag> >("recHitCollections")) {
        recHitTokens_.push_back(
                consumes<HGCRecHitCollection>(recHitCollection));
    }

}

TracksterDNNTrainingData::~TracksterDNNTrainingData() {
    if(tracksterTree_)
        delete tracksterTree_;
    tracksterTree_=0;
}

//
// member functions
//

// ------------ method called for each event  ------------
void TracksterDNNTrainingData::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  auto tracksters = iEvent.get(tracksterToken_);
  auto layerclusters = iEvent.get(layerClusterToken_);
  auto simclusters = iEvent.get(simClusterToken_);
  auto simvertices = iEvent.get(simVertexToken_);
  auto simtracks = iEvent.get(simTrackToken_);

  std::vector<const HGCRecHit *> rechits;
  for(auto& t: recHitTokens_){
      for(const auto& rh: iEvent.get(t))
          rechits.push_back(&rh);
  }

  tracksterTree_->setLayerClusters(layerclusters);
  //first same as inference loop
  for(const auto& t: tracksters){
      tracksterTree_->addTrackster(t);
  }

  tracksterTree_->calculateFeatures();

  tracksterTree_->assignTruth(simclusters,simtracks,simvertices,rechits);

  tracksterTree_->fill();

  tracksterTree_->clear();

}

// ------------ method called once each job just before starting event loop  ------------
void TracksterDNNTrainingData::beginJob() {
    if (!fs_)
        throw edm::Exception(edm::errors::Configuration,
                "TFile Service is not registered in cfg file");

    outTree_ = fs_->make<TTree>("tree", "tree");
    tracksterTree_ = new TiclTrackstersDNNTree(outTree_);
    tracksterTree_->setTruthFractionThreshold(truthFractionThresh_);
}

// ------------ method called once each job just after ending the event loop  ------------
void TracksterDNNTrainingData::endJob() {
  // please remove this method if not needed
    delete tracksterTree_;
    tracksterTree_=0;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TracksterDNNTrainingData::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TracksterDNNTrainingData);
