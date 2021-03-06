#ifndef HIBestVertexProducer_H
#define HIBestVertexProducer_H

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

namespace edm {
  class Event;
  class EventSetup;
  class ConfigurationDescriptions;
}  // namespace edm

class HIBestVertexProducer : public edm::stream::EDProducer<> {
public:
  explicit HIBestVertexProducer(const edm::ParameterSet& ps);
  ~HIBestVertexProducer() override;
  void produce(edm::Event& ev, const edm::EventSetup& es) override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob();
  edm::ParameterSet theConfig;
  edm::EDGetTokenT<reco::BeamSpot> theBeamSpotTag;
  edm::EDGetTokenT<reco::VertexCollection> theMedianVertexCollection;
  edm::EDGetTokenT<reco::VertexCollection> theAdaptiveVertexCollection;
  edm::EDGetTokenT<reco::VertexCollection> theFinalAdaptiveVertexCollection;
  bool theUseFinalAdaptiveVertexCollection;
};
#endif
