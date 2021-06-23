#include "SimDataFormats/PFAnalysis/interface/PFTruthParticle.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <FWCore/MessageLogger/interface/MessageLogger.h>

PFTruthParticle::PFTruthParticle() {
  // No operation
}

PFTruthParticle::~PFTruthParticle() {}

PFTruthParticle::PFTruthParticle(const TrackingParticleRefVector& trackingParticles, const SimClusterRefVector& simClusters) {
    setTrackingParticles(trackingParticles);
    setSimClusters(simClusters);
}

void PFTruthParticle::setTrackingParticles(const TrackingParticleRefVector& refs) { trackingParticles_ = refs; }

void PFTruthParticle::setSimClusters(const SimClusterRefVector& refs) { simClusters_ = refs; }

void PFTruthParticle::addSimCluster(const SimClusterRef& sc) { simClusters_.push_back(sc); }

void PFTruthParticle::addTrackingParticle(const TrackingParticleRef& tp) { trackingParticles_.push_back(tp); }
