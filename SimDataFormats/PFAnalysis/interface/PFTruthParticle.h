#ifndef SimDataFormats_PFTruthParticle_h
#define SimDataFormats_PFTruthParticle_h

#include <vector>
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimClusterFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//
// Forward declarations
//
class TrackingVertex;
class SimTrack;
class EncodedEventId;

class PFTruthParticle {
  friend std::ostream& operator<<(std::ostream& s, PFTruthParticle const& tp);

public:
  typedef math::XYZTLorentzVectorD LorentzVector;           ///< Lorentz vector
  typedef math::XYZPointD Point;                            ///< point in the space
  typedef math::XYZVectorD Vector;                          ///< point in the space

  /** @brief Default constructor. Note that the object will be useless until it is provided
     * with a SimTrack and parent TrackingVertex.
     *
     * Most of the methods assume there is a SimTrack and parent TrackingVertex set, so will either
     * crash or give undefined results if this isn't true. This constructor should only be used to
     * create a placeholder until setParentVertex() and addG4Track() can be called.
     */
  PFTruthParticle();

  PFTruthParticle(const TrackingParticleRefVector& trackingParticles, const SimClusterRefVector& simClusters);

  // destructor
  ~PFTruthParticle();
    void setTrackingParticles(const TrackingParticleRefVector& refs);
    void setSimClusters(const SimClusterRefVector& refs);
    void setPdgId(int pdgId);
    void setCharge(int charge);
    void setP4(LorentzVector p4);
    void addSimCluster(const SimClusterRef sc);
    void addTrackingParticle(const TrackingParticleRef tp);

    SimClusterRefVector& simClusters() { return simClusters_; }
    TrackingParticleRefVector& trackingParticles() { return trackingParticles_; }
    size_t nSimCluster() const { return simClusters_.size(); }
    size_t nTrackingParticle() const { return trackingParticles_.size(); }

  /** @brief PDG ID.
     *
     * Returns the PDG ID of the first associated gen particle. If there are no gen particles associated
     * then it returns type() from the first SimTrack. */
  int pdgId() const {
    return pdgId_;
  }

  void addG4Track(const SimTrack& t);
  
  const std::vector<SimTrack>& g4Tracks() const { return g4Tracks_; }

  /// @brief Electric charge. Note this is taken from the first SimTrack only.
  float charge() const { return charge_; }

  /// @brief Four-momentum Lorentz vector. Note this is taken from the first SimTrack only.
  const LorentzVector& p4() const { return p4_; }

  /// @brief spatial momentum vector
  Vector momentum() const { return p4().Vect(); }

  /// @brief Magnitude of momentum vector. Note this is taken from the first SimTrack only.
  double p() const { return p4().P(); }

  /// @brief Energy. Note this is taken from the first SimTrack only.
  double energy() const { return p4().E(); }

  /// @brief Transverse energy. Note this is taken from the first SimTrack only.
  double et() const { return p4().Et(); }

  /// @brief Mass. Note this is taken from the first SimTrack only.
  double mass() const { return p4().M(); }

  /// @brief Mass squared. Note this is taken from the first SimTrack only.
  double massSqr() const { return pow(mass(), 2); }

  /// @brief Transverse mass. Note this is taken from the first SimTrack only.
  double mt() const { return p4().Mt(); }

  /// @brief Transverse mass squared. Note this is taken from the first SimTrack only.
  double mtSqr() const { return p4().Mt2(); }

  /// @brief x coordinate of momentum vector. Note this is taken from the first SimTrack only.
  double px() const { return p4().Px(); }

  /// @brief y coordinate of momentum vector. Note this is taken from the first SimTrack only.
  double py() const { return p4().Py(); }

  /// @brief z coordinate of momentum vector. Note this is taken from the first SimTrack only.
  double pz() const { return p4().Pz(); }

  /// @brief Transverse momentum. Note this is taken from the first SimTrack only.
  double pt() const { return p4().Pt(); }

  /// @brief Momentum azimuthal angle. Note this is taken from the first SimTrack only.
  double phi() const { return p4().Phi(); }

  /// @brief Momentum polar angle. Note this is taken from the first SimTrack only.
  double theta() const { return p4().Theta(); }

  /// @brief Momentum pseudorapidity. Note this is taken from the first SimTrack only.
  double eta() const { return p4().Eta(); }

  /// @brief Rapidity. Note this is taken from the first SimTrack only.
  double rapidity() const { return p4().Rapidity(); }

  /// @brief Same as rapidity().
  double y() const { return rapidity(); }

private:
  /// references to G4 and reco::GenParticle tracks
  int charge_;
  int pdgId_;
  LorentzVector p4_;
  std::vector<SimTrack> g4Tracks_;
  reco::GenParticleRefVector genParticles_;
  SimClusterRefVector simClusters_;
  TrackingParticleRefVector trackingParticles_;
};

#endif  // SimDataFormats_PFTruthParticle_H

