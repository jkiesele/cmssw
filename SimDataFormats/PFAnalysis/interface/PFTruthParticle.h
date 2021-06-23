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
  void setTrackingParticles(const TrackingParticleRefVector &refs);
  void setSimClusters(const SimClusterRefVector &refs);

  void addSimCluster(const SimClusterRef &sc);
  void addTrackingParticle(const TrackingParticleRef &tp);

  /** @brief PDG ID.
     *
     * Returns the PDG ID of the first associated gen particle. If there are no gen particles associated
     * then it returns type() from the first SimTrack. */
  const int & pdgId() const {
    return pdgid_;
  }

  void setPdgId(int id){
      pdgid_=id;
  }

  void setP4(const math::PtEtaPhiMLorentzVectorF& m){
      p4_=m;
  }
  //just for convenience
  void setP4(const math::PtEtaPhiMLorentzVectorD& m){
      p4_=m;
  }

  void setVertex(const math::XYZTLorentzVectorF& v){
      vertex_=v;
  }
  void setVertex(const math::XYZTLorentzVectorD& v){
      vertex_=v;
  }

  /// @brief Four-momentum Lorentz vector.
  const math::PtEtaPhiMLorentzVectorF& p4() const { return p4_; }

  /// @brief spatial momentum vector
  math::XYZVectorF momentum() const { return p4().Vect(); }

  /// @brief Vector to boost to the particle centre of mass frame.
  math::XYZVectorF boostToCM() const { return p4().BoostToCM(); }

  /// @brief Magnitude of momentum vector.
  float p() const { return p4().P(); }

  /// @brief Energy.
  float energy() const { return p4().E(); }

  /// @brief Transverse energy.
  float et() const { return p4().Et(); }

  /// @brief Mass.
  float mass() const { return p4().M(); }

  /// @brief Mass squared.
  float massSqr() const { return pow(mass(), 2); }

  /// @brief Transverse mass.
  float mt() const { return p4().Mt(); }

  /// @brief Transverse mass squared.
  float mtSqr() const { return p4().Mt2(); }

  /// @brief x coordinate of momentum vector.
  float px() const { return p4().Px(); }

  /// @brief y coordinate of momentum vector.
  float py() const { return p4().Py(); }

  /// @brief z coordinate of momentum vector.
  float pz() const { return p4().Pz(); }

  /// @brief Transverse momentum.
  float pt() const { return p4().Pt(); }

  /// @brief Momentum azimuthal angle.
  float phi() const { return p4().Phi(); }

  /// @brief Momentum polar angle.
  float theta() const { return p4().Theta(); }

  /// @brief Momentum pseudorapidity.
  float eta() const { return p4().Eta(); }

  /// @brief Rapidity.
  float rapidity() const { return p4().Rapidity(); }

  /// @brief Same as rapidity().
  float y() const { return rapidity(); }


private:

  int pdgid_; //this can also be ambiguous

  //this could coincide with the first tracking particle but does not necessarily.
  //add this as members for now - but as float
  math::PtEtaPhiMLorentzVectorF p4_;
  math::XYZTLorentzVectorF vertex_; //position and time

  SimClusterRefVector simClusters_;
  TrackingParticleRefVector trackingParticles_;
};

#endif  // SimDataFormats_PFTruthParticle_H

