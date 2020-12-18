/*! \brief   
 *  \details Here, in the source file, the methods which do depend
 *           on the specific type <T> that can fit the template.
 *
 *  \author Nicola Pozzobon
 *  \date   2013, Jul 19
 *
 */

#include "SimTracker/TrackTriggerAssociation/plugins/TTClusterAssociator.h"

/// Implement the producer
template <>
void TTClusterAssociator<Ref_Phase2TrackerDigi_>::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  /// Exit if real data
  if (iEvent.isRealData())
    return;

  /// Get the PixelDigiSimLink
  iEvent.getByToken(digisimLinkToken, thePixelDigiSimLinkHandle);

  /// Get the TrackingParticles

  iEvent.getByToken(tpToken, trackingParticleHandle);

  //  const TrackerTopology* const tTopo = theTrackerTopology.product();
  const TrackerGeometry* const theTrackerGeom = theTrackerGeometry.product();

  /// Preliminary task: map SimTracks by TrackingParticle
  /// Prepare the map
  std::map<std::pair<unsigned int, EncodedEventId>, TrackingParticlePtr> simTrackUniqueToTPMap;
  simTrackUniqueToTPMap.clear();

  if (not trackingParticleHandle->empty()) {
    /// Loop over TrackingParticles
    std::vector<TrackingParticle>::const_iterator iterTPart;
    for (unsigned int tpCnt = 0; tpCnt < trackingParticleHandle->size(); tpCnt++) {
      /// Make the pointer to the TrackingParticle
      TrackingParticlePtr tempTPPtr(trackingParticleHandle, tpCnt);

      /// Get the EncodedEventId
      EncodedEventId eventId = EncodedEventId(tempTPPtr->eventId());

      /// Loop over SimTracks inside TrackingParticle
      for (const auto& simTrack : tempTPPtr->g4Tracks()) {
        /// Build the unique SimTrack Id (which is SimTrack ID + EncodedEventId)
        std::pair<unsigned int, EncodedEventId> simTrackUniqueId(simTrack.trackId(), eventId);
        simTrackUniqueToTPMap.insert(std::make_pair(simTrackUniqueId, tempTPPtr));
      }
    }  /// End of loop over TrackingParticles
  }

  /// Loop over InputTags to handle multiple collections

  int ncont1 = 0;

  for (const auto& iTag : ttClustersTokens) {
    /// Prepare output
    auto associationMapForOutput = std::make_unique<TTClusterAssociationMap<Ref_Phase2TrackerDigi_>>();

    /// Get the Clusters already stored away
    edm::Handle<TTClusterDetSetVec> TTClusterHandle;

    iEvent.getByToken(iTag, TTClusterHandle);

    /// Prepare the necessary maps
    std::map<TTClusterRef, std::vector<TrackingParticlePtr>> clusterToTrackingParticleVectorMap;
    std::map<TrackingParticlePtr, std::vector<TTClusterRef>> trackingParticleToClusterVectorMap;
    clusterToTrackingParticleVectorMap.clear();
    trackingParticleToClusterVectorMap.clear();

    /// Loop over the input Clusters
    for (const auto& gd : theTrackerGeom->dets()) {
      DetId detid = gd->geographicalId();
      if (detid.subdetId() != StripSubdetector::TOB && detid.subdetId() != StripSubdetector::TID)
        continue;  // only run on OT

      if (TTClusterHandle->find(detid) == TTClusterHandle->end())
        continue;

      /// Get the DetSets of the Clusters
      edmNew::DetSet<TTCluster<Ref_Phase2TrackerDigi_>> clusters = (*TTClusterHandle)[detid];

      for (auto contentIter = clusters.begin(); contentIter != clusters.end(); ++contentIter) {
        /// Make the reference to be put in the map
        TTClusterRef tempCluRef = edmNew::makeRefTo(TTClusterHandle, contentIter);

        /// Prepare the maps wrt TTCluster
        if (clusterToTrackingParticleVectorMap.find(tempCluRef) == clusterToTrackingParticleVectorMap.end()) {
          std::vector<TrackingParticlePtr> tpVector;
          tpVector.clear();
          clusterToTrackingParticleVectorMap.insert(std::make_pair(tempCluRef, tpVector));
        }

        /// Get the PixelDigiSimLink
        /// Safety check added after new digitizer (Oct 2014)
        if (thePixelDigiSimLinkHandle->find(detid) == thePixelDigiSimLinkHandle->end()) {
          /// Sensor is not found in DigiSimLink.
          /// Set MC truth to NULL for all hits in this sensor. Period.

          /// Get the Digis and loop over them
          std::vector<Ref_Phase2TrackerDigi_> theseHits = tempCluRef->getHits();
          for (unsigned int i = 0; i < theseHits.size(); i++) {
            /// No SimLink is found by definition
            /// Then store NULL MC truth for all the digis
            TrackingParticlePtr tempTPPtr;  // = new edm::Ptr< TrackingParticle >();
            clusterToTrackingParticleVectorMap.find(tempCluRef)->second.push_back(tempTPPtr);
          }

          /// Go to the next sensor
          continue;
        }

        edm::DetSet<PixelDigiSimLink> thisDigiSimLink = (*(thePixelDigiSimLinkHandle))[detid];
        edm::DetSet<PixelDigiSimLink>::const_iterator iterSimLink;

        /// Get the Digis and loop over them
        std::vector<Ref_Phase2TrackerDigi_> theseHits = tempCluRef->getHits();
        for (unsigned int i = 0; i < theseHits.size(); i++) {
          /// Loop over PixelDigiSimLink
          for (iterSimLink = thisDigiSimLink.data.begin(); iterSimLink != thisDigiSimLink.data.end(); iterSimLink++) {
            /// Find the link and, if there's not, skip
            if (static_cast<int>(iterSimLink->channel()) != static_cast<int>(theseHits.at(i)->channel()))
              continue;

            /// Get SimTrack Id and type
            unsigned int curSimTrkId = iterSimLink->SimTrackId();
            EncodedEventId curSimEvId = iterSimLink->eventId();

            /// Prepare the SimTrack Unique ID
            std::pair<unsigned int, EncodedEventId> thisUniqueId = std::make_pair(curSimTrkId, curSimEvId);

            /// Get the corresponding TrackingParticle
            if (simTrackUniqueToTPMap.find(thisUniqueId) != simTrackUniqueToTPMap.end()) {
              TrackingParticlePtr thisTrackingParticle = simTrackUniqueToTPMap.find(thisUniqueId)->second;

              /// Store the TrackingParticle
              clusterToTrackingParticleVectorMap.find(tempCluRef)->second.push_back(thisTrackingParticle);

              /// Prepare the maps wrt TrackingParticle
              if (trackingParticleToClusterVectorMap.find(thisTrackingParticle) ==
                  trackingParticleToClusterVectorMap.end()) {
                std::vector<TTClusterRef> clusterVector;
                clusterVector.clear();
                trackingParticleToClusterVectorMap.insert(std::make_pair(thisTrackingParticle, clusterVector));
              }
              trackingParticleToClusterVectorMap.find(thisTrackingParticle)
                  ->second.push_back(tempCluRef);  /// Fill the auxiliary map
            } else {
              /// In case no TrackingParticle is found, store a NULL pointer

              TrackingParticlePtr tempTPPtr;  // = new edm::Ptr< TrackingParticle >();
              clusterToTrackingParticleVectorMap.find(tempCluRef)->second.push_back(tempTPPtr);
            }
          }  /// End of loop over PixelDigiSimLink
        }    /// End of loop over all the hits composing the Cluster

        /// Check that the cluster has a non-NULL TP pointer
        std::vector<TrackingParticlePtr> theseClusterTrackingParticlePtrs =
            clusterToTrackingParticleVectorMap.find(tempCluRef)->second;
        bool allOfThemAreNull = true;
        for (unsigned int tpi = 0; tpi < theseClusterTrackingParticlePtrs.size() && allOfThemAreNull; tpi++) {
          if (theseClusterTrackingParticlePtrs.at(tpi).isNull() == false)
            allOfThemAreNull = false;
        }

        if (allOfThemAreNull) {
          /// In case no TrackingParticle is found at all, drop the map element
          clusterToTrackingParticleVectorMap.erase(tempCluRef);  /// Use "erase by key"
        }
      }
    }  /// End of loop over all the TTClusters of the event

    /// Clean the maps that need cleaning
    /// Prepare the output map wrt TrackingParticle
    std::map<TrackingParticlePtr, std::vector<TTClusterRef>>::iterator iterMapToClean;
    for (iterMapToClean = trackingParticleToClusterVectorMap.begin();
         iterMapToClean != trackingParticleToClusterVectorMap.end();
         ++iterMapToClean) {
      /// Get the vector of references to TTCluster
      std::vector<TTClusterRef> tempVector = iterMapToClean->second;

      /// Sort and remove duplicates
      std::sort(tempVector.begin(), tempVector.end());
      tempVector.erase(std::unique(tempVector.begin(), tempVector.end()), tempVector.end());
      iterMapToClean->second = tempVector;
    }

    /// Put the maps in the association object
    associationMapForOutput->setTTClusterToTrackingParticlesMap(clusterToTrackingParticleVectorMap);
    associationMapForOutput->setTrackingParticleToTTClustersMap(trackingParticleToClusterVectorMap);

    /// Put output in the event
    iEvent.put(std::move(associationMapForOutput), ttClustersInputTags.at(ncont1).instance());

    ++ncont1;

  }  /// End of loop over input tags
}
