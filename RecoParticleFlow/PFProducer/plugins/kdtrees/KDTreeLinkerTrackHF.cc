#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "RecoParticleFlow/PFProducer/interface/KDTreeLinkerBase.h"
#include "CommonTools/RecoAlgos/interface/KDTreeLinkerAlgo.h"

// This class is used to find all links between Tracks and HF clusters
// using a KDTree algorithm.
// It is used in PFBlockAlgo.cc in the function links().
class KDTreeLinkerTrackHF : public KDTreeLinkerBase {
public:
  KDTreeLinkerTrackHF();
  ~KDTreeLinkerTrackHF() override;

  // With this method, we create the list of psCluster that we want to link.
  void insertTargetElt(reco::PFBlockElement* track) override;

  // Here, we create the list of hfCluster that we want to link. From hfCluster
  // and fraction, we will create a second list of rechits that will be used to
  // build the KDTree.
  void insertFieldClusterElt(reco::PFBlockElement* hfCluster) override;

  // The KDTree building from rechits list.
  void buildTree() override;

  // Here we will iterate over all tracks. For each track intersection point with HF,
  // we will search the closest rechits in the KDTree, from rechits we will find the
  // hfClusters and after that we will check the links between the track and
  // all closest hfClusters.
  void searchLinks() override;

  // Here, we will store all PS/HF founded links in the PFBlockElement class
  // of each psCluster in the PFmultilinks field.
  void updatePFBlockEltWithLinks() override;

  // Here we free all allocated structures.
  void clear() override;

private:
  // Data used by the KDTree algorithm : sets of Tracks and HF clusters.
  BlockEltSet targetSet_;
  BlockEltSet fieldClusterSet_;

  // Sets of rechits that compose the HF clusters.
  RecHitSet rechitsSet_;

  // Map of linked Track/HF clusters.
  BlockElt2BlockEltMap cluster2TargetLinks_;

  // Map of the HF clusters associated to a rechit.
  RecHit2BlockEltMap rechit2ClusterLinks_;

  // KD trees
  KDTreeLinkerAlgo<reco::PFRecHit const*> tree_;
};

// the text name is different so that we can easily
// construct it when calling the factory
DEFINE_EDM_PLUGIN(KDTreeLinkerFactory, KDTreeLinkerTrackHF, "KDTreeTrackAndHFLinker");

KDTreeLinkerTrackHF::KDTreeLinkerTrackHF() : KDTreeLinkerBase() {
  cristalPhiEtaMaxSize_ = 0.2;
  phiOffset_ = 0.32;
}

KDTreeLinkerTrackHF::~KDTreeLinkerTrackHF() { clear(); }

void KDTreeLinkerTrackHF::insertTargetElt(reco::PFBlockElement* track) {
  if (track->trackRefPF()->extrapolatedPoint(reco::PFTrajectoryPoint::VFcalEntrance).isValid()) {
    targetSet_.insert(track);
  }
}

void KDTreeLinkerTrackHF::insertFieldClusterElt(reco::PFBlockElement* hfCluster) {
  const reco::PFClusterRef& clusterref = hfCluster->clusterRef();

  // This test is more or less done in PFBlockAlgo.h. In others cases, it should be switch on.
  //   if (!((clusterref->layer() == PFLayer::HF_ENDCAP) ||
  // 	(clusterref->layer() == PFLayer::HF_BARREL1)))
  //     return;

  const std::vector<reco::PFRecHitFraction>& fraction = clusterref->recHitFractions();

  // We create a list of hfCluster
  fieldClusterSet_.insert(hfCluster);
  for (size_t rhit = 0; rhit < fraction.size(); ++rhit) {
    const reco::PFRecHitRef& rh = fraction[rhit].recHitRef();
    double fract = fraction[rhit].fraction();

    if ((rh.isNull()) || (fract < 1E-4))
      continue;

    const reco::PFRecHit& rechit = *rh;

    // We save the links rechit to HFClusters
    rechit2ClusterLinks_[&rechit].insert(hfCluster);

    // We create a liste of rechits
    rechitsSet_.insert(&rechit);
  }
}

void KDTreeLinkerTrackHF::buildTree() {
  // List of pseudo-rechits that will be used to create the KDTree
  std::vector<KDTreeNodeInfo<reco::PFRecHit const*, 2>> eltList;

  // Filling of this list
  for (RecHitSet::const_iterator it = rechitsSet_.begin(); it != rechitsSet_.end(); it++) {
    const reco::PFRecHit::REPPoint& posrep = (*it)->positionREP();

    KDTreeNodeInfo<reco::PFRecHit const*, 2> rh1(*it, posrep.eta(), posrep.phi());
    eltList.push_back(rh1);

    // Here we solve the problem of phi circular set by duplicating some rechits
    // too close to -Pi (or to Pi) and adding (substracting) to them 2 * Pi.
    if (rh1.dims[1] > (M_PI - phiOffset_)) {
      float phi = rh1.dims[1] - 2 * M_PI;
      KDTreeNodeInfo<reco::PFRecHit const*, 2> rh2(*it, float(posrep.eta()), phi);
      eltList.push_back(rh2);
    }

    if (rh1.dims[1] < (M_PI * -1.0 + phiOffset_)) {
      float phi = rh1.dims[1] + 2 * M_PI;
      KDTreeNodeInfo<reco::PFRecHit const*, 2> rh3(*it, float(posrep.eta()), phi);
      eltList.push_back(rh3);
    }
  }

  // Here we define the upper/lower bounds of the 2D space (eta/phi).
  float phimin = -1.0 * M_PI - phiOffset_;
  float phimax = M_PI + phiOffset_;

  // etamin-etamax, phimin-phimax
  KDTreeBox region(-3.0f, 3.0f, phimin, phimax);

  // We may now build the KDTree
  tree_.build(eltList, region);
}

void KDTreeLinkerTrackHF::searchLinks() {
  // Must of the code has been taken from LinkByRecHit.cc

  // We iterate over the tracks.
  for (BlockEltSet::iterator it = targetSet_.begin(); it != targetSet_.end(); it++) {
    reco::PFRecTrackRef trackref = (*it)->trackRefPF();

    const reco::PFTrajectoryPoint& atHF = trackref->extrapolatedPoint(reco::PFTrajectoryPoint::VFcalEntrance);

    // The track didn't reach hf
    if (!atHF.isValid())
      continue;

    double dHeta = 0.;
    float dHphi = 0.;
    if (dHphi > M_PI)
      dHphi = dHphi - 2. * M_PI;
    else if (dHphi < -M_PI)
      dHphi = dHphi + 2. * M_PI;

    float tracketa = atHF.positionREP().eta() + 0.1 * dHeta;
    float trackphi = atHF.positionREP().phi() + 0.1 * dHphi;

    if (trackphi > M_PI)
      trackphi -= 2 * M_PI;
    else if (trackphi < -M_PI)
      trackphi += 2 * M_PI;

    // Estimate the maximal envelope in phi/eta that will be used to find rechit candidates.
    // Same envelope for cap et barrel rechits.
    double inflation = 1.;
    float rangeeta = (cristalPhiEtaMaxSize_ * (1.5 + 0.5) + 0.2 * fabs(dHeta)) * inflation;
    float rangephi = (cristalPhiEtaMaxSize_ * (1.5 + 0.5) + 0.2 * fabs(dHphi)) * inflation;

    // We search for all candidate recHits, ie all recHits contained in the maximal size envelope.
    std::vector<reco::PFRecHit const*> recHits;
    KDTreeBox trackBox(tracketa - rangeeta, tracketa + rangeeta, trackphi - rangephi, trackphi + rangephi);
    tree_.search(trackBox, recHits);

    // Here we check all rechit candidates using the non-approximated method.
    for (auto const& recHit : recHits) {
      const auto& rhrep = recHit->positionREP();
      const auto& corners = recHit->getCornersREP();

      double rhsizeeta = fabs(corners[3].eta() - corners[1].eta());
      double rhsizephi = fabs(corners[3].phi() - corners[1].phi());
      if (rhsizephi > M_PI)
        rhsizephi = 2. * M_PI - rhsizephi;

      double deta = fabs(rhrep.eta() - tracketa);
      double dphi = fabs(rhrep.phi() - trackphi);
      if (dphi > M_PI)
        dphi = 2. * M_PI - dphi;

      // Find all clusters associated to given rechit
      RecHit2BlockEltMap::iterator ret = rechit2ClusterLinks_.find(recHit);

      for (BlockEltSet::iterator clusterIt = ret->second.begin(); clusterIt != ret->second.end(); clusterIt++) {
        const reco::PFClusterRef clusterref = (*clusterIt)->clusterRef();
        int fracsNbr = clusterref->recHitFractions().size();

        double _rhsizeeta = rhsizeeta * (1.5 + 0.5 / fracsNbr) + 0.2 * fabs(dHeta);
        double _rhsizephi = rhsizephi * (1.5 + 0.5 / fracsNbr) + 0.2 * fabs(dHphi);

        // Check if the track and the cluster are linked
        if (deta < (_rhsizeeta / 2.) && dphi < (_rhsizephi / 2.))
          cluster2TargetLinks_[*clusterIt].insert(*it);
      }
    }
  }
}

void KDTreeLinkerTrackHF::updatePFBlockEltWithLinks() {
  //TODO YG : Check if cluster positionREP() is valid ?

  // Here we save in each HF cluster the list of phi/eta values of linked clusters.
  for (BlockElt2BlockEltMap::iterator it = cluster2TargetLinks_.begin(); it != cluster2TargetLinks_.end(); ++it) {
    reco::PFMultiLinksTC multitracks(true);

    for (BlockEltSet::iterator jt = it->second.begin(); jt != it->second.end(); ++jt) {
      reco::PFRecTrackRef trackref = (*jt)->trackRefPF();
      const reco::PFTrajectoryPoint& atHF = trackref->extrapolatedPoint(reco::PFTrajectoryPoint::VFcalEntrance);
      double tracketa = atHF.positionREP().eta();
      double trackphi = atHF.positionREP().phi();

      multitracks.linkedClusters.push_back(std::make_pair(trackphi, tracketa));
    }

    it->first->setMultilinks(multitracks);
  }

  // We set the multilinks flag of the track to true. It will allow us to
  // use in an optimized way our algo results in the recursive linking algo.
  for (BlockEltSet::iterator it = fieldClusterSet_.begin(); it != fieldClusterSet_.end(); ++it)
    (*it)->setIsValidMultilinks(true);
}

void KDTreeLinkerTrackHF::clear() {
  targetSet_.clear();
  fieldClusterSet_.clear();

  rechitsSet_.clear();

  rechit2ClusterLinks_.clear();
  cluster2TargetLinks_.clear();

  tree_.clear();
}
