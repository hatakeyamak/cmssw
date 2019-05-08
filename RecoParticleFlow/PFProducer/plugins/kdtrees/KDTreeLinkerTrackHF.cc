#include "KDTreeLinkerTrackHF.h"

#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "TMath.h"

// the text name is different so that we can easily
// construct it when calling the factory
DEFINE_EDM_PLUGIN(KDTreeLinkerFactory, 
		  KDTreeLinkerTrackHF, 
		  "KDTreeTrackAndHFLinker"); 


KDTreeLinkerTrackHF::KDTreeLinkerTrackHF()
  : KDTreeLinkerBase()
{
  cristalPhiEtaMaxSize_ = 0.2;
  phiOffset_ = 0.32;
}

KDTreeLinkerTrackHF::~KDTreeLinkerTrackHF()
{
  clear();
}

void
KDTreeLinkerTrackHF::insertTargetElt(reco::PFBlockElement	*track)
{
  //KH
  std::cout << "KDTreeLinkerTrackHF::insertTargetElt " << track->trackRef()->pt() << " " << track->trackRef()->eta() << std::endl;
  if (fabs(track->trackRef()->eta())<3.0) return;
  
  if( track->trackRefPF()->extrapolatedPoint( reco::PFTrajectoryPoint::VFcalEntrance ).isValid() ) {
    std::cout << "KDTreeLinkerTrackHF::insertTargetElt - isValid" << std::endl;
    targetSet_.insert(track);
  }

  std::cout << "KDTreeLinkerTrackHF::insertTargetElt ends" << std::endl;

}


void
KDTreeLinkerTrackHF::insertFieldClusterElt(reco::PFBlockElement	*hfCluster)
{
  const reco::PFClusterRef& clusterref = hfCluster->clusterRef();

  //KH
  std::cout << "KDTreeLinkerTrackHF::insertFieldClusterElt "
	    << clusterref->layer() << " "
	    << std::endl;
  
  // This test is more or less done in PFBlockAlgo.h. In others cases, it should be switch on.
  // KH: Keeping it on to be sure
  if (!((clusterref->layer() == PFLayer::HF_HAD) ||
   	(clusterref->layer() == PFLayer::HF_EM)))
    return;

  const std::vector<reco::PFRecHitFraction> &fraction = clusterref->recHitFractions();

  std::cout << "fraction.size(): " << fraction.size() << std::endl;
  
  // We create a list of hfCluster
  fieldClusterSet_.insert(hfCluster);
  for(size_t rhit = 0; rhit < fraction.size(); ++rhit) {
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

void 
KDTreeLinkerTrackHF::buildTree()
{
  // List of pseudo-rechits that will be used to create the KDTree
  std::vector<KDTreeNodeInfo> eltList;

  // Filling of this list
  for(RecHitSet::const_iterator it = rechitsSet_.begin(); 
      it != rechitsSet_.end(); it++) {
    
    const reco::PFRecHit::REPPoint &posrep = (*it)->positionREP();
    
    KDTreeNodeInfo rh1 (*it, posrep.eta(), posrep.phi());
    eltList.push_back(rh1);
    
    // Here we solve the problem of phi circular set by duplicating some rechits
    // too close to -Pi (or to Pi) and adding (substracting) to them 2 * Pi.
    if (rh1.dim2 > (M_PI - phiOffset_)) {
      double phi = rh1.dim2 - 2 * M_PI;
      KDTreeNodeInfo rh2(*it, posrep.eta(), phi); 
      eltList.push_back(rh2);
    }

    if (rh1.dim2 < (M_PI * -1.0 + phiOffset_)) {
      double phi = rh1.dim2 + 2 * M_PI;
      KDTreeNodeInfo rh3(*it, posrep.eta(), phi); 
      eltList.push_back(rh3);
    }
  }

  // Here we define the upper/lower bounds of the 2D space (eta/phi).
  double phimin = -1.0 * M_PI - phiOffset_;
  double phimax = M_PI + phiOffset_;

  // etamin-etamax, phimin-phimax
  KDTreeBox region(-3.0, 3.0, phimin, phimax);

  // We may now build the KDTree
  tree_.build(eltList, region);
}

void
KDTreeLinkerTrackHF::searchLinks()
{
  // Must of the code has been taken from LinkByRecHit.cc

  //KH
  std::cout << "KDTreeLinkerTrackHF::searchLinks" << std::endl;
  
  // We iterate over the tracks.
  for(BlockEltSet::iterator it = targetSet_.begin(); 
      it != targetSet_.end(); it++) {
	
    reco::PFRecTrackRef trackref = (*it)->trackRefPF();

    const reco::PFTrajectoryPoint& atHF = 
      trackref->extrapolatedPoint(reco::PFTrajectoryPoint::VFcalEntrance);
    
    // The track didn't reach hf??
    if( ! atHF.isValid()) continue;
    
    double dHeta = 0.;
    double dHphi = 0.;
    if ( dHphi > M_PI ) dHphi = dHphi - 2. * M_PI;
    else if ( dHphi < -M_PI ) dHphi = dHphi + 2. * M_PI; 
    
    double tracketa = atHF.positionREP().eta() + 0.1 * dHeta;
    double trackphi = atHF.positionREP().phi() + 0.1 * dHphi;
    
    if (trackphi > M_PI) trackphi -= 2 * M_PI;
    else if (trackphi < -M_PI) trackphi += 2 * M_PI;

    // Estimate the maximal envelope in phi/eta that will be used to find rechit candidates.
    // Same envelope for cap et barrel rechits.
    double inflation = 1.;
    double rangeeta = (cristalPhiEtaMaxSize_ * (1.5 + 0.5) + 0.2 * fabs(dHeta)) * inflation; 
    double rangephi = (cristalPhiEtaMaxSize_ * (1.5 + 0.5) + 0.2 * fabs(dHphi)) * inflation; 

    // We search for all candidate recHits, ie all recHits contained in the maximal size envelope.
    std::vector<KDTreeNodeInfo> recHits;
    KDTreeBox trackBox(tracketa - rangeeta, tracketa + rangeeta, 
		       trackphi - rangephi, trackphi + rangephi);
    tree_.search(trackBox, recHits);

    int aaa=0;
    
    // Here we check all rechit candidates using the non-approximated method.
    for(std::vector<KDTreeNodeInfo>::const_iterator rhit = recHits.begin(); 
	rhit != recHits.end(); ++rhit) {

      aaa++;
      
      const auto &rhrep		   = rhit->ptr->positionREP();
      const auto & corners = rhit->ptr->getCornersREP();
      
      double rhsizeeta = fabs(corners[3].eta() - corners[1].eta());
      double rhsizephi = fabs(corners[3].phi() - corners[1].phi());
      if ( rhsizephi > M_PI ) rhsizephi = 2.*M_PI - rhsizephi;
      
      double deta = fabs(rhrep.eta() - tracketa);
      double dphi = fabs(rhrep.phi() - trackphi);
      if ( dphi > M_PI ) dphi = 2.*M_PI - dphi;
      
      // Find all clusters associated to given rechit
      RecHit2BlockEltMap::iterator ret = rechit2ClusterLinks_.find(rhit->ptr);
      
      for(BlockEltSet::iterator clusterIt = ret->second.begin(); 
	  clusterIt != ret->second.end(); clusterIt++) {
	
	const reco::PFClusterRef clusterref = (*clusterIt)->clusterRef();
	int fracsNbr = clusterref->recHitFractions().size();
	
	double _rhsizeeta = rhsizeeta * (1.5 + 0.5 / fracsNbr) + 0.2 * fabs(dHeta);
	double _rhsizephi = rhsizephi * (1.5 + 0.5 / fracsNbr) + 0.2 * fabs(dHphi);
	
	// Check if the track and the cluster are linked
	if(deta < (_rhsizeeta / 2.) && dphi < (_rhsizephi / 2.)){
	  std::cout << "KDTreeLinkerTrackHF::searchLinks" << " linked! "
		    << tracketa << " " << trackphi << " " << aaa << " "
		    << rhrep.eta() << " " << rhrep.phi() << std::endl;
	  cluster2TargetLinks_[*clusterIt].insert(*it);
	}
      }
    }
  }
}

void
KDTreeLinkerTrackHF::updatePFBlockEltWithLinks()
{
  //TODO YG : Check if cluster positionREP() is valid ?

  //KH
  std::cout << "KDTreeLinkerTrackHF::updatePFBlockEltWithLinks" << std::endl;
  
  // Here we save in each HF cluster the list of phi/eta values of linked clusters.
  for (BlockElt2BlockEltMap::iterator it = cluster2TargetLinks_.begin();
       it != cluster2TargetLinks_.end(); ++it) {
    reco::PFMultiLinksTC multitracks(true);

    for (BlockEltSet::iterator jt = it->second.begin();
	 jt != it->second.end(); ++jt) {

      reco::PFRecTrackRef trackref = (*jt)->trackRefPF();
      const reco::PFTrajectoryPoint& atHF = 
	trackref->extrapolatedPoint(reco::PFTrajectoryPoint::VFcalEntrance);
      double tracketa = atHF.positionREP().eta();
      double trackphi = atHF.positionREP().phi();
      
      multitracks.linkedClusters.push_back(std::make_pair(trackphi, tracketa));
    }

    it->first->setMultilinks(multitracks);
  }

  // We set the multilinks flag of the track to true. It will allow us to 
  // use in an optimized way our algo results in the recursive linking algo.
  for (BlockEltSet::iterator it = fieldClusterSet_.begin();
       it != fieldClusterSet_.end(); ++it)
    (*it)->setIsValidMultilinks(true);

}

void
KDTreeLinkerTrackHF::clear()
{
  targetSet_.clear();
  fieldClusterSet_.clear();

  rechitsSet_.clear();

  rechit2ClusterLinks_.clear();
  cluster2TargetLinks_.clear();

  tree_.clear();
}
