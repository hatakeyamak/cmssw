#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoParticleFlow/PFProducer/interface/BlockElementImporterBase.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "RecoParticleFlow/PFProducer/interface/PFMuonAlgo.h"
#include "RecoParticleFlow/PFTracking/interface/PFTrackAlgoTools.h"

class GeneralTracksImporter : public BlockElementImporterBase {
public:
  GeneralTracksImporter(const edm::ParameterSet& conf, edm::ConsumesCollector& sumes)
      : BlockElementImporterBase(conf, sumes),
        src_(sumes.consumes<reco::PFRecTrackCollection>(conf.getParameter<edm::InputTag>("source"))),
        vetoEndcap_(conf.getParameter<bool>("vetoEndcap")),
        muons_(sumes.consumes<reco::MuonCollection>(conf.getParameter<edm::InputTag>("muonSrc"))),
        trackQuality_(reco::TrackBase::qualityByName(conf.getParameter<std::string>("trackQuality"))),
        DPtovPtCut_(conf.getParameter<std::vector<double> >("DPtOverPtCuts_byTrackAlgo")),
        NHitCut_(conf.getParameter<std::vector<unsigned> >("NHitCuts_byTrackAlgo")),
        useIterTracking_(conf.getParameter<bool>("useIterativeTracking")),
        cleanBadConvBrems_(conf.getParameter<bool>("cleanBadConvertedBrems")),
        muonMaxDPtOPt_(conf.getParameter<double>("muonMaxDPtOPt")) {
    if (vetoEndcap_) {
      vetoEndcapSrc_ = sumes.consumes<reco::PFRecTrackCollection>(conf.getParameter<edm::InputTag>("vetoEndcapSource"));
      excludeMuonRefFromVeto_ = conf.getParameter<bool>("excludeMuonRefFromVeto");
    }
  }

  void importToBlock(const edm::Event&, ElementList&) const override;

private:
  const edm::EDGetTokenT<reco::PFRecTrackCollection> src_;
  const bool vetoEndcap_;
  const edm::EDGetTokenT<reco::MuonCollection> muons_;
  const reco::TrackBase::TrackQuality trackQuality_;
  const std::vector<double> DPtovPtCut_;
  const std::vector<unsigned> NHitCut_;
  const bool useIterTracking_, cleanBadConvBrems_;
  const double muonMaxDPtOPt_;
  edm::EDGetTokenT<reco::PFRecTrackCollection> vetoEndcapSrc_;
  bool excludeMuonRefFromVeto_;
};

DEFINE_EDM_PLUGIN(BlockElementImporterFactory, GeneralTracksImporter, "GeneralTracksImporter");

void GeneralTracksImporter::importToBlock(const edm::Event& e, BlockElementImporterBase::ElementList& elems) const {
  typedef BlockElementImporterBase::ElementList::value_type ElementType;
  auto tracks = e.getHandle(src_);
  auto vetosH = e.getHandle(vetoEndcapSrc_);
  const auto& vetos = *vetosH;
  std::unordered_set<unsigned> vetoed;
  if (vetoEndcap_) {
    for (unsigned i = 0; i < vetos.size(); ++i)
      vetoed.insert(vetos[i].trackRef().key());
  }
  auto muons = e.getHandle(muons_);
  elems.reserve(elems.size() + tracks->size());
  std::vector<bool> mask(tracks->size(), true);
  reco::MuonRef muonref;

  // remove converted brems with bad pT resolution if requested
  // this reproduces the old behavior of PFBlockAlgo
  if (cleanBadConvBrems_) {
    auto itr = elems.begin();
    while (itr != elems.end()) {
      if ((*itr)->type() == reco::PFBlockElement::TRACK) {
        const reco::PFBlockElementTrack* trkel = static_cast<reco::PFBlockElementTrack*>(itr->get());
        const reco::ConversionRefVector& cRef = trkel->convRefs();
        const reco::PFDisplacedTrackerVertexRef& dvRef = trkel->displacedVertexRef(reco::PFBlockElement::T_FROM_DISP);
        const reco::VertexCompositeCandidateRef& v0Ref = trkel->V0Ref();
        // if there is no displaced vertex reference  and it is marked
        // as a conversion it's gotta be a converted brem
        if (trkel->trackType(reco::PFBlockElement::T_FROM_GAMMACONV) && cRef.empty() && dvRef.isNull() &&
            v0Ref.isNull()) {
          // if the Pt resolution is bad we kill this element
          if (!PFTrackAlgoTools::goodPtResolution(
                  trkel->trackRef(), DPtovPtCut_, NHitCut_, useIterTracking_, trackQuality_)) {
            itr = elems.erase(itr);
            continue;
          }
        }
      }
      ++itr;
    }  // loop on existing elements
  }
  // preprocess existing tracks in the element list and create a mask
  // so that we do not import tracks twice, tag muons we find
  // in this collection
  auto TKs_end = std::partition(
      elems.begin(), elems.end(), [](const ElementType& a) { return a->type() == reco::PFBlockElement::TRACK; });
  auto btk_elems = elems.begin();
  auto btrack = tracks->cbegin();
  auto etrack = tracks->cend();
  for (auto track = btrack; track != etrack; ++track) {
    auto tk_elem =
        std::find_if(btk_elems, TKs_end, [&](const ElementType& a) { return (a->trackRef() == track->trackRef()); });
    if (tk_elem != TKs_end) {
      mask[std::distance(tracks->cbegin(), track)] = false;
      // check and update if this track is a muon
      const int muId = PFMuonAlgo::muAssocToTrack((*tk_elem)->trackRef(), muons);
      if (muId != -1) {
        muonref = reco::MuonRef(muons, muId);
        if (PFMuonAlgo::isLooseMuon(muonref) || PFMuonAlgo::isMuon(muonref)) {
          static_cast<reco::PFBlockElementTrack*>(tk_elem->get())->setMuonRef(muonref);
        }
      }
    }
  }
  // now we actually insert tracks, again tagging muons along the way
  reco::PFRecTrackRef pftrackref;
  reco::PFBlockElementTrack* trkElem = nullptr;
  for (auto track = btrack; track != etrack; ++track) {
    const unsigned idx = std::distance(btrack, track);
    // since we already set muon refs in the previously imported tracks,
    // here we can skip everything that is already imported
    if (!mask[idx])
      continue;
    muonref = reco::MuonRef();
    pftrackref = reco::PFRecTrackRef(tracks, idx);
    // Get the eventual muon associated to this track
    const int muId = PFMuonAlgo::muAssocToTrack(pftrackref->trackRef(), muons);
    bool thisIsAPotentialMuon = false;
    if (muId != -1) {
      muonref = reco::MuonRef(muons, muId);
      thisIsAPotentialMuon =
          ((PFMuonAlgo::hasValidTrack(muonref, true, muonMaxDPtOPt_) && PFMuonAlgo::isLooseMuon(muonref)) ||
           (PFMuonAlgo::hasValidTrack(muonref, false, muonMaxDPtOPt_) && PFMuonAlgo::isMuon(muonref)));
    }
    if (thisIsAPotentialMuon || PFTrackAlgoTools::goodPtResolution(
                                    pftrackref->trackRef(), DPtovPtCut_, NHitCut_, useIterTracking_, trackQuality_)) {
      trkElem = new reco::PFBlockElementTrack(pftrackref);
      if (thisIsAPotentialMuon) {
        LogDebug("GeneralTracksImporter")
            << "Potential Muon P " << pftrackref->trackRef()->p() << " pt " << pftrackref->trackRef()->p() << std::endl;
      }
      if (muId != -1)
        trkElem->setMuonRef(muonref);
      // if _vetoEndcap is false, add this trk automatically.
      // if _vetoEndcap is true, veto against the hgcal region tracks unless there is muonref is null.
      // when ticl is used, we want to reconstruct tracks including those with muonref in ticl and exclude them here.
      if (!vetoEndcap_ || vetoed.count(pftrackref->trackRef().key()) == 0 ||
          (excludeMuonRefFromVeto_ && muonref.isNonnull())) {
        elems.emplace_back(trkElem);
      } else
        delete trkElem;
    }
  }
  elems.shrink_to_fit();
}
