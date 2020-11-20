// This producer converts a list of TICLCandidates to a list of PFCandidates.

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/View.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "DataFormats/HGCalReco/interface/TICLCandidate.h"

#include "RecoParticleFlow/PFProducer/interface/PFMuonAlgo.h"

class PFTICLProducer : public edm::global::EDProducer<> {
public:
  PFTICLProducer(const edm::ParameterSet&);
  ~PFTICLProducer() override {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

private:
  // inputs
  const edm::EDGetTokenT<edm::View<TICLCandidate>> ticl_candidates_;
  const edm::EDGetTokenT<reco::MuonCollection> muons_;
  // For PFMuonAlgo
  std::unique_ptr<PFMuonAlgo> pfmu_;
};

DEFINE_FWK_MODULE(PFTICLProducer);

PFTICLProducer::PFTICLProducer(const edm::ParameterSet& conf)
    : ticl_candidates_(consumes<edm::View<TICLCandidate>>(conf.getParameter<edm::InputTag>("ticlCandidateSrc"))),
      muons_(consumes<reco::MuonCollection>(conf.getParameter<edm::InputTag>("muonSrc"))) {
  produces<reco::PFCandidateCollection>();
  // For PFMuonAlgo
  const edm::ParameterSet pfMuonAlgoParams = conf.getParameter<edm::ParameterSet>("PFMuonAlgoParameters");
  const bool postMuonCleaning = false;
  pfmu_ = std::make_unique<PFMuonAlgo>(pfMuonAlgoParams, postMuonCleaning);
}

void PFTICLProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("ticlCandidateSrc", edm::InputTag("ticlTrackstersMerge"));
  // For PFMuonAlgo
  desc.add<edm::InputTag>("muonSrc", edm::InputTag("muons1stStep"));
  edm::ParameterSetDescription psd_PFMuonAlgo;
  PFMuonAlgo::fillPSetDescription(psd_PFMuonAlgo);
  desc.add<edm::ParameterSetDescription>("PFMuonAlgoParameters", psd_PFMuonAlgo);
  //
  descriptions.add("pfTICLProducer", desc);
}

void PFTICLProducer::produce(edm::StreamID, edm::Event& evt, const edm::EventSetup& es) const {
  //get TICLCandidates
  edm::Handle<edm::View<TICLCandidate>> ticl_cand_h;
  evt.getByToken(ticl_candidates_, ticl_cand_h);
  const auto ticl_candidates = *ticl_cand_h;
  const auto muons = evt.getHandle(muons_);

  auto candidates = std::make_unique<reco::PFCandidateCollection>();

  for (const auto& ticl_cand : ticl_candidates) {
    const auto abs_pdg_id = std::abs(ticl_cand.pdgId());
    const auto charge = ticl_cand.charge();
    const auto& four_mom = ticl_cand.p4();
    double ecal_energy = 0.;

    for (const auto& t : ticl_cand.tracksters()) {
      double ecal_energy_fraction = t->raw_em_pt() / t->raw_pt();
      ecal_energy += t->raw_energy() * ecal_energy_fraction;
    }
    double hcal_energy = ticl_cand.rawEnergy() - ecal_energy;
    // fix for floating point rounding could go slightly below 0
    hcal_energy = hcal_energy < 0 ? 0 : hcal_energy;

    reco::PFCandidate::ParticleType part_type;
    switch (abs_pdg_id) {
      case 11:
        part_type = reco::PFCandidate::e;
        break;
      case 13:
        part_type = reco::PFCandidate::mu;
        break;
      case 22:
        part_type = reco::PFCandidate::gamma;
        break;
      case 130:
        part_type = reco::PFCandidate::h0;
        break;
      case 211:
        part_type = reco::PFCandidate::h;
        break;
      // default also handles neutral pions (111) for the time being (not yet foreseen in PFCandidate)
      default:
        part_type = reco::PFCandidate::X;
    }

    candidates->emplace_back(charge, four_mom, part_type);

    auto& candidate = candidates->back();
    candidate.setEcalEnergy(ecal_energy, ecal_energy);
    candidate.setHcalEnergy(hcal_energy, hcal_energy);
    if (candidate.charge()) {  // otherwise PFCandidate throws
      // Construct edm::Ref from edm::Ptr. As of now, assumes type to be reco::Track. To be extended (either via
      // dynamic type checking or configuration) if additional track types are needed.
      reco::TrackRef trackref(ticl_cand.trackPtr().id(), int(ticl_cand.trackPtr().key()), &evt.productGetter());
      candidate.setTrackRef(trackref);
      std::cout << "PFTICL trackkey " << ticl_cand.trackPtr().key() << std::endl;
      //
      // Utilize PFMuonAlgo
      const int muId = PFMuonAlgo::muAssocToTrack(trackref, muons);
      if (muId != -1) {
        const reco::MuonRef muonref = reco::MuonRef(muons, muId);
        const bool allowLoose = part_type == reco::PFCandidate::mu ? true : false;
        // Redefine pfmuon candidate kinematics and add muonref
        pfmu_->reconstructMuon(candidate, muonref, allowLoose);
      }
    }

    candidate.setTime(ticl_cand.time(), ticl_cand.timeError());
  }

  evt.put(std::move(candidates));
}
