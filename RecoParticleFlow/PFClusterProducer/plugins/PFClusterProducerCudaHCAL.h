#ifndef __newpf_PFClusterProducerCudaHCAL_H__
#define __newpf_PFClusterProducerCudaHCAL_H__

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "RecoParticleFlow/PFClusterProducer/interface/RecHitTopologicalCleanerBase.h"

#include "RecoParticleFlow/PFClusterProducer/interface/SeedFinderBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/InitialClusteringStepBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFClusterBuilderBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFCPositionCalculatorBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFClusterEnergyCorrectorBase.h"

#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <memory>

class PFClusterProducerCudaHCAL : public edm::stream::EDProducer<> {
  typedef RecHitTopologicalCleanerBase RHCB;
  typedef InitialClusteringStepBase ICSB;
  typedef PFClusterBuilderBase PFCBB;
  typedef PFCPositionCalculatorBase PosCalc;

public:
  PFClusterProducerCudaHCAL(const edm::ParameterSet&);
  ~PFClusterProducerCudaHCAL();

  void beginLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  //void endJob();
  //  void beginStream(edm::StreamID);
 
  
  
  // inputs
  std::vector< std::vector<double> > theThresh;
  edm::EDGetTokenT<reco::PFRecHitCollection> _rechitsLabel;
  edm::EDGetTokenT< std::vector<double> > _trash;
  
  // options
  const bool _prodInitClusters;
  // the actual algorithm
  std::vector<std::unique_ptr<RecHitTopologicalCleanerBase> > _cleaners;
  std::unique_ptr<SeedFinderBase> _seedFinder;
  std::unique_ptr<InitialClusteringStepBase> _initialClustering;
  std::unique_ptr<PFClusterBuilderBase> _pfClusterBuilder;
  std::unique_ptr<PFCPositionCalculatorBase> _positionReCalc;
  std::unique_ptr<PFCPositionCalculatorBase> _allCellsPosCalc;
  std::unique_ptr<PFCPositionCalculatorBase> _allCellsPosCalcCuda;
  std::unique_ptr<PFClusterEnergyCorrectorBase> _energyCorrector;

  TFile *MyFile = new TFile("EventHCAL.root","recreate");
  TH1F *nTopo_CPU = new TH1F("nTopo_CPU","nTopo_CPU",500,0,500);
  TH1F *nTopo_GPU = new TH1F("nTopo_GPU","nTopo_GPU",500,0,500);

  TH1F *sumSeed_CPU = new TH1F("sumSeed_CPU", "sumSeed_CPU",201, -0.5, 200.5);
  TH1F *sumSeed_GPU = new TH1F("sumSeed_GPU", "sumSeed_GPU",201, -0.5, 200.5);

  TH1F *topoEn_CPU = new TH1F("topoEn_CPU", "topoEn_CPU", 500, 0, 500);
  TH1F *topoEn_GPU = new TH1F("topoEn_GPU", "topoEn_GPU", 500, 0, 500);

  TH1F *topoEta_CPU = new TH1F("topoEta_CPU", "topoEta_CPU", 100, -3, 3);
  TH1F *topoEta_GPU = new TH1F("topoEta_GPU", "topoEta_GPU", 100, -3, 3);

  TH1F *topoPhi_CPU = new TH1F("topoPhi_CPU", "topoPhi_CPU", 100, -3.1415926, 3.1415926);
  TH1F *topoPhi_GPU = new TH1F("topoPhi_GPU", "topoPhi_GPU", 100, -3.1415926, 3.1415926);

  TH1F *nPFCluster_CPU = new TH1F("nPFCluster_CPU","nPFCluster_CPU",1000,0,1000);
  TH1F *nPFCluster_GPU = new TH1F("nPFCluster_GPU","nPFCluster_GPU",1000,0,1000);

  TH1F *enPFCluster_CPU = new TH1F("enPFCluster_CPU","enPFCluster_CPU",500,0,500);
  TH1F *enPFCluster_GPU = new TH1F("enPFCluster_GPU","enPFCluster_GPU",500,0,500);

  TH1F *nRH_perPFCluster_CPU = new TH1F("nRH_perPFCluster_CPU","nRH_perPFCluster_CPU",100,0,100);
  TH1F *nRH_perPFCluster_GPU = new TH1F("nRH_perPFCluster_GPU","nRH_perPFCluster_GPU",100,0,100);

  TH2F *nRh_CPUvsGPU = new TH2F("nRh_CPUvsGPU","nRh_CPUvsGPU",100,0,100,100,0,100);
  TH2F *enPFCluster_CPUvsGPU = new TH2F("enPFCluster_CPUvsGPU","enPFCluster_CPUvsGPU",5000,0,500,5000,0,500);
  TH1F *enPFCluster_CPUvsGPU_1d = new TH1F("enPFCluster_CPUvsGPU_1d","enPFCluster_CPUvsGPU_1d",400,-2,2);

  bool doComparison=true;

  TH1F *deltaSumSeed  = new TH1F("deltaSumSeed", "sumSeed_{GPU} - sumSeed_{CPU}", 201, -100.5, 100.5);
  TH1F *deltaRH  = new TH1F("deltaRH", "nRH_{GPU} - nRH_{CPU}", 41, -20.5, 20.5);
  TH1F *deltaEn  = new TH1F("deltaEn", "E_{GPU} - E_{CPU}", 200, -10, 10);
  TH1F *deltaEta = new TH1F("deltaEta", "#eta_{GPU} - #eta_{CPU}", 200, -0.2, 0.2);
  TH1F *deltaPhi = new TH1F("deltaPhi", "#phi_{GPU} - #phi_{CPU}", 200, -0.2, 0.2);

  TH2F *coordinate = new TH2F("coordinate","coordinate",100,-3,3,100,-3.1415926,3.14159);
  TH1F *layer = new TH1F("layer","layer",7,0,7);
};

DEFINE_FWK_MODULE(PFClusterProducerCudaHCAL);

#endif
