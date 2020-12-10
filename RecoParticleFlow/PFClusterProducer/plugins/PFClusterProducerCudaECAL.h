#ifndef __newpf_PFClusterProducerCudaECAL_H__
#define __newpf_PFClusterProducerCudaECAL_H__

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

class PFClusterProducerCudaECAL : public edm::stream::EDProducer<> {
  typedef RecHitTopologicalCleanerBase RHCB;
  typedef InitialClusteringStepBase ICSB;
  typedef PFClusterBuilderBase PFCBB;
  typedef PFCPositionCalculatorBase PosCalc;

public:
  PFClusterProducerCudaECAL(const edm::ParameterSet&);
  ~PFClusterProducerCudaECAL();

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
  std::unique_ptr<PFClusterEnergyCorrectorBase> _energyCorrector;

TFile *MyFile = new TFile("EventECAL.root","recreate");
  TH1F *nTopo_CPU = new TH1F("nTopo_CPU","nTopo_CPU",500,0,500);
  TH1F *nTopo_GPU = new TH1F("nTopo_GPU","nTopo_GPU",500,0,500);

  TH1F *nPFCluster_CPU = new TH1F("nPFCluster_CPU","nPFCluster_CPU",500,0,500);
  TH1F *nPFCluster_GPU = new TH1F("nPFCluster_GPU","nPFCluster_GPU",500,0,500);

  TH1F *enPFCluster_CPU = new TH1F("enPFCluster_CPU","enPFCluster_CPU",500,0,500);
  TH1F *enPFCluster_GPU = new TH1F("enPFCluster_GPU","enPFCluster_GPU",500,0,500);

  TH1F *nRH_perPFCluster_CPU = new TH1F("nRH_perPFCluster_CPU","nRH_perPFCluster_CPU",100,0,100);
  TH1F *nRH_perPFCluster_GPU = new TH1F("nRH_perPFCluster_GPU","nRH_perPFCluster_GPU",100,0,100);

  TH2F *nRh_CPUvsGPU = new TH2F("nRh_CPUvsGPU","nRh_CPUvsGPU",100,0,100,100,0,100);
  TH2F *enPFCluster_CPUvsGPU = new TH2F("enPFCluster_CPUvsGPU","enPFCluster_CPUvsGPU",500,0,500,500,0,500);

  bool doComparison=false;

};

DEFINE_FWK_MODULE(PFClusterProducerCudaECAL);

#endif
