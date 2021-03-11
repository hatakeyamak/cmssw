#include "PFClusterProducerCudaECAL.h"
#include "RecoParticleFlow/PFClusterProducer/plugins/PFClusterCudaECAL.h"
#include <TFile.h>
#include <TH1F.h>
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "HeterogeneousCore/CUDAUtilities/interface/device_unique_ptr.h"
#include "HeterogeneousCore/CUDAUtilities/interface/host_unique_ptr.h"


#ifdef PFLOW_DEBUG
#define LOGVERB(x) edm::LogVerbatim(x)
#define LOGWARN(x) edm::LogWarning(x)
#define LOGERR(x) edm::LogError(x)
#define LOGDRESSED(x) edm::LogInfo(x)
#else
#define LOGVERB(x) LogTrace(x)
#define LOGWARN(x) edm::LogWarning(x)
#define LOGERR(x) edm::LogError(x)
#define LOGDRESSED(x) LogDebug(x)
#endif

PFClusterProducerCudaECAL::PFClusterProducerCudaECAL(const edm::ParameterSet& conf)
  : 

  _prodInitClusters(conf.getUntrackedParameter<bool>("prodInitialClusters", false)) {
  _rechitsLabel = consumes<reco::PFRecHitCollection>(conf.getParameter<edm::InputTag>("recHitsSource"));

  //setup TTree
  clusterTree->Branch("initialClusters", "PFClusterCollection", &__initialClusters);
  clusterTree->Branch("pfClusters", "PFClusterCollection", &__pfClusters);
  clusterTree->Branch("pfClustersFromCuda", "PFClusterCollection", &__pfClustersFromCuda);

  //setup rechit cleaners
  const edm::VParameterSet& cleanerConfs = conf.getParameterSetVector("recHitCleaners");

  for (const auto& conf : cleanerConfs) {
    const std::string& cleanerName = conf.getParameter<std::string>("algoName");
    _cleaners.emplace_back(RecHitTopologicalCleanerFactory::get()->create(cleanerName, conf));
  }

  edm::ConsumesCollector sumes = consumesCollector();

  // setup seed finding
  const edm::ParameterSet& sfConf = conf.getParameterSet("seedFinder");
  const std::string& sfName = sfConf.getParameter<std::string>("algoName");
  _seedFinder = SeedFinderFactory::get()->create(sfName, sfConf);

  const edm::VParameterSet& seedFinderConfs = sfConf.getParameterSetVector("thresholdsByDetector"); 
  

  if (conf.exists("allCellsPositionCalc")) {
    const edm::ParameterSet& acConf = conf.getParameterSet("allCellsPositionCalc");
    const std::string& algoac = acConf.getParameter<std::string>("algoName");
    _allCellsPosCalc = PFCPositionCalculatorFactory::get()->create(algoac, acConf);
  }
 
  //setup topo cluster builder
  const edm::ParameterSet& initConf = conf.getParameterSet("initialClusteringStep");
  const std::string& initName = initConf.getParameter<std::string>("algoName");
  _initialClustering = InitialClusteringStepFactory::get()->create(initName, initConf, sumes);
  //setup pf cluster builder if requested
  const edm::ParameterSet& pfcConf = conf.getParameterSet("pfClusterBuilder");
  if (!pfcConf.empty()) {
    const std::string& pfcName = pfcConf.getParameter<std::string>("algoName");
    _pfClusterBuilder = PFClusterBuilderFactory::get()->create(pfcName, pfcConf);
  }
  //setup (possible) recalcuation of positions
  const edm::ParameterSet& pConf = conf.getParameterSet("positionReCalc");
  if (!pConf.empty()) {
    const std::string& pName = pConf.getParameter<std::string>("algoName");
    _positionReCalc = PFCPositionCalculatorFactory::get()->create(pName, pConf);
  }
  // see if new need to apply corrections, setup if there.
  const edm::ParameterSet& cConf = conf.getParameterSet("energyCorrector");
  if (!cConf.empty()) {
    const std::string& cName = cConf.getParameter<std::string>("algoName");
    _energyCorrector = PFClusterEnergyCorrectorFactory::get()->create(cName, cConf);
  }
  
  //bool doComparison = false;

  





  if (_prodInitClusters) {
    produces<reco::PFClusterCollection>("initialClusters");
  }
  produces<reco::PFClusterCollection>();
 
  
}

PFClusterProducerCudaECAL::~PFClusterProducerCudaECAL()
{
  MyFile->cd();
  clusterTree->Write();
  nTopo_CPU->Write();
  nTopo_GPU->Write();
  sumSeed_CPU->Write();
  sumSeed_GPU->Write();
  topoEn_CPU->Write();
  topoEn_GPU->Write();
  topoEta_CPU->Write();
  topoEta_GPU->Write();
  topoPhi_CPU->Write();
  topoPhi_GPU->Write();
  nPFCluster_CPU->Write();
  nPFCluster_GPU->Write();
  enPFCluster_CPU->Write();
  enPFCluster_GPU->Write();
  pfcEta_CPU->Write();
  pfcEta_GPU->Write();
  pfcPhi_CPU->Write();
  pfcPhi_GPU->Write();
  nRH_perPFCluster_CPU->Write();
  nRH_perPFCluster_GPU->Write();
  matched_pfcRh_CPU->Write();
  matched_pfcRh_GPU->Write();
  matched_pfcEn_CPU->Write();
  matched_pfcEn_GPU->Write();
  matched_pfcEta_CPU->Write();
  matched_pfcEta_GPU->Write();
  matched_pfcPhi_CPU->Write();
  matched_pfcPhi_GPU->Write();
  nRh_CPUvsGPU->Write();
  enPFCluster_CPUvsGPU->Write();
  coordinate->Write();
  layer->Write();
  deltaSumSeed->Write();
  deltaRH->Write();
  deltaEn->Write();
  deltaEta->Write();
  deltaPhi->Write();

  // MyFile->Close();
  delete MyFile;
}

void PFClusterProducerCudaECAL::beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup& es) {
  _initialClustering->update(es);
  if (_pfClusterBuilder)
    _pfClusterBuilder->update(es);
  if (_positionReCalc)
    _positionReCalc->update(es);
}

void PFClusterProducerCudaECAL::produce(edm::Event& e, const edm::EventSetup& es) {
  _initialClustering->reset();
  if (_pfClusterBuilder)
    _pfClusterBuilder->reset();

  edm::Handle<reco::PFRecHitCollection> rechits;
  e.getByToken(_rechitsLabel, rechits);
  
  _initialClustering->updateEvent(e);

  std::vector<bool> mask(rechits->size(), true);
  for (const auto& cleaner : _cleaners) {
    cleaner->clean(rechits, mask);
  }

  size_t rh_size = 2000;
  int numbytes_float = rh_size*sizeof(float);
  int numbytes_int = rh_size*sizeof(int);

  std::vector<float>                                    h_cuda_pfRhFrac=std::vector<float>(rh_size*50,-1.);
  std::vector<float>                                    h_cuda_pcRhFrac=std::vector<float>(rh_size*50,-1.);
  std::vector<int>                                      h_cuda_pfRhFracInd=std::vector<int>(rh_size*50,-1);
  std::vector<int>                                      h_cuda_pfNeighEightInd=std::vector<int>(rh_size*8,0);
  std::vector<int>                                      h_cuda_pcRhFracInd=std::vector<int>(rh_size*50,-1);
  std::vector<float>                                    h_cuda_fracsum=std::vector<float>(rh_size,0);
  std::vector<int>                                      h_cuda_rhcount=std::vector<int>(rh_size,1);
  std::vector<float>                                    h_cuda_pfrh_x=std::vector<float>(rh_size,0);
  std::vector<float>                                    h_cuda_pfrh_y=std::vector<float>(rh_size,0);
  std::vector<float>                                    h_cuda_pfrh_z=std::vector<float>(rh_size,0);
  std::vector<float>                                    h_cuda_pfrh_energy=std::vector<float>(rh_size,0);
  std::vector<float>                                    h_cuda_pfrh_pt2=std::vector<float>(rh_size,0);
  std::vector<int>                                      h_cuda_pfrh_topoId=std::vector<int>(rh_size,-1);
  std::vector<int>                                      h_cuda_pfrh_isSeed=std::vector<int>(rh_size,0);
  std::vector<int>                                      h_cuda_pfrh_layer=std::vector<int>(rh_size,-999);

auto d_cuda_pfrh_x = cms::cuda::make_device_unique<float[]>(numbytes_float, nullptr);
auto d_cuda_pfrh_y = cms::cuda::make_device_unique<float[]>(numbytes_float, nullptr);
auto d_cuda_pfrh_z = cms::cuda::make_device_unique<float[]>(numbytes_float, nullptr);
auto d_cuda_pfrh_energy = cms::cuda::make_device_unique<float[]>(numbytes_float, nullptr);
auto d_cuda_pfrh_pt2 = cms::cuda::make_device_unique<float[]>(numbytes_float, nullptr);
auto d_cuda_pfRhFrac = cms::cuda::make_device_unique<float[]>(numbytes_float*50, nullptr);
auto d_cuda_pcRhFrac = cms::cuda::make_device_unique<float[]>(numbytes_float*50, nullptr);
auto d_cuda_rhcount = cms::cuda::make_device_unique<int[]>(numbytes_int, nullptr);
auto d_cuda_fracsum = cms::cuda::make_device_unique<float[]>(numbytes_float, nullptr);
auto d_cuda_pfrh_topoId = cms::cuda::make_device_unique<int[]>(numbytes_int, nullptr);
auto d_cuda_pfrh_isSeed = cms::cuda::make_device_unique<int[]>(numbytes_int, nullptr);
auto d_cuda_pfrh_layer = cms::cuda::make_device_unique<int[]>(numbytes_int, nullptr);
auto d_cuda_pfNeighEightInd = cms::cuda::make_device_unique<int[]>(numbytes_int*8, nullptr);
auto d_cuda_pfRhFracInd = cms::cuda::make_device_unique<int[]>(numbytes_int*50, nullptr);
auto d_cuda_pcRhFracInd = cms::cuda::make_device_unique<int[]>(numbytes_int*50, nullptr);

  int p=0; 
  for (auto rh: *rechits){
    h_cuda_pfrh_x[p]=rh.position().x();
    h_cuda_pfrh_y[p]=rh.position().y();
    h_cuda_pfrh_z[p]=rh.position().z();
    h_cuda_pfrh_energy[p]=rh.energy();
    h_cuda_pfrh_pt2[p]=rh.pt2();
    h_cuda_pfrh_layer[p]=(int)rh.layer();
    h_cuda_pfrh_topoId[p]=p;
    auto theneighboursEight = rh.neighbours8();
    int z = 0;
    for(auto nh: theneighboursEight)
      {
	h_cuda_pfNeighEightInd[8*p+z] = nh;
	z++;
      }
    for(int l=z; l<8; l++)
      {
	h_cuda_pfNeighEightInd[8*p+l] = -1;
      }
    p++;
  }//end of rechit loop  
  //std::cout<<"p: "<<p<<std::endl;
  
  cudaCheck(cudaMemcpy(d_cuda_fracsum.get(), h_cuda_fracsum.data(), numbytes_float, cudaMemcpyHostToDevice));
  cudaCheck(cudaMemcpy(d_cuda_rhcount.get(), h_cuda_rhcount.data(), numbytes_int, cudaMemcpyHostToDevice));
  cudaCheck(cudaMemcpy(d_cuda_pfrh_x.get(), h_cuda_pfrh_x.data(), numbytes_float, cudaMemcpyHostToDevice));
  cudaCheck(cudaMemcpy(d_cuda_pfrh_y.get(), h_cuda_pfrh_y.data(), numbytes_float, cudaMemcpyHostToDevice));
  cudaCheck(cudaMemcpy(d_cuda_pfrh_z.get(), h_cuda_pfrh_z.data(), numbytes_float, cudaMemcpyHostToDevice));
  cudaCheck(cudaMemcpy(d_cuda_pfrh_energy.get(), h_cuda_pfrh_energy.data(), numbytes_float, cudaMemcpyHostToDevice));
  cudaCheck(cudaMemcpy(d_cuda_pfrh_pt2.get(), h_cuda_pfrh_pt2.data(), numbytes_float, cudaMemcpyHostToDevice));
  cudaCheck(cudaMemcpy(d_cuda_pfrh_topoId.get(), h_cuda_pfrh_topoId.data(), numbytes_int, cudaMemcpyHostToDevice));
  cudaCheck(cudaMemcpy(d_cuda_pfrh_isSeed.get(), h_cuda_pfrh_isSeed.data(), numbytes_int, cudaMemcpyHostToDevice));
  cudaCheck(cudaMemcpy(d_cuda_pfrh_layer.get(), h_cuda_pfrh_layer.data(), numbytes_int, cudaMemcpyHostToDevice));
  cudaCheck(cudaMemcpy(d_cuda_pfNeighEightInd.get(), h_cuda_pfNeighEightInd.data(), numbytes_int*8, cudaMemcpyHostToDevice));  
  cudaCheck(cudaMemcpy(d_cuda_pfRhFrac.get(), h_cuda_pfRhFrac.data(), numbytes_float*50, cudaMemcpyHostToDevice));
  cudaCheck(cudaMemcpy(d_cuda_pcRhFrac.get(), h_cuda_pcRhFrac.data(), numbytes_float*50, cudaMemcpyHostToDevice));
  cudaCheck(cudaMemcpy(d_cuda_pfRhFracInd.get(), h_cuda_pfRhFracInd.data(), numbytes_int*50, cudaMemcpyHostToDevice));
  cudaCheck(cudaMemcpy(d_cuda_pcRhFracInd.get(), h_cuda_pcRhFracInd.data(), numbytes_int*50, cudaMemcpyHostToDevice));
 
  
  //PFClusterCudaECAL::PFRechitToPFCluster_ECAL_serialize(rechits->size(), 
  //PFClusterCudaECAL::PFRechitToPFCluster_ECAL_serialize_topoParallel(rechits->size(), 
  //PFClusterCudaECAL::PFRechitToPFCluster_ECAL_serialize_seedingParallel(rechits->size(), 
  //PFClusterCudaECAL::PFRechitToPFCluster_ECAL_serialize_step1Parallel(rechits->size(), 
  //PFClusterCudaECAL::PFRechitToPFCluster_ECAL_serialize_step2Parallel(rechits->size(), 
  PFClusterCudaECAL::PFRechitToPFCluster_ECALV2(rechits->size(), 
					      d_cuda_pfrh_x.get(),  
					      d_cuda_pfrh_y.get(),  
					      d_cuda_pfrh_z.get(), 
					      d_cuda_pfrh_energy.get(), 
					      d_cuda_pfrh_pt2.get(), 	
					      d_cuda_pfrh_isSeed.get(),
					      d_cuda_pfrh_topoId.get(),
					      d_cuda_pfrh_layer.get(), 
					      d_cuda_pfNeighEightInd.get(), 
					      d_cuda_pfRhFrac.get(), 
					      d_cuda_pfRhFracInd.get(), 
					      d_cuda_pcRhFracInd.get(),
					      d_cuda_pcRhFrac.get(),
					      d_cuda_fracsum.get(),
					      d_cuda_rhcount.get()
					      );

  /*
  PFClusterCudaECAL::PFRechitToPFCluster_ECALV2(rechits->size(), 
					      d_cuda_pfrh_x.get(),  
					      d_cuda_pfrh_y.get(),  
					      d_cuda_pfrh_z.get(), 
					      d_cuda_pfrh_energy.get(), 
					      d_cuda_pfrh_pt2.get(), 	
					      d_cuda_pfrh_isSeed.get(),
					      d_cuda_pfrh_topoId.get(),
					      d_cuda_pfrh_layer.get(), 
					      d_cuda_pfNeighEightInd.get(), 
					      d_cuda_pfRhFrac.get(), 
					      d_cuda_pfRhFracInd.get(), 
					      d_cuda_pcRhFracInd.get(),
					      d_cuda_pcRhFrac.get(),
					      d_cuda_fracsum.get(),
					      d_cuda_rhcount.get()
					      );
  
  
  PFClusterCudaECAL::PFRechitToPFCluster_ECALV1(rh_size, 
					      d_cuda_pfrh_x.get(),  
					      d_cuda_pfrh_y.get(),  
					      d_cuda_pfrh_z.get(), 
					      d_cuda_pfrh_energy.get(), 
					      d_cuda_pfrh_pt2.get(), 	
					      d_cuda_pfrh_isSeed.get(),
					      d_cuda_pfrh_topoId.get(),
					      d_cuda_pfrh_layer.get(), 
					      d_cuda_pfNeighEightInd.get(), 
					      d_cuda_pfRhFrac.get(), 
					      d_cuda_pfRhFracInd.get(), 
					      d_cuda_pcRhFracInd.get(),
					      d_cuda_pcRhFrac.get()
					      );
  */
     
    
  cudaMemcpy(h_cuda_pcRhFracInd.data()    , d_cuda_pcRhFracInd.get()  , numbytes_int*50 , cudaMemcpyDeviceToHost);  
  cudaMemcpy(h_cuda_pcRhFrac.data()       , d_cuda_pcRhFrac.get()  , numbytes_float*50 , cudaMemcpyDeviceToHost);  
  cudaMemcpy(h_cuda_pfrh_isSeed.data()    , d_cuda_pfrh_isSeed.get()  , numbytes_int , cudaMemcpyDeviceToHost);  
  
  auto pfClustersFromCuda = std::make_unique<reco::PFClusterCollection>();
  pfClustersFromCuda.reset(new reco::PFClusterCollection);
  for(int n=0; n<(int)rh_size; n++){
    if(h_cuda_pfrh_isSeed[n]==1){
      reco::PFCluster temp;
      temp.setSeed((*rechits)[n].detId()); 
      for(int k=0;k<50;k++){
	if(h_cuda_pcRhFracInd[n*50+k] > -1){
	  const reco::PFRecHitRef& refhit = reco::PFRecHitRef(rechits,h_cuda_pcRhFracInd[n*50+k]);
	  temp.addRecHitFraction( reco::PFRecHitFraction(refhit, h_cuda_pcRhFrac[n*50+k]) );
	}
	if(h_cuda_pcRhFracInd[n*50+k] < 0.) break;
      }    
      pfClustersFromCuda->push_back(temp);
    }   
  }
  _positionReCalc->calculateAndSetPositions(*pfClustersFromCuda);



  if(doComparison)
  {
    std::vector<bool> seedable(rechits->size(), false);
    _seedFinder->findSeeds(rechits, mask, seedable);
    auto initialClusters = std::make_unique<reco::PFClusterCollection>();
    _initialClustering->buildClusters(rechits, mask, seedable, *initialClusters);
    LOGVERB("PFClusterProducer::produce()") << *_initialClustering;
    __initialClusters = *initialClusters;  // For TTree

    int topoRhCount=0;
    for(auto pfc : *initialClusters)
      {
        nTopo_CPU->Fill(pfc.recHitFractions().size());
        topoEn_CPU->Fill(pfc.energy());
        topoEta_CPU->Fill(pfc.eta());
        topoPhi_CPU->Fill(pfc.phi());
        
        topoRhCount=topoRhCount+pfc.recHitFractions().size();
      }

    nPFCluster_CPU->Fill(initialClusters->size());
    std::sort (h_cuda_pfrh_topoId.begin(), h_cuda_pfrh_topoId.end());
    
    int topoCount=1;
    int intTopoCount=0;
    for(int l=1; l<(int)h_cuda_pfrh_topoId.size();l++){
      if((h_cuda_pfrh_topoId[l]==h_cuda_pfrh_topoId[l+1]) && h_cuda_pfrh_topoId[l]>-1.) topoCount++;
      else if(h_cuda_pfrh_topoId[l]>-1.){
    nTopo_GPU->Fill(topoCount);
    topoCount=1;
    intTopoCount++;
      }
    }
    /*
    std::cout<<"ECAL:"<<std::endl;
    std::cout<<"sum rechits          : "<<rh_size<<std::endl;
    std::cout<<"sum rechits in topo  : "<<topoRhCount<<std::endl;
    */
    nPFCluster_GPU->Fill(intTopoCount);

    int seedSumCPU=0;
    int seedSumGPU=0;
    int maskSize = 0;
    for (int j=0;j<(int)seedable.size(); j++) seedSumCPU=seedSumCPU+seedable[j];
    for (int j=0;j<(int)h_cuda_pfrh_isSeed.size(); j++) seedSumGPU=seedSumGPU +h_cuda_pfrh_isSeed[j];
    for (int j=0;j<(int)mask.size(); j++) maskSize=maskSize +mask[j];

    sumSeed_CPU->Fill(seedSumCPU);
    sumSeed_GPU->Fill(seedSumGPU);
    deltaSumSeed->Fill(seedSumGPU - seedSumCPU);

    auto pfClusters = std::make_unique<reco::PFClusterCollection>();
    pfClusters.reset(new reco::PFClusterCollection);
    if (_pfClusterBuilder) {  // if we've defined a re-clustering step execute it
      _pfClusterBuilder->buildClusters(*initialClusters, seedable, *pfClusters);
    LOGVERB("PFClusterProducer::produce()") << *_pfClusterBuilder;
    } else {
      pfClusters->insert(pfClusters->end(), initialClusters->begin(), initialClusters->end());
    }
    
    __pfClusters = *pfClusters;  // For TTree  
    for(auto pfc : *pfClusters)
    {
      nRH_perPFCluster_CPU->Fill(pfc.recHitFractions().size());
	  enPFCluster_CPU->Fill(pfc.energy());
      pfcEta_CPU->Fill(pfc.eta());
      pfcPhi_CPU->Fill(pfc.phi());
	for(auto pfcx : *pfClustersFromCuda)
	  {
	    if(pfc.seed()==pfcx.seed()){
	      matched_pfcRh_CPU->Fill(pfc.recHitFractions().size());
          matched_pfcRh_GPU->Fill(pfcx.recHitFractions().size());
          matched_pfcEn_CPU->Fill(pfc.energy());
          matched_pfcEn_GPU->Fill(pfcx.energy());
          matched_pfcEta_CPU->Fill(pfc.eta());
          matched_pfcEta_GPU->Fill(pfcx.eta());
          matched_pfcPhi_CPU->Fill(pfc.phi());
          matched_pfcPhi_GPU->Fill(pfcx.phi());
          
          nRh_CPUvsGPU->Fill(pfcx.recHitFractions().size(),pfc.recHitFractions().size());
	      enPFCluster_CPUvsGPU->Fill(pfcx.energy(),pfc.energy());

          if(abs((pfcx.energy()-pfc.energy())/pfc.energy())>0.05){

            coordinate->Fill(pfcx.eta(),pfcx.phi());
            deltaRH->Fill((int)pfcx.recHitFractions().size() - (int)pfc.recHitFractions().size());
            deltaEn->Fill(pfcx.energy() - pfc.energy());
            deltaEta->Fill(pfcx.eta() - pfc.eta());
            deltaPhi->Fill(pfcx.phi() - pfc.phi());

            for(auto rhf: pfc.recHitFractions()){
              if(rhf.fraction()==1)layer->Fill(rhf.recHitRef()->depth());
             }
            }
	      /*if(pfcx.recHitFractions().size()>30){
		std::cout<<"fractions"<<std::endl;
		for(auto rhf: pfcx.recHitFractions()) std::cout<<rhf.fraction()<<"  ";
		std::cout<<std::endl;
		for(auto rhf: pfc.recHitFractions()) std::cout<<rhf.fraction()<<"  ";
		std::cout<<std::endl;*/
	      //}
	    }
	  }
    }
    
    __pfClustersFromCuda = *pfClustersFromCuda;      // For TTree
    for(auto pfc : *pfClustersFromCuda)
    {
  	    nRH_perPFCluster_GPU->Fill(pfc.recHitFractions().size());
	    enPFCluster_GPU->Fill(pfc.energy());
        pfcEta_GPU->Fill(pfc.eta());
        pfcPhi_GPU->Fill(pfc.phi());
    }
    
  }

  clusterTree->Fill();
  if (_prodInitClusters)
    e.put(std::move(pfClustersFromCuda), "initialClusters");
  e.put(std::move(pfClustersFromCuda));
}



