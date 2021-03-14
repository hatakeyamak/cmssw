#include "PFClusterProducerCudaHCAL.h"
#include "RecoParticleFlow/PFClusterProducer/plugins/PFClusterCudaHCAL.h"
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

PFClusterProducerCudaHCAL::PFClusterProducerCudaHCAL(const edm::ParameterSet& conf)
  :

  _prodInitClusters(conf.getUntrackedParameter<bool>("prodInitialClusters", false)) {
  _rechitsLabel = consumes<reco::PFRecHitCollection>(conf.getParameter<edm::InputTag>("recHitsSource"));

  //setup TTree
  clusterTree->Branch("Event", &numEvents);
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




  //setup topo cluster builder
  const edm::ParameterSet& initConf = conf.getParameterSet("initialClusteringStep");
  const std::string& initName = initConf.getParameter<std::string>("algoName");
  _initialClustering = InitialClusteringStepFactory::get()->create(initName, initConf, sumes);
  //setup pf cluster builder if requested
  const edm::ParameterSet& pfcConf = conf.getParameterSet("pfClusterBuilder");
  if (!pfcConf.empty()) {
    const std::string& pfcName = pfcConf.getParameter<std::string>("algoName");
    _pfClusterBuilder = PFClusterBuilderFactory::get()->create(pfcName, pfcConf);
    /*if (pfcConf.exists("allCellsPositionCalc")) {
    const edm::ParameterSet& acConf = pfcConf.getParameterSet("allCellsPositionCalc");
    const std::string& algoac = acConf.getParameter<std::string>("algoName");
    _allCellsPosCalcCuda = PFCPositionCalculatorFactory::get()->create(algoac, acConf);*/

    if (pfcConf.exists("positionCalc")) {
    const edm::ParameterSet& acConf = pfcConf.getParameterSet("positionCalc");
    const std::string& algoac = acConf.getParameter<std::string>("algoName");
    _allCellsPosCalcCuda = PFCPositionCalculatorFactory::get()->create(algoac, acConf);

  }
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

  if (_prodInitClusters) {
    produces<reco::PFClusterCollection>("initialClusters");
  }
  produces<reco::PFClusterCollection>();

  // KH:
  // Histograms
  // for eta binning
  double etaBins[83] = {-5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489, -3.314, -3.139, -2.964, -2.853, -2.65,  -2.5, -2.322, -2.172, -2.043, -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0., 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191, 4.363, 4.538, 4.716, 4.889, 5.191};
  double nPhiBins = 73;
  double phiBinWidth = M_PI / (nPhiBins - 1) * 2.;
  for (int jevt=0; jevt<100; jevt++){
    for (int idepth=0; idepth<7; idepth++){
      std::string str_gpu = "topoMap_GPU_d"+std::to_string(idepth+1)+"_evt"+std::to_string(jevt);
      std::string str_cpu = "topoMap_CPU_d"+std::to_string(idepth+1)+"_evt"+std::to_string(jevt);
      topoMap_GPU[idepth][jevt] = new TH2F(str_gpu.c_str(),str_gpu.c_str(),82,etaBins,nPhiBins,-M_PI - 0.25 * phiBinWidth, +M_PI + 0.75 * phiBinWidth);
      topoMap_CPU[idepth][jevt] = new TH2F(str_cpu.c_str(),str_cpu.c_str(),82,etaBins,nPhiBins,-M_PI - 0.25 * phiBinWidth, +M_PI + 0.75 * phiBinWidth);
      //
      str_gpu = "rhMap_GPU_d"+std::to_string(idepth+1)+"_evt"+std::to_string(jevt);
      str_cpu = "rhMap_CPU_d"+std::to_string(idepth+1)+"_evt"+std::to_string(jevt);
      rhMap_GPU[idepth][jevt] = new TH2F(str_gpu.c_str(),str_gpu.c_str(),82,etaBins,nPhiBins,-M_PI - 0.25 * phiBinWidth, +M_PI + 0.75 * phiBinWidth);
      rhMap_CPU[idepth][jevt] = new TH2F(str_cpu.c_str(),str_cpu.c_str(),82,etaBins,nPhiBins,-M_PI - 0.25 * phiBinWidth, +M_PI + 0.75 * phiBinWidth);

    }
  }
  ievt=-1; // initialize

}

PFClusterProducerCudaHCAL::~PFClusterProducerCudaHCAL()
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
  enPFCluster_CPUvsGPU_1d->Write();
  coordinate->Write();
  layer->Write();
  deltaSumSeed->Write();
  deltaRH->Write();
  deltaEn->Write();
  deltaEta->Write();
  deltaPhi->Write();
  // MyFile->Close();

  // KH:
  for (int jevt=0; jevt<100; jevt++){
    for (int idepth=0; idepth<7; idepth++){
      topoMap_GPU[idepth][jevt]->Write();
      topoMap_CPU[idepth][jevt]->Write();
      rhMap_GPU[idepth][jevt]->Write();
      rhMap_CPU[idepth][jevt]->Write();
    }
  }

  delete MyFile;
}

void PFClusterProducerCudaHCAL::beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup& es) {
  _initialClustering->update(es);
  if (_pfClusterBuilder)
    _pfClusterBuilder->update(es);
  if (_positionReCalc)
    _positionReCalc->update(es);
}

void PFClusterProducerCudaHCAL::produce(edm::Event& e, const edm::EventSetup& es) {
  _initialClustering->reset();
  if (_pfClusterBuilder)
    _pfClusterBuilder->reset();

  ievt++; // KH: ad-hoc event counter

  edm::Handle<reco::PFRecHitCollection> rechits;
  e.getByToken(_rechitsLabel, rechits);

  _initialClustering->updateEvent(e);

  std::cout << "Input rechit size: " << rechits->size() << std::endl;
  std::vector<bool> mask(rechits->size(), true);
  for (const auto& cleaner : _cleaners) {
    cleaner->clean(rechits, mask);
  }
  std::cout<<std::endl;
  for(int l=0;l<(int)mask.size();l++) {
    if (!mask[l]) std::cout<<"mask: "<<l<<" "<<mask[l]<<std::endl;
  }

  size_t rh_size = rechits->size();
  //std::cout<<rh_size<<std::endl;

  std::vector<float>                                    h_cuda_fracsum=std::vector<float>(rh_size,0);
  std::vector<int>                                      h_cuda_rhcount=std::vector<int>(rh_size,1);

  std::vector<float>                                    h_cuda_pfRhFrac(rechits->size()*100,-1.);
  std::vector<float>                                    h_cuda_pcRhFrac(rechits->size()*100,-1.);
  std::vector<int>                                      h_cuda_pfRhFracInd(rechits->size()*100,-1);
  std::vector<int>                                      h_cuda_pfNeighEightInd(rechits->size()*8,-1);
  std::vector<int>                                      h_cuda_pfNeighFourInd(rechits->size()*4,-1);
  std::vector<int>                                      h_cuda_pcRhFracInd(rechits->size()*100,-1);

  std::vector<float>                                    h_cuda_pfrh_x(rechits->size(),0);
  std::vector<float>                                    h_cuda_pfrh_y(rechits->size(),0);
  std::vector<float>                                    h_cuda_pfrh_z(rechits->size(),0);
  std::vector<double>                                    h_cuda_pfrh_energy(rechits->size(),0);
  std::vector<double>                                    h_cuda_pfrh_pt2(rechits->size(),0);
  std::vector<int>                                      h_cuda_pfrh_topoId(rechits->size(),-1);
  std::vector<int>                                      h_cuda_pfrh_isSeed(rechits->size(),0);
  std::vector<int>                                      h_cuda_pfrh_layer(rechits->size(),-999);
  std::vector<int>                                      h_cuda_pfrh_depth(rechits->size(),-999);



  int numbytes_float = rh_size*sizeof(float);
  int numbytes_double = rh_size*sizeof(double);
  int numbytes_int = rh_size*sizeof(int);

  auto d_cuda_rhcount = cms::cuda::make_device_unique<int[]>(numbytes_int, nullptr);
  auto d_cuda_fracsum = cms::cuda::make_device_unique<float[]>(numbytes_float, nullptr);

  float*                                    d_cuda_pfrh_x;
  cudaCheck(cudaMalloc(&d_cuda_pfrh_x, numbytes_float));
  float*                                    d_cuda_pfrh_y;
  cudaCheck(cudaMalloc(&d_cuda_pfrh_y, numbytes_float));
  float*                                    d_cuda_pfrh_z;
  cudaCheck(cudaMalloc(&d_cuda_pfrh_z, numbytes_float));
  double*                                    d_cuda_pfrh_energy;
  cudaCheck(cudaMalloc(&d_cuda_pfrh_energy, numbytes_double));
  double*                                    d_cuda_pfrh_pt2;
  cudaCheck(cudaMalloc(&d_cuda_pfrh_pt2, numbytes_double));
  int*                                      d_cuda_pfrh_topoId;
  cudaCheck(cudaMalloc(&d_cuda_pfrh_topoId, numbytes_int));
  int*                                      d_cuda_pfrh_isSeed;
  cudaCheck(cudaMalloc(&d_cuda_pfrh_isSeed, numbytes_int));
  int*                                      d_cuda_pfrh_layer;
  cudaCheck(cudaMalloc(&d_cuda_pfrh_layer, numbytes_int));
  int*                                      d_cuda_pfrh_depth;
  cudaCheck(cudaMalloc(&d_cuda_pfrh_depth, numbytes_int));
  int*                                      d_cuda_pfNeighEightInd;
  cudaCheck(cudaMalloc(&d_cuda_pfNeighEightInd, numbytes_int*8));
  int*                                      d_cuda_pfNeighFourInd;
  cudaCheck(cudaMalloc(&d_cuda_pfNeighFourInd, numbytes_int*4));

  int *d_cuda_pfRhFracInd;
  cudaCheck(cudaMalloc(&d_cuda_pfRhFracInd, numbytes_int*100));
  int *d_cuda_pcRhFracInd;
  cudaCheck(cudaMalloc(&d_cuda_pcRhFracInd, numbytes_int*100));
  float *d_cuda_pfRhFrac;
  cudaCheck(cudaMalloc(&d_cuda_pfRhFrac, numbytes_float*100));
  float *d_cuda_pcRhFrac;
  cudaCheck(cudaMalloc(&d_cuda_pcRhFrac, numbytes_float*100));

  int p=0;
  for (auto rh: *rechits){

    h_cuda_pfrh_x[p]=rh.position().x();
    h_cuda_pfrh_y[p]=rh.position().y();
    h_cuda_pfrh_z[p]=rh.position().z();
    h_cuda_pfrh_energy[p]=rh.energy();
    h_cuda_pfrh_pt2[p]=rh.pt2();
    h_cuda_pfrh_layer[p]=(int)rh.layer();
    h_cuda_pfrh_depth[p]=(int)rh.depth();
    h_cuda_pfrh_topoId[p]=p;
    // std::cout<<"depth  "<<h_cuda_pfrh_depth[p]<<std::endl;
    //std::cout<<"layer  "<<h_cuda_pfrh_layer[p]<<std::endl;

    auto theneighboursEight = rh.neighbours8();
    int z = 0;
    // h_cuda_pfNeighEightInd[9*p] = p;
    for(auto nh: theneighboursEight)
      {
	h_cuda_pfNeighEightInd[8*p+z] = nh;
	z++;
      }

    auto theneighboursFour = rh.neighbours4();
    int y = 0;
    for(auto nh: theneighboursFour)
      {
	h_cuda_pfNeighFourInd[4*p+y] = nh;
	y++;
      }


    p++;
  }//end of rechit loop

  cudaCheck(cudaMemcpy(d_cuda_fracsum.get(), h_cuda_fracsum.data(), numbytes_float, cudaMemcpyHostToDevice));
  cudaCheck(cudaMemcpy(d_cuda_rhcount.get(), h_cuda_rhcount.data(), numbytes_int, cudaMemcpyHostToDevice));

  cudaCheck(cudaMemcpy(d_cuda_pfrh_x, h_cuda_pfrh_x.data(), numbytes_float, cudaMemcpyHostToDevice));
  cudaCheck(cudaMemcpy(d_cuda_pfrh_y, h_cuda_pfrh_y.data(), numbytes_float, cudaMemcpyHostToDevice));
  cudaCheck(cudaMemcpy(d_cuda_pfrh_z, h_cuda_pfrh_z.data(), numbytes_float, cudaMemcpyHostToDevice));
  cudaCheck(cudaMemcpy(d_cuda_pfrh_energy, h_cuda_pfrh_energy.data(), numbytes_double, cudaMemcpyHostToDevice));
  cudaCheck(cudaMemcpy(d_cuda_pfrh_pt2, h_cuda_pfrh_pt2.data(), numbytes_double, cudaMemcpyHostToDevice));
  cudaCheck(cudaMemcpy(d_cuda_pfrh_topoId, h_cuda_pfrh_topoId.data(), numbytes_int, cudaMemcpyHostToDevice));
  cudaCheck(cudaMemcpy(d_cuda_pfrh_isSeed, h_cuda_pfrh_isSeed.data(), numbytes_int, cudaMemcpyHostToDevice));
  cudaCheck(cudaMemcpy(d_cuda_pfrh_layer, h_cuda_pfrh_layer.data(), numbytes_int, cudaMemcpyHostToDevice));
  cudaCheck(cudaMemcpy(d_cuda_pfrh_depth, h_cuda_pfrh_depth.data(), numbytes_int, cudaMemcpyHostToDevice));
  cudaCheck(cudaMemcpy(d_cuda_pfNeighEightInd, h_cuda_pfNeighEightInd.data(), numbytes_int*8, cudaMemcpyHostToDevice));
  cudaCheck(cudaMemcpy(d_cuda_pfNeighFourInd, h_cuda_pfNeighFourInd.data(), numbytes_int*4, cudaMemcpyHostToDevice));

  cudaCheck(cudaMemcpy(d_cuda_pfRhFrac, h_cuda_pfRhFrac.data(), numbytes_float*100, cudaMemcpyHostToDevice));
  cudaCheck(cudaMemcpy(d_cuda_pcRhFrac, h_cuda_pcRhFrac.data(), numbytes_float*100, cudaMemcpyHostToDevice));
  cudaCheck(cudaMemcpy(d_cuda_pfRhFracInd, h_cuda_pfRhFracInd.data(), numbytes_int*100, cudaMemcpyHostToDevice));
  cudaCheck(cudaMemcpy(d_cuda_pcRhFracInd, h_cuda_pcRhFracInd.data(), numbytes_int*100, cudaMemcpyHostToDevice));

  /*  PFClusterCudaHCAL::PFRechitToPFCluster_HCALV1(rh_size,
					      d_cuda_pfrh_x,
					      d_cuda_pfrh_y,
					      d_cuda_pfrh_z,
					      d_cuda_pfrh_energy,
					      d_cuda_pfrh_pt2,
					      d_cuda_pfrh_isSeed,
					      d_cuda_pfrh_topoId,
					      d_cuda_pfrh_layer,
					      d_cuda_pfrh_depth,
					      d_cuda_pfNeighEightInd,
					      d_cuda_pfNeighFourInd,
					      d_cuda_pfRhFrac,
					      d_cuda_pfRhFracInd,
					      d_cuda_pcRhFracInd,
					      d_cuda_pcRhFrac
					      );*/
     /*
     PFClusterCudaHCAL::PFRechitToPFCluster_HCALV2(rh_size, 
					      d_cuda_pfrh_x,  
					      d_cuda_pfrh_y,  
					      d_cuda_pfrh_z, 
					      d_cuda_pfrh_energy, 
					      d_cuda_pfrh_pt2, 	
					      d_cuda_pfrh_isSeed,
					      d_cuda_pfrh_topoId,
					      d_cuda_pfrh_layer,
					      d_cuda_pfrh_depth,
					      d_cuda_pfNeighEightInd,
					      d_cuda_pfNeighFourInd,

					      d_cuda_pcRhFracInd,
					      d_cuda_pcRhFrac,
					      d_cuda_fracsum.get(),
					      d_cuda_rhcount.get()
					      );
     */
     
     //PFClusterCudaHCAL::PFRechitToPFCluster_HCAL_serialize(rh_size, 
     //PFClusterCudaHCAL::PFRechitToPFCluster_HCAL_serialize_topoParallel(rh_size, 
     //PFClusterCudaHCAL::PFRechitToPFCluster_HCAL_serialize_seedingParallel(rh_size, 
     //PFClusterCudaHCAL::PFRechitToPFCluster_HCAL_serialize_step1Parallel(rh_size, 
     //PFClusterCudaHCAL::PFRechitToPFCluster_HCAL_serialize_step2Parallel(rh_size, 
     PFClusterCudaHCAL::PFRechitToPFCluster_HCAL_serialize(rh_size, 
					      d_cuda_pfrh_x,  
					      d_cuda_pfrh_y,  
					      d_cuda_pfrh_z, 
					      d_cuda_pfrh_energy, 
					      d_cuda_pfrh_pt2, 	
					      d_cuda_pfrh_isSeed,
					      d_cuda_pfrh_topoId,
					      d_cuda_pfrh_layer, 
					      d_cuda_pfrh_depth, 
					      d_cuda_pfNeighEightInd, 
					      d_cuda_pfNeighFourInd, 
					      
					      d_cuda_pcRhFracInd,
					      d_cuda_pcRhFrac,
					      d_cuda_fracsum.get(),
					      d_cuda_rhcount.get()
					      );
					

  cudaMemcpy(h_cuda_pcRhFracInd.data()    , d_cuda_pcRhFracInd  , numbytes_int*100 , cudaMemcpyDeviceToHost);  
  cudaMemcpy(h_cuda_pcRhFrac.data()       , d_cuda_pcRhFrac  , numbytes_float*100 , cudaMemcpyDeviceToHost);  
  cudaMemcpy(h_cuda_pfrh_isSeed.data()    , d_cuda_pfrh_isSeed  , numbytes_int , cudaMemcpyDeviceToHost);  
  cudaMemcpy(h_cuda_pfrh_topoId.data()    , d_cuda_pfrh_topoId  , numbytes_int , cudaMemcpyDeviceToHost);  
  cudaMemcpy(h_cuda_pfNeighEightInd.data()    , d_cuda_pfNeighEightInd  , numbytes_int*8 , cudaMemcpyDeviceToHost);  
  
  if(doComparison){ 
    for(unsigned int i=0;i<rh_size;i++){
      int topoIda=h_cuda_pfrh_topoId[i];
      //
      //
      for(unsigned int j=0;j<8;j++){
	if(h_cuda_pfNeighEightInd[i*8+j]>-1 && h_cuda_pfrh_topoId[h_cuda_pfNeighEightInd[i*8+j]]!=topoIda) std::cout<<"DIFFERENT TOPOID "<<i<<"  "<<j<<"  "<<topoIda<<"  "<<h_cuda_pfrh_topoId[h_cuda_pfNeighEightInd[i*8+j]]<<std::endl;
      }
      //
      //KH
      int z = 0;
      int p = i;
      auto rh = (*rechits)[i];
      auto theneighboursEight = rh.neighbours8();
      GlobalPoint gp1(h_cuda_pfrh_x[p],h_cuda_pfrh_y[p],h_cuda_pfrh_z[p]);
      for(auto nh: theneighboursEight){
	GlobalPoint gp2(h_cuda_pfrh_x[nh],h_cuda_pfrh_y[nh],h_cuda_pfrh_z[nh]);
	// std::cout << "p,z,nh: " 
	// 	  << p  << " " << h_cuda_pfrh_depth[p]  << " " << gp1.eta() << " " << gp1.phi() << " " 
	// 	  << h_cuda_pfrh_topoId[p] << " "
	// 	  << nh << " " 
	// 	  << z  << " " << h_cuda_pfrh_depth[nh] << " " << gp2.eta() << " " << gp2.phi() << " "
	// 	  << h_cuda_pfrh_topoId[nh]
	// 	  << std::endl;
	if (h_cuda_pfrh_topoId[p]!=h_cuda_pfrh_topoId[nh]) std::cout 
							     << "KH: topo mismatch" << std::endl;
	z++;
      }
    }
    //KH
    // for(int l=0; l<(int)h_cuda_pfrh_topoId.size();l++){
    //   GlobalPoint gp(h_cuda_pfrh_x[l],h_cuda_pfrh_y[l],h_cuda_pfrh_z[l]);
    //   std::cout << "KHFill2: " << " "
    // 		<< l << " " 
    // 		<< gp.eta() << " " << gp.phi() << " "
    // 		<< h_cuda_pfrh_depth[l] << " " << h_cuda_pfrh_topoId[l]
    // 		<< std::endl;
    // }
    // KH:
    std::vector< std::pair<int, std::pair<double, double> > > myvec; // store list of pf rechits
    myvec.clear();
    std::vector<int> topoid;
    std::vector<int> topoidd[7];
    for(int l=0; l<(int)h_cuda_pfrh_topoId.size();l++){
    if (h_cuda_pfrh_topoId[l]>-1.) { 
      topoid.push_back(h_cuda_pfrh_topoId[l]);
      topoidd[h_cuda_pfrh_depth[l]-1].push_back(h_cuda_pfrh_topoId[l]);
      GlobalPoint gp(h_cuda_pfrh_x[l],h_cuda_pfrh_y[l],h_cuda_pfrh_z[l]);
      myvec.push_back(std::make_pair(h_cuda_pfrh_depth[l],
				     std::make_pair(gp.eta(),gp.phi())));
      
      if (ievt<100){ // check only first 100 events
	int ieta=topoMap_GPU[h_cuda_pfrh_depth[l]-1][ievt]->GetXaxis()->FindBin(gp.eta());
	int iphi=topoMap_GPU[h_cuda_pfrh_depth[l]-1][ievt]->GetYaxis()->FindBin(gp.phi());
	if (topoMap_GPU[h_cuda_pfrh_depth[l]-1][ievt]->GetBinContent(ieta,iphi)!=0.){
	  std::cout << "KH already filled" << " "
		    << l << " " << ieta << " " << iphi << " " 
		    << h_cuda_pfrh_depth[l] << " " << h_cuda_pfrh_topoId[l]
		    << std::endl; 
	} else {
	  topoMap_GPU[h_cuda_pfrh_depth[l]-1][ievt]->SetBinContent(ieta,iphi,h_cuda_pfrh_topoId[l]);
	  rhMap_GPU[h_cuda_pfrh_depth[l]-1][ievt]->SetBinContent(ieta,iphi,1.);
	  // std::cout << "KHFill: " << " "
	  // 	    << l << " " << ieta << " " << iphi << " " 
	  // 	    << gp.eta() << " " << gp.phi() << " "
	  // 	    << h_cuda_pfrh_depth[l] << " " << h_cuda_pfrh_topoId[l]
	  // 	    << std::endl; 
	}
      }
    } // if
    } // for
    std::cout << "number of rechits (GPU): " << topoid.size() << std::endl;
    sort(topoid.begin(), topoid.end());
    std::vector<int>::iterator it;
    it = unique(topoid.begin(), topoid.end()); topoid.resize(distance(topoid.begin(),it)); 
    std::cout << "number of topo    (GPU): " << topoid.size() << std::endl;
    //
    std::vector<int> topoidd2[7];
    for (int i=0; i<7; i++) 
      std::copy(topoidd[i].begin(), topoidd[i].end(), back_inserter(topoidd2[i]));
    for (int i=0; i<7; i++) {
      sort(topoidd[i].begin(), topoidd[i].end());
      it = unique(topoidd[i].begin(), topoidd[i].end()); topoidd[i].resize(distance(topoidd[i].begin(),it)); 
      std::cout << "number of topo (rechits) d" << i+1 << " (GPU): " << topoidd[i].size() 
		<< " " << topoidd2[i].size() << std::endl;
    }
    std::cout << "number of topo    (GPU): " << topoidd[0].size()+topoidd[1].size()+topoidd[2].size()+topoidd[3].size()+topoidd[4].size()+topoidd[5].size()+topoidd[6].size() << std::endl;
    std::cout << "number of hits    (GPU): " << topoidd2[0].size()+topoidd2[1].size()+topoidd2[2].size()+topoidd2[3].size()+topoidd2[4].size()+topoidd2[5].size()+topoidd2[6].size() << std::endl;
  }
  

  //free up
  cudaFree(d_cuda_pfrh_x);
  cudaFree(d_cuda_pfrh_y);
  cudaFree(d_cuda_pfrh_z);
  cudaFree(d_cuda_pfrh_energy);
  cudaFree(d_cuda_pfrh_layer);
  cudaFree(d_cuda_pfrh_depth);
  cudaFree(d_cuda_pfrh_isSeed);
  cudaFree(d_cuda_pfrh_topoId);
  cudaFree(d_cuda_pfrh_pt2);
  cudaFree(d_cuda_pfNeighEightInd);
  cudaFree(d_cuda_pfNeighFourInd);
  cudaFree(d_cuda_pfRhFracInd);
  cudaFree(d_cuda_pcRhFracInd);
  cudaFree(d_cuda_pfRhFrac);
  cudaFree(d_cuda_pcRhFrac);

  auto pfClustersFromCuda = std::make_unique<reco::PFClusterCollection>();
  pfClustersFromCuda.reset(new reco::PFClusterCollection);
  for(int n=0; n<(int)rh_size; n++){
    if(h_cuda_pfrh_isSeed[n]==1){
      reco::PFCluster temp;
      temp.setSeed((*rechits)[n].detId());
      //if((*rechits)[n]==nullptr) std::cout<<"null det seed: "<<n<<std::endl;
      for(int k=0;k<100;k++){
	if(h_cuda_pcRhFracInd[n*100+k] > -1){
	  const reco::PFRecHitRef& refhit = reco::PFRecHitRef(rechits,h_cuda_pcRhFracInd[n*100+k]);
	  temp.addRecHitFraction( reco::PFRecHitFraction(refhit, h_cuda_pcRhFrac[n*100+k]) );
	}
	if(h_cuda_pcRhFracInd[n*100+k] < 0.) break;
      }
      pfClustersFromCuda->push_back(temp);
    }
  }
  //_positionReCalc->calculateAndSetPositions(*pfClustersFromCuda);
  _allCellsPosCalcCuda->calculateAndSetPositions(*pfClustersFromCuda);

  //if (_energyCorrector) {
  //  _energyCorrector->correctEnergies(*pfClustersFromCuda);
  //}

  float sumEn_CPU = 0.f;
  if(doComparison)
  {

    //
    // KH:
    // CPU version topo clustering and checks
    //
    std::vector<bool> seedable(rechits->size(), false);
    _seedFinder->findSeeds(rechits, mask, seedable);
    auto initialClusters = std::make_unique<reco::PFClusterCollection>();
    _initialClustering->buildClusters(rechits, mask, seedable, *initialClusters);

    __initialClusters = *initialClusters;  // For TTree
    
    int topoRhCount=0;
    int topoIndex=0;
    int rhIndex=0;
    std::vector<int> topoSizeCPU; // vector of topo cluster size
    std::vector< std::pair<int, std::pair<double, double> > > myvec; // store list of pf rechits
    std::vector<int> topoidd[7];
    std::vector<int> rhidd[7];

    for(auto pfc : *initialClusters)
      {
  	nTopo_CPU->Fill(pfc.recHitFractions().size());
	topoSizeCPU.push_back(pfc.recHitFractions().size());
	topoIndex++;
	topoRhCount=topoRhCount+pfc.recHitFractions().size();
	//
	// check individual topo clusters
	myvec.clear();
	int depth=0;
	for (const reco::PFRecHitFraction& rhf : pfc.recHitFractions()) {
	  const reco::PFRecHitRef& refhit = rhf.recHitRef();
	  if (depth==0) depth = refhit->depth();
	  else {
	    if (depth!=refhit->depth()) std::cout << "topo include multiple depths?" << std::endl;
	  }
	  myvec.push_back(std::make_pair(refhit->depth(),
					 std::make_pair(refhit->position().eta(),refhit->position().phi())));
	  if (depth>0 && depth<=7) {
	    if (ievt<100){ // check only first 100 events
	      topoMap_CPU[depth-1][ievt]->Fill(refhit->position().eta(),refhit->position().phi(),topoIndex);
	      rhMap_CPU[depth-1][ievt]->Fill(refhit->position().eta(),refhit->position().phi(),1.);
	    }
	    rhidd[depth-1].push_back(rhIndex);
	  }
	  rhIndex++;
	}
	if (depth>0 && depth<=7) topoidd[depth-1].push_back(topoIndex);
	/*
	  std::sort(myvec.begin(), myvec.end());
	  for (unsigned int i=0; i<myvec.size(); i++)
	    std::cout << i << " " << myvec[i].first  << " "
	    << myvec[i].second.first << " "
	    << myvec[i].second.second << std::endl;
	    std::cout << std::endl;
	*/
	topoEn_CPU->Fill(pfc.energy());
	topoEta_CPU->Fill(pfc.eta());
	topoPhi_CPU->Fill(pfc.phi());
      }
    for (int i=0; i<7; i++) {
      std::vector<int>::iterator it
	= unique(topoidd[i].begin(), topoidd[i].end()); topoidd[i].resize(distance(topoidd[i].begin(),it)); 
      std::cout << "number of topo (rechits) d" << i+1 << " (CPU): " << topoidd[i].size() << " "
		<< rhidd[i].size() << std::endl;
    }
    std::cout << "number of topo    (CPU): " << topoidd[0].size()+topoidd[1].size()+topoidd[2].size()+topoidd[3].size()+topoidd[4].size()+topoidd[5].size()+topoidd[6].size() << std::endl;
    std::cout << "number of hits    (CPU): " << rhidd[0].size()+rhidd[1].size()+rhidd[2].size()+rhidd[3].size()+rhidd[4].size()+rhidd[5].size()+rhidd[6].size() << std::endl;

    nPFCluster_CPU->Fill(initialClusters->size());
    std::sort (h_cuda_pfrh_topoId.begin(), h_cuda_pfrh_topoId.end());

    //
    // KH:
    // GPU topo cluster checks
    //
    int topoCount=1;
    int intTopoCount=0;
    std::vector<int> topoSizeGPU; // size of each topo cluster
    std::vector<int> topoIDGPU; // ID of each topo cluster
    for(int l=0; l<(int)h_cuda_pfrh_topoId.size();l++){
      if((h_cuda_pfrh_topoId[l]==h_cuda_pfrh_topoId[l+1]) && h_cuda_pfrh_topoId[l]>-1.) topoCount++;
      else if(h_cuda_pfrh_topoId[l]>-1.){
	topoSizeGPU.push_back(topoCount);
	if (topoCount>0) topoIDGPU.push_back(h_cuda_pfrh_topoId[l]); // store topoID info
	nTopo_GPU->Fill(topoCount);
	topoCount=1;
	intTopoCount++;
      }
    }
    nPFCluster_GPU->Fill(intTopoCount);
    LOGVERB("PFClusterProducer::produce()") << *_initialClustering;

    int seedSumCPU=0;
    int seedSumGPU=0;
    int maskSize = 0;
    for (int j=0;j<(int)seedable.size(); j++) seedSumCPU=seedSumCPU+seedable[j];
    for (int j=0;j<(int)h_cuda_pfrh_isSeed.size(); j++) seedSumGPU=seedSumGPU +h_cuda_pfrh_isSeed[j];
    for (int j=0;j<(int)mask.size(); j++) maskSize=maskSize +mask[j];

    /*
    for (int j=0;j<(int)seedable.size(); j++){
      if(seedable[j]!=h_cuda_pfrh_isSeed[j]){
	std::cout<<j<<" "<<seedable[j]<<"  "<<h_cuda_pfrh_isSeed[j]<<", depth:  "<<(*rechits)[j].depth()<<", layer: "<<(*rechits)[j].layer()<<std::endl;
	std::cout<<"pt2: "<<(*rechits)[j].pt2()<<std::endl;
	std::cout<<"energy: "<<(*rechits)[j].energy()<<std::endl;
	auto theneighboursFour = (*rechits)[j].neighbours4();
	for(auto nh: theneighboursFour)
	  {
	    std::cout<<"neigh: "<<(*rechits)[nh].energy()<<std::endl;
	  }
      }
    }
    */

    std::cout<<"sum CPU seeds: "<<seedSumCPU<<std::endl;
    std::cout<<"sum GPU seeds: "<<seedSumGPU<<std::endl;
    std::cout<<"sum mask  : "<<maskSize<<std::endl;

    std::cout<<"sum rechits               : "<<rh_size<<std::endl;
    std::cout<<"sum rechits in topo  (CPU): "<<topoRhCount<<std::endl;

    std::cout<<"# of topo clusters (CPU): " << initialClusters->size() <<std::endl;
    std::cout<<"# of topo clusters (GPU): " << intTopoCount <<std::endl;
    
    /*
    std::cout<<"sum CPU seeds: "<<seedSumCPU<<std::endl;
    std::cout<<"sum GPU seeds: "<<seedSumGPU<<std::endl;
    //std::cout<<"sum rechits  : "<<rh_size<<std::endl;
    std::cout<<"sum mask  : "<<maskSize<<std::endl;
    */

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

    std::cout<<"# of PF clusters (CPU): " << pfClusters->size() <<std::endl;
    std::cout<<"# of PF clusters (GPU): " << pfClustersFromCuda->size() <<std::endl;

    //
    // KH:
    // Check CPU version PF clusters vs GPU version PF clusters
    //
    //std::cout<<"HCAL pfClusters->size() = "<<pfClusters->size()<<std::endl; 
    __pfClusters = *pfClusters;  // For TTree
    for(auto pfc : *pfClusters)
    {
      nRH_perPFCluster_CPU->Fill(pfc.recHitFractions().size());
	  enPFCluster_CPU->Fill(pfc.energy());
      pfcEta_CPU->Fill(pfc.eta());
      pfcPhi_CPU->Fill(pfc.phi());
    sumEn_CPU += pfc.energy();
    //if (numEvents < 1) std::cout<<pfc.energy()<<std::endl;	
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
	      enPFCluster_CPUvsGPU_1d->Fill((pfcx.energy()-pfc.energy())/pfc.energy());
	      if(abs((pfcx.energy()-pfc.energy())/pfc.energy())>0.05){

		coordinate->Fill(pfcx.eta(),pfcx.phi());
		deltaRH->Fill((int)pfcx.recHitFractions().size() - (int)pfc.recHitFractions().size());
		deltaEn->Fill(pfcx.energy() - pfc.energy());
		deltaEta->Fill(pfcx.eta() - pfc.eta());
		deltaPhi->Fill(pfcx.phi() - pfc.phi());

		for(auto rhf: pfc.recHitFractions()){
		  if(rhf.fraction()==1)layer->Fill(rhf.recHitRef()->depth());
		}
		std::cout<<std::endl;
		std::cout<<"fractions"<<std::endl;
		for(auto rhf: pfcx.recHitFractions()) std::cout<<rhf.fraction()<<", eta:"<<rhf.recHitRef()->positionREP().eta()<<", phi:"<< rhf.recHitRef()->positionREP().phi()<<"  ";
		std::cout<<std::endl;
		for(auto rhf: pfc.recHitFractions()) std::cout<<rhf.fraction()<<", eta:"<<rhf.recHitRef()->positionREP().eta()<<", phi:"<< rhf.recHitRef()->positionREP().phi()<<"  ";
		std::cout<<std::endl;
	      }
	      /*if(abs((int)(pfcx.recHitFractions().size() - pfc.recHitFractions().size() ))>30){
		std::cout<<"fractions"<<std::endl;
		for(auto rhf: pfcx.recHitFractions()) std::cout<<rhf.fraction()<<"  ";
		std::cout<<std::endl;
		for(auto rhf: pfc.recHitFractions()) std::cout<<rhf.fraction()<<"  ";
		std::cout<<std::endl;
		}*/
	    }
	  }
    }

    //
    // KH:
    // Check GPU version PF clustering
    //
    __pfClustersFromCuda = *pfClustersFromCuda;      // For TTree
    for(auto pfc : *pfClustersFromCuda)
    {
        nRH_perPFCluster_GPU->Fill(pfc.recHitFractions().size());
        enPFCluster_GPU->Fill(pfc.energy());
        pfcEta_GPU->Fill(pfc.eta());
        pfcPhi_GPU->Fill(pfc.phi());
    }
  }

  numEvents++;
  //std::cout<<"Sum En CPU = "<<sumEn_CPU<<std::endl;
  //std::cout<<"***** Filling event "<<numEvents<<std::endl;
  clusterTree->Fill();
  if (_prodInitClusters)
    e.put(std::move(pfClustersFromCuda), "initialClusters");
  e.put(std::move(pfClustersFromCuda));
}
