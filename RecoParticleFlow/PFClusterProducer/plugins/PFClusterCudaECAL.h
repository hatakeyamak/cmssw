#ifndef PFClusterCudaECAL_h
#define PFClusterCudaECAL_h
//#include <thrust/host_vector.h>
//#include <thrust/device_vector.h>
#include <Eigen/Dense>

namespace PFClusterCudaECAL {

  


  void PFRechitToPFCluster_ECALV2(size_t size, 
				float* pfrh_x,
				float* pfrh_y, 
				float* pfrh_z,	
				float* pfrh_energy,	
				float* pfrh_pt2, 
				int* pfrh_isSeed, 
				int* pfrh_topoId, 
				int* pfrh_layer, 
				int* neigh8_Ind, 				
				float* pfrhfrac, 
				int* pfrhind, 
				int* pcrhind,
				float* pcrhfracind,
				float* fracSum,
				int* rhCount
				);

 void PFRechitToPFCluster_ECALV1(size_t size, 
				float* pfrh_x,
				float* pfrh_y, 
				float* pfrh_z,	
				float* pfrh_energy,	
				float* pfrh_pt2, 
				int* pfrh_isSeed, 
				int* pfrh_topoId, 
				int* pfrh_layer, 
				int* neigh8_Ind, 				
				float* pfrhfrac, 
				int* pfrhind, 
				int* pcrhind,
				float* pcrhfracind);

  
  void PFRechitToPFCluster_ECAL_serialize(size_t size, 
				float* pfrh_x,
				float* pfrh_y, 
				float* pfrh_z,	
				float* pfrh_energy,	
				float* pfrh_pt2, 
				int* pfrh_isSeed, 
				int* pfrh_topoId, 
				int* pfrh_layer, 
				int* neigh8_Ind, 				
				float* pfrhfrac, 
				int* pfrhind, 
				int* pcrhind,
				float* pcrhfracind,
				float* fracSum,
				int* rhCount
				);

  void PFRechitToPFCluster_ECAL_serialize_seedingParallel(size_t size, 
				float* pfrh_x,
				float* pfrh_y, 
				float* pfrh_z,	
				float* pfrh_energy,	
				float* pfrh_pt2, 
				int* pfrh_isSeed, 
				int* pfrh_topoId, 
				int* pfrh_layer, 
				int* neigh8_Ind, 				
				float* pfrhfrac, 
				int* pfrhind, 
				int* pcrhind,
				float* pcrhfracind,
				float* fracSum,
				int* rhCount
				);

  void PFRechitToPFCluster_ECAL_serialize_topoParallel(size_t size, 
				float* pfrh_x,
				float* pfrh_y, 
				float* pfrh_z,	
				float* pfrh_energy,	
				float* pfrh_pt2, 
				int* pfrh_isSeed, 
				int* pfrh_topoId, 
				int* pfrh_layer, 
				int* neigh8_Ind, 				
				float* pfrhfrac, 
				int* pfrhind, 
				int* pcrhind,
				float* pcrhfracind,
				float* fracSum,
				int* rhCount
				);

  void PFRechitToPFCluster_ECAL_serialize_step1Parallel(size_t size, 
				float* pfrh_x,
				float* pfrh_y, 
				float* pfrh_z,	
				float* pfrh_energy,	
				float* pfrh_pt2, 
				int* pfrh_isSeed, 
				int* pfrh_topoId, 
				int* pfrh_layer, 
				int* neigh8_Ind, 				
				float* pfrhfrac, 
				int* pfrhind, 
				int* pcrhind,
				float* pcrhfracind,
				float* fracSum,
				int* rhCount
				);

  void PFRechitToPFCluster_ECAL_serialize_step2Parallel(size_t size, 
				float* pfrh_x,
				float* pfrh_y, 
				float* pfrh_z,	
				float* pfrh_energy,	
				float* pfrh_pt2, 
				int* pfrh_isSeed, 
				int* pfrh_topoId, 
				int* pfrh_layer, 
				int* neigh8_Ind, 				
				float* pfrhfrac, 
				int* pfrhind, 
				int* pcrhind,
				float* pcrhfracind,
				float* fracSum,
				int* rhCount
				);


}  // namespace cudavectors

#endif  // ClusterCudaECAL_h
