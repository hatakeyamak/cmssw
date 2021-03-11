#ifndef PFClusterCudaHCAL_h
#define PFClusterCudaHCAL_h
//#include <thrust/host_vector.h>
//#include <thrust/device_vector.h>
#include <Eigen/Dense>

namespace PFClusterCudaHCAL { 

  void PFRechitToPFCluster_HCALV1(size_t size, 
				float* pfrh_x,
				float* pfrh_y, 
				float* pfrh_z,	
				float* pfrh_energy,	
				float* pfrh_pt2, 
				int* pfrh_isSeed, 
				int* pfrh_topoId, 
				int* pfrh_layer, 
				int* pfrh_depth, 
				int* neigh8_Ind, 				
				int* neigh4_Ind, 				
				float* pfrhfrac, 
				int* pfrhind, 
				int* pcrhind,
				float* pcrhfracind
				);

void PFRechitToPFCluster_HCALV2(size_t size, 
				float* pfrh_x,
				float* pfrh_y, 
				float* pfrh_z,	
				double* pfrh_energy,	
				double* pfrh_pt2, 
				int* pfrh_isSeed, 
				int* pfrh_topoId, 
				int* pfrh_layer, 
				int* pfrh_depth, 
				int* neigh8_Ind, 				
				int* neigh4_Ind, 				
				int* pcrhind,
				float* pcrhfracind,
				float* fracSum,
				int* rhCount
				);

void PFRechitToPFCluster_HCAL_serialize(size_t size, 
				float* pfrh_x,
				float* pfrh_y, 
				float* pfrh_z,	
				double* pfrh_energy,	
				double* pfrh_pt2, 
				int* pfrh_isSeed, 
				int* pfrh_topoId, 
				int* pfrh_layer, 
				int* pfrh_depth, 
				int* neigh8_Ind, 				
				int* neigh4_Ind, 				
				int* pcrhind,
				float* pcrhfracind,
				float* fracSum,
				int* rhCount
				);

void PFRechitToPFCluster_HCAL_serialize_seedingParallel(size_t size, 
				float* pfrh_x,
				float* pfrh_y, 
				float* pfrh_z,	
				double* pfrh_energy,	
				double* pfrh_pt2, 
				int* pfrh_isSeed, 
				int* pfrh_topoId, 
				int* pfrh_layer, 
				int* pfrh_depth, 
				int* neigh8_Ind, 				
				int* neigh4_Ind, 				
				int* pcrhind,
				float* pcrhfracind,
				float* fracSum,
				int* rhCount
				);

void PFRechitToPFCluster_HCAL_serialize_topoParallel(size_t size, 
				float* pfrh_x,
				float* pfrh_y, 
				float* pfrh_z,	
				double* pfrh_energy,	
				double* pfrh_pt2, 
				int* pfrh_isSeed, 
				int* pfrh_topoId, 
				int* pfrh_layer, 
				int* pfrh_depth, 
				int* neigh8_Ind, 				
				int* neigh4_Ind, 				
				int* pcrhind,
				float* pcrhfracind,
				float* fracSum,
				int* rhCount
				);

void PFRechitToPFCluster_HCAL_serialize_step1Parallel(size_t size, 
				float* pfrh_x,
				float* pfrh_y, 
				float* pfrh_z,	
				double* pfrh_energy,	
				double* pfrh_pt2, 
				int* pfrh_isSeed, 
				int* pfrh_topoId, 
				int* pfrh_layer, 
				int* pfrh_depth, 
				int* neigh8_Ind, 				
				int* neigh4_Ind, 				
				int* pcrhind,
				float* pcrhfracind,
				float* fracSum,
				int* rhCount
				);

void PFRechitToPFCluster_HCAL_serialize_step2Parallel(size_t size, 
				float* pfrh_x,
				float* pfrh_y, 
				float* pfrh_z,	
				double* pfrh_energy,	
				double* pfrh_pt2, 
				int* pfrh_isSeed, 
				int* pfrh_topoId, 
				int* pfrh_layer, 
				int* pfrh_depth, 
				int* neigh8_Ind, 				
				int* neigh4_Ind, 				
				int* pcrhind,
				float* pcrhfracind,
				float* fracSum,
				int* rhCount
				);


}  // namespace cudavectors

#endif  // ClusterCudaECAL_h
