
#include <cmath>
#include <iostream>

// CUDA include files
#include <cuda_runtime.h>

// CMSSW include files
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"
#include "RecoParticleFlow/PFClusterProducer/plugins/PFClusterCudaECAL.h"
#include <Eigen/Dense>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>

namespace PFClusterCudaECAL {

  // THE ART OF HARDCODING
  // these numbers should be copied over during initialization
  __constant__ float showerSigma = 1.5;
  __constant__ float recHitEnergyNormEB = 0.08;
  __constant__ float recHitEnergyNormEE = 0.3;
  __constant__ float minFracToKeep = 0.0000001;

  __constant__ float seedEThresholdEB = 0.23;
  __constant__ float seedEThresholdEE = 0.6;
  __constant__ float seedPt2ThresholdEB = 0.0;
  __constant__ float seedPt2hresholdEE = 0.0225;

  __constant__ float topoEThresholdEB = 0.08;
  __constant__ float topoEThresholdEE = 0.3;

  __constant__ int nNeigh = 8;
  __constant__ int maxSize = 50;

 
   
 __global__ void seedingKernel_ECAL(
     				    size_t size, 
				    float* pfrh_energy,
				    float* pfrh_pt2,
				    int*   pfrh_isSeed,
				    int*   pfrh_topoId,
				    int*   pfrh_layer,
				    int*   neigh8_Ind
				    ) {

   int i = threadIdx.x+blockIdx.x*blockDim.x;

   if(i<size) {        
     if( ( pfrh_layer[i] == -1 && pfrh_energy[i]>seedEThresholdEB && pfrh_pt2[i]>seedPt2ThresholdEB) || ( (pfrh_layer[i] == -2) && pfrh_energy[i]>seedEThresholdEE && pfrh_pt2[i]>seedPt2hresholdEE) )
       {
	 pfrh_isSeed[i]=1;		 
	 for(int k=0; k<nNeigh; k++){
	   if(neigh8_Ind[nNeigh*i+k]<0) continue; 
	   if(pfrh_energy[i]<pfrh_energy[neigh8_Ind[nNeigh*i+k]]){
	     pfrh_isSeed[i]=0;
	     pfrh_topoId[i]=-1;	     
	     break;
	   }
	 }		
       }
     else{ 
       pfrh_topoId[i]=-1;
       pfrh_isSeed[i]=0;
       	    
     }     
   }
 }
   
 __global__ void seedingKernel_ECAL_serialize(
     				    size_t size, 
				    float* pfrh_energy,
				    float* pfrh_pt2,
				    int*   pfrh_isSeed,
				    int*   pfrh_topoId,
				    int*   pfrh_layer,
				    int*   neigh8_Ind
				    ) {

   //int i = threadIdx.x+blockIdx.x*blockDim.x;

   //if(i<size) {        
   for (int i=0; i<size; i++) {
     if( ( pfrh_layer[i] == -1 && pfrh_energy[i]>seedEThresholdEB && pfrh_pt2[i]>seedPt2ThresholdEB) || ( (pfrh_layer[i] == -2) && pfrh_energy[i]>seedEThresholdEE && pfrh_pt2[i]>seedPt2hresholdEE) )
       {
	 pfrh_isSeed[i]=1;		 
	 for(int k=0; k<nNeigh; k++){
	   if(neigh8_Ind[nNeigh*i+k]<0) continue; 
	   if(pfrh_energy[i]<pfrh_energy[neigh8_Ind[nNeigh*i+k]]){
	     pfrh_isSeed[i]=0;
	     pfrh_topoId[i]=-1;	     
	     break;
	   }
	 }		
       }
     else{ 
       pfrh_topoId[i]=-1;
       pfrh_isSeed[i]=0;
       	    
     }     
   }
 }
  
  __global__ void topoKernel_ECAL( 
				  size_t size,
				  float* pfrh_energy,
				  int* pfrh_topoId,
				  int* pfrh_layer,
				  int* neigh8_Ind
				  ) {
    for(int j=0;j<16;j++){
      int l = threadIdx.x+blockIdx.x*blockDim.x;
      if(l<size) {
	//printf("layer: %d",pfrh_layer[l]);
	for(int k=0; k<nNeigh; k++){
	  if( neigh8_Ind[nNeigh*l+k] > -1 && 
	      pfrh_topoId[l] < pfrh_topoId[neigh8_Ind[nNeigh*l+k]] && 
	      ( (pfrh_layer[l] == -2 && pfrh_energy[l]>topoEThresholdEE) || 
		(pfrh_layer[l] == -1 && pfrh_energy[l]>topoEThresholdEB) ) )
	      {
		pfrh_topoId[l]=pfrh_topoId[neigh8_Ind[nNeigh*l+k]];
	      }
	}				       
      }//loop over neighbours end
      
    }//loop over neumann neighbourhood clustering end
  }
  
  
  __global__ void topoKernel_ECAL_serialize( 
				  size_t size,
				  float* pfrh_energy,
				  int* pfrh_topoId,
				  int* pfrh_layer,
				  int* neigh8_Ind
				  ) {
    for(int j=0;j<16;j++){
      //int l = threadIdx.x+blockIdx.x*blockDim.x;
      //if(l<size) {
      for (int l = 0; l < size; l++) {
        
    //printf("layer: %d",pfrh_layer[l]);
	for(int k=0; k<nNeigh; k++){
	  if( neigh8_Ind[nNeigh*l+k] > -1 && 
	      pfrh_topoId[l] < pfrh_topoId[neigh8_Ind[nNeigh*l+k]] && 
	      ( (pfrh_layer[l] == -2 && pfrh_energy[l]>topoEThresholdEE) || 
		(pfrh_layer[l] == -1 && pfrh_energy[l]>topoEThresholdEB) ) )
	      {
		pfrh_topoId[l]=pfrh_topoId[neigh8_Ind[nNeigh*l+k]];
	      }
	}				       
      }//loop over neighbours end
      
    }//loop over neumann neighbourhood clustering end
  }
  



__global__ void fastCluster_step1( size_t size,
					     float* pfrh_x,
					     float* pfrh_y,
					     float* pfrh_z,
					     float* pfrh_energy,
					     int* pfrh_topoId,
					     int* pfrh_isSeed,
					     int* pfrh_layer,
					     float* pcrhfrac, 
					     int* pcrhfracind,
					     float* fracSum,
					     int* rhCount
					     ) {

    int i = threadIdx.x+blockIdx.x*blockDim.x;
    int j = threadIdx.y+blockIdx.y*blockDim.y;
    //make sure topoID, Layer is the same, i is seed and j is not seed
    if( i<size && j<size){
      
      if( pfrh_topoId[i] == pfrh_topoId[j] && pfrh_isSeed[i]==1 ){

      float dist2 = 
	       (pfrh_x[i] - pfrh_x[j])*(pfrh_x[i] - pfrh_x[j])
	      +(pfrh_y[i] - pfrh_y[j])*(pfrh_y[i] - pfrh_y[j])
	      +(pfrh_z[i] - pfrh_z[j])*(pfrh_z[i] - pfrh_z[j]);	

      float d2 = dist2 / (showerSigma*showerSigma);
      float fraction = -1.;

      if(pfrh_layer[j] == -1) { fraction = pfrh_energy[i] / recHitEnergyNormEB * expf(-0.5 * d2); }
      if(pfrh_layer[j] == -2) { fraction = pfrh_energy[i] / recHitEnergyNormEE * expf(-0.5 * d2); }
      if(fraction==-1.) printf("FRACTION is NEGATIVE!!!");
      
      if( pfrh_isSeed[j]!=1 && d2<100)
	{ 
	  atomicAdd(&fracSum[j],fraction);	  
	}
      }
    }
  }


__global__ void fastCluster_step2( size_t size,
					     float* pfrh_x,
					     float* pfrh_y,
					     float* pfrh_z,
					     float* pfrh_energy,
					     int* pfrh_topoId,
					     int* pfrh_isSeed,
					     int* pfrh_layer,
					     float* pcrhfrac, 
					     int* pcrhfracind,
					     float* fracSum,
					     int* rhCount
					     ) {

    int i = threadIdx.x+blockIdx.x*blockDim.x;
    int j = threadIdx.y+blockIdx.y*blockDim.y;
    //make sure topoID, Layer is the same, i is seed and j is not seed
    if( i<size && j<size){
      if( pfrh_topoId[i] == pfrh_topoId[j] && pfrh_isSeed[i]==1 ){
      if(i==j) 
	{
	  pcrhfrac[i*maxSize]    = 1.;
	  pcrhfracind[i*maxSize] = j;
	}
      if( pfrh_isSeed[j]!=1 ){
	float dist2 = 
	   (pfrh_x[i] - pfrh_x[j])*(pfrh_x[i] - pfrh_x[j])
	  +(pfrh_y[i] - pfrh_y[j])*(pfrh_y[i] - pfrh_y[j])
	  +(pfrh_z[i] - pfrh_z[j])*(pfrh_z[i] - pfrh_z[j]);	

	float d2 = dist2 / (showerSigma*showerSigma);
	float fraction = -1.;

	if(pfrh_layer[j] == -1) { fraction = pfrh_energy[i] / recHitEnergyNormEB * expf(-0.5 * d2); }
	if(pfrh_layer[j] == -2) { fraction = pfrh_energy[i] / recHitEnergyNormEE * expf(-0.5 * d2); }
	if(fraction==-1.) printf("FRACTION is NEGATIVE!!!");
	if(d2 < 100. )
	  { 
	    if ((fraction/fracSum[j])>minFracToKeep){
	      int k = atomicAdd(&rhCount[i],1);
	      pcrhfrac[i*maxSize+k] = fraction/fracSum[j];
	      pcrhfracind[i*maxSize+k] = j;
	    }
	  }
      }
      }
    }        
}

__global__ void fastCluster_step1_serialize( size_t size,
					     float* pfrh_x,
					     float* pfrh_y,
					     float* pfrh_z,
					     float* pfrh_energy,
					     int* pfrh_topoId,
					     int* pfrh_isSeed,
					     int* pfrh_layer,
					     float* pcrhfrac, 
					     int* pcrhfracind,
					     float* fracSum,
					     int* rhCount
					     ) {

    //int i = threadIdx.x+blockIdx.x*blockDim.x;
    //int j = threadIdx.y+blockIdx.y*blockDim.y;
    //make sure topoID, Layer is the same, i is seed and j is not seed
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if( i<size && j<size){
              
              if( pfrh_topoId[i] == pfrh_topoId[j] && pfrh_isSeed[i]==1 ){

              float dist2 = 
                   (pfrh_x[i] - pfrh_x[j])*(pfrh_x[i] - pfrh_x[j])
                  +(pfrh_y[i] - pfrh_y[j])*(pfrh_y[i] - pfrh_y[j])
                  +(pfrh_z[i] - pfrh_z[j])*(pfrh_z[i] - pfrh_z[j]);	

              float d2 = dist2 / (showerSigma*showerSigma);
              float fraction = -1.;

              if(pfrh_layer[j] == -1) { fraction = pfrh_energy[i] / recHitEnergyNormEB * expf(-0.5 * d2); }
              if(pfrh_layer[j] == -2) { fraction = pfrh_energy[i] / recHitEnergyNormEE * expf(-0.5 * d2); }
              if(fraction==-1.) printf("FRACTION is NEGATIVE!!!");
              
              if( pfrh_isSeed[j]!=1 && d2<100)
            { 
              atomicAdd(&fracSum[j],fraction);	  
            }
              }
            }
        }
    }
  }


__global__ void fastCluster_step2_serialize( size_t size,
					     float* pfrh_x,
					     float* pfrh_y,
					     float* pfrh_z,
					     float* pfrh_energy,
					     int* pfrh_topoId,
					     int* pfrh_isSeed,
					     int* pfrh_layer,
					     float* pcrhfrac, 
					     int* pcrhfracind,
					     float* fracSum,
					     int* rhCount
					     ) {

    //int i = threadIdx.x+blockIdx.x*blockDim.x;
    //int j = threadIdx.y+blockIdx.y*blockDim.y;
    //make sure topoID, Layer is the same, i is seed and j is not seed
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
    
            if( i<size && j<size){
              if( pfrh_topoId[i] == pfrh_topoId[j] && pfrh_isSeed[i]==1 ){
              if(i==j) 
            {
              pcrhfrac[i*maxSize]    = 1.;
              pcrhfracind[i*maxSize] = j;
            }
              if( pfrh_isSeed[j]!=1 ){
            float dist2 = 
               (pfrh_x[i] - pfrh_x[j])*(pfrh_x[i] - pfrh_x[j])
              +(pfrh_y[i] - pfrh_y[j])*(pfrh_y[i] - pfrh_y[j])
              +(pfrh_z[i] - pfrh_z[j])*(pfrh_z[i] - pfrh_z[j]);	

            float d2 = dist2 / (showerSigma*showerSigma);
            float fraction = -1.;

            if(pfrh_layer[j] == -1) { fraction = pfrh_energy[i] / recHitEnergyNormEB * expf(-0.5 * d2); }
            if(pfrh_layer[j] == -2) { fraction = pfrh_energy[i] / recHitEnergyNormEE * expf(-0.5 * d2); }
            if(fraction==-1.) printf("FRACTION is NEGATIVE!!!");
            if(d2 < 100. )
              { 
                if ((fraction/fracSum[j])>minFracToKeep){
                  int k = atomicAdd(&rhCount[i],1);
                  pcrhfrac[i*maxSize+k] = fraction/fracSum[j];
                  pcrhfracind[i*maxSize+k] = j;
                }
              }
              }
              }
            }        
        }
    }
}
  

   

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
				int* pfrhfracind,
				int* pcrhfracind,
				float* pcrhfrac,
				float* fracSum,
				int* rhCount
				)
  { 
    //seeding
    if(size>0) seedingKernel_ECAL<<<(size+512-1)/512, 512>>>( size,  pfrh_energy,   pfrh_pt2,   pfrh_isSeed,  pfrh_topoId,  pfrh_layer,  neigh8_Ind);
    //cudaDeviceSynchronize();
      
    // for(int a=0;a<16;a++){
    if(size>0) topoKernel_ECAL<<<(size+512-1)/512, 512>>>( size, pfrh_energy,  pfrh_topoId,  pfrh_layer, neigh8_Ind);
    //}	    
    cudaDeviceSynchronize();

    dim3 grid( (size+32-1)/32, (size+32-1)/32 );
    dim3 block( 32, 32);

    //if(size>0) std::cout<<std::endl<<"NEW EVENT !!"<<std::endl<<std::endl;

     if(size>0) fastCluster_step1<<<grid, block>>>( size, pfrh_x,  pfrh_y,  pfrh_z,  pfrh_energy, pfrh_topoId,  pfrh_isSeed,  pfrh_layer, pcrhfrac, pcrhfracind, fracSum, rhCount);
     cudaDeviceSynchronize();

    if(size>0) fastCluster_step2<<<grid, block>>>( size, pfrh_x,  pfrh_y,  pfrh_z,  pfrh_energy, pfrh_topoId,  pfrh_isSeed,  pfrh_layer, pcrhfrac, pcrhfracind, fracSum, rhCount);  
    cudaDeviceSynchronize();
    
   
  }

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
				int* pfrhfracind,
				int* pcrhfracind,
				float* pcrhfrac,
				float* fracSum,
				int* rhCount
				)
  { 
    //seeding
    if(size>0) seedingKernel_ECAL_serialize<<<1,1>>>( size,  pfrh_energy,   pfrh_pt2,   pfrh_isSeed,  pfrh_topoId,  pfrh_layer,  neigh8_Ind);
    //cudaDeviceSynchronize();
      
    // for(int a=0;a<16;a++){
    if(size>0) topoKernel_ECAL_serialize<<<1,1>>>( size, pfrh_energy,  pfrh_topoId,  pfrh_layer, neigh8_Ind);
    //}	    
    cudaDeviceSynchronize();

    dim3 grid( (size+32-1)/32, (size+32-1)/32 );
    dim3 block( 32, 32);

    //if(size>0) std::cout<<std::endl<<"NEW EVENT !!"<<std::endl<<std::endl;

     if(size>0) fastCluster_step1_serialize<<<1,1>>>( size, pfrh_x,  pfrh_y,  pfrh_z,  pfrh_energy, pfrh_topoId,  pfrh_isSeed,  pfrh_layer, pcrhfrac, pcrhfracind, fracSum, rhCount);
     cudaDeviceSynchronize();

    if(size>0) fastCluster_step2_serialize<<<1,1>>>( size, pfrh_x,  pfrh_y,  pfrh_z,  pfrh_energy, pfrh_topoId,  pfrh_isSeed,  pfrh_layer, pcrhfrac, pcrhfracind, fracSum, rhCount);  
    cudaDeviceSynchronize();
    
   
  }
  
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
				int* pfrhfracind,
				int* pcrhfracind,
				float* pcrhfrac,
				float* fracSum,
				int* rhCount
				)
  { 
    //seeding
    if(size>0) seedingKernel_ECAL<<<(size+512-1)/512, 512>>>( size,  pfrh_energy,   pfrh_pt2,   pfrh_isSeed,  pfrh_topoId,  pfrh_layer,  neigh8_Ind);
    //cudaDeviceSynchronize();
      
    // for(int a=0;a<16;a++){
    if(size>0) topoKernel_ECAL_serialize<<<1,1>>>( size, pfrh_energy,  pfrh_topoId,  pfrh_layer, neigh8_Ind);
    //}	    
    cudaDeviceSynchronize();

    dim3 grid( (size+32-1)/32, (size+32-1)/32 );
    dim3 block( 32, 32);

    //if(size>0) std::cout<<std::endl<<"NEW EVENT !!"<<std::endl<<std::endl;

     if(size>0) fastCluster_step1_serialize<<<1,1>>>( size, pfrh_x,  pfrh_y,  pfrh_z,  pfrh_energy, pfrh_topoId,  pfrh_isSeed,  pfrh_layer, pcrhfrac, pcrhfracind, fracSum, rhCount);
     cudaDeviceSynchronize();

    if(size>0) fastCluster_step2_serialize<<<1,1>>>( size, pfrh_x,  pfrh_y,  pfrh_z,  pfrh_energy, pfrh_topoId,  pfrh_isSeed,  pfrh_layer, pcrhfrac, pcrhfracind, fracSum, rhCount);  
    cudaDeviceSynchronize();
    
   
  }
  
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
				int* pfrhfracind,
				int* pcrhfracind,
				float* pcrhfrac,
				float* fracSum,
				int* rhCount
				)
  { 
    //seeding
    if(size>0) seedingKernel_ECAL_serialize<<<1,1>>>( size,  pfrh_energy,   pfrh_pt2,   pfrh_isSeed,  pfrh_topoId,  pfrh_layer,  neigh8_Ind);
    //cudaDeviceSynchronize();
      
    // for(int a=0;a<16;a++){
    if(size>0) topoKernel_ECAL<<<(size+512-1)/512, 512>>>( size, pfrh_energy,  pfrh_topoId,  pfrh_layer, neigh8_Ind);
    //}	    
    cudaDeviceSynchronize();

    dim3 grid( (size+32-1)/32, (size+32-1)/32 );
    dim3 block( 32, 32);

    //if(size>0) std::cout<<std::endl<<"NEW EVENT !!"<<std::endl<<std::endl;

     if(size>0) fastCluster_step1_serialize<<<1,1>>>( size, pfrh_x,  pfrh_y,  pfrh_z,  pfrh_energy, pfrh_topoId,  pfrh_isSeed,  pfrh_layer, pcrhfrac, pcrhfracind, fracSum, rhCount);
     cudaDeviceSynchronize();

    if(size>0) fastCluster_step2_serialize<<<1,1>>>( size, pfrh_x,  pfrh_y,  pfrh_z,  pfrh_energy, pfrh_topoId,  pfrh_isSeed,  pfrh_layer, pcrhfrac, pcrhfracind, fracSum, rhCount);  
    cudaDeviceSynchronize();
    
   
  }
  
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
				int* pfrhfracind,
				int* pcrhfracind,
				float* pcrhfrac,
				float* fracSum,
				int* rhCount
				)
  { 
    //seeding
    if(size>0) seedingKernel_ECAL_serialize<<<1,1>>>( size,  pfrh_energy,   pfrh_pt2,   pfrh_isSeed,  pfrh_topoId,  pfrh_layer,  neigh8_Ind);
    //cudaDeviceSynchronize();
      
    // for(int a=0;a<16;a++){
    if(size>0) topoKernel_ECAL_serialize<<<1,1>>>( size, pfrh_energy,  pfrh_topoId,  pfrh_layer, neigh8_Ind); 
    //}	    
    cudaDeviceSynchronize();

    dim3 grid( (size+32-1)/32, (size+32-1)/32 );
    dim3 block( 32, 32);

    //if(size>0) std::cout<<std::endl<<"NEW EVENT !!"<<std::endl<<std::endl;
     if(size>0) fastCluster_step1<<<grid, block>>>( size, pfrh_x,  pfrh_y,  pfrh_z,  pfrh_energy, pfrh_topoId,  pfrh_isSeed,  pfrh_layer, pcrhfrac, pcrhfracind, fracSum, rhCount);
     cudaDeviceSynchronize();

     if(size>0) fastCluster_step2_serialize<<<1,1>>>( size, pfrh_x,  pfrh_y,  pfrh_z,  pfrh_energy, pfrh_topoId,  pfrh_isSeed,  pfrh_layer, pcrhfrac, pcrhfracind, fracSum, rhCount);  
     cudaDeviceSynchronize();
    
   
  }
  
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
				int* pfrhfracind,
				int* pcrhfracind,
				float* pcrhfrac,
				float* fracSum,
				int* rhCount
				)
  { 
    //seeding
    if(size>0) seedingKernel_ECAL_serialize<<<1,1>>>( size,  pfrh_energy,   pfrh_pt2,   pfrh_isSeed,  pfrh_topoId,  pfrh_layer,  neigh8_Ind);
    //cudaDeviceSynchronize();
      
    // for(int a=0;a<16;a++){
    if(size>0) topoKernel_ECAL_serialize<<<1,1>>>( size, pfrh_energy,  pfrh_topoId,  pfrh_layer, neigh8_Ind); 
    //}	    
    cudaDeviceSynchronize();

    dim3 grid( (size+32-1)/32, (size+32-1)/32 );
    dim3 block( 32, 32);

    //if(size>0) std::cout<<std::endl<<"NEW EVENT !!"<<std::endl<<std::endl;

     if(size>0) fastCluster_step1_serialize<<<1,1>>>( size, pfrh_x,  pfrh_y,  pfrh_z,  pfrh_energy, pfrh_topoId,  pfrh_isSeed,  pfrh_layer, pcrhfrac, pcrhfracind, fracSum, rhCount);
     cudaDeviceSynchronize();
     
     if(size>0) fastCluster_step2<<<grid, block>>>( size, pfrh_x,  pfrh_y,  pfrh_z,  pfrh_energy, pfrh_topoId,  pfrh_isSeed,  pfrh_layer, pcrhfrac, pcrhfracind, fracSum, rhCount);
     cudaDeviceSynchronize();
    
   
  }
}  // namespace cudavectors
