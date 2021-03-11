
#include <cmath>
#include <iostream>

// CUDA include files
#include <cuda_runtime.h>

// CMSSW include files
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"
#include "RecoParticleFlow/PFClusterProducer/plugins/PFClusterCudaHCAL.h"
#include <Eigen/Dense>

namespace PFClusterCudaHCAL {

  // THE ART OF HARDCODING
  // these numbers should be copied over during initialization
  __constant__ float showerSigma = 10;
  
  __constant__ float recHitEnergyNormEB_1 = 0.1;
  __constant__ float recHitEnergyNormEB_2 = 0.2;
  __constant__ float recHitEnergyNormEB_3 = 0.3;
  __constant__ float recHitEnergyNormEB_4 = 0.3;
  __constant__ float recHitEnergyNormEE_1 = 0.1;
  __constant__ float recHitEnergyNormEE_2_7 = 0.2;
  
  __constant__ float minFracToKeep = 0.0000001;

  __constant__ float seedEThresholdEB_1 = 0.125;
  __constant__ float seedEThresholdEB_2 = 0.25;
  __constant__ float seedEThresholdEB_3 = 0.35;
  __constant__ float seedEThresholdEB_4 = 0.35;
  __constant__ float seedEThresholdEE_1 = 0.1375;
  __constant__ float seedEThresholdEE_2_7 = 0.275;

  __constant__ float seedPt2ThresholdEB = 0.;
  __constant__ float seedPt2hresholdEE = 0.;

  __constant__ float topoEThresholdEB_1 = 0.1;
  __constant__ float topoEThresholdEB_2 = 0.2;
  __constant__ float topoEThresholdEB_3 = 0.3;
  __constant__ float topoEThresholdEB_4 = 0.3;
  __constant__ float topoEThresholdEE_1 = 0.1;
  __constant__ float topoEThresholdEE_2_7 = 0.2;
  __constant__ float topoEThresholdEB_vec[4] = {0.1,0.2,0.3,0.3};
  __constant__ float topoEThresholdEE_vec[7] = {0.1,0.2,0.2,0.2,0.2,0.2,0.2};
  

  __constant__ int nNeighTopo = 8;
  __constant__ int nNT = 8;
  __constant__ int nNeigh = 4;
  __constant__ int maxSize = 100;
  
  //int nTopoLoops = 100; // 35;
  int nTopoLoops = 35;


 __global__ void seedingKernel_HCAL(
     				    size_t size, 
				    double* pfrh_energy,
				    double* pfrh_pt2,
				    int*   pfrh_isSeed,
				    int*   pfrh_topoId,
				    int*   pfrh_layer,
				    int*   pfrh_depth,
				    int*   neigh4_Ind
				    ) {

   int i = threadIdx.x+blockIdx.x*blockDim.x;

   if(i<size) {        
     if( ( pfrh_layer[i] == 1 && 
	   pfrh_depth[i] == 1 &&
	   pfrh_energy[i]>seedEThresholdEB_1 && 
	   pfrh_pt2[i]>seedPt2ThresholdEB ) ||
	 ( pfrh_layer[i] == 1 && 
	   pfrh_depth[i] == 2 &&
	   pfrh_energy[i]>seedEThresholdEB_2 && 
	   pfrh_pt2[i]>seedPt2ThresholdEB ) ||
	 ( pfrh_layer[i] == 1 && 
	   pfrh_depth[i] == 3 &&
	   pfrh_energy[i]>seedEThresholdEB_3 && 
	   pfrh_pt2[i]>seedPt2ThresholdEB ) ||
	 ( pfrh_layer[i] == 1 && 
	   pfrh_depth[i] == 4 &&
	   pfrh_energy[i]>seedEThresholdEB_4 && 
	   pfrh_pt2[i]>seedPt2ThresholdEB ) ||
	 ( pfrh_layer[i] == 3  && 
	   pfrh_depth[i] == 1  &&
	   pfrh_energy[i]>seedEThresholdEE_1 && 
	   pfrh_pt2[i]>seedPt2hresholdEE)   ||
	 ( pfrh_layer[i] == 3  && 
	   pfrh_depth[i] > 1   &&
	   pfrh_energy[i]>seedEThresholdEE_2_7 && 
	   pfrh_pt2[i]>seedPt2hresholdEE))
       {
	 pfrh_isSeed[i]=1;		 
	 for(int k=0; k<nNeigh; k++){
	   if(neigh4_Ind[nNeigh*i+k]<0) continue; 
	   if(pfrh_energy[i]<pfrh_energy[neigh4_Ind[nNeigh*i+k]]){
	     pfrh_isSeed[i]=0;
	     //pfrh_topoId[i]=-1;	     
	     break;
	   }
	 }		
       }
     else{ 
       // pfrh_topoId[i]=-1;
       pfrh_isSeed[i]=0;
       	    
     }     
   }
 }
  
  
 __global__ void seedingKernel_HCAL_serialize(
     				    size_t size, 
				    double* pfrh_energy,
				    double* pfrh_pt2,
				    int*   pfrh_isSeed,
				    int*   pfrh_topoId,
				    int*   pfrh_layer,
				    int*   pfrh_depth,
				    int*   neigh4_Ind
				    ) {

   //int i = threadIdx.x+blockIdx.x*blockDim.x;
   for (int i = 0; i < size; i++) {
       if(i<size) {        
         if( ( pfrh_layer[i] == 1 && 
           pfrh_depth[i] == 1 &&
           pfrh_energy[i]>seedEThresholdEB_1 && 
           pfrh_pt2[i]>seedPt2ThresholdEB ) ||
         ( pfrh_layer[i] == 1 && 
           pfrh_depth[i] == 2 &&
           pfrh_energy[i]>seedEThresholdEB_2 && 
           pfrh_pt2[i]>seedPt2ThresholdEB ) ||
         ( pfrh_layer[i] == 1 && 
           pfrh_depth[i] == 3 &&
           pfrh_energy[i]>seedEThresholdEB_3 && 
           pfrh_pt2[i]>seedPt2ThresholdEB ) ||
         ( pfrh_layer[i] == 1 && 
           pfrh_depth[i] == 4 &&
           pfrh_energy[i]>seedEThresholdEB_4 && 
           pfrh_pt2[i]>seedPt2ThresholdEB ) ||
         ( pfrh_layer[i] == 3  && 
           pfrh_depth[i] == 1  &&
           pfrh_energy[i]>seedEThresholdEE_1 && 
           pfrh_pt2[i]>seedPt2hresholdEE)   ||
         ( pfrh_layer[i] == 3  && 
           pfrh_depth[i] > 1   &&
           pfrh_energy[i]>seedEThresholdEE_2_7 && 
           pfrh_pt2[i]>seedPt2hresholdEE))
           {
         pfrh_isSeed[i]=1;		 
         for(int k=0; k<nNeigh; k++){
           if(neigh4_Ind[nNeigh*i+k]<0) continue; 
           if(pfrh_energy[i]<pfrh_energy[neigh4_Ind[nNeigh*i+k]]){
             pfrh_isSeed[i]=0;
             //pfrh_topoId[i]=-1;	     
             break;
           }
         }		
           }
         else{ 
           // pfrh_topoId[i]=-1;
           pfrh_isSeed[i]=0;
                
         }     
       }
    }
 }
  

   __global__ void topoKernel_HCALV2( 
				  size_t size,
				  double* pfrh_energy,
				  int* pfrh_topoId,
				  int* pfrh_layer,
				  int* pfrh_depth,
				  int* neigh8_Ind
				  ) {
     
     int l = threadIdx.x+blockIdx.x*blockDim.x;
     //int k = threadIdx.y+blockIdx.y*blockDim.y;
     int k = (threadIdx.y+blockIdx.y*blockDim.y) % nNT;
           
      if(l<size && k<nNT) {

	while( neigh8_Ind[nNT*l+k] > -1 && pfrh_topoId[l] != pfrh_topoId[neigh8_Ind[nNT*l+k]] && 
	    ( (pfrh_layer[neigh8_Ind[nNT*l+k]] == 3 && pfrh_energy[neigh8_Ind[nNT*l+k]]>topoEThresholdEE_vec[pfrh_depth[neigh8_Ind[nNT*l+k]]-1]) ||
	      (pfrh_layer[neigh8_Ind[nNT*l+k]] == 1 && pfrh_energy[neigh8_Ind[nNT*l+k]]>topoEThresholdEB_vec[pfrh_depth[neigh8_Ind[nNT*l+k]]-1]) ) &&
	    ( (pfrh_layer[l] == 3 && pfrh_energy[l]>topoEThresholdEE_vec[pfrh_depth[l]-1]) ||
	      (pfrh_layer[l] == 1 && pfrh_energy[l]>topoEThresholdEB_vec[pfrh_depth[l]-1]))
	    )
	    {
	      if(pfrh_topoId[l] > pfrh_topoId[neigh8_Ind[nNT*l+k]]){
		atomicMax(&pfrh_topoId[neigh8_Ind[nNT*l+k]],pfrh_topoId[l]);
	      }
	      if(pfrh_topoId[l] < pfrh_topoId[neigh8_Ind[nNT*l+k]]){
		atomicMax(&pfrh_topoId[l], pfrh_topoId[neigh8_Ind[nNT*l+k]]);
	      }
	    }	
      }
   }
 

   __global__ void topoKernel_HCAL_serialize( 
				  size_t size,
				  double* pfrh_energy,
				  int* pfrh_topoId,
				  int* pfrh_layer,
				  int* pfrh_depth,
				  int* neigh8_Ind
				  ) {
     
     //int l = threadIdx.x+blockIdx.x*blockDim.x;
     //int k = threadIdx.y+blockIdx.y*blockDim.y;
     
     for (int l = 0; l < size; l++) {
        //for (int k = 0; k < size; k++) {
        for (int k = 0; k < 8; k++) {
           
            if(l<size && k<nNT) {

            while( neigh8_Ind[nNT*l+k] > -1 && pfrh_topoId[l] != pfrh_topoId[neigh8_Ind[nNT*l+k]] && 
                ( (pfrh_layer[neigh8_Ind[nNT*l+k]] == 3 && pfrh_energy[neigh8_Ind[nNT*l+k]]>topoEThresholdEE_vec[pfrh_depth[neigh8_Ind[nNT*l+k]]-1]) ||
                  (pfrh_layer[neigh8_Ind[nNT*l+k]] == 1 && pfrh_energy[neigh8_Ind[nNT*l+k]]>topoEThresholdEB_vec[pfrh_depth[neigh8_Ind[nNT*l+k]]-1]) ) &&
                ( (pfrh_layer[l] == 3 && pfrh_energy[l]>topoEThresholdEE_vec[pfrh_depth[l]-1]) ||
                  (pfrh_layer[l] == 1 && pfrh_energy[l]>topoEThresholdEB_vec[pfrh_depth[l]-1]))
                )
                {
                  if(pfrh_topoId[l] > pfrh_topoId[neigh8_Ind[nNT*l+k]]){
                atomicMax(&pfrh_topoId[neigh8_Ind[nNT*l+k]],pfrh_topoId[l]);
                  }
                  if(pfrh_topoId[l] < pfrh_topoId[neigh8_Ind[nNT*l+k]]){
                atomicMax(&pfrh_topoId[l], pfrh_topoId[neigh8_Ind[nNT*l+k]]);
                  }
                }	
              }
         }
      }
   }
 
__global__ void hcalFastCluster_step1( size_t size,
					     float* pfrh_x,
					     float* pfrh_y,
					     float* pfrh_z,
					     double* pfrh_energy,
					     int* pfrh_topoId,
					     int* pfrh_isSeed,
					     int* pfrh_layer,
				             int* pfrh_depth,
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

      if(pfrh_layer[j] == 1 && pfrh_depth[j] == 1) { fraction = pfrh_energy[i] / recHitEnergyNormEB_1 * expf(-0.5 * d2); }
      if(pfrh_layer[j] == 1 && pfrh_depth[j] == 2) { fraction = pfrh_energy[i] / recHitEnergyNormEB_2 * expf(-0.5 * d2); }
      if(pfrh_layer[j] == 1 && pfrh_depth[j] == 3) { fraction = pfrh_energy[i] / recHitEnergyNormEB_3 * expf(-0.5 * d2); }
      if(pfrh_layer[j] == 1 && pfrh_depth[j] == 4) { fraction = pfrh_energy[i] / recHitEnergyNormEB_4 * expf(-0.5 * d2); }
      if(pfrh_layer[j] == 3 && pfrh_depth[j] == 1) { fraction = pfrh_energy[i] / recHitEnergyNormEE_1 * expf(-0.5 * d2); }
      if(pfrh_layer[j] == 3 && pfrh_depth[j] > 1 ) { fraction = pfrh_energy[i] / recHitEnergyNormEE_2_7 * expf(-0.5 * d2); }
	  
      if(fraction==-1.) printf("FRACTION is NEGATIVE!!!");

      if( pfrh_isSeed[j]!=1 && d2<100.)
	{
	  atomicAdd(&fracSum[j],fraction);
	}
      }
    }
  }

 
__global__ void hcalFastCluster_step1_serialize( size_t size,
					     float* pfrh_x,
					     float* pfrh_y,
					     float* pfrh_z,
					     double* pfrh_energy,
					     int* pfrh_topoId,
					     int* pfrh_isSeed,
					     int* pfrh_layer,
				             int* pfrh_depth,
					     float* pcrhfrac,
					     int* pcrhfracind,
					     float* fracSum,
					     int* rhCount
					     ) {

    //int i = threadIdx.x+blockIdx.x*blockDim.x;
    //int j = threadIdx.y+blockIdx.y*blockDim.y;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            //make sure topoID, Layer is the same, i is seed and j is not seed
            if( i<size && j<size){

              if( pfrh_topoId[i] == pfrh_topoId[j] && pfrh_isSeed[i]==1 ){

              float dist2 =
                   (pfrh_x[i] - pfrh_x[j])*(pfrh_x[i] - pfrh_x[j])
                  +(pfrh_y[i] - pfrh_y[j])*(pfrh_y[i] - pfrh_y[j])
                  +(pfrh_z[i] - pfrh_z[j])*(pfrh_z[i] - pfrh_z[j]);

              float d2 = dist2 / (showerSigma*showerSigma);
              float fraction = -1.;

              if(pfrh_layer[j] == 1 && pfrh_depth[j] == 1) { fraction = pfrh_energy[i] / recHitEnergyNormEB_1 * expf(-0.5 * d2); }
              if(pfrh_layer[j] == 1 && pfrh_depth[j] == 2) { fraction = pfrh_energy[i] / recHitEnergyNormEB_2 * expf(-0.5 * d2); }
              if(pfrh_layer[j] == 1 && pfrh_depth[j] == 3) { fraction = pfrh_energy[i] / recHitEnergyNormEB_3 * expf(-0.5 * d2); }
              if(pfrh_layer[j] == 1 && pfrh_depth[j] == 4) { fraction = pfrh_energy[i] / recHitEnergyNormEB_4 * expf(-0.5 * d2); }
              if(pfrh_layer[j] == 3 && pfrh_depth[j] == 1) { fraction = pfrh_energy[i] / recHitEnergyNormEE_1 * expf(-0.5 * d2); }
              if(pfrh_layer[j] == 3 && pfrh_depth[j] > 1 ) { fraction = pfrh_energy[i] / recHitEnergyNormEE_2_7 * expf(-0.5 * d2); }
              
              if(fraction==-1.) printf("FRACTION is NEGATIVE!!!");

              if( pfrh_isSeed[j]!=1 && d2<100.)
            {
              atomicAdd(&fracSum[j],fraction);
            }
              }
            }
        }
    }
  }

__global__ void hcalFastCluster_step2( size_t size,
					     float* pfrh_x,
					     float* pfrh_y,
					     float* pfrh_z,
					     double* pfrh_energy,
					     int* pfrh_topoId,
					     int* pfrh_isSeed,
					     int* pfrh_layer,
				             int* pfrh_depth,
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

	if(pfrh_layer[j] == 1 && pfrh_depth[j] == 1) { fraction = pfrh_energy[i] / recHitEnergyNormEB_1 * expf(-0.5 * d2); }
	if(pfrh_layer[j] == 1 && pfrh_depth[j] == 2) { fraction = pfrh_energy[i] / recHitEnergyNormEB_2 * expf(-0.5 * d2); }
	if(pfrh_layer[j] == 1 && pfrh_depth[j] == 3) { fraction = pfrh_energy[i] / recHitEnergyNormEB_3 * expf(-0.5 * d2); }
	if(pfrh_layer[j] == 1 && pfrh_depth[j] == 4) { fraction = pfrh_energy[i] / recHitEnergyNormEB_4 * expf(-0.5 * d2); }
	if(pfrh_layer[j] == 3 && pfrh_depth[j] == 1) { fraction = pfrh_energy[i] / recHitEnergyNormEE_1 * expf(-0.5 * d2); }
	if(pfrh_layer[j] == 3 && pfrh_depth[j] > 1 ) { fraction = pfrh_energy[i] / recHitEnergyNormEE_2_7 * expf(-0.5 * d2); }
	  
	
	if(fraction==-1.) printf("FRACTION is NEGATIVE!!!");
	if(d2 < 100. )
	  {
	    if ((fraction/fracSum[j])>minFracToKeep){
	      int k = atomicAdd(&rhCount[i],1);
	      pcrhfrac[i*maxSize+k] = fraction/fracSum[j];
	      pcrhfracind[i*maxSize+k] = j;
	      //printf("(i,j)=(%d,%d), rhCount=%d, fraction=%f, fracsum=%f \n",i,j,rhCount[i], fraction, fracSum[j]);
	    }
	  }
      }
      }
    }
}


__global__ void hcalFastCluster_step2_serialize( size_t size,
					     float* pfrh_x,
					     float* pfrh_y,
					     float* pfrh_z,
					     double* pfrh_energy,
					     int* pfrh_topoId,
					     int* pfrh_isSeed,
					     int* pfrh_layer,
				             int* pfrh_depth,
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

            if(pfrh_layer[j] == 1 && pfrh_depth[j] == 1) { fraction = pfrh_energy[i] / recHitEnergyNormEB_1 * expf(-0.5 * d2); }
            if(pfrh_layer[j] == 1 && pfrh_depth[j] == 2) { fraction = pfrh_energy[i] / recHitEnergyNormEB_2 * expf(-0.5 * d2); }
            if(pfrh_layer[j] == 1 && pfrh_depth[j] == 3) { fraction = pfrh_energy[i] / recHitEnergyNormEB_3 * expf(-0.5 * d2); }
            if(pfrh_layer[j] == 1 && pfrh_depth[j] == 4) { fraction = pfrh_energy[i] / recHitEnergyNormEB_4 * expf(-0.5 * d2); }
            if(pfrh_layer[j] == 3 && pfrh_depth[j] == 1) { fraction = pfrh_energy[i] / recHitEnergyNormEE_1 * expf(-0.5 * d2); }
            if(pfrh_layer[j] == 3 && pfrh_depth[j] > 1 ) { fraction = pfrh_energy[i] / recHitEnergyNormEE_2_7 * expf(-0.5 * d2); }
              
            
            if(fraction==-1.) printf("FRACTION is NEGATIVE!!!");
            if(d2 < 100. )
              {
                if ((fraction/fracSum[j])>minFracToKeep){
                  int k = atomicAdd(&rhCount[i],1);
                  pcrhfrac[i*maxSize+k] = fraction/fracSum[j];
                  pcrhfracind[i*maxSize+k] = j;
                  //printf("(i,j)=(%d,%d), rhCount=%d, fraction=%f, fracsum=%f \n",i,j,rhCount[i], fraction, fracSum[j]);
                }
              }
              }
              }
            }
        }
    }
}



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
				
				int* pcrhfracind,
				float* pcrhfrac,
				float* fracSum,
				int* rhCount
				)
  { 
    //seeding
    if(size>0) seedingKernel_HCAL<<<(size+512-1)/512, 512>>>( size,  pfrh_energy,   pfrh_pt2,   pfrh_isSeed,  pfrh_topoId,  pfrh_layer,pfrh_depth,  neigh4_Ind);
    //cudaDeviceSynchronize();
    
    //topoclustering 
      
      //dim3 gridT( (size+64-1)/64, 1 );
      //dim3 blockT( 64, 8);
      dim3 gridT( (size+64-1)/64, 8 );
      dim3 blockT( 64, 16); // 16 threads in a half-warp
      for(int h=0;h<nTopoLoops; h++){
    
      if(size>0) topoKernel_HCALV2<<<gridT, blockT>>>( size, pfrh_energy,  pfrh_topoId,  pfrh_layer, pfrh_depth, neigh8_Ind);	     
      }
     cudaDeviceSynchronize();
    
    

      dim3 grid( (size+32-1)/32, (size+32-1)/32 );
      dim3 block( 32, 32);

      if(size>0) hcalFastCluster_step1<<<grid, block>>>( size, pfrh_x,  pfrh_y,  pfrh_z,  pfrh_energy, pfrh_topoId,  pfrh_isSeed,  pfrh_layer, pfrh_depth, pcrhfrac, pcrhfracind, fracSum, rhCount);
     //cudaDeviceSynchronize();

      if(size>0) hcalFastCluster_step2<<<grid, block>>>( size, pfrh_x,  pfrh_y,  pfrh_z,  pfrh_energy, pfrh_topoId,  pfrh_isSeed,  pfrh_layer, pfrh_depth, pcrhfrac, pcrhfracind, fracSum, rhCount);


  }

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
				
				int* pcrhfracind,
				float* pcrhfrac,
				float* fracSum,
				int* rhCount
				)
  { 
    //seeding
    if(size>0) seedingKernel_HCAL_serialize<<<1,1>>>( size,  pfrh_energy,   pfrh_pt2,   pfrh_isSeed,  pfrh_topoId,  pfrh_layer,pfrh_depth,  neigh4_Ind);
    //cudaDeviceSynchronize();
    
    //topoclustering 
      
      //dim3 gridT( (size+64-1)/64, 1 );
      //dim3 blockT( 64, 8);
      for(int h=0;h<nTopoLoops; h++){
    
      if(size>0) topoKernel_HCAL_serialize<<<1,1>>>( size, pfrh_energy,  pfrh_topoId,  pfrh_layer, pfrh_depth, neigh8_Ind);	     
      }
     cudaDeviceSynchronize();
    
    

      dim3 grid( (size+32-1)/32, (size+32-1)/32 );
      dim3 block( 32, 32);

      if(size>0) hcalFastCluster_step1_serialize<<<1,1>>>( size, pfrh_x,  pfrh_y,  pfrh_z,  pfrh_energy, pfrh_topoId,  pfrh_isSeed,  pfrh_layer, pfrh_depth, pcrhfrac, pcrhfracind, fracSum, rhCount);
     //cudaDeviceSynchronize();

      if(size>0) hcalFastCluster_step2_serialize<<<1,1>>>( size, pfrh_x,  pfrh_y,  pfrh_z,  pfrh_energy, pfrh_topoId,  pfrh_isSeed,  pfrh_layer, pfrh_depth, pcrhfrac, pcrhfracind, fracSum, rhCount);


  }

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
				
				int* pcrhfracind,
				float* pcrhfrac,
				float* fracSum,
				int* rhCount
				)
  { 
    //seeding
    if(size>0) seedingKernel_HCAL<<<(size+512-1)/512, 512>>>( size,  pfrh_energy,   pfrh_pt2,   pfrh_isSeed,  pfrh_topoId,  pfrh_layer,pfrh_depth,  neigh4_Ind);
    //cudaDeviceSynchronize();
    
    //topoclustering 
      
      dim3 gridT( (size+64-1)/64, 1 );
      dim3 blockT( 64, 8);
      for(int h=0;h<nTopoLoops; h++){
      if(size>0) topoKernel_HCAL_serialize<<<1,1>>>( size, pfrh_energy,  pfrh_topoId,  pfrh_layer, pfrh_depth, neigh8_Ind);
      }
     cudaDeviceSynchronize();
    
    

      dim3 grid( (size+32-1)/32, (size+32-1)/32 );
      dim3 block( 32, 32);

      if(size>0) hcalFastCluster_step1_serialize<<<1,1>>>( size, pfrh_x,  pfrh_y,  pfrh_z,  pfrh_energy, pfrh_topoId,  pfrh_isSeed,  pfrh_layer, pfrh_depth, pcrhfrac, pcrhfracind, fracSum, rhCount);
     //cudaDeviceSynchronize();

      if(size>0) hcalFastCluster_step2_serialize<<<1,1>>>( size, pfrh_x,  pfrh_y,  pfrh_z,  pfrh_energy, pfrh_topoId,  pfrh_isSeed,  pfrh_layer, pfrh_depth, pcrhfrac, pcrhfracind, fracSum, rhCount);


  }

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
				
				int* pcrhfracind,
				float* pcrhfrac,
				float* fracSum,
				int* rhCount
				)
  { 
    //seeding
    if(size>0) seedingKernel_HCAL_serialize<<<1,1>>>( size,  pfrh_energy,   pfrh_pt2,   pfrh_isSeed,  pfrh_topoId,  pfrh_layer,pfrh_depth,  neigh4_Ind);
    //cudaDeviceSynchronize();
    
    //topoclustering 
      
      dim3 gridT( (size+64-1)/64, 1 );
      dim3 blockT( 64, 8);
      for(int h=0;h<nTopoLoops; h++){
   
      if(size>0) topoKernel_HCALV2<<<gridT, blockT>>>( size, pfrh_energy,  pfrh_topoId,  pfrh_layer, pfrh_depth, neigh8_Ind);
      }
     cudaDeviceSynchronize();
    
    

      dim3 grid( (size+32-1)/32, (size+32-1)/32 );
      dim3 block( 32, 32);

      if(size>0) hcalFastCluster_step1_serialize<<<1,1>>>( size, pfrh_x,  pfrh_y,  pfrh_z,  pfrh_energy, pfrh_topoId,  pfrh_isSeed,  pfrh_layer, pfrh_depth, pcrhfrac, pcrhfracind, fracSum, rhCount);
     //cudaDeviceSynchronize();

      if(size>0) hcalFastCluster_step2_serialize<<<1,1>>>( size, pfrh_x,  pfrh_y,  pfrh_z,  pfrh_energy, pfrh_topoId,  pfrh_isSeed,  pfrh_layer, pfrh_depth, pcrhfrac, pcrhfracind, fracSum, rhCount);


  }

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
				
				int* pcrhfracind,
				float* pcrhfrac,
				float* fracSum,
				int* rhCount
				)
  { 
    //seeding
    if(size>0) seedingKernel_HCAL_serialize<<<1,1>>>( size,  pfrh_energy,   pfrh_pt2,   pfrh_isSeed,  pfrh_topoId,  pfrh_layer,pfrh_depth,  neigh4_Ind);
    //cudaDeviceSynchronize();
    
    //topoclustering 
      
      dim3 gridT( (size+64-1)/64, 1 );
      dim3 blockT( 64, 8);
      for(int h=0;h<nTopoLoops; h++){
      if(size>0) topoKernel_HCAL_serialize<<<1,1>>>( size, pfrh_energy,  pfrh_topoId,  pfrh_layer, pfrh_depth, neigh8_Ind); 
      }
     cudaDeviceSynchronize();
    
    

      dim3 grid( (size+32-1)/32, (size+32-1)/32 );
      dim3 block( 32, 32);

      if(size>0) hcalFastCluster_step1<<<grid, block>>>( size, pfrh_x,  pfrh_y,  pfrh_z,  pfrh_energy, pfrh_topoId,  pfrh_isSeed,  pfrh_layer, pfrh_depth, pcrhfrac, pcrhfracind, fracSum, rhCount);
     //cudaDeviceSynchronize();

      if(size>0) hcalFastCluster_step2_serialize<<<1,1>>>( size, pfrh_x,  pfrh_y,  pfrh_z,  pfrh_energy, pfrh_topoId,  pfrh_isSeed,  pfrh_layer, pfrh_depth, pcrhfrac, pcrhfracind, fracSum, rhCount);


  }

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
				
				int* pcrhfracind,
				float* pcrhfrac,
				float* fracSum,
				int* rhCount
				)
  { 
    //seeding
    if(size>0) seedingKernel_HCAL_serialize<<<1,1>>>( size,  pfrh_energy,   pfrh_pt2,   pfrh_isSeed,  pfrh_topoId,  pfrh_layer,pfrh_depth,  neigh4_Ind);
    //cudaDeviceSynchronize();
    
    //topoclustering 
      
      dim3 gridT( (size+64-1)/64, 1 );
      dim3 blockT( 64, 8);
      for(int h=0;h<nTopoLoops; h++){
      if(size>0) topoKernel_HCAL_serialize<<<1,1>>>( size, pfrh_energy,  pfrh_topoId,  pfrh_layer, pfrh_depth, neigh8_Ind); 
      }
     cudaDeviceSynchronize();
    
    

      dim3 grid( (size+32-1)/32, (size+32-1)/32 );
      dim3 block( 32, 32);

      if(size>0) hcalFastCluster_step1_serialize<<<1,1>>>( size, pfrh_x,  pfrh_y,  pfrh_z,  pfrh_energy, pfrh_topoId,  pfrh_isSeed,  pfrh_layer, pfrh_depth, pcrhfrac, pcrhfracind, fracSum, rhCount);
     //cudaDeviceSynchronize();
      
      if(size>0) hcalFastCluster_step2<<<grid, block>>>( size, pfrh_x,  pfrh_y,  pfrh_z,  pfrh_energy, pfrh_topoId,  pfrh_isSeed,  pfrh_layer, pfrh_depth, pcrhfrac, pcrhfracind, fracSum, rhCount);

  }

} // namespace cudavectors
