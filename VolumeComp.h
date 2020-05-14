#ifndef VOLUMECOMP_H_
#define VOLUMECOMP_H_ 

#include "SystemStructures.h"

void ComputeVolume(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    LJInfoVecs& ljInfoVecs);
    
struct VolumeCompFunctor {
    
    double spring_constant;
    double* locXAddr;
    double* locYAddr;
    double* locZAddr;

	__host__ __device__ VolumeCompFunctor(
        double& _spring_constant,
        double* _locXAddr,
        double* _locYAddr,
        double* _locZAddr):

        spring_constant(_spring_constant),
        locXAddr(_locXAddr),
        locYAddr(_locYAddr),
        locZAddr(_locZAddr){}

	//hand in counting iterator and id of two edges and preferred length
	__device__ double operator()(const Tuuuu &u4) {
        		
        //counter ranges from 0 to num_edges. 
        int counter = thrust::get<0>(u4);
        int r1 = thrust::get<1>(u4);
        int r2 = thrust::get<2>(u4);
        int r3 = thrust::get<3>(u4);

        if (r1 != INT_MAX && r2 != INT_MAX && r3 != INT_MAX){
        double r1x = locXAddr[r1];
        double r1y = locYAddr[r1];
        double r1z = locZAddr[r1];
        double r2x = locXAddr[r2];
        double r2y = locYAddr[r2];
        double r2z = locZAddr[r2];
        double r3x = locXAddr[r3];
        double r3y = locYAddr[r3];
        double r3z = locZAddr[r3];

        //(v1x + v2x + v3x)*((v2y-v1y)*(v3z-v1z) - (v3y-v2y)*(v2z-v1z))+
        //(v1y + v2y + v3y)*(-(v2x-v1x)*(v3z-v2z) + (v3x-v2x)*(v2z-v1z))+
        //(v1z + v2z + v3z)*((v2x-v1x)*(v3y-v2y) - (v3x-v2x)*(v2y-v1y));

        double N1 = (r2y - r1y)*(r3z - r1z) - (r3y - r1y)*(r2z - r1z);
        double N2 = -(r2x - r1x)*(r3z - r1z) + (r3x - r1x)*(r2z - r1z);
        double N3 = (r2x - r1x)*(r3y - r1y) - (r3x - r1x)*(r2y - r1y);
        double normN = sqrt(N1*N1 + N2*N2 + N3*N3);
        N1 = N1/normN;
        N2 = N2/normN;
        N3 = N3/normN;
        
        double r1_dot_N = r1x*N1 + r1y*N2 + r1z*N3;
        double r1cr2x = r1y*r2z - r2y*r1z;
        double r1cr2y = -r1x*r2z + r2x*r1z;
        double r1cr2z = r1x*r2y - r2x*r1y;
        double r2cr3x = r2y*r3z - r3y*r2z;
        double r2cr3y = -r2x*r3z + r3x*r2z;
        double r2cr3z = r2x*r3y - r3x*r2y;
        double r3cr1x = r3y*r1z - r1y*r3z;
        double r3cr1y = -r3x*r1z + r1x*r3z;
        double r3cr1z = r3x*r1y - r1x*r3y;

        double NN = N1*(r1cr2x + r2cr3x + r3cr1x) + N2*(r1cr2y + r2cr3y + r3cr1y) + N3*(r1cr2z + r2cr3z + r3cr1z);
        double volume = r1_dot_N*sqrt(NN*NN);
        //if (normN != 0.0){
         //   volume = 1.0;   
        //}

        return volume;
        
        }
        else{
            double volume = 0.0;
            
            return volume;
        }


    }
};

#endif