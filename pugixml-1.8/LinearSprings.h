#ifndef LINEARSPRINGS_H_
#define LINEARSPRINGS_H_ 

#include "SystemStructures.h"

void ComputeLinearSprings(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    LJInfoVecs& ljInfoVecs);
    
struct LinearSpringFunctor {
    
    double spring_constant;
    double* locXAddr;
    double* locYAddr;
    double* locZAddr;

    unsigned* idKey;
    double* forceXAddr;
    double* forceYAddr;
    double* forceZAddr;
    
	__host__ __device__ LinearSpringFunctor(
        double& _spring_constant,
        double* _locXAddr,
        double* _locYAddr,
        double* _locZAddr,
        
        unsigned* _idKey,
        double* _forceXAddr,
        double* _forceYAddr,
        double* _forceZAddr) :
        spring_constant(_spring_constant),
        locXAddr(_locXAddr),
        locYAddr(_locYAddr),
        locZAddr(_locZAddr),
        idKey(_idKey),
        forceXAddr(_forceXAddr),
        forceYAddr(_forceYAddr),
        forceZAddr(_forceZAddr) {}

	//hand in counting iterator and id of two edges and preferred length
	__device__ double operator()(const Tuuud &u3d) {
        		
        //counter ranges from 0 to num_edges. 
        unsigned counter = thrust::get<0>(u3d);
		unsigned place = 2 * counter;//represents location in write to vector.

        unsigned edgeL = thrust::get<1>(u3d);
        unsigned edgeR = thrust::get<2>(u3d);
        double length_zero = thrust::get<3>(u3d);

        // compute forces.
        double xLoc_LR = locXAddr[edgeL] - locXAddr[edgeR];
        double yLoc_LR = locYAddr[edgeL] - locYAddr[edgeR];
        double zLoc_LR = locZAddr[edgeL] - locZAddr[edgeR];
   

        double length_current = sqrt( (xLoc_LR) * (xLoc_LR) + 
                                    (yLoc_LR) * (yLoc_LR)  + 
                                    (zLoc_LR) * (zLoc_LR) );
        

        double magnitude = -spring_constant * (length_current - length_zero);
        
        
            idKey[place] = edgeL;

            
        //issue here writing to vectors with force????
            forceXAddr[place] = magnitude * (xLoc_LR/length_current);
            forceYAddr[place] = magnitude * (yLoc_LR/length_current);
            forceZAddr[place] = magnitude * (zLoc_LR/length_current);

            idKey[place + 1] = edgeR;
            forceXAddr[place + 1] = -magnitude * (xLoc_LR/length_current);
            forceYAddr[place + 1] = -magnitude * (yLoc_LR/length_current);
            forceZAddr[place + 1] = -magnitude * (zLoc_LR/length_current);
        
        double energy = (spring_constant/2.0) * (length_current - length_zero) * (length_current - length_zero);
        return energy;


    }
};

#endif