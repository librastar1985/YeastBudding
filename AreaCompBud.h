#ifndef AREACOMPBUD_H_
#define AREACOMPBUD_H_ 

#include "SystemStructures.h"

void ComputeAreaBud(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs//,
    //LinearSpringInfoVecs& linearSpringInfoVecs,
    //LJInfoVecs& ljInfoVecs
    );
    
struct AreaBudCompFunctor {
    
    //double spring_constant;
    double* locXAddr;
    double* locYAddr;
    double* locZAddr;
    int* triangles_in_upperhem;
	__host__ __device__ AreaBudCompFunctor(
      //  double& _spring_constant,
        double* _locXAddr,
        double* _locYAddr,
        double* _locZAddr,
        int* _triangles_in_upperhem
        ):

        //spring_constant(_spring_constant),
        locXAddr(_locXAddr),
        locYAddr(_locYAddr),
        locZAddr(_locZAddr),
        triangles_in_upperhem(_triangles_in_upperhem){}

	//hand in counting iterator and id of two edges and preferred length
	__device__ double operator()(const Tuuuu &u4) {
        		
        //counter ranges from 0 to num_edges. 
        int counter = thrust::get<0>(u4);
        int r1 = thrust::get<1>(u4);
        int r2 = thrust::get<2>(u4);
        int r3 = thrust::get<3>(u4);

        if (triangles_in_upperhem[counter] == 1 && r1 != INT_MAX && r2 != INT_MAX && r3 != INT_MAX){
        double r1x = locXAddr[r1];
        double r1y = locYAddr[r1];
        double r1z = locZAddr[r1];
        double r2x = locXAddr[r2];
        double r2y = locYAddr[r2];
        double r2z = locZAddr[r2];
        double r3x = locXAddr[r3];
        double r3y = locYAddr[r3];
        double r3z = locZAddr[r3];
        double norm_r1r2 = sqrt((r2x-r1x)*(r2x-r1x) + (r2y-r1y)*(r2y-r1y) + (r2z-r1z)*(r2z-r1z));
        double norm_r2r3 = sqrt((r3x-r2x)*(r3x-r2x) + (r3y-r2y)*(r3y-r2y) + (r3z-r2z)*(r3z-r2z));
        double norm_r3r1 = sqrt((r3x-r1x)*(r3x-r1x) + (r3y-r1y)*(r3y-r1y) + (r3z-r1z)*(r3z-r1z));
        double s = (norm_r1r2 + norm_r2r3 + norm_r3r1)/2.0;
        double area = sqrt(s*(s-norm_r1r2)*(s-norm_r2r3)*(s-norm_r3r1));

        return area;
        
        }
        else{
            double area = 0.0;
            
            return area;
        }


    }
};

#endif