#ifndef TURGORFORCE_H_
#define TURGORFORCE_H_ 

#include "SystemStructures.h"

void ComputeTurgorSprings(
  
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    AreaTriangleInfoVecs& areaTriangleInfoVecs);


struct TurgorSpringFunctor {
    double spring_constant;
    double* locXAddr;
    double* locYAddr;
    double* locZAddr;

    int* idKey;
    double* forceXAddr;
    double* forceYAddr;
    double* forceZAddr;
    
	__host__ __device__ TurgorSpringFunctor(
        double& _spring_constant,
        double* _locXAddr,
        double* _locYAddr,
        double* _locZAddr,
        int* _idKey,
        double* _forceXAddr,
        double* _forceYAddr,
        double* _forceZAddr):
        spring_constant(_spring_constant),
        locXAddr(_locXAddr),
        locYAddr(_locYAddr),
        locZAddr(_locZAddr),
        idKey(_idKey),
        forceXAddr(_forceXAddr),
        forceYAddr(_forceYAddr),
        forceZAddr(_forceZAddr) {}

	//hand in counting iterator and id of triangle
	__device__ double operator()(const Tuuuuuuu &u7) {
        //test placing the ids of the nodes and then get positions. 
        //double scaling_pow = 4.0;
		int counter = thrust::get<0>(u7);
		int place = 3 * counter;//represents location in write to vector.

        int id_i = thrust::get<1>(u7);
        int id_j = thrust::get<2>(u7);
        int id_k = thrust::get<3>(u7);
        int e_id_i = thrust::get<4>(u7);
        int e_id_j = thrust::get<5>(u7);
        int e_id_k = thrust::get<6>(u7);
        
        //if (id_i != INT_MAX && id_j != INT_MAX && id_k != INT_MAX){
        if ((id_i != INT_MAX || id_i < 0) && (id_j != INT_MAX || id_j < 0) && (id_k != INT_MAX || id_k < 0)){           
            CVec3 ri = thrust::make_tuple<double>(locXAddr[id_i], locYAddr[id_i], locZAddr[id_i]);
            CVec3 rj = thrust::make_tuple<double>(locXAddr[id_j], locYAddr[id_j], locZAddr[id_j]);
            CVec3 rk = thrust::make_tuple<double>(locXAddr[id_k], locYAddr[id_k], locZAddr[id_k]);

            CVec3 rkj = CVec3_subtract(rk, rj);
            CVec3 rij = CVec3_subtract(ri, rj);

            double area_current = sqrt( CVec3_dot( CVec3_cross(rkj, rij), CVec3_cross(rkj, rij) ) )/2.0;

            //computes derivative wrt to area
            CVec3 A = CVec3_cross(rkj, rij);//rkj must come first
            double A1 = thrust::get<0>(A);
            double A2 = thrust::get<1>(A);
            double A3 = thrust::get<2>(A);
            double magnitude = sqrt(A1*A1 + A2*A2 + A3*A3);
            A1 = A1/magnitude;
            A2 = A2/magnitude;
            A3 = A3/magnitude;
            CVec3 rj_force = thrust::make_tuple<double>((1.0/3.0)*area_current*A1*spring_constant,
                                                        (1.0/3.0)*area_current*A2*spring_constant,
                                                        (1.0/3.0)*area_current*A3*spring_constant);
                
            CVec3 rk_force = thrust::make_tuple<double>((1.0/3.0)*area_current*A1*spring_constant,
                                                        (1.0/3.0)*area_current*A2*spring_constant,
                                                        (1.0/3.0)*area_current*A3*spring_constant);

            CVec3 ri_force = thrust::make_tuple<double>((1.0/3.0)*area_current*A1*spring_constant,
                                                        (1.0/3.0)*area_current*A2*spring_constant,
                                                        (1.0/3.0)*area_current*A3*spring_constant);

            idKey[place] = id_i;
            forceXAddr[place] = thrust::get<0>( ri_force );
            forceYAddr[place] = thrust::get<1>( ri_force );
            forceZAddr[place] = thrust::get<2>( ri_force );

            idKey[place+1] = id_j;
            forceXAddr[place+1] = thrust::get<0>( rj_force );
            forceYAddr[place+1] = thrust::get<1>( rj_force );
            forceZAddr[place+1] = thrust::get<2>( rj_force );
            
            idKey[place+2] = id_k;
            forceXAddr[place+2] = thrust::get<0>( rk_force );
            forceYAddr[place+2] = thrust::get<1>( rk_force );
            forceZAddr[place+2] = thrust::get<2>( rk_force );

            //double energy =  what_spring_constant/(2.0) * (area_current - area_0) * (area_current - area_0) / area_0;
            //return energy;
        }
        else{
            //double energy = 0.0;
            //return energy;
        }

    };
};
#endif //AREATRIANGLES_H_
