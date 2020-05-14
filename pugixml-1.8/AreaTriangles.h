#ifndef AREATRIANGLES_H_
#define AREATRIANGLES_H_ 

#include "SystemStructures.h"

void ComputeAreaTriangleSprings(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    AreaTriangleInfoVecs& areaTriangleInfoVecs);


struct AreaSpringFunctor {
	double area_0;
    double spring_constant;
    double* locXAddr;
    double* locYAddr;
    double* locZAddr;

    unsigned* idKey;
    double* forceXAddr;
    double* forceYAddr;
    double* forceZAddr;
    
	__host__ __device__ AreaSpringFunctor(
        double& _area_0,
        double& _spring_constant,
        double* _locXAddr,
        double* _locYAddr,
        double* _locZAddr,
        unsigned* _idKey,
        double* _forceXAddr,
        double* _forceYAddr,
        double* _forceZAddr):
        area_0(_area_0),
        spring_constant(_spring_constant),
        locXAddr(_locXAddr),
        locYAddr(_locYAddr),
        locZAddr(_locZAddr),
        idKey(_idKey),
        forceXAddr(_forceXAddr),
        forceYAddr(_forceYAddr),
        forceZAddr(_forceZAddr) {}

	//hand in counting iterator and id of triangle
	__device__ double operator()(const Tuuuu &u4) {
        //test placing the ids of the nodes and then get positions. 
		unsigned counter = thrust::get<0>(u4);
		unsigned place = 3 * counter;//represents location in write to vector.

        unsigned id_i = thrust::get<1>(u4);
        unsigned id_j = thrust::get<2>(u4);
        unsigned id_k = thrust::get<3>(u4);

		
		CVec3 ri = thrust::make_tuple<double>(locXAddr[id_i], locYAddr[id_i], locZAddr[id_i]);
		CVec3 rj = thrust::make_tuple<double>(locXAddr[id_j], locYAddr[id_j], locZAddr[id_j]);
		CVec3 rk = thrust::make_tuple<double>(locXAddr[id_k], locYAddr[id_k], locZAddr[id_k]);

		CVec3 rkj = CVec3_subtract(rk, rj);
		CVec3 rij = CVec3_subtract(ri, rj);

        double area_current = sqrt( CVec3_dot( CVec3_cross(rkj, rij), CVec3_cross(rkj, rij) ) )/2;

        //computes derivative wrt to area
        CVec3 A = CVec3_cross(rkj, rij);//rkj must come first
        double A1 = thrust::get<0>(A);
        double A2 = thrust::get<1>(A);
        double A3 = thrust::get<2>(A);
        
        CVec3 A1Rj = thrust::make_tuple<double>(
            0.0, -thrust::get<2>(rij) + thrust::get<2>(rkj), -thrust::get<1>(rkj) + thrust::get<1>(rij));
            //[0, -Rij(3)+Rkj(3), -Rkj(2)+Rij(2)];
        CVec3 A2Rj = thrust::make_tuple<double>(    
            thrust::get<2>(rij) - thrust::get<2>(rkj), 0.0, thrust::get<0>(rkj) - thrust::get<0>(rij));
            //[Rij(3)-Rkj(3), 0, Rkj(1)-Rij(1)];
        CVec3 A3Rj = thrust::make_tuple<double>(
            -thrust::get<1>(rij) + thrust::get<1>(rkj), -thrust::get<0>(rkj) + thrust::get<0>(rij), 0.0);
            //[-Rij(2)+Rkj(2), -Rkj(1)+Rij(1), 0];

        CVec3 A1Rk = thrust::make_tuple<double>(
            0.0, thrust::get<2>(rij), -thrust::get<1>(rij));
            //[0, Rij(3), -Rij(2)];
        //Derivative of A1 with respect to [Rkx, Rky, Rkz].
        CVec3 A2Rk = thrust::make_tuple<double>(
            -thrust::get<2>(rij), 0.0, thrust::get<0>(rij));
            //[-Rij(3), 0, Rij(1)];
        CVec3 A3Rk = thrust::make_tuple<double>(
            thrust::get<1>(rij), -thrust::get<0>(rij), 0.0);
            //[Rij(2), -Rij(1), 0];
        CVec3 A1Ri = thrust::make_tuple<double>(
            0.0, -thrust::get<2>(rkj), thrust::get<1>(rkj));
            //[0, -Rkj(3), Rkj(2)];
        CVec3 A2Ri = thrust::make_tuple<double>(
            thrust::get<2>(rkj), 0.0, -thrust::get<0>(rkj));
            //[Rkj(3), 0, -Rkj(1)];
        CVec3 A3Ri = thrust::make_tuple<double>(
            -thrust::get<1>(rkj), thrust::get<0>(rkj), 0.0);
            //[-Rkj(2), Rkj(1), 0];
    
        double magnitude = -((spring_constant) * (area_current - area_0)/area_0) / (2 * area_current) ;
        CVec3 rj_force = CVec3_scalermult(magnitude, CVec3_plus( 
            CVec3_scalermult(A1, A1Rj), CVec3_scalermult(A2, A2Rj), CVec3_scalermult(A3, A3Rj)));
            // -(k/A0)*(AREA(i)-A0) * (A1*A1Rj(1) + A2*A2Rj(1) + A3*A3Rj(1))/(2*AREA(i));
            // -(k/A0)*(AREA(i)-A0) * (A1*A1Rj(2) + A2*A2Rj(2) + A3*A3Rj(2))/(2*AREA(i));
            // -(k/A0)*(AREA(i)-A0) * (A1*A1Rj(3) + A2*A2Rj(3) + A3*A3Rj(3))/(2*AREA(i));
            
        CVec3 rk_force = CVec3_scalermult(magnitude, CVec3_plus( 
            CVec3_scalermult(A1, A1Rk), CVec3_scalermult(A2, A2Rk), CVec3_scalermult(A3, A3Rk)));
            //-(k/A0)*(AREA(i)-A0)*(A1*A1Rk(1) + A2*A2Rk(1) + A3*A3Rk(1))/(2*AREA(i));
            //-(k/A0)*(AREA(i)-A0)*(A1*A1Rk(2) + A2*A2Rk(2) + A3*A3Rk(2))/(2*AREA(i));
            //-(k/A0)*(AREA(i)-A0)*(A1*A1Rk(3) + A2*A2Rk(3) + A3*A3Rk(3))/(2*AREA(i));

        CVec3 ri_force = CVec3_scalermult(magnitude, CVec3_plus( 
            CVec3_scalermult(A1, A1Ri), CVec3_scalermult(A2, A2Ri), CVec3_scalermult(A3, A3Ri)));
            //-(k/A0)*(AREA(i)-A0)*(A1*A1Ri(1) + A2*A2Ri(1) + A3*A3Ri(1))/(2*AREA(i));
            //-(k/A0)*(AREA(i)-A0)*(A1*A1Ri(2) + A2*A2Ri(2) + A3*A3Ri(2))/(2*AREA(i));
            //-(k/A0)*(AREA(i)-A0)*(A1*A1Ri(3) + A2*A2Ri(3) + A3*A3Ri(3))/(2*AREA(i));
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

        double energy =  spring_constant/2 * (area_current - area_0) * (area_current - area_0) / area_0;
        return energy;

    };
};
#endif //AREATRIANGLES_H_