
#ifndef LJSPRINGS_H_
#define LJSPRINGS_H_ 

#include "SystemStructures.h"

void ComputeLJSprings(
    CoordInfoVecs& coordInfoVecs,
    LJInfoVecs& ljInfoVecs,
    GeneralParams& generalParams);

void AdvanceLJParticle(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LJInfoVecs& ljInfoVecs);

struct LJSpringFunctor {
    double Rcutoff;
    double Rmin;
    double epsilon;

    double LJLocX;
    double LJLocY;
    double LJLocZ;

    double* locXAddr;
    double* locYAddr;
    double* locZAddr;

    double* forceXAddr;
    double* forceYAddr;
    double* forceZAddr;
    
	__host__ __device__ LJSpringFunctor(
        double& _Rcutoff,
        double& _Rmin,
        double& _epsilon,

        double& _LJLocX,
        double& _LJLocY,
        double& _LJLocZ,

        double* _locXAddr,
        double* _locYAddr,
        double* _locZAddr,

        double* _forceXAddr,
        double* _forceYAddr,
        double* _forceZAddr) :

        Rcutoff(_Rcutoff),
        Rmin(_Rmin),
        epsilon(_epsilon),

        LJLocX(_LJLocX),
        LJLocY(_LJLocY),
        LJLocZ(_LJLocZ),

        locXAddr(_locXAddr),
        locYAddr(_locYAddr),
        locZAddr(_locZAddr),

        forceXAddr(_forceXAddr),
        forceYAddr(_forceYAddr),
        forceZAddr(_forceZAddr) {}

	//hand in counting iterator and id of two edges and preferred length
	__device__ CVec4 operator()(const unsigned & id) {
        		
        // compute forces
        double xLoc = locXAddr[id];
        double yLoc = locYAddr[id];
        double zLoc = locZAddr[id];
   

        double R = sqrt( 
            (LJLocX - xLoc) * (LJLocX - xLoc) + 
            (LJLocY - yLoc) * (LJLocY - yLoc) + 
            (LJLocZ - zLoc) * (LJLocZ - zLoc) );
        
        double forceX=0.0;
        double forceY=0.0;
        double forceZ=0.0;
        if (R < Rcutoff) {
            //double mag = 50*(R-4.0);
            forceX = //mag*((LJLocX - xLoc));
            epsilon * (-12.0 * pow(Rmin, 12.0) * (LJLocX - xLoc)/pow(R, 14.0) + 
                            12 * pow(Rmin, 6.0)  * (LJLocX - xLoc)/pow(R, 8.0) );

            forceY = //mag*(LJLocY - yLoc);
            epsilon * (-12.0 * pow(Rmin, 12.0) * (LJLocY - yLoc)/pow(R, 14.0) + 
                            12.0 * pow(Rmin, 6.0) * (LJLocY - yLoc)/pow(R, 8.0) );

            forceZ = //mag*(LJLocZ - zLoc);
            epsilon * (-12.0 * pow(Rmin, 12)  * (LJLocZ - zLoc)/pow(R, 14.0) + 
                            12.0 * pow(Rmin, 6.0) * (LJLocZ - zLoc)/pow(R, 8.0) );

            //WARNING: since this function modifies nodeForceX etc, you cannot rewrite entries. 
            forceXAddr[id] += forceX;
            forceYAddr[id] += forceY;
            forceZAddr[id] += forceZ;
        }

        

        double energy = epsilon * (pow( Rmin/R, 12.0 ) - 2.0 * pow( Rmin / R, 6.0));

        return thrust::make_tuple(-forceX, -forceY, -forceZ, energy);


    }
};

#endif