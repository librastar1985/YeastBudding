
#ifndef LJENERGY_H_
#define LJENERGY_H_ 

#include "SystemStructures.h"

double ComputeLJEnergy(
    CoordInfoVecs& coordInfoVecs,
    LJInfoVecs& ljInfoVecs,
    GeneralParams& generalParams);



struct LJEnergyFunctor {
    double Rcutoff;
    double Rmin;
    double epsilon;

    double LJLocX;
    double LJLocY;
    double LJLocZ;

    double* locXAddr;
    double* locYAddr;
    double* locZAddr;

    
	__host__ __device__ LJEnergyFunctor(
        double& _Rcutoff,
        double& _Rmin,
        double& _epsilon,

        double& _LJLocX,
        double& _LJLocY,
        double& _LJLocZ,

        double* _locXAddr,
        double* _locYAddr,
        double* _locZAddr

        ):

        Rcutoff(_Rcutoff),
        Rmin(_Rmin),
        epsilon(_epsilon),

        LJLocX(_LJLocX),
        LJLocY(_LJLocY),
        LJLocZ(_LJLocZ),

        locXAddr(_locXAddr),
        locYAddr(_locYAddr),
        locZAddr(_locZAddr){}

	//hand in counting iterator and id of two edges and preferred length
	__device__ double operator()(const int & id) {
        		
        // compute forces
        double xLoc = locXAddr[id];
        double yLoc = locYAddr[id];
        double zLoc = locZAddr[id];
   

        double R = sqrt( 
            (LJLocX - xLoc) * (LJLocX - xLoc) + 
            (LJLocY - yLoc) * (LJLocY - yLoc) + 
            (LJLocZ - zLoc) * (LJLocZ - zLoc) );
        

		double energy = 0.0;
        if (R < Rcutoff) {
            //energy = epsilon + (epsilon * (pow( Rmin/R, 12.0 ) - 2.0 * pow( Rmin / R, 6.0)));
            energy = (epsilon * (pow( Rmin/R, 12.0 ) - 2.0 * pow( Rmin / R, 6.0)));
        }

        

        //double energy = epsilon * (pow( Rmin/R, 12.0 ) - 2.0 * pow( Rmin / R, 6.0));

        return energy;

    }
};

#endif