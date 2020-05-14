#ifndef LJSPRINGS_LJ_H_
#define LJSPRINGS_LJ_H_ 

#include "SystemStructures.h"

void ComputeLJSprings_LJ(
    CoordInfoVecs& coordInfoVecs,
    LJInfoVecs& ljInfoVecs,
    GeneralParams& generalParams);

struct LJSpringFunctor_LJ {
    double Rcutoff;
    double Rmin;
    double epsilon_rep1;
    double epsilon_rep2;

    double LJLocX;
    double LJLocY;
    double LJLocZ;

    double* LJLocX_all;
    double* LJLocY_all;
    double* LJLocZ_all;

    double* LJforceXAddr;
    double* LJforceYAddr;
    double* LJforceZAddr;
    
	__host__ __device__ LJSpringFunctor_LJ(
        double& _Rcutoff,
        double& _Rmin,
        double& _epsilon_rep1,
        double& _epsilon_rep2,

        double& _LJLocX,
        double& _LJLocY,
        double& _LJLocZ,

        double* _LJLocX_all,
        double* _LJLocY_all,
        double* _LJLocZ_all,

        double* _LJforceXAddr,
        double* _LJforceYAddr,
        double* _LJforceZAddr) :

        Rcutoff(_Rcutoff),
        Rmin(_Rmin),
        epsilon_rep1(_epsilon_rep1),
        epsilon_rep2(_epsilon_rep2),

        LJLocX(_LJLocX),
        LJLocY(_LJLocY),
        LJLocZ(_LJLocZ),

        LJLocX_all(_LJLocX_all),
        LJLocY_all(_LJLocY_all),
        LJLocZ_all(_LJLocZ_all),

        LJforceXAddr(_LJforceXAddr),
        LJforceYAddr(_LJforceYAddr),
        LJforceZAddr(_LJforceZAddr) {}

	//hand in counting iterator and id of two edges and preferred length
	__device__ CVec4 operator()(const int & id) {
        		
        // compute forces
        double xLoc = LJLocX_all[id];
        double yLoc = LJLocY_all[id];
        double zLoc = LJLocZ_all[id];
   

        double R = sqrt( 
            (LJLocX - xLoc) * (LJLocX - xLoc) + 
            (LJLocY - yLoc) * (LJLocY - yLoc) + 
            (LJLocZ - zLoc) * (LJLocZ - zLoc) );
        
        double forceX=0.0;
        double forceY=0.0;
        double forceZ=0.0;
        double energy = 0.0;
        if (R < Rcutoff && R != 0.0) {
            
            double magnitude = 2*epsilon_rep1*
                                    (1-exp(-epsilon_rep2*(R-Rmin)))*
                                    (-exp(-epsilon_rep2*(R-Rmin)))*
                                    (epsilon_rep2/R);

            forceX = -magnitude*(LJLocX - xLoc);//xLoc_LR;
            forceY = -magnitude*(LJLocY - yLoc);//yLoc_LR;
            forceZ = -magnitude*(LJLocZ - zLoc);//zLoc_LR;

           // forceX = epsilon * (-4.0 * pow(Rmin, 4.0) * (LJLocX - xLoc)/pow(R, 6.0) + 
           //                 4.0 * pow(Rmin, 2.0)  * (LJLocX - xLoc)/pow(R, 4.0) );

           // forceY = epsilon * (-4.0 * pow(Rmin, 4.0) * (LJLocY - yLoc)/pow(R, 6.0) + 
           //                 4.0 * pow(Rmin, 2.0) * (LJLocY - yLoc)/pow(R, 4.0) );

           // forceZ = epsilon * (-4.0 * pow(Rmin, 4.0)  * (LJLocZ - zLoc)/pow(R, 6.0) + 
           //                 4.0 * pow(Rmin, 2.0) * (LJLocZ - zLoc)/pow(R, 4.0) );

            //WARNING: since this function modifies nodeForceX etc, you cannot rewrite entries. 
            //LJforceXAddr[id] += forceX;
            //LJforceYAddr[id] += forceY;
            //LJforceZAddr[id] += forceZ;
            //energy = epsilon + epsilon * (pow( Rmin/R, 12.0 ) - 2.0 * pow( Rmin / R, 6.0));//force positivity
            energy = epsilon_rep1 * (1-exp(-epsilon_rep2*(R - Rmin))) * (1-exp(-epsilon_rep2*(R - Rmin)));
        }

        

        //double energy = epsilon + epsilon * (pow( Rmin/R, 12.0 ) - 2.0 * pow( Rmin / R, 6.0));//force positivity

        return thrust::make_tuple(-forceX, -forceY, -forceZ, energy);


    }
};

#endif