#ifndef LINEARSPRINGS_H_
#define LINEARSPRINGS_H_ 

#include "SystemStructures.h"

void ComputeLinearSprings(
  
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    LJInfoVecs& ljInfoVecs);
    
struct LinearSpringFunctor {
    int SCALE_TYPE;
    double scaling_pow;
    double gausssigma;
    double hilleqnconst;
    double hilleqnpow;
    double* scaling_per_edge;
    double spring_constant;
    double spring_constant_weak;
    double length_zero;
    //double length_zero_growth;
    int* edges_in_upperhem;
    int* edges_in_tip;
    //int* boundaries_in_upperhem;

    double* locXAddr;
    double* locYAddr;
    double* locZAddr;

    int* idKey;
    double* forceXAddr;
    double* forceYAddr;
    double* forceZAddr;
    
	__host__ __device__ LinearSpringFunctor(
        int& _SCALE_TYPE,
        double& _scaling_pow,
        double& _gausssigma,
        double& _hilleqnconst,
        double& _hilleqnpow,
        double* _scaling_per_edge,
        double& _spring_constant,
        double& _spring_constant_weak,
        double& _length_zero,
       // double& _length_zero_growth,
        int* _edges_in_upperhem,
        int* _edges_in_tip,
       // int* _boundaries_in_upperhem,
        double* _locXAddr,
        double* _locYAddr,
        double* _locZAddr,
        
        int* _idKey,
        double* _forceXAddr,
        double* _forceYAddr,
        double* _forceZAddr) :
        SCALE_TYPE(_SCALE_TYPE),
        scaling_pow(_scaling_pow),
        gausssigma(_gausssigma),
        hilleqnconst(_hilleqnconst),
        hilleqnpow(_hilleqnpow),
        scaling_per_edge(_scaling_per_edge),
        spring_constant(_spring_constant),
        spring_constant_weak(_spring_constant_weak),
        length_zero(_length_zero),
       // length_zero_growth(_length_zero_growth),
        edges_in_upperhem(_edges_in_upperhem),
        edges_in_tip(_edges_in_tip),
       // boundaries_in_upperhem(_boundaries_in_upperhem),
        locXAddr(_locXAddr),
        locYAddr(_locYAddr),
        locZAddr(_locZAddr),
        idKey(_idKey),
        forceXAddr(_forceXAddr),
        forceYAddr(_forceYAddr),
        forceZAddr(_forceZAddr) {}

	//hand in counting iterator and id of two edges and preferred length
	__device__ double operator()(const Tuuu &u3d) {
        		
             //   double scaling_pow = 4.0;
        //counter ranges from 0 to num_edges. 
        int counter = thrust::get<0>(u3d);
		int place = 2 * counter;//represents location in write to vector.

        int edgeL = thrust::get<1>(u3d);
        int edgeR = thrust::get<2>(u3d);

        //if (edgeL != INT_MAX && edgeL >= 0 && edgeR != INT_MAX && edgeR >= 0){
        if (edgeL < (INT_MAX-100) && edgeL >= 0 && edgeR < (INT_MAX-100) && edgeR >= 0){
            double what_spring_constant;
            if (SCALE_TYPE == 0){
                what_spring_constant = spring_constant*(1.0 - ((1.0/sqrt(2*3.14159*gausssigma))*exp(-(scaling_per_edge[counter]*scaling_per_edge[counter])/gausssigma)));
                if (what_spring_constant < spring_constant_weak){what_spring_constant = spring_constant;}
            }
            else if (SCALE_TYPE == 1){
                if (edges_in_upperhem[counter] == 1 || edges_in_upperhem[counter] == 0){
                // what_spring_constant = (spring_constant_weak*2.0)*pow(scaling_per_edge[counter],scaling_pow) +
                //      spring_constant_weak*(1-pow(scaling_per_edge[counter], scaling_pow));
                // }
                // else{
                //     what_spring_constant = spring_constant;
                // }
                what_spring_constant = spring_constant*pow(scaling_per_edge[counter],scaling_pow) +
                     spring_constant_weak*(1-pow(scaling_per_edge[counter], scaling_pow));
                }
            }
            else if (SCALE_TYPE == 2){
                what_spring_constant = spring_constant - (spring_constant - spring_constant_weak)*scaling_per_edge[counter];
            }
            else if (SCALE_TYPE == 3){
                if (edges_in_upperhem[counter] == 1 && edges_in_tip[counter]==1){
                    what_spring_constant = spring_constant_weak;
                    //length_zero = length_zero_growth;
                }
                else if (edges_in_upperhem[counter] == 1 && edges_in_tip[counter] != 1){
                    what_spring_constant = spring_constant_weak*2.0;
                    //what_spring_constant = (spring_constant_weak + spring_constant)/2.0;
                }
                else{
                    what_spring_constant = spring_constant;
                }
            }
            else if (SCALE_TYPE == 4){
                //double scaling = 0.0;//spring_constant_weak/spring_constant;
                double spectrum = spring_constant - spring_constant_weak;
			    //what_spring_constant = spring_constant*((1.0/(1.0+pow(hilleqnconst/scaling_per_edge[counter], hilleqnpow)))*(1-scaling) + scaling);
			    what_spring_constant = spring_constant_weak + ((1.0/(1.0+pow(hilleqnconst/scaling_per_edge[counter], hilleqnpow)))*spectrum);
                if (what_spring_constant < spring_constant_weak){what_spring_constant = spring_constant_weak;}
		    }

            
            //double length_zero = thrust::get<3>(u3d);

            // compute forces.
            double xLoc_LR = locXAddr[edgeL] - locXAddr[edgeR];
            double yLoc_LR = locYAddr[edgeL] - locYAddr[edgeR];
            double zLoc_LR = locZAddr[edgeL] - locZAddr[edgeR];
    

            double length_current = sqrt( (xLoc_LR) * (xLoc_LR) + 
                                        (yLoc_LR) * (yLoc_LR)  + 
                                        (zLoc_LR) * (zLoc_LR) );

        double energy = 0.0;
        if (length_current != length_zero){
            //double magnitude = -(what_spring_constant/(length_zero*length_zero)) * (length_current - length_zero);
            double magnitude = -(what_spring_constant) * (length_current - length_zero);
            
                idKey[place] = edgeL;

                
            //issue here writing to vectors with force????
                forceXAddr[place] = magnitude * (xLoc_LR/length_current);
                forceYAddr[place] = magnitude * (yLoc_LR/length_current);
                forceZAddr[place] = magnitude * (zLoc_LR/length_current);

                idKey[place + 1] = edgeR;
                forceXAddr[place + 1] = -magnitude * (xLoc_LR/length_current);
                forceYAddr[place + 1] = -magnitude * (yLoc_LR/length_current);
                forceZAddr[place + 1] = -magnitude * (zLoc_LR/length_current);
            
            //energy = (what_spring_constant/(2.0*length_zero*length_zero)) * (length_current - length_zero) * (length_current - length_zero);
             energy = (what_spring_constant/(2.0)) * (length_current - length_zero) * (length_current - length_zero);
        }
            return energy;
        }

        else{
            double energy = 0.0;
            return energy;
        }
    }
    
};

#endif
