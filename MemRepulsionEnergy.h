#ifndef MEMREPULSIONENERGY_H_
#define MEMREPULSIONENERGY_H_ 

#include "SystemStructures.h"

double ComputeMemRepulsionEnergy(
    CoordInfoVecs& coordInfoVecs,
	LinearSpringInfoVecs& linearSpringInfoVecs, 
	CapsidInfoVecs& capsidInfoVecs,
    GeneralParams& generalParams,
    AuxVecs& auxVecs);

struct MemRepulsionEnergyFunctor : public thrust::unary_function< UCVec3, double> {
    double Rmin;
    double abs_Rmin;
    double epsilon_rep1;
    double epsilon_rep2;
    int membraneMaxNode;

    int* nndata1;
    int* nndata2;
    int* nndata3;
    int* nndata4;
    int* nndata5;
    int* nndata6;
    int* nndata7;
    int* nndata8;
    int* nndata9;
    //int* nndata10;
    //int* nndata11;
    //int* nndata12;

    double* membraneNodeXAddr;
    double* membraneNodeYAddr;
    double* membraneNodeZAddr;

    int* id_value_expanded;
    int* keyBegin;
    int* keyEnd;
    
	__host__ __device__ 
    MemRepulsionEnergyFunctor(
        double& _Rmin,
        double& _abs_Rmin,
        double& _epsilon_rep1,
        double& _epsilon_rep2,
        int& _membraneMaxNode,

        int* _nndata1,
        int* _nndata2,
        int* _nndata3,
        int* _nndata4,
        int* _nndata5,
        int* _nndata6,
        int* _nndata7,
        int* _nndata8,
        int* _nndata9,
     //   int* _nndata10,
     //   int* _nndata11,
      //  int* _nndata12,
     
        double* _membraneNodeXAddr,
        double* _membraneNodeYAddr,
        double* _membraneNodeZAddr,

        int* _id_value_expanded,
        int* _keyBegin,
        int* _keyEnd):

        Rmin(_Rmin),
        abs_Rmin(_abs_Rmin),
        epsilon_rep1(_epsilon_rep1),
        epsilon_rep2(_epsilon_rep2),
        membraneMaxNode(_membraneMaxNode),

        nndata1(_nndata1),
        nndata2(_nndata2),
        nndata3(_nndata3),
        nndata4(_nndata4),
        nndata5(_nndata5),
        nndata6(_nndata6),
        nndata7(_nndata7),
        nndata8(_nndata8),
        nndata9(_nndata9),
      //  nndata10(_nndata10),
      //  nndata11(_nndata11),
      //  nndata12(_nndata12),
    
        membraneNodeXAddr(_membraneNodeXAddr),
        membraneNodeYAddr(_membraneNodeYAddr),
        membraneNodeZAddr(_membraneNodeZAddr),

        id_value_expanded(_id_value_expanded),
        keyBegin(_keyBegin),
        keyEnd(_keyEnd) {}

	//hand in counting iterator id for membrane nodes
	__device__ 
    double operator()(const UCVec3& u1d3) {

        int node_id = thrust::get<0>(u1d3);
 
		double posX = thrust::get<1>(u1d3);
        double posY = thrust::get<2>(u1d3);
        double posZ = thrust::get<3>(u1d3);

        double xLoc_LR;
        double yLoc_LR;
        double zLoc_LR;

        
        double epsilon_rep, energy = 0.0;
       for (unsigned memId_count = 0; memId_count < membraneMaxNode; memId_count++ ){
            
            unsigned memId = memId_count;//id_value_expanded[ memId_count ];
            if (nndata1[node_id] == memId ||
                 nndata2[node_id] == memId ||
                 nndata3[node_id] == memId ||
                 nndata4[node_id] == memId ||
                 nndata5[node_id] == memId ||
                 nndata6[node_id] == memId ||
                 nndata7[node_id] == memId ||
                 nndata8[node_id] == memId ||
                 nndata9[node_id] == memId ){
                     continue;
                 }


            if ((memId < membraneMaxNode) && (memId != node_id)) {
                //calculate distance
                // if (nodes_in_upperhem[node_id] == 1 && nodes_in_upperhem[memId]==1){
                //     epsilon_rep = epsilon_rep1;
                // }
                // else if (nodes_in_upperhem[node_id] != 1 || nodes_in_upperhem[memId] != 1){
                //     epsilon_rep = epsilon_rep2;
                // }

                xLoc_LR = -( posX - membraneNodeXAddr[memId] ); //This is like Rk - Rj, where Rj is the point we are trying to update
                yLoc_LR = -( posY - membraneNodeYAddr[memId] );
                zLoc_LR = -( posZ - membraneNodeZAddr[memId] );

                double R = sqrt( 
                            (xLoc_LR) * (xLoc_LR) + 
                            (yLoc_LR) * (yLoc_LR) + 
                            (zLoc_LR) * (zLoc_LR) );

                if (R < abs_Rmin) {
                     energy += epsilon_rep1*(1.0-exp(-epsilon_rep2*(R - abs_Rmin)))*(1.0-exp(-epsilon_rep2*(R - abs_Rmin)));
                }
                // if( R < abs_Rmin){
                //   current_forceX += (epsilon_rep)*(R - abs_Rmin)*(xLoc_LR)/R;
				// 	current_forceY += (epsilon_rep)*(R - abs_Rmin)*(yLoc_LR)/R;
				// 	current_forceZ += (epsilon_rep)*(R - abs_Rmin)*(zLoc_LR)/R;
                // }
            }
        }

       // energy += energy1 + energy2 + energy3 + energy4 + energy5 + energy6 + energy7 + energy8 + energy9;// + energy10 + energy11 + energy12; 
       //}
    
        return energy;
    }
};

#endif

