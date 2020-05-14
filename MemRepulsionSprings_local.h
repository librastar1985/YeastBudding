
#ifndef MEMREPULSIONSPRINGS_LOCAL_H_
#define MEMREPULSIONSPRINGS_LOCAL_H_ 

#include "SystemStructures.h"

void ComputeMemRepulsionSprings_local(
    CoordInfoVecs& coordInfoVecs,
	LinearSpringInfoVecs& linearSpringInfoVecs, 
	CapsidInfoVecs& capsidInfoVecs,
    GeneralParams& generalParams,
    AuxVecs& auxVecs);

struct MemRepulsionSpringFunctor_local : public thrust::unary_function<U2CVec3,CVec3> {
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
    int* nodes_in_upperhem;
    double* membraneNodeXAddr;
    double* membraneNodeYAddr;
    double* membraneNodeZAddr;

    int* id_value_expanded;
    int* keyBegin;
    int* keyEnd;
    
	__host__ __device__ 
    MemRepulsionSpringFunctor_local(
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
        //int* _nndata10,
       // int* _nndata11,
       // int* _nndata12,
        int* _nodes_in_upperhem,
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
        //nndata10(_nndata10),
        //nndata11(_nndata11),
        //nndata12(_nndata12),
        nodes_in_upperhem(_nodes_in_upperhem),
        membraneNodeXAddr(_membraneNodeXAddr),
        membraneNodeYAddr(_membraneNodeYAddr),
        membraneNodeZAddr(_membraneNodeZAddr),

        id_value_expanded(_id_value_expanded),
        keyBegin(_keyBegin),
        keyEnd(_keyEnd) {}

	//hand in counting iterator id for membrane nodes
	__device__ 
    CVec3 operator()(const U2CVec3& u2d3) {

        int node_id = thrust::get<0>(u2d3);
        
        int bucketId = thrust::get<1>(u2d3);//bucket containing nodeId
	
		//beginning and end of attempted attachment id's in id_value_expanded
        //these indices are membrane
	
        double xLoc_LR;
        double yLoc_LR;
        double zLoc_LR;

        //extract current force
        double current_forceX = thrust::get<2>(u2d3);
        double current_forceY = thrust::get<3>(u2d3);
        double current_forceZ = thrust::get<4>(u2d3);

		double posX = membraneNodeXAddr[node_id];
        double posY = membraneNodeYAddr[node_id];
        double posZ = membraneNodeZAddr[node_id];

        //for now iterate through all membrane id's
       // int begin = keyBegin[bucketId];
       // int end = keyEnd[bucketId];

       //First, create a new vector detailing all the neighboring nodes of node i//
       
       int neighbor1 = nndata1[node_id];
        int neighbor2 = nndata2[node_id];
        int neighbor3 = nndata3[node_id];
        int neighbor4 = nndata4[node_id];
        int neighbor5 = nndata5[node_id];
        int neighbor6 = nndata6[node_id];
        int neighbor7 = nndata7[node_id];
        int neighbor8 = nndata8[node_id];
        int neighbor9 = nndata9[node_id];
        //int neighbor10 = nndata10[node_id];
        //int neighbor11 = nndata11[node_id];
        //int neighbor12 = nndata12[node_id];

       //NN tells us the neighboring nodes that we will search for its neighbor again
        double current_forceX1, current_forceX2, current_forceX3, current_forceX4, 
                current_forceX5, current_forceX6, current_forceX7, current_forceX8, 
                current_forceX9;//, current_forceX10, current_forceX11, current_forceX12;
        double current_forceY1, current_forceY2, current_forceY3, current_forceY4, 
                current_forceY5, current_forceY6, current_forceY7, current_forceY8, 
                current_forceY9;//, current_forceY10, current_forceY11, current_forceY12;
        double current_forceZ1, current_forceZ2, current_forceZ3, current_forceZ4, 
                current_forceZ5, current_forceZ6, current_forceZ7, current_forceZ8, 
                current_forceZ9;//, current_forceZ10, current_forceZ11, current_forceZ12;
        double magnitude;
        int neighbor;
    
    for (int j = 0; j < 12; j ++){
        if (j == 0 && (nndata1[node_id] >= 0)){neighbor = nndata1[node_id];}
        //else{continue;}
        else if (j == 1 && (nndata2[node_id] >= 0)){neighbor = nndata2[node_id];}
        //else{continue;}
        else if (j == 2 && (nndata3[node_id] >= 0)){neighbor = nndata3[node_id];}
        //else{continue;}
        else if (j == 3 && (nndata4[node_id] >= 0)){neighbor = nndata4[node_id];}
        //else{continue;}
        else if (j == 4 && (nndata5[node_id] >= 0)){neighbor = nndata5[node_id];}
        //else{continue;}
        else if (j == 5 && (nndata6[node_id] >= 0)){neighbor = nndata6[node_id];}
        //else{continue;}
        else if (j == 6 && (nndata7[node_id] >= 0)){neighbor = nndata7[node_id];}
        //else{continue;}
        else if (j == 7 && (nndata8[node_id] >= 0)){neighbor = nndata8[node_id];}
        //else{continue;}
        else if (j == 8 && (nndata9[node_id] >= 0)){neighbor = nndata9[node_id];}
        //else{continue;}
       // else if (j == 9 && (nndata10[node_id] >= 0)){neighbor = nndata10[node_id];}
        //else{continue;}
       // else if (j == 10 && (nndata11[node_id] >= 0)){neighbor = nndata11[node_id];}
        //else{continue;}
        //else if (j == 11 && (nndata12[node_id] >= 0)){neighbor = nndata12[node_id];}
        else{continue;}
       //for (int i = 0; i < 12; i++){
           if (nndata1[neighbor] >= 0 &&
                 nndata1[neighbor] != node_id &&
                 nndata1[neighbor] != neighbor1 &&
                 nndata1[neighbor] != neighbor2 &&
                 nndata1[neighbor] != neighbor3 &&
                 nndata1[neighbor] != neighbor4 &&
                 nndata1[neighbor] != neighbor5 &&
                 nndata1[neighbor] != neighbor6 &&
                 nndata1[neighbor] != neighbor7 &&
                 nndata1[neighbor] != neighbor8 &&
                 nndata1[neighbor] != neighbor9 /*&&
                 nndata1[neighbor] != neighbor10 &&
                 nndata1[neighbor] != neighbor11 &&
                 nndata1[neighbor] != neighbor12*/ ){
                xLoc_LR = -( posX - membraneNodeXAddr[nndata1[neighbor]] ); //This is like Rk - Rj, where Rj is the point we are trying to update
                yLoc_LR = -( posY - membraneNodeYAddr[nndata1[neighbor]] );
                zLoc_LR = -( posZ - membraneNodeZAddr[nndata1[neighbor]] );
                double R = sqrt( 
                            (xLoc_LR) * (xLoc_LR) + 
                            (yLoc_LR) * (yLoc_LR) + 
                            (zLoc_LR) * (zLoc_LR) );
                if (R < abs_Rmin){
                    
                    // magnitude = 2*epsilon_rep1*
                    //                         (1-exp(-epsilon_rep2*(R-abs_Rmin)))*
                    //                         (-exp(-epsilon_rep2*(R-abs_Rmin)))*
                    //                         (epsilon_rep2/R);

                    //     current_forceX1 = -magnitude*xLoc_LR;
                    //     current_forceY1 = -magnitude*yLoc_LR;
                    //     current_forceZ1 = -magnitude*zLoc_LR;
                    double epsilon_rep;
                    if (nodes_in_upperhem[node_id] == 1 && nodes_in_upperhem[nndata1[neighbor]]==1){
                        epsilon_rep = epsilon_rep1;
                    }
                    else if (nodes_in_upperhem[node_id] != 1 || nodes_in_upperhem[nndata1[neighbor]] != 1){
                        epsilon_rep = epsilon_rep2;
                    }
                    current_forceX1 = (epsilon_rep)*(R - abs_Rmin)*(xLoc_LR)/R;
					current_forceY1 = (epsilon_rep)*(R - abs_Rmin)*(yLoc_LR)/R;
					current_forceZ1 = (epsilon_rep)*(R - abs_Rmin)*(zLoc_LR)/R;
                }
                else{
                    current_forceX1 = 0.0;
                    current_forceY1 = 0.0;
                    current_forceZ1 = 0.0;
                }
           }
           else{
               current_forceX1 = 0.0;
                current_forceY1 = 0.0;
                current_forceZ1 = 0.0;
           }
           
           if (nndata2[neighbor] >= 0 &&
                 nndata2[neighbor] != node_id &&
                 nndata2[neighbor] != neighbor1 &&
                 nndata2[neighbor] != neighbor2 &&
                 nndata2[neighbor] != neighbor3 &&
                 nndata2[neighbor] != neighbor4 &&
                 nndata2[neighbor] != neighbor5 &&
                 nndata2[neighbor] != neighbor6 &&
                 nndata2[neighbor] != neighbor7 &&
                 nndata2[neighbor] != neighbor8 &&
                 nndata2[neighbor] != neighbor9/* &&
                 nndata2[neighbor] != neighbor10 &&
                 nndata2[neighbor] != neighbor11 &&
                 nndata2[neighbor] != neighbor12*/){
                xLoc_LR = -( posX - membraneNodeXAddr[nndata2[neighbor]] ); //This is like Rk - Rj, where Rj is the point we are trying to update
                yLoc_LR = -( posY - membraneNodeYAddr[nndata2[neighbor]] );
                zLoc_LR = -( posZ - membraneNodeZAddr[nndata2[neighbor]] );
                double R = sqrt( 
                            (xLoc_LR) * (xLoc_LR) + 
                            (yLoc_LR) * (yLoc_LR) + 
                            (zLoc_LR) * (zLoc_LR) );
                if (R < abs_Rmin){
                    // magnitude = 2*epsilon_rep1*
                    //                         (1-exp(-epsilon_rep2*(R-abs_Rmin)))*
                    //                         (-exp(-epsilon_rep2*(R-abs_Rmin)))*
                    //                         (epsilon_rep2/R);

                    //     current_forceX2 = -magnitude*xLoc_LR;
                    //     current_forceY2 = -magnitude*yLoc_LR;
                    //     current_forceZ2 = -magnitude*zLoc_LR;
                    double epsilon_rep;
                    if (nodes_in_upperhem[node_id] == 1 && nodes_in_upperhem[nndata2[neighbor]]==1){
                        epsilon_rep = epsilon_rep1;
                    }
                    else if (nodes_in_upperhem[node_id] != 1 || nodes_in_upperhem[nndata2[neighbor]] != 1){
                        epsilon_rep = epsilon_rep2;
                    }
                    current_forceX2 = (epsilon_rep)*(R - abs_Rmin)*(xLoc_LR)/R;
					current_forceY2 = (epsilon_rep)*(R - abs_Rmin)*(yLoc_LR)/R;
					current_forceZ2 = (epsilon_rep)*(R - abs_Rmin)*(zLoc_LR)/R;
                }
                else{
                    current_forceX2 = 0.0;
                    current_forceY2 = 0.0;
                    current_forceZ2 = 0.0;
                }
           }
           else{
               current_forceX2 = 0.0;
                current_forceY2 = 0.0;
                current_forceZ2 = 0.0;
           }
           
           if (nndata3[neighbor] >= 0 &&
                nndata3[neighbor] != node_id &&
                nndata3[neighbor] != neighbor1 &&
                nndata3[neighbor] != neighbor2 &&
                nndata3[neighbor] != neighbor3 &&
                nndata3[neighbor] != neighbor4 &&
                nndata3[neighbor] != neighbor5 &&
                nndata3[neighbor] != neighbor6 &&
                nndata3[neighbor] != neighbor7 &&
                nndata3[neighbor] != neighbor8 &&
                nndata3[neighbor] != neighbor9 /*&&
                nndata3[neighbor] != neighbor10 &&
                nndata3[neighbor] != neighbor11 &&
                nndata3[neighbor] != neighbor12*/){
                xLoc_LR = -( posX - membraneNodeXAddr[nndata3[neighbor]] ); //This is like Rk - Rj, where Rj is the point we are trying to update
                yLoc_LR = -( posY - membraneNodeYAddr[nndata3[neighbor]] );
                zLoc_LR = -( posZ - membraneNodeZAddr[nndata3[neighbor]] );
                double R = sqrt( 
                            (xLoc_LR) * (xLoc_LR) + 
                            (yLoc_LR) * (yLoc_LR) + 
                            (zLoc_LR) * (zLoc_LR) );
                if (R < abs_Rmin){
                    // magnitude = 2*epsilon_rep1*
                    //                         (1-exp(-epsilon_rep2*(R-abs_Rmin)))*
                    //                         (-exp(-epsilon_rep2*(R-abs_Rmin)))*
                    //                         (epsilon_rep2/R);

                    //     current_forceX3 = -magnitude*xLoc_LR;
                    //     current_forceY3 = -magnitude*yLoc_LR;
                    //     current_forceZ3 = -magnitude*zLoc_LR;
                    double epsilon_rep;
                    if (nodes_in_upperhem[node_id] == 1 && nodes_in_upperhem[nndata3[neighbor]]==1){
                        epsilon_rep = epsilon_rep1;
                    }
                    else if (nodes_in_upperhem[node_id] != 1 || nodes_in_upperhem[nndata3[neighbor]] != 1){
                        epsilon_rep = epsilon_rep2;
                    }
                    current_forceX3 = (epsilon_rep)*(R - abs_Rmin)*(xLoc_LR)/R;
					current_forceY3 = (epsilon_rep)*(R - abs_Rmin)*(yLoc_LR)/R;
					current_forceZ3 = (epsilon_rep)*(R - abs_Rmin)*(zLoc_LR)/R;
                    }
                    else{
                        current_forceX3 = 0.0;
                        current_forceY3 = 0.0;
                        current_forceZ3 = 0.0;
                    }
            }
           else{
               current_forceX3 = 0.0;
                current_forceY3 = 0.0;
                current_forceZ3 = 0.0;
           }
           
           if (nndata4[neighbor] >= 0 &&
                 nndata4[neighbor] != node_id &&
                 nndata4[neighbor] != neighbor1 &&
                 nndata4[neighbor] != neighbor2 &&
                 nndata4[neighbor] != neighbor3 &&
                 nndata4[neighbor] != neighbor4 &&
                 nndata4[neighbor] != neighbor5 &&
                 nndata4[neighbor] != neighbor6 &&
                 nndata4[neighbor] != neighbor7 &&
                 nndata4[neighbor] != neighbor8 &&
                 nndata4[neighbor] != neighbor9 /*&&
                 nndata4[neighbor] != neighbor10 &&
                 nndata4[neighbor] != neighbor11 &&
                 nndata4[neighbor] != neighbor12*/){
                xLoc_LR = -( posX - membraneNodeXAddr[nndata4[neighbor]] ); //This is like Rk - Rj, where Rj is the point we are trying to update
                yLoc_LR = -( posY - membraneNodeYAddr[nndata4[neighbor]] );
                zLoc_LR = -( posZ - membraneNodeZAddr[nndata4[neighbor]] );
                double R = sqrt( 
                            (xLoc_LR) * (xLoc_LR) + 
                            (yLoc_LR) * (yLoc_LR) + 
                            (zLoc_LR) * (zLoc_LR) );
            if (R < abs_Rmin){
                // magnitude = 2*epsilon_rep1*
                //                             (1-exp(-epsilon_rep2*(R-abs_Rmin)))*
                //                             (-exp(-epsilon_rep2*(R-abs_Rmin)))*
                //                             (epsilon_rep2/R);

                //         current_forceX4 = -magnitude*xLoc_LR;
                //         current_forceY4 = -magnitude*yLoc_LR;
                //         current_forceZ4 = -magnitude*zLoc_LR;
                double epsilon_rep;
                    if (nodes_in_upperhem[node_id] == 1 && nodes_in_upperhem[nndata4[neighbor]]==1){
                        epsilon_rep = epsilon_rep1;
                    }
                    else if (nodes_in_upperhem[node_id] != 1 || nodes_in_upperhem[nndata4[neighbor]] != 1){
                        epsilon_rep = epsilon_rep2;
                    }
                    current_forceX4 = (epsilon_rep)*(R - abs_Rmin)*(xLoc_LR)/R;
					current_forceY4 = (epsilon_rep)*(R - abs_Rmin)*(yLoc_LR)/R;
					current_forceZ4 = (epsilon_rep)*(R - abs_Rmin)*(zLoc_LR)/R;
                    }
                    else{
                        current_forceX4 = 0.0;
                        current_forceY4 = 0.0;
                        current_forceZ4 = 0.0;
                    }
            }
           else{
               current_forceX4 = 0.0;
                current_forceY4 = 0.0;
                current_forceZ4 = 0.0;
           }
           
           if (nndata5[neighbor] >= 0 &&
                 nndata5[neighbor] != node_id &&
                 nndata5[neighbor] != neighbor1 &&
                 nndata5[neighbor] != neighbor2 &&
                 nndata5[neighbor] != neighbor3 &&
                 nndata5[neighbor] != neighbor4 &&
                 nndata5[neighbor] != neighbor5 &&
                 nndata5[neighbor] != neighbor6 &&
                 nndata5[neighbor] != neighbor7 &&
                 nndata5[neighbor] != neighbor8 &&
                 nndata5[neighbor] != neighbor9 /*&&
                 nndata5[neighbor] != neighbor10 &&
                 nndata5[neighbor] != neighbor11 &&
                 nndata5[neighbor] != neighbor12*/){
                xLoc_LR = -( posX - membraneNodeXAddr[nndata5[neighbor]] ); //This is like Rk - Rj, where Rj is the point we are trying to update
                yLoc_LR = -( posY - membraneNodeYAddr[nndata5[neighbor]] );
                zLoc_LR = -( posZ - membraneNodeZAddr[nndata5[neighbor]] );
                double R = sqrt( 
                            (xLoc_LR) * (xLoc_LR) + 
                            (yLoc_LR) * (yLoc_LR) + 
                            (zLoc_LR) * (zLoc_LR) );
            if (R < abs_Rmin){
                // magnitude = 2*epsilon_rep1*
                //                             (1-exp(-epsilon_rep2*(R-abs_Rmin)))*
                //                             (-exp(-epsilon_rep2*(R-abs_Rmin)))*
                //                             (epsilon_rep2/R);

                //         current_forceX5 = -magnitude*xLoc_LR;
                //         current_forceY5 = -magnitude*yLoc_LR;
                //         current_forceZ5 = -magnitude*zLoc_LR;
                double epsilon_rep;
                    if (nodes_in_upperhem[node_id] == 1 && nodes_in_upperhem[nndata5[neighbor]]==1){
                        epsilon_rep = epsilon_rep1;
                    }
                    else if (nodes_in_upperhem[node_id] != 1 || nodes_in_upperhem[nndata5[neighbor]] != 1){
                        epsilon_rep = epsilon_rep2;
                    }
                    current_forceX5 = (epsilon_rep)*(R - abs_Rmin)*(xLoc_LR)/R;
					current_forceY5 = (epsilon_rep)*(R - abs_Rmin)*(yLoc_LR)/R;
					current_forceZ5 = (epsilon_rep)*(R - abs_Rmin)*(zLoc_LR)/R;
                    }
                    else{
                        current_forceX5 = 0.0;
                        current_forceY5 = 0.0;
                        current_forceZ5 = 0.0;
                    }
            }
           else{
               current_forceX5 = 0.0;
                current_forceY5 = 0.0;
                current_forceZ5 = 0.0;
           }
          
           if (nndata6[neighbor] >= 0 && 
                nndata6[neighbor] != node_id &&
                nndata6[neighbor] != neighbor1 &&
                 nndata6[neighbor] != neighbor2 &&
                 nndata6[neighbor] != neighbor3 &&
                 nndata6[neighbor] != neighbor4 &&
                 nndata6[neighbor] != neighbor5 &&
                 nndata6[neighbor] != neighbor6 &&
                 nndata6[neighbor] != neighbor7 &&
                 nndata6[neighbor] != neighbor8 &&
                 nndata6[neighbor] != neighbor9 /*&&
                 nndata6[neighbor] != neighbor10 &&
                 nndata6[neighbor] != neighbor11 &&
                 nndata6[neighbor] != neighbor12*/){
                xLoc_LR = -( posX - membraneNodeXAddr[nndata6[neighbor]] ); //This is like Rk - Rj, where Rj is the point we are trying to update
                yLoc_LR = -( posY - membraneNodeYAddr[nndata6[neighbor]] );
                zLoc_LR = -( posZ - membraneNodeZAddr[nndata6[neighbor]] );
                double R = sqrt( 
                            (xLoc_LR) * (xLoc_LR) + 
                            (yLoc_LR) * (yLoc_LR) + 
                            (zLoc_LR) * (zLoc_LR) );
            if (R < abs_Rmin){
                // magnitude = 2*epsilon_rep1*
                //                                 (1-exp(-epsilon_rep2*(R-abs_Rmin)))*
                //                                 (-exp(-epsilon_rep2*(R-abs_Rmin)))*
                //                                 (epsilon_rep2/R);

                //             current_forceX6 = -magnitude*xLoc_LR;
                //             current_forceY6 = -magnitude*yLoc_LR;
                //             current_forceZ6 = -magnitude*zLoc_LR;
                double epsilon_rep;
                    if (nodes_in_upperhem[node_id] == 1 && nodes_in_upperhem[nndata6[neighbor]]==1){
                        epsilon_rep = epsilon_rep1;
                    }
                    else if (nodes_in_upperhem[node_id] != 1 || nodes_in_upperhem[nndata6[neighbor]] != 1){
                        epsilon_rep = epsilon_rep2;
                    }
                    current_forceX6 = (epsilon_rep)*(R - abs_Rmin)*(xLoc_LR)/R;
					current_forceY6 = (epsilon_rep)*(R - abs_Rmin)*(yLoc_LR)/R;
					current_forceZ6 = (epsilon_rep)*(R - abs_Rmin)*(zLoc_LR)/R;
                    }
                    else{
                        current_forceX6 = 0.0;
                        current_forceY6 = 0.0;
                        current_forceZ6 = 0.0;
                    }
            }
           else{
               current_forceX6 = 0.0;
                current_forceY6 = 0.0;
                current_forceZ6 = 0.0;
           }
           
           if (nndata7[neighbor] >= 0 &&
                nndata7[neighbor] != node_id &&
                nndata7[neighbor] != neighbor1 &&
                 nndata7[neighbor] != neighbor2 &&
                 nndata7[neighbor] != neighbor3 &&
                 nndata7[neighbor] != neighbor4 &&
                 nndata7[neighbor] != neighbor5 &&
                 nndata7[neighbor] != neighbor6 &&
                 nndata7[neighbor] != neighbor7 &&
                 nndata7[neighbor] != neighbor8 &&
                 nndata7[neighbor] != neighbor9 /*&&
                 nndata7[neighbor] != neighbor10 &&
                 nndata7[neighbor] != neighbor11 &&
                 nndata7[neighbor] != neighbor12*/){
                xLoc_LR = -( posX - membraneNodeXAddr[nndata7[neighbor]] ); //This is like Rk - Rj, where Rj is the point we are trying to update
                yLoc_LR = -( posY - membraneNodeYAddr[nndata7[neighbor]] );
                zLoc_LR = -( posZ - membraneNodeZAddr[nndata7[neighbor]] );
                double R = sqrt( 
                            (xLoc_LR) * (xLoc_LR) + 
                            (yLoc_LR) * (yLoc_LR) + 
                            (zLoc_LR) * (zLoc_LR) );
            if (R < abs_Rmin){
                // magnitude = 2*epsilon_rep1*
                //                                 (1-exp(-epsilon_rep2*(R-abs_Rmin)))*
                //                                 (-exp(-epsilon_rep2*(R-abs_Rmin)))*
                //                                 (epsilon_rep2/R);

                //             current_forceX7 = -magnitude*xLoc_LR;
                //             current_forceY7 = -magnitude*yLoc_LR;
                //             current_forceZ7 = -magnitude*zLoc_LR;
                double epsilon_rep;
                    if (nodes_in_upperhem[node_id] == 1 && nodes_in_upperhem[nndata7[neighbor]]==1){
                        epsilon_rep = epsilon_rep1;
                    }
                    else if (nodes_in_upperhem[node_id] != 1 || nodes_in_upperhem[nndata7[neighbor]] != 1){
                        epsilon_rep = epsilon_rep2;
                    }
                    current_forceX7 = (epsilon_rep)*(R - abs_Rmin)*(xLoc_LR)/R;
					current_forceY7 = (epsilon_rep)*(R - abs_Rmin)*(yLoc_LR)/R;
					current_forceZ7 = (epsilon_rep)*(R - abs_Rmin)*(zLoc_LR)/R;
                    }
                    else{
                        current_forceX7 = 0.0;
                        current_forceY7 = 0.0;
                        current_forceZ7 = 0.0;
                    }
            }
           else{
               current_forceX7 = 0.0;
                current_forceY7 = 0.0;
                current_forceZ7 = 0.0;
           }
           
           if (nndata8[neighbor] >= 0 && 
                nndata8[neighbor] != node_id &&
                nndata8[neighbor] != neighbor1 &&
                nndata8[neighbor] != neighbor2 &&
                nndata8[neighbor] != neighbor3 &&
                nndata8[neighbor] != neighbor4 &&
                nndata8[neighbor] != neighbor5 &&
                nndata8[neighbor] != neighbor6 &&
                nndata8[neighbor] != neighbor7 &&
                nndata8[neighbor] != neighbor8 &&
                nndata8[neighbor] != neighbor9 /*&&
                nndata8[neighbor] != neighbor10 &&
                nndata8[neighbor] != neighbor11 &&
                nndata8[neighbor] != neighbor12*/){
                xLoc_LR = -( posX - membraneNodeXAddr[nndata8[neighbor]] ); //This is like Rk - Rj, where Rj is the point we are trying to update
                yLoc_LR = -( posY - membraneNodeYAddr[nndata8[neighbor]] );
                zLoc_LR = -( posZ - membraneNodeZAddr[nndata8[neighbor]] );
                double R = sqrt( 
                            (xLoc_LR) * (xLoc_LR) + 
                            (yLoc_LR) * (yLoc_LR) + 
                            (zLoc_LR) * (zLoc_LR) );
            if (R < abs_Rmin){
                // magnitude = 2*epsilon_rep1*
                //                                 (1-exp(-epsilon_rep2*(R-abs_Rmin)))*
                //                                 (-exp(-epsilon_rep2*(R-abs_Rmin)))*
                //                                 (epsilon_rep2/R);

                //             current_forceX8 = -magnitude*xLoc_LR;
                //             current_forceY8 = -magnitude*yLoc_LR;
                //             current_forceZ8 = -magnitude*zLoc_LR;
                double epsilon_rep;
                    if (nodes_in_upperhem[node_id] == 1 && nodes_in_upperhem[nndata8[neighbor]]==1){
                        epsilon_rep = epsilon_rep1;
                    }
                    else if (nodes_in_upperhem[node_id] != 1 || nodes_in_upperhem[nndata8[neighbor]] != 1){
                        epsilon_rep = epsilon_rep2;
                    }
                    current_forceX8 = (epsilon_rep)*(R - abs_Rmin)*(xLoc_LR)/R;
					current_forceY8 = (epsilon_rep)*(R - abs_Rmin)*(yLoc_LR)/R;
					current_forceZ8 = (epsilon_rep)*(R - abs_Rmin)*(zLoc_LR)/R;
                    }
                    else{
                        current_forceX8 = 0.0;
                        current_forceY8 = 0.0;
                        current_forceZ8 = 0.0;
                    }
            }
           else{
               current_forceX8 = 0.0;
                current_forceY8 = 0.0;
                current_forceZ8 = 0.0;
           }
           
           if (nndata9[neighbor] >= 0 && 
                nndata9[neighbor] != node_id &&
                nndata9[neighbor] != neighbor1 &&
                 nndata9[neighbor] != neighbor2 &&
                 nndata9[neighbor] != neighbor3 &&
                 nndata9[neighbor] != neighbor4 &&
                 nndata9[neighbor] != neighbor5 &&
                 nndata9[neighbor] != neighbor6 &&
                 nndata9[neighbor] != neighbor7 &&
                 nndata9[neighbor] != neighbor8 &&
                 nndata9[neighbor] != neighbor9 /*&&
                 nndata9[neighbor] != neighbor10 &&
                 nndata9[neighbor] != neighbor11 &&
                 nndata9[neighbor] != neighbor12*/ ){
                xLoc_LR = -( posX - membraneNodeXAddr[nndata9[neighbor]] ); //This is like Rk - Rj, where Rj is the point we are trying to update
                yLoc_LR = -( posY - membraneNodeYAddr[nndata9[neighbor]] );
                zLoc_LR = -( posZ - membraneNodeZAddr[nndata9[neighbor]] );
                double R = sqrt( 
                            (xLoc_LR) * (xLoc_LR) + 
                            (yLoc_LR) * (yLoc_LR) + 
                            (zLoc_LR) * (zLoc_LR) );
            if (R < abs_Rmin){
                // magnitude = 2*epsilon_rep1*
                //                                 (1-exp(-epsilon_rep2*(R-abs_Rmin)))*
                //                                 (-exp(-epsilon_rep2*(R-abs_Rmin)))*
                //                                 (epsilon_rep2/R);

                //             current_forceX9 = -magnitude*xLoc_LR;
                //             current_forceY9 = -magnitude*yLoc_LR;
                //             current_forceZ9 = -magnitude*zLoc_LR;
                double epsilon_rep;
                    if (nodes_in_upperhem[node_id] == 1 && nodes_in_upperhem[nndata9[neighbor]]==1){
                        epsilon_rep = epsilon_rep1;
                    }
                    else if (nodes_in_upperhem[node_id] != 1 || nodes_in_upperhem[nndata9[neighbor]] != 1){
                        epsilon_rep = epsilon_rep2;
                    }
                    current_forceX9 = (epsilon_rep)*(R - abs_Rmin)*(xLoc_LR)/R;
					current_forceY9 = (epsilon_rep)*(R - abs_Rmin)*(yLoc_LR)/R;
					current_forceZ9 = (epsilon_rep)*(R - abs_Rmin)*(zLoc_LR)/R;
                    }
                    else{
                        current_forceX9 = 0.0;
                        current_forceY9 = 0.0;
                        current_forceZ9 = 0.0;
                    }
            }
           else{
               current_forceX9 = 0.0;
                current_forceY9 = 0.0;
                current_forceZ9 = 0.0;
           }
           
           /*if (nndata10[neighbor] >= 0 && 
                nndata10[neighbor] != node_id &&
                nndata10[neighbor] != neighbor1 &&
                nndata10[neighbor] != neighbor2 &&
                nndata10[neighbor] != neighbor3 &&
                nndata10[neighbor] != neighbor4 &&
                nndata10[neighbor] != neighbor5 &&
                nndata10[neighbor] != neighbor6 &&
                nndata10[neighbor] != neighbor7 &&
                nndata10[neighbor] != neighbor8 &&
                nndata10[neighbor] != neighbor9 &&
                nndata10[neighbor] != neighbor10 &&
                nndata10[neighbor] != neighbor11 &&
                nndata10[neighbor] != neighbor12){
                xLoc_LR = -( posX - membraneNodeXAddr[nndata10[neighbor]] ); //This is like Rk - Rj, where Rj is the point we are trying to update
                yLoc_LR = -( posY - membraneNodeYAddr[nndata10[neighbor]] );
                zLoc_LR = -( posZ - membraneNodeZAddr[nndata10[neighbor]] );
                double R = sqrt( 
                            (xLoc_LR) * (xLoc_LR) + 
                            (yLoc_LR) * (yLoc_LR) + 
                            (zLoc_LR) * (zLoc_LR) );
            if (R < abs_Rmin){
                magnitude = 2*epsilon_rep1*
                                                (1-exp(-epsilon_rep2*(R-abs_Rmin)))*
                                                (-exp(-epsilon_rep2*(R-abs_Rmin)))*
                                                (epsilon_rep2/R);

                            current_forceX10 = -magnitude*xLoc_LR;
                            current_forceY10 = -magnitude*yLoc_LR;
                            current_forceZ10 = -magnitude*zLoc_LR;
                    }
                    else{
                        current_forceX10 = 0.0;
                        current_forceY10 = 0.0;
                        current_forceZ10 = 0.0;
                    }
            }
           else{
               current_forceX10 = 0.0;
                current_forceY10 = 0.0;
                current_forceZ10 = 0.0;
           }*/
          
           /* if (nndata11[neighbor] >= 0 && 
                nndata11[neighbor] != node_id &&
                nndata11[neighbor] != neighbor1 &&
                nndata11[neighbor] != neighbor2 &&
                nndata11[neighbor] != neighbor3 &&
                nndata11[neighbor] != neighbor4 &&
                nndata11[neighbor] != neighbor5 &&
                nndata11[neighbor] != neighbor6 &&
                nndata11[neighbor] != neighbor7 &&
                nndata11[neighbor] != neighbor8 &&
                nndata11[neighbor] != neighbor9 &&
                nndata11[neighbor] != neighbor10 &&
                nndata11[neighbor] != neighbor11 &&
                nndata11[neighbor] != neighbor12){
                xLoc_LR = -( posX - membraneNodeXAddr[nndata11[neighbor]] ); //This is like Rk - Rj, where Rj is the point we are trying to update
                yLoc_LR = -( posY - membraneNodeYAddr[nndata11[neighbor]] );
                zLoc_LR = -( posZ - membraneNodeZAddr[nndata11[neighbor]] );
                double R = sqrt( 
                            (xLoc_LR) * (xLoc_LR) + 
                            (yLoc_LR) * (yLoc_LR) + 
                            (zLoc_LR) * (zLoc_LR) );
            if (R < abs_Rmin){
                magnitude = 2*epsilon_rep1*
                                                (1-exp(-epsilon_rep2*(R-abs_Rmin)))*
                                                (-exp(-epsilon_rep2*(R-abs_Rmin)))*
                                                (epsilon_rep2/R);

                            current_forceX11 = -magnitude*xLoc_LR;
                            current_forceY11 = -magnitude*yLoc_LR;
                            current_forceZ11 = -magnitude*zLoc_LR;
                    }
                    else{
                        current_forceX11 = 0.0;
                        current_forceY11 = 0.0;
                        current_forceZ11 = 0.0;
                    }
            }
           else{
               current_forceX11 = 0.0;
                current_forceY11 = 0.0;
                current_forceZ11 = 0.0;
           } */
           
           /* if (nndata12[neighbor] >= 0 && 
                nndata12[neighbor] != node_id &&
                nndata12[neighbor] != neighbor1 &&
                nndata12[neighbor] != neighbor2 &&
                nndata12[neighbor] != neighbor3 &&
                nndata12[neighbor] != neighbor4 &&
                nndata12[neighbor] != neighbor5 &&
                nndata12[neighbor] != neighbor6 &&
                nndata12[neighbor] != neighbor7 &&
                nndata12[neighbor] != neighbor8 &&
                nndata12[neighbor] != neighbor9 &&
                nndata12[neighbor] != neighbor10 &&
                nndata12[neighbor] != neighbor11 &&
                nndata12[neighbor] != neighbor12){
                xLoc_LR = -( posX - membraneNodeXAddr[nndata12[neighbor]] ); //This is like Rk - Rj, where Rj is the point we are trying to update
                yLoc_LR = -( posY - membraneNodeYAddr[nndata12[neighbor]] );
                zLoc_LR = -( posZ - membraneNodeZAddr[nndata12[neighbor]] );
                double R = sqrt( 
                            (xLoc_LR) * (xLoc_LR) + 
                            (yLoc_LR) * (yLoc_LR) + 
                            (zLoc_LR) * (zLoc_LR) );
            if (R < abs_Rmin){
                magnitude = 2*epsilon_rep1*
                                                (1-exp(-epsilon_rep2*(R-abs_Rmin)))*
                                                (-exp(-epsilon_rep2*(R-abs_Rmin)))*
                                                (epsilon_rep2/R);

                            current_forceX12 = -magnitude*xLoc_LR;
                            current_forceY12 = -magnitude*yLoc_LR;
                            current_forceZ12 = -magnitude*zLoc_LR;
                    }
                    else{
                        current_forceX12 = 0.0;
                        current_forceY12 = 0.0;
                        current_forceZ12 = 0.0;
                    }
            } 
           else{
               current_forceX12 = 0.0;
                current_forceY12 = 0.0;
                current_forceZ12 = 0.0;
           }*/

        current_forceX += current_forceX1 + current_forceX2 + current_forceX3 + current_forceX4 +
                        current_forceX5 + current_forceX6 + current_forceX7 + current_forceX8 +
                        current_forceX9;// + current_forceX10 ;//+ current_forceX11 + current_forceX12 ;
        current_forceY += current_forceY1 + current_forceY2 + current_forceY3 + current_forceY4 +
                        current_forceY5 + current_forceY6 + current_forceY7 + current_forceY8 +
                        current_forceY9;// + current_forceY10 ;//+ current_forceY11 + current_forceY12 ;
        current_forceZ += current_forceZ1 + current_forceZ2 + current_forceZ3 + current_forceZ4 +
                        current_forceZ5 + current_forceZ6 + current_forceZ7 + current_forceZ8 +
                        current_forceZ9;// + current_forceZ10 ;//+ current_forceZ11 + current_forceZ12 ;    
       //}
    }

        return thrust::make_tuple(current_forceX, current_forceY, current_forceZ);
    }
};

#endif