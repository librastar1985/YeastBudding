#include "System.h"
#include <random>
#include "Edgeswap_test.h"

Edgeswap::Edgeswap(CoordInfoVecs& coordInfoVecs) {
            
    unsigned nnode = coordInfoVecs.nodeLocX.size();
    std::vector<bool> boundary_node_temp(nnode,false);
    for (unsigned i = 0; i < nnode; i++){
        if (coordInfoVecs.edges2Nodes_1[i] == coordInfoVecs.edges2Nodes_2[i]){
            boundary_node_temp[coordInfoVecs.edges2Nodes_1[i]] = true;
            boundary_node_temp[coordInfoVecs.edges2Nodes_2[i]] = true;
        }
    }
    
    //This creates a unsigned vector whose length equals to number of nodes.
    //The initial mesh has every node paired with 6 neighboring nodes. 
    //During the simulation, the number will be updated accordingly. Therefore this has to be moved
    //to another location to avoid re-initialization every time Edgeswap is called.

    std::vector<unsigned> nndata_temp(nnode, 6);
    
    boundary_node = boundary_node_temp;
    nndata = nndata_temp;
};




//The goal is to perform the swap without explicitly creating a copy of the whole data structure.
//This is achieved by extracting the full info of a smaller affected system.
void Edgeswap::edge_swap(
    unsigned iedge, 
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    BendingTriangleInfoVecs& bendingTriangleInfoVecs) {
        
    std::cout<<"here1"<<std::endl;
    unsigned HEAD,TAIL;
    unsigned H0, T0,H1,H2,T1,T2;
    unsigned edge_start, edge_end;
    unsigned a1, b1, c1, a2, b2, c2;
    double temp_bend = 0.0;
        
    if ( coordInfoVecs.edges2Triangles_1[iedge] != coordInfoVecs.edges2Triangles_2[iedge]){
        H0 = coordInfoVecs.edges2Triangles_1[iedge];//index of the 1st triangle to i-th edge
        T0 = coordInfoVecs.edges2Triangles_2[iedge];//index of the 2nd triangle to i-th edge
        edge_start = coordInfoVecs.edges2Nodes_1[iedge];//index of the 1st node of i-th edge
        edge_end = coordInfoVecs.edges2Nodes_2[iedge];//index of the 2nd node of i-th edge

        a1 = coordInfoVecs.triangles2Edges_1[H0];//index of the 1st node of triangle H0
        b1 = coordInfoVecs.triangles2Edges_2[H0];//index of the 2nd node of triangle H0
        c1 = coordInfoVecs.triangles2Edges_3[H0];//index of the 3rd node of triangle H0
        
        a2 = coordInfoVecs.triangles2Edges_1[T0];//index of the 1st node of triangle T0
        b2 = coordInfoVecs.triangles2Edges_2[T0];
        c2 = coordInfoVecs.triangles2Edges_3[T0];
        
    std::cout<<"here2"<<std::endl;
        //Now we identify the edge indices associated with the small subsystem.
        //This gives us the indices for H1, H2, T1, T2 (see the figure below).
        if (a1 != iedge && coordInfoVecs.edges2Nodes_1[a1] == edge_start){H1 = a1;}
        else if (a1 != iedge && coordInfoVecs.edges2Nodes_2[a1] == edge_start){H1 = a1;}
        else if (a1 != iedge && coordInfoVecs.edges2Nodes_1[a1] == edge_end){H2 = a1;}
        else if (a1 != iedge && coordInfoVecs.edges2Nodes_2[a1] == edge_end){H2 = a1;}

        if (b1 != iedge && coordInfoVecs.edges2Nodes_1[b1] == edge_start){H1 = b1;}
        else if (b1 != iedge && coordInfoVecs.edges2Nodes_2[b1] == edge_start){H1 = b1;}
        else if (b1 != iedge && coordInfoVecs.edges2Nodes_1[b1] == edge_end){H2 = b1;}
        else if (b1 != iedge && coordInfoVecs.edges2Nodes_2[b1] == edge_end){H2 = b1;}

        if (c1 != iedge && coordInfoVecs.edges2Nodes_1[c1] == edge_start){H1 = c1;}
        else if (c1 != iedge && coordInfoVecs.edges2Nodes_2[c1] == edge_start){H1 = c1;}
        else if (c1 != iedge && coordInfoVecs.edges2Nodes_1[c1] == edge_end){H2 = c1;}
        else if (c1 != iedge && coordInfoVecs.edges2Nodes_2[c1] == edge_end){H2 = c1;}
        
        if (a2 != iedge && coordInfoVecs.edges2Nodes_1[a2] == edge_start){T1 = a2;}
        else if (a2 != iedge && coordInfoVecs.edges2Nodes_2[a2] == edge_start){T1 = a2;}
        else if (a2 != iedge && coordInfoVecs.edges2Nodes_1[a2] == edge_end){T2 = a2;}
        else if (a2 != iedge && coordInfoVecs.edges2Nodes_2[a2] == edge_end){T2 = a2;}

        if (b2 != iedge && coordInfoVecs.edges2Nodes_1[b2] == edge_start){T1 = b2;}
        else if (b2 != iedge && coordInfoVecs.edges2Nodes_2[b2] == edge_start){T1 = b2;}
        else if (b2 != iedge && coordInfoVecs.edges2Nodes_1[b2] == edge_end){T2 = b2;}
        else if (b2 != iedge && coordInfoVecs.edges2Nodes_2[b2] == edge_end){T2 = b2;}

        if (c2 != iedge && coordInfoVecs.edges2Nodes_1[c2] == edge_start){T1 = c2;}
        else if (c2 != iedge && coordInfoVecs.edges2Nodes_2[c2] == edge_start){T1 = c2;}
        else if (c2 != iedge && coordInfoVecs.edges2Nodes_1[c2] == edge_end){T2 = c2;}
        else if (c2 != iedge && coordInfoVecs.edges2Nodes_2[c2] == edge_end){T2 = c2;}

        //Now search for the associated 

        unsigned CANDIDATE1_1 = coordInfoVecs.triangles2Nodes_1[H0];
        unsigned CANDIDATE1_2 = coordInfoVecs.triangles2Nodes_2[H0];
        unsigned CANDIDATE1_3 = coordInfoVecs.triangles2Nodes_3[H0];
        unsigned CANDIDATE2_1 = coordInfoVecs.triangles2Nodes_1[T0];
        unsigned CANDIDATE2_2 = coordInfoVecs.triangles2Nodes_2[T0];
        unsigned CANDIDATE2_3 = coordInfoVecs.triangles2Nodes_3[T0];
        std::cout<< CANDIDATE1_1<< " "<< CANDIDATE1_2 << " "<<CANDIDATE1_3 << " "<< CANDIDATE2_1 << " "<< CANDIDATE2_2 << " "<< CANDIDATE2_3 << std::endl;
        if ((CANDIDATE1_1 != edge_start) 
            && (CANDIDATE1_1 != edge_end)) {
            HEAD = CANDIDATE1_1;
        }
        else if ((CANDIDATE1_2 != edge_start) && (CANDIDATE1_2 != edge_end)){HEAD = CANDIDATE1_2;}
        else if (CANDIDATE1_3 != edge_start && CANDIDATE1_3 != edge_end){HEAD = CANDIDATE1_3;}
        else {std::cout<<"head not set" <<std::endl;}

        if (CANDIDATE2_1 != edge_start && CANDIDATE2_1 != edge_end){TAIL = CANDIDATE2_1;}
        else if (CANDIDATE2_2 != edge_start && CANDIDATE2_2 != edge_end){TAIL = CANDIDATE2_2;}
        else if (CANDIDATE2_3 != edge_start && CANDIDATE2_3 != edge_end){TAIL = CANDIDATE2_3;}
        else {std::cout<<"tail not set" <<std::endl;}


        unsigned temp_edges2Nodes_2 = HEAD;
        std::cout<<"head tail in loop: "<< HEAD << " "<< TAIL <<std::endl;
        //The small subsystem we will be working with is
        //          
        //           edge_start
        //    T10    *   |    *     H10
        //         T1    |     H1
        //        *      |       *
        //    TAIL   T0  |  H0    HEAD
        //        *      |       *
        //         T2    |     H2
        //    T20    *   v    *     H20
        //            edge_end
        //
        //H10 is the triangle sharing the same edge H1 with triangle H0.

        //energy E_0 calculation
        //Since nodes are NOT moved and area is not changed, we only need to calculate 
        //linear spring energy and bending energy.
        //Furthermore, since linear spring energy will only be nontrivial for the edge swapped,
        //we can condense the linear spring energy computation to only one edge.
        //Bending energy is more complicated due to the need of unit normals.
        
    std::cout<<"here3"<<std::endl;
        std::vector<unsigned> edges_iteration(5);
        edges_iteration[0] = iedge;
        edges_iteration[1] = H1;
        edges_iteration[2] = H2;
        edges_iteration[3] = T1;
        edges_iteration[4] = T2;
            for (unsigned j = 0; j < 5; j++){
               
                unsigned Tri1 = coordInfoVecs.edges2Triangles_1[edges_iteration[j]];//index of the 1st triangle
                unsigned Tri2 = coordInfoVecs.edges2Triangles_2[edges_iteration[j]];
                //unsigned id_k = coordInfoVecs.edges2Nodes_1[edges_iteration[j]];
                //unsigned id_i = coordInfoVecs.edges2Nodes_2[edges_iteration[j]];

                double vec1x, vec1y, vec1z, vec2x, vec2y, vec2z;
                if (Tri1 != Tri2) {
                    unsigned Tri1_n1 = coordInfoVecs.triangles2Nodes_1[Tri1];
                    unsigned Tri1_n2 = coordInfoVecs.triangles2Nodes_2[Tri1];
                    unsigned Tri1_n3 = coordInfoVecs.triangles2Nodes_3[Tri1];
                    vec1x = coordInfoVecs.nodeLocX[Tri1_n2] - coordInfoVecs.nodeLocX[Tri1_n1];
                    vec1y = coordInfoVecs.nodeLocY[Tri1_n2] - coordInfoVecs.nodeLocY[Tri1_n1];
                    vec1z = coordInfoVecs.nodeLocZ[Tri1_n2] - coordInfoVecs.nodeLocZ[Tri1_n1];
                    vec2x = coordInfoVecs.nodeLocX[Tri1_n3] - coordInfoVecs.nodeLocX[Tri1_n1];
                    vec2y = coordInfoVecs.nodeLocY[Tri1_n3] - coordInfoVecs.nodeLocY[Tri1_n1];
                    vec2z = coordInfoVecs.nodeLocZ[Tri1_n3] - coordInfoVecs.nodeLocZ[Tri1_n1];
                    std::vector<double> N1(3);
                    N1[0] = vec1y*vec2z - vec2y*vec1z;
                    N1[1] = -(vec1x*vec2z - vec2x*vec1z);
                    N1[2] = vec1x*vec2y - vec2x*vec1y;
                    double nN1 = sqrt(pow(N1[0],2)+pow(N1[1],2)+pow(N1[2],2));

                    unsigned Tri2_n1 = coordInfoVecs.triangles2Nodes_1[Tri2];
                    unsigned Tri2_n2 = coordInfoVecs.triangles2Nodes_2[Tri2];
                    unsigned Tri2_n3 = coordInfoVecs.triangles2Nodes_3[Tri2];
                    vec1x = coordInfoVecs.nodeLocX[Tri2_n2] - coordInfoVecs.nodeLocX[Tri2_n1];
                    vec1y = coordInfoVecs.nodeLocY[Tri2_n2] - coordInfoVecs.nodeLocY[Tri2_n1];
                    vec1z = coordInfoVecs.nodeLocZ[Tri2_n2] - coordInfoVecs.nodeLocZ[Tri2_n1];
                    vec2x = coordInfoVecs.nodeLocX[Tri2_n3] - coordInfoVecs.nodeLocX[Tri2_n1];
                    vec2y = coordInfoVecs.nodeLocY[Tri2_n3] - coordInfoVecs.nodeLocY[Tri2_n1];
                    vec2z = coordInfoVecs.nodeLocZ[Tri2_n3] - coordInfoVecs.nodeLocZ[Tri2_n1];
                    std::vector<double> N2(3);
                    N2[0] = vec1y*vec2z - vec2y*vec1z;
                    N2[1] = -(vec1x*vec2z - vec2x*vec1z);
                    N2[2] = vec1x*vec2y - vec2x*vec1y; 
                    double nN2 = sqrt(pow(N2[0],2)+pow(N2[1],2)+pow(N2[2],2));;

                    double cosAngle = (N1[0]*N2[0] + N1[1]*N2[1] + N1[2]*N2[2])/(nN1*nN2);
                   
    std::cout<<"here4"<<std::endl; 
    
    std::cout<<"head tailin loop 2: "<< HEAD << " "<< TAIL<<std::endl;
                    if (cosAngle > 1.0) {
                        cosAngle = 1.0;
                    }
                    else if (cosAngle < -1.0){
                        cosAngle = -1.0;
                    }

                    double theta_current = acos( cosAngle );
                    
                    double local_energy = bendingTriangleInfoVecs.spring_constant * (1 - cos(theta_current - 0.0) );
                    temp_bend = temp_bend + local_energy;
                }
            }
        }       

        double bend_0 = temp_bend;
        //
        double linear_0 = linearSpringInfoVecs.spring_constant*(sqrt(
            pow(coordInfoVecs.nodeLocX[edge_end] - coordInfoVecs.nodeLocX[edge_start], 2.0) + 
            pow(coordInfoVecs.nodeLocY[edge_end] - coordInfoVecs.nodeLocY[edge_start], 2.0) + 
            pow(coordInfoVecs.nodeLocZ[edge_end] - coordInfoVecs.nodeLocZ[edge_start], 2.0)) - linearSpringInfoVecs.edge_initial_length[0]);
        
        double E_0 = linear_0 + bend_0;

        
        
    std::cout<<"here5"<<std::endl;
        //Flip the edge, build the data structure for the smaller system.
        bool BAD_CHOICE;
        unsigned temp_edges2Nodes_1 = TAIL;
        unsigned temp_edges2Nodes_2 = HEAD;
    std::cout<<"head tail: "<< HEAD << " "<< TAIL<<std::endl;
    
    std::cout<<"head tailnndata: "<< nndata[HEAD] << " "<< nndata[TAIL]<<std::endl;

        unsigned temp_nndata_HEAD = nndata[HEAD] + 1;
        unsigned temp_nndata_TAIL = nndata[TAIL] + 1;
        unsigned temp_nndata_edge_start = nndata[edge_start] - 1;
        unsigned temp_nndata_edge_end = nndata[edge_end] - 1;
        if (boundary_node[HEAD] == false && temp_nndata_HEAD < 3){
            BAD_CHOICE = true;
        }
        else if (boundary_node[TAIL] == false && temp_nndata_TAIL < 3){
            BAD_CHOICE = true;
        }
        else if (boundary_node[edge_start] == false && temp_nndata_edge_start < 3){
            BAD_CHOICE = true;
        }
        else if (boundary_node[edge_end] == false && temp_nndata_edge_end < 3){
            BAD_CHOICE = true;
        }
        else {
            BAD_CHOICE = false;
        }

        std::cout<<"badchocie:" << BAD_CHOICE<<std::endl;

        if (BAD_CHOICE == false) {
            /*temp_edges2Nodes_1[iedge] = TAIL;
            temp_edges2Nodes_2[iedge] = HEAD;
            temp_nndata[HEAD] = temp_nndata[HEAD] + 1;
            temp_nndata[TAIL] = temp_nndata[TAIL] + 1;
            temp_nndata[edge_start] = temp_nndata[edge_start] - 1;
            temp_nndata[edge_end] = temp_nndata[edge_end] - 1;*/

            //The labeling of neighboring edge will as follows after swap:
            //          
            //           edge_start
            //           *        *
            //         T1    H0     H1
            //        *              *
            //      TAIL ----------> HEAD
            //        *              *
            //         T2    T0    H2
            //           *        *
            //            edge_end
            //
            //Now we will update the temporary data structure to accomodate the edgeswap
            
            //Update the new triangles2Nodes information
            /*temp_triangles2Nodes_1[H0] = HEAD;
            temp_triangles2Nodes_2[H0] = edge_start;
            temp_triangles2Nodes_3[H0] = TAIL;
            temp_triangles2Nodes_1[T0] = HEAD;
            temp_triangles2Nodes_2[T0] = TAIL;
            temp_triangles2Nodes_3[T0] = edge_end;*/

            std::cout<<"h1: "<< H1 << " " << H2 << std::endl;
            //Creating vectors to compute the normal vectors under the swapped configuration.
            unsigned H1t1 = coordInfoVecs.edges2Triangles_1[H1];
            unsigned H1t2 = coordInfoVecs.edges2Triangles_2[H1]; //These are the associated triangles to edge H1
            //For the following if statement, we identify the triangles that are affected by the edge-swap.
            //Since we do not know the exact index of the affected triangle, we use the if statement to consider possible cases.
            //This gives us the vectors necessary to compute unit normal vectors required for bending energy.
            

    std::cout<<"here6"<<std::endl;
            double H1t1_vec1x;
            double H1t1_vec1y;
            double H1t1_vec1z;
            double H1t1_vec2x;
            double H1t1_vec2y;
            double H1t1_vec2z;
            double H1t2_vec1x;
            double H1t2_vec1y;
            double H1t2_vec1z;
            double H1t2_vec2x;
            double H1t2_vec2y;
            double H1t2_vec2z;

            double H2t1_vec1x;
            double H2t1_vec1y;
            double H2t1_vec1z;
            double H2t1_vec2x;
            double H2t1_vec2y;
            double H2t1_vec2z;
            double H2t2_vec1x;
            double H2t2_vec1y;
            double H2t2_vec1z;
            double H2t2_vec2x;
            double H2t2_vec2y;
            double H2t2_vec2z;

            double T1t2_vec1x;
            double T1t2_vec1y;
            double T1t2_vec1z;
            double T1t2_vec2x;
            double T1t2_vec2y;
            double T1t2_vec2z;
            double T1t1_vec1x;
            double T1t1_vec1y;
            double T1t1_vec1z;
            double T1t1_vec2x;
            double T1t1_vec2y;
            double T1t1_vec2z;

            double T2t2_vec1x;
            double T2t2_vec1y;
            double T2t2_vec1z;
            double T2t2_vec2x;
            double T2t2_vec2y;
            double T2t2_vec2z;
            double T2t1_vec1x;
            double T2t1_vec1y;
            double T2t1_vec1z;
            double T2t1_vec2x;
            double T2t1_vec2y;
            double T2t1_vec2z;
            
            if (H1t1 == H0){H1t1_vec1x = coordInfoVecs.nodeLocX[edge_start] - coordInfoVecs.nodeLocX[HEAD];
                            H1t1_vec1y = coordInfoVecs.nodeLocY[edge_start] - coordInfoVecs.nodeLocY[HEAD];
                            H1t1_vec1z = coordInfoVecs.nodeLocZ[edge_start] - coordInfoVecs.nodeLocZ[HEAD];
                            H1t1_vec2x = coordInfoVecs.nodeLocX[TAIL] - coordInfoVecs.nodeLocX[HEAD];
                            H1t1_vec2y = coordInfoVecs.nodeLocY[TAIL] - coordInfoVecs.nodeLocY[HEAD];
                            H1t1_vec2z = coordInfoVecs.nodeLocZ[TAIL] - coordInfoVecs.nodeLocZ[HEAD];
                            H1t2_vec1x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_2[H1t2]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[H1t2]];
                            H1t2_vec1y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_2[H1t2]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[H1t2]];
                            H1t2_vec1z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[H1t2]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[H1t2]];
                            H1t2_vec2x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_3[H1t2]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[H1t2]];
                            H1t2_vec2y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_3[H1t2]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[H1t2]];
                            H1t2_vec2z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[H1t2]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[H1t2]];
                            }
            else if (H1t2 == H0){H1t2_vec1x = coordInfoVecs.nodeLocX[edge_start] - coordInfoVecs.nodeLocX[HEAD];
                            H1t2_vec1y = coordInfoVecs.nodeLocY[edge_start] - coordInfoVecs.nodeLocY[HEAD];
                            H1t2_vec1z = coordInfoVecs.nodeLocZ[edge_start] - coordInfoVecs.nodeLocZ[HEAD];
                            H1t2_vec2x = coordInfoVecs.nodeLocX[TAIL] - coordInfoVecs.nodeLocX[HEAD];
                            H1t2_vec2y = coordInfoVecs.nodeLocY[TAIL] - coordInfoVecs.nodeLocY[HEAD];
                            H1t2_vec2z = coordInfoVecs.nodeLocZ[TAIL] - coordInfoVecs.nodeLocZ[HEAD];
                            H1t1_vec1x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_2[H1t1]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[H1t1]];
                            H1t1_vec1y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_2[H1t1]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[H1t1]];
                            H1t1_vec1z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[H1t1]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[H1t1]];
                            H1t1_vec2x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_3[H1t1]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[H1t1]];
                            H1t1_vec2y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_3[H1t1]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[H1t1]];
                            H1t1_vec2z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[H1t1]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[H1t1]];
                            }
            unsigned H2t1 = coordInfoVecs.edges2Triangles_1[H2];
            unsigned H2t2 = coordInfoVecs.edges2Triangles_2[H2]; //These are the associated triangles to edge H2        
            if (H2t1 == H0){//In this case H2t1 turns into T0.
                            H2t1_vec1x = coordInfoVecs.nodeLocX[TAIL] - coordInfoVecs.nodeLocX[HEAD];
                            H2t1_vec1y = coordInfoVecs.nodeLocY[TAIL] - coordInfoVecs.nodeLocY[HEAD];
                            H2t1_vec1z = coordInfoVecs.nodeLocZ[TAIL] - coordInfoVecs.nodeLocZ[HEAD];
                            H2t1_vec2x = coordInfoVecs.nodeLocX[edge_end] - coordInfoVecs.nodeLocX[HEAD];
                            H2t1_vec2y = coordInfoVecs.nodeLocY[edge_end] - coordInfoVecs.nodeLocY[HEAD];
                            H2t1_vec2z = coordInfoVecs.nodeLocZ[edge_end] - coordInfoVecs.nodeLocZ[HEAD];
                            H2t2_vec1x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_2[H2t2]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[H2t2]];
                            H2t2_vec1y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_2[H2t2]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[H2t2]];
                            H2t2_vec1z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[H2t2]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[H2t2]];
                            H2t2_vec2x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_3[H2t2]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[H2t2]];
                            H2t2_vec2y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_3[H2t2]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[H2t2]];
                            H2t2_vec2z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[H2t2]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[H2t2]];
                            }
            else if (H2t2 == H0){//In this case H2t2 tunrs into T0
                            H2t2_vec1x = coordInfoVecs.nodeLocX[TAIL] - coordInfoVecs.nodeLocX[HEAD];
                            H2t2_vec1y = coordInfoVecs.nodeLocY[TAIL] - coordInfoVecs.nodeLocY[HEAD];
                            H2t2_vec1z = coordInfoVecs.nodeLocZ[TAIL] - coordInfoVecs.nodeLocZ[HEAD];
                            H2t2_vec2x = coordInfoVecs.nodeLocX[edge_end] - coordInfoVecs.nodeLocX[HEAD];
                            H2t2_vec2y = coordInfoVecs.nodeLocY[edge_end] - coordInfoVecs.nodeLocY[HEAD];
                            H2t2_vec2z = coordInfoVecs.nodeLocZ[edge_end] - coordInfoVecs.nodeLocZ[HEAD];
                            H2t1_vec1x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_2[H2t1]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[H2t1]];
                            H2t1_vec1y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_2[H2t1]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[H2t1]];
                            H2t1_vec1z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[H2t1]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[H2t1]];
                            H2t1_vec2x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_3[H2t1]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[H2t1]];
                            H2t1_vec2y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_3[H2t1]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[H2t1]];
                            H2t1_vec2z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[H2t1]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[H2t1]];
                            }
            unsigned T1t1 = coordInfoVecs.edges2Triangles_1[T1];
            unsigned T1t2 = coordInfoVecs.edges2Triangles_2[T1];
            if (T1t1 == T0){//In this case T1t1 turns into H0.
                            T1t1_vec1x = coordInfoVecs.nodeLocX[HEAD] - coordInfoVecs.nodeLocX[TAIL];
                            T1t1_vec1y = coordInfoVecs.nodeLocY[HEAD] - coordInfoVecs.nodeLocY[TAIL];
                            T1t1_vec1z = coordInfoVecs.nodeLocZ[HEAD] - coordInfoVecs.nodeLocZ[TAIL];
                            T1t1_vec2x = coordInfoVecs.nodeLocX[edge_start] - coordInfoVecs.nodeLocX[TAIL];
                            T1t1_vec2y = coordInfoVecs.nodeLocY[edge_start] - coordInfoVecs.nodeLocY[TAIL];
                            T1t1_vec2z = coordInfoVecs.nodeLocZ[edge_start] - coordInfoVecs.nodeLocZ[TAIL];
                            T1t2_vec1x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_2[T1t2]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[T1t2]];
                            T1t2_vec1y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_2[T1t2]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[T1t2]];
                            T1t2_vec1z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[T1t2]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[T1t2]];
                            T1t2_vec2x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_3[T1t2]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[T1t2]];
                            T1t2_vec2y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_3[T1t2]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[T1t2]];
                            T1t2_vec2z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[T1t2]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[T1t2]];
                            }
            else if (T1t2 == T0){//In this case T1t2 turns into H0.
                            T1t2_vec1x = coordInfoVecs.nodeLocX[HEAD] - coordInfoVecs.nodeLocX[TAIL];
                            T1t2_vec1y = coordInfoVecs.nodeLocY[HEAD] - coordInfoVecs.nodeLocY[TAIL];
                            T1t2_vec1z = coordInfoVecs.nodeLocZ[HEAD] - coordInfoVecs.nodeLocZ[TAIL];
                            T1t2_vec2x = coordInfoVecs.nodeLocX[edge_start] - coordInfoVecs.nodeLocX[TAIL];
                            T1t2_vec2y = coordInfoVecs.nodeLocY[edge_start] - coordInfoVecs.nodeLocY[TAIL];
                            T1t2_vec2z = coordInfoVecs.nodeLocZ[edge_start] - coordInfoVecs.nodeLocZ[TAIL];
                            T1t1_vec1x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_2[T1t1]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[T1t1]];
                            T1t1_vec1y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_2[T1t1]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[T1t1]];
                            T1t1_vec1z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[T1t1]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[T1t1]];
                            T1t1_vec2x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_3[T1t1]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[T1t1]];
                            T1t1_vec2y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_3[T1t1]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[T1t1]];
                            T1t1_vec2z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[T1t1]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[T1t1]];
                            }
            unsigned T2t1 = coordInfoVecs.edges2Triangles_1[T2];
            unsigned T2t2 = coordInfoVecs.edges2Triangles_2[T2];
            if (T2t1 == T0){T2t1_vec1x = coordInfoVecs.nodeLocX[edge_end] - coordInfoVecs.nodeLocX[TAIL];
                            T2t1_vec1y = coordInfoVecs.nodeLocY[edge_end] - coordInfoVecs.nodeLocY[TAIL];
                            T2t1_vec1z = coordInfoVecs.nodeLocZ[edge_end] - coordInfoVecs.nodeLocZ[TAIL];
                            T2t1_vec2x = coordInfoVecs.nodeLocX[HEAD] - coordInfoVecs.nodeLocX[TAIL];
                            T2t1_vec2y = coordInfoVecs.nodeLocY[HEAD] - coordInfoVecs.nodeLocY[TAIL];
                            T2t1_vec2z = coordInfoVecs.nodeLocZ[HEAD] - coordInfoVecs.nodeLocZ[TAIL];
                            T2t2_vec1x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_2[T2t2]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[T2t2]];
                            T2t2_vec1y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_2[T2t2]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[T2t2]];
                            T2t2_vec1z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[T2t2]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[T2t2]];
                            T2t2_vec2x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_3[T2t2]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[T2t2]];
                            T2t2_vec2y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_3[T2t2]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[T2t2]];
                            T2t2_vec2z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[T2t2]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[T2t2]];
                            }
            else if (T2t2 == T0){
                            T2t2_vec1x = coordInfoVecs.nodeLocX[edge_end] - coordInfoVecs.nodeLocX[TAIL];
                            T2t2_vec1y = coordInfoVecs.nodeLocY[edge_end] - coordInfoVecs.nodeLocY[TAIL];
                            T2t2_vec1z = coordInfoVecs.nodeLocZ[edge_end] - coordInfoVecs.nodeLocZ[TAIL];
                            T2t2_vec2x = coordInfoVecs.nodeLocX[HEAD] - coordInfoVecs.nodeLocX[TAIL];
                            T2t2_vec2y = coordInfoVecs.nodeLocY[HEAD] - coordInfoVecs.nodeLocY[TAIL];
                            T2t2_vec2z = coordInfoVecs.nodeLocZ[HEAD] - coordInfoVecs.nodeLocZ[TAIL];
                            T2t1_vec1x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_2[T2t1]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[T2t1]];
                            T2t1_vec1y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_2[T2t1]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[T2t1]];
                            T2t1_vec1z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[T2t1]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[T2t1]];
                            T2t1_vec2x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_3[T2t1]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[T2t1]];
                            T2t1_vec2y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_3[T2t1]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[T2t1]];
                            T2t1_vec2z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[T2t1]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[T2t1]];
                            }
            
            //First calculate the linear spring energy due to edge-swap.
            double linear_1 = linearSpringInfoVecs.spring_constant*(sqrt(
                pow(coordInfoVecs.nodeLocX[HEAD] - coordInfoVecs.nodeLocX[TAIL], 2) + 
                pow(coordInfoVecs.nodeLocY[HEAD] - coordInfoVecs.nodeLocY[TAIL], 2) + 
                pow(coordInfoVecs.nodeLocZ[HEAD] - coordInfoVecs.nodeLocZ[TAIL], 2)) - linearSpringInfoVecs.edge_initial_length[0]);
            
            //WARNING: RESET BENDING COUNTER
            temp_bend = 0.0;


    std::cout<<"here7"<<std::endl;
            double N1_vec1x, N1_vec1y, N1_vec1z, N1_vec2x, N1_vec2y, N1_vec2z, N2_vec1x, N2_vec1y, N2_vec1z, N2_vec2x, N2_vec2y, N2_vec2z;
            for (unsigned j = 0; j < 5; j++){
                    if (j == 0){
                        N1_vec1x = coordInfoVecs.nodeLocX[HEAD] - coordInfoVecs.nodeLocX[TAIL];//x component of the 1st vector to calculate N1
                        N1_vec1y = coordInfoVecs.nodeLocY[HEAD] - coordInfoVecs.nodeLocY[TAIL];
                        N1_vec1z = coordInfoVecs.nodeLocZ[HEAD] - coordInfoVecs.nodeLocZ[TAIL];
                        N1_vec2x = coordInfoVecs.nodeLocX[edge_start] - coordInfoVecs.nodeLocX[TAIL];//x component of the 2nd vector to calculate N1
                        N1_vec2y = coordInfoVecs.nodeLocY[edge_start] - coordInfoVecs.nodeLocY[TAIL];
                        N1_vec2z = coordInfoVecs.nodeLocZ[edge_start] - coordInfoVecs.nodeLocZ[TAIL];
                        N2_vec1x = coordInfoVecs.nodeLocX[edge_end] - coordInfoVecs.nodeLocX[TAIL];
                        N2_vec1y = coordInfoVecs.nodeLocY[edge_end] - coordInfoVecs.nodeLocY[TAIL];
                        N2_vec1z = coordInfoVecs.nodeLocZ[edge_end] - coordInfoVecs.nodeLocZ[TAIL];
                        N2_vec1x = coordInfoVecs.nodeLocX[HEAD] - coordInfoVecs.nodeLocX[TAIL];
                        N2_vec1y = coordInfoVecs.nodeLocY[HEAD] - coordInfoVecs.nodeLocY[TAIL];
                        N2_vec1z = coordInfoVecs.nodeLocZ[HEAD] - coordInfoVecs.nodeLocZ[TAIL];
                    }
                    else if (j == 1){
                        N1_vec1x = H1t1_vec1x;
                        N1_vec1y = H1t1_vec1y;
                        N1_vec1z = H1t1_vec1z;
                        N1_vec2x = H1t1_vec2x;
                        N1_vec2y = H1t1_vec2y;
                        N1_vec2z = H1t1_vec2z;
                        N2_vec1x = H1t2_vec1x;
                        N2_vec1y = H1t2_vec1y;
                        N2_vec1z = H1t2_vec1z;
                        N2_vec2x = H1t2_vec2x;
                        N2_vec2y = H1t2_vec2y;
                        N2_vec2z = H1t2_vec2z;
                    }
                    else if (j == 2){
                        N1_vec1x = H2t1_vec1x;
                        N1_vec1y = H2t1_vec1y;
                        N1_vec1z = H2t1_vec1z;
                        N1_vec2x = H2t1_vec2x;
                        N1_vec2y = H2t1_vec2y;
                        N1_vec2z = H2t1_vec2z;
                        N2_vec1x = H2t2_vec1x;
                        N2_vec1y = H2t2_vec1y;
                        N2_vec1z = H2t2_vec1z;
                        N2_vec2x = H2t2_vec2x;
                        N2_vec2y = H2t2_vec2y;
                        N2_vec2z = H2t2_vec2z;
                    }
                    else if (j == 3){
                        N1_vec1x = T1t1_vec1x;
                        N1_vec1y = T1t1_vec1y;
                        N1_vec1z = T1t1_vec1z;
                        N1_vec2x = T1t1_vec2x;
                        N1_vec2y = T1t1_vec2y;
                        N1_vec2z = T1t1_vec2z;
                        N2_vec1x = T1t2_vec1x;
                        N2_vec1y = T1t2_vec1y;
                        N2_vec1z = T1t2_vec1z;
                        N2_vec2x = T1t2_vec2x;
                        N2_vec2y = T1t2_vec2y;
                        N2_vec2z = T1t2_vec2z;
                    }
                    else if (j == 4){
                        N1_vec1x = T2t1_vec1x;
                        N1_vec1y = T2t1_vec1y;
                        N1_vec1z = T2t1_vec1z;
                        N1_vec2x = T2t1_vec2x;
                        N1_vec2y = T2t1_vec2y;
                        N1_vec2z = T2t1_vec2z;
                        N2_vec1x = T2t2_vec1x;
                        N2_vec1y = T2t2_vec1y;
                        N2_vec1z = T2t2_vec1z;
                        N2_vec2x = T2t2_vec2x;
                        N2_vec2y = T2t2_vec2y;
                        N2_vec2z = T2t2_vec2z;
                    }                        
                    std::vector<double> N1(3);
                    N1[0] = N1_vec1y*N1_vec2z - N1_vec2y*N1_vec1z;
                    N1[1] = -(N1_vec1x*N1_vec2z - N1_vec2x*N1_vec1z);
                    N1[2] = N1_vec1x*N1_vec2y - N1_vec2x*N1_vec1y;
                    double nN1 = sqrt(pow(N1[0],2)+pow(N1[1],2)+pow(N1[2],2));

    std::cout<<"here8"<<std::endl;
                    std::vector<double> N2(3);
                    N2[0] = N2_vec1y*N2_vec2z - N2_vec2y*N2_vec1z;
                    N2[1] = -(N2_vec1x*N2_vec2z - N2_vec2x*N2_vec1z);
                    N2[2] = N2_vec1x*N2_vec2y - N2_vec2x*N2_vec1y; 
                    double nN2 = sqrt(pow(N2[0],2)+pow(N2[1],2)+pow(N2[2],2));;

                    double cosAngle = (N1[0]*N2[0] + N1[1]*N2[1] + N1[2]*N2[2])/(nN1*nN2);
                    
                    if (cosAngle > 1.0) {
                        cosAngle = 1.0;
                    }
                    else if (cosAngle < -1.0){
                        cosAngle = -1.0;
                    }

                    double theta_current = acos( cosAngle );
                    
                    double local_energy = bendingTriangleInfoVecs.spring_constant * (1 - cos(theta_current - 0.0) );
                    temp_bend = temp_bend + local_energy;
                }
            double bend_1 = temp_bend;
            double E_1 = linear_1 + bend_1;
        
        //Now compute the Boltzmann factor to determine if a swap occurs.
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> dis(0.0, 1.0);
        double random_number = dis(gen);
        double Edif = (E_1 - E_0);
        
//////////////////////////////////////////////////////////////
        //WARNING< CHANGE BACK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//////////////////////////////////////////////////////////////
        double prob = 1.0;//exp(-Edif/generalParams.kT);

        bool ACCEPT2;

    std::cout<<"here9"<<std::endl;
        if (prob >= 1){ACCEPT2 = true;}
        else if (prob < 1 && random_number <= prob){ACCEPT2 = true;}
        else if (prob < 1 && random_number > prob){ACCEPT2 = false;}

        //Perform real update
        if (ACCEPT2 == true){
            coordInfoVecs.triangles2Nodes_1[H0] = HEAD;
            coordInfoVecs.triangles2Nodes_2[H0] = edge_start;
            coordInfoVecs.triangles2Nodes_3[H0] = TAIL;
            coordInfoVecs.triangles2Nodes_1[T0] = HEAD;
            coordInfoVecs.triangles2Nodes_2[T0] = TAIL;
            coordInfoVecs.triangles2Nodes_3[T0] = edge_end;
            coordInfoVecs.edges2Nodes_1[iedge] = TAIL;
            coordInfoVecs.edges2Nodes_2[iedge] = HEAD;
            H1t1 = coordInfoVecs.edges2Triangles_1[H1];
            H1t2 = coordInfoVecs.edges2Triangles_2[H1]; //These are the associated triangles to edge H1
            if (H1t1 == H0){coordInfoVecs.edges2Triangles_1[H1] = H0;}
            if (H1t2 == H0){coordInfoVecs.edges2Triangles_2[H1] = H0;}
            H2t1 = coordInfoVecs.edges2Triangles_1[H2];
            H2t2 = coordInfoVecs.edges2Triangles_2[H2]; //These are the associated triangles to edge H2        
            if (H2t1 == H0){coordInfoVecs.edges2Triangles_1[H2] = T0;}
            if (H2t2 == H0){coordInfoVecs.edges2Triangles_2[H2] = T0;}
            T1t1 = coordInfoVecs.edges2Triangles_1[T1];
            T1t2 = coordInfoVecs.edges2Triangles_2[T1];
            if (T1t1 == T0){coordInfoVecs.edges2Triangles_1[T1] = H0;}
            if (T1t2 == T0){coordInfoVecs.edges2Triangles_2[T1] = H0;}
            T2t1 = coordInfoVecs.edges2Triangles_1[T2];
            T2t2 = coordInfoVecs.edges2Triangles_2[T2];
            if (T2t1 == T0){coordInfoVecs.edges2Triangles_1[T2] = T0;}
            if (T2t2 == T0){coordInfoVecs.edges2Triangles_2[T2] = T0;}
            nndata[HEAD] += 1;
            nndata[TAIL] += 1;
            nndata[edge_start] -= 1;
            nndata[edge_end] -= 1;
            
    std::cout<<"here10"<<std::endl;
        }
    } 
    
};        
//This completes the update (if necessary) of the following data structures: triangles2Nodes, edges2Nodes, edges2Triangles.