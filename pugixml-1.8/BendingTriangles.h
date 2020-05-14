#ifndef BENDINGTRIANGLES_H_
#define BENDINGTRIANGLES_H_

#include "SystemStructures.h"

void ComputeCosTriangleSprings(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    BendingTriangleInfoVecs& bendingTriangleInfoVecs);
/*
for i = 1:size(edge2elem,1)
        %edge2elem is a matrix data structure that records the indices of the two
        %triangles sharing edge "i".
        if (edge2elem(i,1) == edge2elem(i,2))
            continue;
            %The data structure would show two identitical triangles if edge
            %"i" is a boundary edge. In this case, there is no bending
            %force possible. "Continue" skips the computation for this
            %particular edge.
        end
        CANDIDATE = elem(edge2elem(i,1),:);
        CANDIDATE2 = elem(edge2elem(i,2),:);
        HEAD = setdiff(CANDIDATE, edgeinfo(i,:));
        TAIL = setdiff(CANDIDATE2, edgeinfo(i,:));
        %This identifies that vertices in the two-triangle system that are
        %not on the edge shared.
        rj = node(HEAD,:);
        rl = node(TAIL,:);
        rk = node(edgeinfo(i,1),:);
        ri = node(edgeinfo(i,2),:);*/
/*
Given a pair of triangle, one can identify the vertices in the following
fashion:
       rl
      /  \
     /    \
   ri ---- rk
     \    /
      \  /
       rj


rj, rk, ri, rl must first be identified in order to carry out the
computation. rjk implies (rk - rj).
Unit normals are N1 = cross(rjk, rji) and N2 = cross(rli, rlk).
The cosine bending energy is in the form of E = k(1-cos(theta - theta0)).
The cosine bending energy can be rewritten using dot product and cross
product into the form:
E = k(1 - dot(N1,N2)/(norm(N1)*norm(N2))*cos(theta0) 
     - dot(cross(N1,N2),(rki/norm(rki)))/(norm(N1)*norm(N2))*sin(theta0)).
Let COSE represents the term associated with cos(theta0) and SINE
represents the term associated with sin(theta0).

Now rjx represents the derivative with respect to the x-component of
vector rj.
*/



//we want to take 
//output forces on six nodes
struct CosBendingFunctor {
	double spring_constant;
	double angle_0;
	double* locXAddr;
    double* locYAddr;
    double* locZAddr;

    unsigned* idKey;
    double* forceXAddr;
    double* forceYAddr;
    double* forceZAddr;

	unsigned* triangle2Nodes_1Addr;
	unsigned* triangle2Nodes_2Addr;
	unsigned* triangle2Nodes_3Addr;

	__host__ __device__ CosBendingFunctor(
		double& _spring_constant,
		double& _angle_0,
		double* _locXAddr,
		double* _locYAddr,
		double* _locZAddr,

		unsigned* _idKey,
		double* _forceXAddr,
		double* _forceYAddr,
		double* _forceZAddr,

		unsigned* _triangle2Nodes_1Addr,
		unsigned* _triangle2Nodes_2Addr,
		unsigned* _triangle2Nodes_3Addr) :
		


		spring_constant(_spring_constant),
		angle_0(_angle_0),
		locXAddr(_locXAddr),
		locYAddr(_locYAddr),
		locZAddr(_locZAddr),
		
		idKey(_idKey),
		forceXAddr(_forceXAddr),
		forceYAddr(_forceYAddr),
		forceZAddr(_forceZAddr),

		triangle2Nodes_1Addr(_triangle2Nodes_1Addr),
		triangle2Nodes_2Addr(_triangle2Nodes_2Addr),
		triangle2Nodes_3Addr(_triangle2Nodes_3Addr) {}

	//hand in counting iterator and id's of two triangles, and two nodes involved. 
	__device__ double operator()(const Tuuuuu &u5) {
		double energy = 0.0;
		unsigned counter = thrust::get<0>(u5);
		unsigned place = 4 * counter;//represents location in write to vector.

		unsigned id_l, id_j;
		//Id's of elements on sides of a given edge
		//these are the same if we are at the edge of the membrane, so do not compute those.
		unsigned T1 = thrust::get<1>(u5);
		unsigned T2 = thrust::get<2>(u5);

		//these id's are accurate
		unsigned id_k = thrust::get<3>(u5);
		unsigned id_i = thrust::get<4>(u5);
		if (T1 != T2) {
			//we need to compute rl and rj from the two involved triangles. 

			//j is from the element on edge2elem_1
			//l is from the element node in edge2elem_1, not i or k
			unsigned n1T1 = triangle2Nodes_1Addr[T1];
			unsigned n2T1 = triangle2Nodes_2Addr[T1];
			unsigned n3T1 = triangle2Nodes_3Addr[T1];
			if ((n1T1 != id_i) && (n1T1 != id_k)) {
				id_j = n1T1;
			}
			else if ((n2T1 != id_i) && (n2T1 != id_k)) {
				id_j = n2T1;
			}
			else if ((n3T1 != id_i) && (n3T1 != id_k)) {
				id_j = n3T1;
			}



			//l is from the element node in edge2elem_2, not i or k
			//one of these is l, find it
			unsigned n1T2 = triangle2Nodes_1Addr[T2];
			unsigned n2T2 = triangle2Nodes_2Addr[T2];
			unsigned n3T2 = triangle2Nodes_3Addr[T2];
			if ((n1T2 != id_i) && (n1T2 != id_k)) {
				id_l = n1T2;
			}
			else if ((n2T2 != id_i) && (n2T2 != id_k)) {
				id_l = n2T2;
			}
			else if ((n3T2 != id_i) && (n3T2 != id_k)) {
				id_l = n3T2;
			}

			CVec3 ri = thrust::make_tuple<double>(locXAddr[id_i], locYAddr[id_i], locZAddr[id_i]);
			CVec3 rj = thrust::make_tuple<double>(locXAddr[id_j], locYAddr[id_j], locZAddr[id_j]);
			CVec3 rk = thrust::make_tuple<double>(locXAddr[id_k], locYAddr[id_k], locZAddr[id_k]);
			CVec3 rl = thrust::make_tuple<double>(locXAddr[id_l], locYAddr[id_l], locZAddr[id_l]);

			CVec3 rjk = CVec3_minus(rk, rj);
			CVec3 rji = CVec3_minus(ri, rj);
			CVec3 rli = CVec3_minus(ri, rl);		
			CVec3 rlk = CVec3_minus(rk, rl);		
			CVec3 rki = CVec3_minus(ri, rk);

			double nrki = sqrt(CVec3_dot(rki, rki));
			//nrki = sqrt(sum(rki.^2));

			CVec3 unitDir =  thrust::make_tuple<double>(thrust::get<0>(rki)/nrki,
														thrust::get<1>(rki)/nrki,
														thrust::get<2>(rki)/nrki);

			//UD is the unit direction we use to check if the cross product is pointing
			//in the right direction.
			double inv_nrki_sq = 1.0/ (CVec3_dot(rki, rki)); //CHANGE(9/13): removing sqrt, this corresponds to norm^2
			CVec3 zero_vec = thrust::make_tuple<double>(0.0,0.0,0.0);
			CVec3 unitX = thrust::make_tuple<double>(1.0,0.0,0.0);
			CVec3 unitY = thrust::make_tuple<double>(0.0,1.0,0.0);
			CVec3 unitZ = thrust::make_tuple<double>(0.0,0.0,1.0);
			
			Mat_3x3 dUD_rj = thrust::make_tuple<CVec3>(
				CVec3_scalermult( inv_nrki_sq,
					CVec3_plus( 
						CVec3_scalermult( nrki, zero_vec) , 
						CVec3_scalermult( ((0.0 + 0.0 + 0.0) * (-1.0/nrki)), rki) ) ),
				CVec3_scalermult( inv_nrki_sq,
					CVec3_plus( 
						CVec3_scalermult( nrki, zero_vec) , 
						CVec3_scalermult( ((0.0 + 0.0 + 0.0) * (-1.0/nrki)), rki)) ),
				CVec3_scalermult( inv_nrki_sq,
					CVec3_plus( 
						CVec3_scalermult( nrki, zero_vec) , 
						CVec3_scalermult( ((0.0 + 0.0 + 0.0) * (-1.0/nrki)), rki)) ) );

				//CVec3_minus( CVec3_scalermult( inv_nrki_sq, zero_vec), CVec3_scalermult( (0 + 0 + 0) * (1/nrki) ,rki )),
				//CVec3_minus( CVec3_scalermult( inv_nrki_sq, zero_vec), CVec3_scalermult( (0 + 0 + 0) * (1/nrki) ,rki )),
				//CVec3_minus( CVec3_scalermult( inv_nrki_sq, zero_vec), CVec3_scalermult( (0 + 0 + 0) * (1/nrki) ,rki ))); 
				//CHANGE(9/13): rewriting the computation to match the original matlab version

				//(1/norm(rki)^2)*[
				//	nrki*[0,0,0] - rki*(1/nrki)*(0+0+0);...
				//  nrki*[0,0,0] - rki*(1/nrki)*(0+0+0);...
				//  nrki*[0,0,0] - rki*(1/nrki)*(0+0+0)];

			//While dUD_rj can be listed as zero vectors directly since it has no
			//dependence on rj, it is written out fully for double-checking.
			Mat_3x3 dUD_rk = thrust::make_tuple<CVec3>(
				CVec3_scalermult( inv_nrki_sq,
					CVec3_plus( 
						CVec3_scalermult( -nrki, unitX) , 
						CVec3_scalermult(-1.0*(thrust::get<0>(rki) * (-1.0) + 0.0 + 0.0) * (1.0/nrki), rki) )),
				CVec3_scalermult( inv_nrki_sq,
					CVec3_plus( 
						CVec3_scalermult( -nrki, unitY) , 
						CVec3_scalermult(-1.0*(0.0 + thrust::get<1>(rki) * (-1.0) + 0.0) * (1.0/nrki), rki) )),
				CVec3_scalermult( inv_nrki_sq,
					CVec3_plus( 
						CVec3_scalermult( -nrki, unitZ) , 
						CVec3_scalermult(-1.0*(0.0 + 0.0 + thrust::get<2>(rki) * (-1.0)) * (1.0/nrki), rki) )) );

				/*CVec3_minus( CVec3_scalermult(-inv_nrki_sq, CVec3_scalermult(nrki, unitX)), 
					CVec3_scalermult(inv_nrki_sq, CVec3_scalermult((thrust::get<0>(rki)*(-1.0) + 0.0 + 0.0) * (1.0/nrki), rki))),
				CVec3_minus( CVec3_scalermult(-inv_nrki_sq, CVec3_scalermult(nrki, unitY)), 
					CVec3_scalermult(inv_nrki_sq, CVec3_scalermult((0.0 + thrust::get<1>(rki)*(-1.0) + 0.0) * (1.0/nrki), rki))),
				CVec3_minus( CVec3_scalermult(-inv_nrki_sq, CVec3_scalermult(nrki, unitZ)), 
					CVec3_scalermult(inv_nrki_sq, CVec3_scalermult((0.0 + 0.0 + thrust::get<2>(rki)*(-1.0)) * (1.0/nrki), rki))));*/
					//CHANGE(9/13): rewriting the computation to match the original matlab version

				//(1/nrki^2)*[
				//	nrki*[-1,0,0] - rki*(1/nrki)*(rki(1)*-1+0+0);...
				//  nrki*[0,-1,0] - rki*(1/nrki)*(0+rki(2)*-1+0);...
				//  nrki*[0,0,-1] - rki*(1/nrki)*(0+0+rki(3)*-1)];
				
			
			
			Mat_3x3 dUD_ri = thrust::make_tuple<CVec3>(
				CVec3_scalermult(inv_nrki_sq, 
					CVec3_plus(
						CVec3_scalermult(nrki, unitX),  
						CVec3_scalermult(-1.0*(thrust::get<0>(rki)*(1.0) + 0.0 + 0.0) * (1.0/nrki), rki) )),
				CVec3_scalermult(inv_nrki_sq, 
					CVec3_plus(
						CVec3_scalermult(nrki, unitY),  
						CVec3_scalermult(-1.0*(0.0 + thrust::get<1>(rki) * (1.0) + 0.0) * (1.0/nrki), rki) )),
				CVec3_scalermult(inv_nrki_sq, 
					CVec3_plus(
						CVec3_scalermult(nrki, unitZ),  
						CVec3_scalermult(-1.0*(0.0 + 0.0 + thrust::get<2>(rki) * (1.0)) * (1.0/nrki), rki) )) );
			
			
				//CVec3_minus( CVec3_scalermult( inv_nrki_sq, unitX), CVec3_scalermult( (thrust::get<0>(rki)*(1) + 0 + 0) * (1/nrki) ,rki )),
				//CVec3_minus( CVec3_scalermult( inv_nrki_sq, unitY), CVec3_scalermult( (0 + thrust::get<1>(rki)*(1) + 0) * (1/nrki) ,rki )),
				//CVec3_minus( CVec3_scalermult( inv_nrki_sq, unitZ), CVec3_scalermult( (0 + 0 + thrust::get<2>(rki)*(1)) * (1/nrki) ,rki ))); 
				//CHANGE(9/13): rewriting the computation to match the original matlab version

				//(1/nrki^2)*[
				//	nrki*[1,0,0] - rki*(1/nrki)*(rki(1)*1+0+0);...
				//	nrki*[0,1,0] - rki*(1/nrki)*(0+rki(2)*1+0);...
				//	nrki*[0,0,1] - rki*(1/nrki)*(0+0+rki(3)*1)];

			Mat_3x3 dUD_rl = thrust::make_tuple<CVec3>(
				CVec3_scalermult(inv_nrki_sq, 
					CVec3_plus(
						CVec3_scalermult(nrki, zero_vec),  
						CVec3_scalermult(-1.0*(0.0 + 0.0 + 0.0) * (1.0/nrki), rki) )),
				CVec3_scalermult(inv_nrki_sq, 
					CVec3_plus(
						CVec3_scalermult(nrki, zero_vec),  
						CVec3_scalermult(-1.0*(0.0 + 0.0 + 0.0) * (1.0/nrki), rki) )),
				CVec3_scalermult(inv_nrki_sq, 
					CVec3_plus(
						CVec3_scalermult(nrki, zero_vec),  
						CVec3_scalermult(-1.0*(0.0 + 0.0 + 0.0) * (1.0/nrki), rki) )) );



				/*CVec3_minus( CVec3_scalermult(inv_nrki_sq, CVec3_scalermult(nrki, zero_vec)), 
					CVec3_scalermult(inv_nrki_sq, CVec3_scalermult((0.0 + 0.0 + 0.0) * (1.0/nrki), rki)) ),
				CVec3_minus( CVec3_scalermult(inv_nrki_sq, CVec3_scalermult(nrki, zero_vec)), 
					CVec3_scalermult(inv_nrki_sq, CVec3_scalermult((0.0 + 0.0 + 0.0) * (1.0/nrki), rki)) ),
				CVec3_minus( CVec3_scalermult(inv_nrki_sq, CVec3_scalermult(nrki, zero_vec)), 
					CVec3_scalermult(inv_nrki_sq, CVec3_scalermult((0.0 + 0.0 + 0.0) * (1.0/nrki), rki)) )
				);*/
				
				//CHANGE(9/13): rewriting the computation to match the original matlab version

				//(1/nrki^2)*[
				//	nrki*[0,0,0] - rki*(1/nrki)*(0+0+0);...
				//	nrki*[0,0,0] - rki*(1/nrki)*(0+0+0);...
				//	nrki*[0,0,0] - rki*(1/nrki)*(0+0+0)];
		 
			
			CVec3 N1 = CVec3_cross(rjk,rji);
			CVec3 N2 = CVec3_cross(rli, rlk);
			double nN1 = sqrt(CVec3_dot(N1,N1));
			double nN2 = sqrt(CVec3_dot(N2,N2));

			//N1 = cross(rjk, rji);
			//N2 = cross(rli, rlk);
			//nN1 = sqrt(sum(N1.^2)); %norm of N1
			//nN2 = sqrt(sum(N2.^2)); %norm of N2

			double A1 = thrust::get<0>(N1);
			double B1 = thrust::get<1>(N1);
			double C1 = thrust::get<2>(N1);
			double A2 = thrust::get<0>(N2);
			double B2 = thrust::get<1>(N2);
			double C2 = thrust::get<2>(N2);

			
			

			//Derivative of 1st component in N1 with respect to rj
			CVec3 A1_rj = thrust::make_tuple<double>(0.0, -thrust::get<2>(rji) + thrust::get<2>(rjk), -thrust::get<1>(rjk) +thrust::get<1>(rji));
			//A1_rj = [0 , -rji(3)+rjk(3) , -rjk(2)+rji(2)];
			CVec3 A1_rk = thrust::make_tuple<double>(0.0, thrust::get<2>(rji), -thrust::get<1>(rji) );
			//A1_rk = [0 , rji(3) , -rji(2)];
			CVec3 A1_ri = thrust::make_tuple<double>(0.0, -thrust::get<2>(rjk), thrust::get<1>(rjk) );
			//A1_ri = [0 , -rjk(3) , rjk(2)];
			CVec3 A1_rl = thrust::make_tuple<double>(0.0, 0.0, 0.0 );
			//A1_rl = [0 , 0 , 0];

			CVec3 B1_rj = thrust::make_tuple<double>(thrust::get<2>(rji) - thrust::get<2>(rjk), 0.0, thrust::get<0>(rjk) - thrust::get<0>(rji));
			//B1_rj = [rji(3)-rjk(3) , 0 , rjk(1)-rji(1)];
			CVec3 B1_rk = thrust::make_tuple<double>(-thrust::get<2>(rji), 0.0, thrust::get<0>(rji));
			//B1_rk = [-rji(3) , 0 , rji(1)];
			CVec3 B1_ri = thrust::make_tuple<double>(thrust::get<2>(rjk), 0.0, -thrust::get<0>(rjk));
			//B1_ri = [rjk(3) , 0 , -rjk(1)];
			CVec3 B1_rl = thrust::make_tuple<double>(0.0,0.0,0.0);
			//B1_rl = [0 , 0 , 0];

			CVec3 C1_rj = thrust::make_tuple<double>(-thrust::get<1>(rji) + thrust::get<1>(rjk), -thrust::get<0>(rjk) + thrust::get<0>(rji), 0.0);
			//C1_rj = [-rji(2)+rjk(2), -rjk(1)+rji(1) , 0];
			CVec3 C1_rk = thrust::make_tuple<double>(thrust::get<1>(rji), -thrust::get<0>(rji), 0.0);
			//C1_rk = [rji(2), -rji(1) , 0];
            CVec3 C1_ri = thrust::make_tuple<double>(-thrust::get<1>(rjk), thrust::get<0>(rjk), 0.0);
            //C1_ri = [-rjk(2), rjk(1) , 0];
			CVec3 C1_rl = thrust::make_tuple<double>(0.0,0.0,0.0);
			//C1_rl = [0 , 0 , 0];

			CVec3 A2_rj = thrust::make_tuple<double>(0.0,0.0,0.0);
			//A2_rj = [0 , 0 , 0];
			CVec3 A2_rk = thrust::make_tuple<double>( 0.0, -thrust::get<2>(rli), thrust::get<1>(rli) );
			//A2_rk = [0 , -rli(3) , rli(2)];
			CVec3 A2_ri = thrust::make_tuple<double>( 0.0, thrust::get<2>(rlk), -thrust::get<1>(rlk) );
			//A2_ri = [0 , rlk(3) , -rlk(2)];
			CVec3 A2_rl = thrust::make_tuple<double>( 0.0, -thrust::get<2>(rlk) + thrust::get<2>(rli), -thrust::get<1>(rli) + thrust::get<1>(rlk) );
			//A2_rl = [0 , -rlk(3)+rli(3) , -rli(2)+rlk(2)];

			CVec3 B2_rj = thrust::make_tuple<double>(0.0,0.0,0.0);
			//B2_rj = [0 , 0 , 0];
			CVec3 B2_rk = thrust::make_tuple<double>(thrust::get<2>(rli), 0.0, -thrust::get<0>(rli));
			//B2_rk = [rli(3) , 0 , -rli(1)];
			CVec3 B2_ri = thrust::make_tuple<double>(-thrust::get<2>(rlk), 0.0, thrust::get<0>(rlk));
			//B2_ri = [-rlk(3) , 0 , rlk(1)];
			CVec3 B2_rl = thrust::make_tuple<double>(thrust::get<2>(rlk) - thrust::get<2>(rli), 0.0, thrust::get<0>(rli) - thrust::get<0>(rlk));
			//B2_rl = [rlk(3)-rli(3) , 0 , rli(1)-rlk(1)];

			CVec3 C2_rj = thrust::make_tuple<double>(0.0,0.0,0.0);
			//C2_rj = [0 , 0 , 0];
			CVec3 C2_rk = thrust::make_tuple<double>(-thrust::get<1>(rli), thrust::get<0>(rli), 0.0);
			//C2_rk = [-rli(2) , rli(1) , 0];
			CVec3 C2_ri = thrust::make_tuple<double>(thrust::get<1>(rlk), -thrust::get<0>(rlk), 0.0);
			//C2_ri = [rlk(2) , -rlk(1) , 0];
			CVec3 C2_rl = thrust::make_tuple<double>(-thrust::get<1>(rlk) + thrust::get<1>(rli), -thrust::get<0>(rli) + thrust::get<0>(rlk), 0.0);

			//Derivative of the dot product of normal vectors

			CVec3 DN1N2_rj = 
				CVec3_plus(
					CVec3_scalermult( A1 , A2_rj ),
					CVec3_scalermult( A2 , A1_rj ),
					CVec3_scalermult( B1 , B2_rj ),
					CVec3_scalermult( B2 , B1_rj ),
					CVec3_scalermult( C1 , C2_rj ),
					CVec3_scalermult( C2 , C1_rj ) );
				//A1*A2_rj + A2*A1_rj + B1*B2_rj + B2*B1_rj + C1*C2_rj + C2*C1_rj;

			//DN1N2 := dot product dot(N1,N2), and "_rj" represents the partial 
			//derivative with respect to 1st, 2nd, and 3rd component of rj;
			//i.e. rj(1), rj(2), rj(3) being the x,y,z component of rj.

			CVec3	DN1N2_rk = 
				CVec3_plus(
					CVec3_scalermult( A1 , A2_rk ),
					CVec3_scalermult( A2 , A1_rk ),
					CVec3_scalermult( B1 , B2_rk ),
					CVec3_scalermult( B2 , B1_rk ),
					CVec3_scalermult( C1 , C2_rk ),
					CVec3_scalermult( C2 , C1_rk ) );
				//A1*A2_rk + A2*A1_rk + B1*B2_rk + B2*B1_rk + C1*C2_rk + C2*C1_rk;

			CVec3 DN1N2_ri = 
				CVec3_plus(
					CVec3_scalermult( A1 , A2_ri ),
					CVec3_scalermult( A2 , A1_ri ),
					CVec3_scalermult( B1 , B2_ri ),
					CVec3_scalermult( B2 , B1_ri ),
					CVec3_scalermult( C1 , C2_ri ),
					CVec3_scalermult( C2 , C1_ri ) );
				//A1*A2_ri + A2*A1_ri + B1*B2_ri + B2*B1_ri + C1*C2_ri + C2*C1_ri;

			CVec3 DN1N2_rl =
				CVec3_plus(
					CVec3_scalermult( A1 , A2_rl ),
					CVec3_scalermult( A2 , A1_rl ),
					CVec3_scalermult( B1 , B2_rl ),
					CVec3_scalermult( B2 , B1_rl ),
					CVec3_scalermult( C1 , C2_rl ),
					CVec3_scalermult( C2 , C1_rl ) );
				//A1*A2_rl + A2*A1_rl + B1*B2_rl + B2*B1_rl + C1*C2_rl + C2*C1_rl;


//WRONG HERE!!
			// Derivative of the product of norms of normal vectors
			CVec3 PnN1nN2_rj = CVec3_plus(
				CVec3_scalermult(nN1/nN2, 
					CVec3_plus(
						CVec3_scalermult( A2 ,A2_rj ), 
						CVec3_scalermult( B2, B2_rj), 
						CVec3_scalermult( C2, C2_rj) ) ),
				CVec3_scalermult(nN2/nN1, 
					CVec3_plus(
						CVec3_scalermult( A1 ,A1_rj ), 
						CVec3_scalermult( B1, B1_rj), 
						CVec3_scalermult( C1, C1_rj) ) ) );
				//PnN1nN2_rj = [
				//nN1*(1/nN2)*(A2*A2_rj(1)+B2*B2_rj(1)+C2*C2_rj(1))+nN2*(1/nN1)*(A1*A1_rj(1)+B1*B1_rj(1)+C1*C1_rj(1));...
    			//nN1*(1/nN2)*(A2*A2_rj(2)+B2*B2_rj(2)+C2*C2_rj(2))+nN2*(1/nN1)*(A1*A1_rj(2)+B1*B1_rj(2)+C1*C1_rj(2));...
    			//nN1*(1/nN2)*(A2*A2_rj(3)+B2*B2_rj(3)+C2*C2_rj(3))+nN2*(1/nN1)*(A1*A1_rj(3)+B1*B1_rj(3)+C1*C1_rj(3))]; 

			CVec3 PnN1nN2_rk = CVec3_plus(
				CVec3_scalermult(nN1/nN2, 
					CVec3_plus(
						CVec3_scalermult( A2, A2_rk ), 
						CVec3_scalermult( B2, B2_rk),
						CVec3_scalermult( C2, C2_rk) ) ),
				CVec3_scalermult(nN2/nN1, 
					CVec3_plus(
						CVec3_scalermult( A1, A1_rk ), 
						CVec3_scalermult( B1, B1_rk), 
						CVec3_scalermult( C1, C1_rk) ) ) );
				//PnN1nN2_rk = [
				//nN1*(1/nN2)*(A2*A2_rk(1)+B2*B2_rk(1)+C2*C2_rk(1))+nN2*(1/nN1)*(A1*A1_rk(1)+B1*B1_rk(1)+C1*C1_rk(1));...
    			//nN1*(1/nN2)*(A2*A2_rk(2)+B2*B2_rk(2)+C2*C2_rk(2))+nN2*(1/nN1)*(A1*A1_rk(2)+B1*B1_rk(2)+C1*C1_rk(2));...
    			//nN1*(1/nN2)*(A2*A2_rk(3)+B2*B2_rk(3)+C2*C2_rk(3))+nN2*(1/nN1)*(A1*A1_rk(3)+B1*B1_rk(3)+C1*C1_rk(3))]; 


			CVec3 PnN1nN2_ri = CVec3_plus(
				CVec3_scalermult(nN1/nN2, 
					CVec3_plus(
						CVec3_scalermult( A2 ,A2_ri ), 
						CVec3_scalermult( B2, B2_ri), 
						CVec3_scalermult( C2, C2_ri) ) ),
				CVec3_scalermult(nN2/nN1, 
					CVec3_plus(
						CVec3_scalermult( A1, A1_ri ), 
						CVec3_scalermult( B1, B1_ri), 
						CVec3_scalermult( C1, C1_ri) ) ) );
				//PnN1nN2_ri=[
				//nN1*(1/nN2)*(A2*A2_ri(1)+B2*B2_ri(1)+C2*C2_ri(1))+nN2*(1/nN1)*(A1*A1_ri(1)+B1*B1_ri(1)+C1*C1_ri(1));...
    			//nN1*(1/nN2)*(A2*A2_ri(2)+B2*B2_ri(2)+C2*C2_ri(2))+nN2*(1/nN1)*(A1*A1_ri(2)+B1*B1_ri(2)+C1*C1_ri(2));...
    			//nN1*(1/nN2)*(A2*A2_ri(3)+B2*B2_ri(3)+C2*C2_ri(3))+nN2*(1/nN1)*(A1*A1_ri(3)+B1*B1_ri(3)+C1*C1_ri(3))]; 

			CVec3 PnN1nN2_rl = CVec3_plus(
				CVec3_scalermult(nN1/nN2, 
					CVec3_plus(
						CVec3_scalermult( A2, A2_rl ), 
						CVec3_scalermult( B2, B2_rl), 
						CVec3_scalermult( C2, C2_rl) ) ),
				CVec3_scalermult(nN2/nN1, 
					CVec3_plus(
						CVec3_scalermult( A1, A1_rl ), 
						CVec3_scalermult( B1, B1_rl), 
						CVec3_scalermult( C1, C1_rl) ) ) );
				
				//PnN1nN2_rl = [
				//nN1*(1/nN2)*(A2*A2_rl(1)+B2*B2_rl(1)+C2*C2_rl(1))+nN2*(1/nN1)*(A1*A1_rl(1)+B1*B1_rl(1)+C1*C1_rl(1));...
   				//nN1*(1/nN2)*(A2*A2_rl(2)+B2*B2_rl(2)+C2*C2_rl(2))+nN2*(1/nN1)*(A1*A1_rl(2)+B1*B1_rl(2)+C1*C1_rl(2));...
   				//nN1*(1/nN2)*(A2*A2_rl(3)+B2*B2_rl(3)+C2*C2_rl(3))+nN2*(1/nN1)*(A1*A1_rl(3)+B1*B1_rl(3)+C1*C1_rl(3))];

			///////////////////////////////////////////////////////////////////////
			// Derivative of the cross product of normal vectors 
			//This cross product can be written as {B1C2 - B2C1, -A1C2 + A2C1, A1B2 -
			//A2B1}. dcN1N2_rj is read as "d"erivative of the "c"ross product of "N1"
			//and "N2" with respect to "rj".
			////////////////////////////////////////////////////////////////////////
		CVec3 dcN1N2_rj_1 = thrust::make_tuple<double>(
			(B1 * thrust::get<0>(C2_rj) + C2 * thrust::get<0>(B1_rj)) -(B2 * thrust::get<0>(C1_rj) + C1 * thrust::get<0>(B2_rj)), 
			//B1*C2_rj(1) + C2*B1_rj(1))-(B2*C1_rj(1) + C1*B2_rj(1)),
			-(A1 * thrust::get<0>(C2_rj) + C2 * thrust::get<0>(A1_rj) ) + (A2 * thrust::get<0>(C1_rj) + C1 * thrust::get<0>(A2_rj)) , //CHANGE(9/19): misplaced parenthesis
			// -(A1*C2_rj(1)+C2*A1_rj(1))+(A2*C1_rj(1)+C1*A2_rj(1)), 
			(A1 * thrust::get<0>(B2_rj) + B2 * thrust::get<0>(A1_rj) - (A1 * thrust::get<0>(B1_rj) + B1 * thrust::get<0>(A1_rj) ) ) );
			//(A1*B2_rj(1)+B2*A1_rj(1))-(A1*B1_rj(1)+B1*A1_rj(1));

		CVec3 dcN1N2_rj_2 = thrust::make_tuple<double>(
			(B2 * thrust::get<1>(C2_rj) + C2 * thrust::get<1>(B1_rj)) - (B2 * thrust::get<1>(C1_rj) + C1 * thrust::get<1>(B2_rj)),
			//(B1*C2_rj(2) + C2*B1_rj(2))-(B2*C1_rj(2) + C1*B2_rj(2)),
			-(A1 * thrust::get<1>(C2_rj)  + C2 * thrust::get<1>(A1_rj))  + ( A2 * thrust::get<1>(C1_rj) + C1 * thrust::get<1>(A2_rj)) , //CHANGE(9/19): misplaced parenthesis
			// -(A1*C2_rj(2)+C2*A1_rj(2))+(A2*C1_rj(2)+C1*A2_rj(2)), 
	 		(A1 * thrust::get<1>(B2_rj) + B2 * thrust::get<1>(A1_rj) - ( A1 * thrust::get<1>(B1_rj) + B1 * thrust::get<1>(A1_rj) ) ) );
			 //(A1*B2_rj(2)+B2*A1_rj(2))-(A1*B1_rj(2)+B1*A1_rj(2));

		CVec3 dcN1N2_rj_3 = thrust::make_tuple<double>(
			(B1 * thrust::get<2>(C2_rj) + C2 * thrust::get<2>(B1_rj)) - (B2 * thrust::get<2>(C1_rj) + C1 * thrust::get<2>(B2_rj)),
			// (B1*C2_rj(3) + C2*B1_rj(3))-(B2*C1_rj(3) + C1*B2_rj(3)), 
			-(A1 * thrust::get<2>(C2_rj)  + C2 * thrust::get<2>(A1_rj))  + ( A2 * thrust::get<2>(C1_rj) + C1 * thrust::get<2>(A2_rj)) , //CHANGE(9/19):misplaced parenthesis 
			//-(A1*C2_rj(3)+C2*A1_rj(3))+(A2*C1_rj(3)+C1*A2_rj(3)),
	 		(A1 * thrust::get<2>(B2_rj) + B2 * thrust::get<2>(A1_rj) - ( A1 * thrust::get<2>(B1_rj) + B1 * thrust::get<2>(A1_rj) ) ) );
			// (A1*B2_rj(3)+B2*A1_rj(3))-(A1*B1_rj(3)+B1*A1_rj(3))];

		Mat_3x3	dcN1N2_rj = thrust::make_tuple<CVec3>(dcN1N2_rj_1, dcN1N2_rj_2, dcN1N2_rj_3);



		CVec3 dcN1N2_rk_1 = thrust::make_tuple<double>(
			(B1 * thrust::get<0>(C2_rk) + C2 * thrust::get<0>(B1_rk)) - (B2 * thrust::get<0>(C1_rk) + C1 * thrust::get<0>(B2_rk)),
			//(B1*C2_rk(1) + C2*B1_rk(1))-(B2*C1_rk(1) + C1*B2_rk(1)), 
			-(A1 * thrust::get<0>(C2_rk) + C2 * thrust::get<0>(A1_rk)) + (A2 * thrust::get<0>(C1_rk) + C1 * thrust::get<0>(A2_rk)) , //CHANGE(9/19): misplaced parenthesis
			//-(A1*C2_rk(1)+C2*A1_rk(1))+(A2*C1_rk(1)+C1*A2_rk(1)), 
			(A1 * thrust::get<0>(B2_rk) + B2 * thrust::get<0>(A1_rk)) - (A1 * thrust::get<0>(B1_rk) + B1 * thrust::get<0>(A1_rk)));
			//(A1*B2_rk(1)+B2*A1_rk(1))-(A1*B1_rk(1)+B1*A1_rk(1));...
		
		CVec3 dcN1N2_rk_2 = thrust::make_tuple<double>(
			(B1 * thrust::get<1>(C2_rk) + C2 * thrust::get<1>(B1_rk)) - (B2 * thrust::get<1>(C1_rk) + C1 * thrust::get<1>(B2_rk)),
			//  (B1*C2_rk(2) + C2*B1_rk(2))-(B2*C1_rk(2) + C1*B2_rk(2)),
			-(A1 * thrust::get<1>(C2_rk) + C2 * thrust::get<1>(A1_rk)) + (A2 * thrust::get<1>(C1_rk) + C1 * thrust::get<1>(A2_rk)) , //CHANGE(9/19): misplaced parenthesis
			//-(A1*C2_rk(2)+C2*A1_rk(2))+(A2*C1_rk(2)+C1*A2_rk(2)),
			(A1 * thrust::get<1>(B2_rk) + B2 * thrust::get<1>(A1_rk)) - (A1 * thrust::get<1>(B1_rk) + B1 * thrust::get<1>(A1_rk)));	
			// (A1*B2_rk(2)+B2*A1_rk(2))-(A1*B1_rk(2)+B1*A1_rk(2));...
		
		CVec3 dcN1N2_rk_3 = thrust::make_tuple<double>(
			(B1 * thrust::get<2>(C2_rk) + C2 * thrust::get<2>(B1_rk)) - (B2 * thrust::get<2>(C1_rk) + C1 * thrust::get<2>(B2_rk)),
			//  (B1*C2_rk(3) + C2*B1_rk(3))-(B2*C1_rk(3) + C1*B2_rk(3)),
			-(A1 * thrust::get<2>(C2_rk) + C2 * thrust::get<2>(A1_rk)) + (A2 * thrust::get<2>(C1_rk) + C1 * thrust::get<2>(A2_rk)) , //CHANGE(9/19): misplaced parenthesis
			//-(A1*C2_rk(3)+C2*A1_rk(3))+(A2*C1_rk(3)+C1*A2_rk(3)),
			(A1 * thrust::get<2>(B2_rk) + B2 * thrust::get<2>(A1_rk)) - (A1 * thrust::get<2>(B1_rk) + B1 * thrust::get<2>(A1_rk)));	
			// (A1*B2_rk(3)+B2*A1_rk(3))-(A1*B1_rk(3)+B1*A1_rk(3))];

		Mat_3x3	dcN1N2_rk = thrust::make_tuple<CVec3>(dcN1N2_rk_1, dcN1N2_rk_2, dcN1N2_rk_3);



		CVec3 dcN1N2_ri_1 = thrust::make_tuple<double>(
			(B1 * thrust::get<0>(C2_ri) + C2 * thrust::get<0>(B1_ri)) - (B2 * thrust::get<0>(C1_ri) + C1 * thrust::get<0>(B2_ri)),	
			//[(B1*C2_ri(1) + C2*B1_ri(1))-(B2*C1_ri(1) + C1*B2_ri(1)),
			-(A1 * thrust::get<0>(C2_ri) + C2 * thrust::get<0>(A1_ri)) + (A2 * thrust::get<0>(C1_ri) + C1 * thrust::get<0>(A2_ri)) , //CHANGE(9/19): misplaced parenthesis
			//-(A1*C2_ri(1)+C2*A1_ri(1))+(A2*C1_ri(1)+C1*A2_ri(1)),
			(A1 * thrust::get<0>(B2_ri) + B2 * thrust::get<0>(A1_ri)) - (A1 * thrust::get<0>(B1_ri) + B1 * thrust::get<0>(A1_ri)));
			//(A1*B2_ri(1)+B2*A1_ri(1))-(A1*B1_ri(1)+B1*A1_ri(1));...
			
		CVec3 dcN1N2_ri_2 = thrust::make_tuple<double>(
			(B1 * thrust::get<1>(C2_ri) + C2 * thrust::get<1>(B1_ri)) - (B2 * thrust::get<1>(C1_ri) + C1 * thrust::get<1>(B2_ri)),	
			//(B1*C2_ri(2) + C2*B1_ri(2))-(B2*C1_ri(2) + C1*B2_ri(2)),
			-(A1 * thrust::get<1>(C2_ri) + C2 * thrust::get<1>(A1_ri)) + (A2 * thrust::get<1>(C1_ri) + C1 * thrust::get<1>(A2_ri)), //CHANGE(9/19): misplaced parenthesis
			// -(A1*C2_ri(2)+C2*A1_ri(2))+(A2*C1_ri(2)+C1*A2_ri(2)),
			(A1 * thrust::get<1>(B2_ri) + B2 * thrust::get<1>(A1_ri)) - (A1 * thrust::get<1>(B1_ri) + B1 * thrust::get<1>(A1_ri)));
			//  (A1*B2_ri(2)+B2*A1_ri(2))-(A1*B1_ri(2)+B1*A1_ri(2));...
			
		CVec3 dcN1N2_ri_3 = thrust::make_tuple<double>(
			(B1 * thrust::get<2>(C2_ri) + C2 * thrust::get<2>(B1_ri)) - (B2 * thrust::get<2>(C1_ri) + C1 * thrust::get<2>(B2_ri)),	
			//(B1*C2_ri(3) + C2*B1_ri(3))-(B2*C1_ri(3) + C1*B2_ri(3)),
			-(A1 * thrust::get<2>(C2_ri) + C2 * thrust::get<2>(A1_ri)) + (A2 * thrust::get<2>(C1_ri) + C1 * thrust::get<2>(A2_ri)) , //CHANGE(9/19): misplaced parenthesis
			// -(A1*C2_ri(3)+C2*A1_ri(3))+(A2*C1_ri(3)+C1*A2_ri(3)),
			(A1 * thrust::get<2>(B2_ri) + B2 * thrust::get<2>(A1_ri)) - (A1 * thrust::get<2>(B1_ri) + B1 * thrust::get<2>(A1_ri)));
			//(A1*B2_ri(3)+B2*A1_ri(3))-(A1*B1_ri(3)+B1*A1_ri(3))];

		Mat_3x3	dcN1N2_ri = thrust::make_tuple<CVec3>(dcN1N2_ri_1, dcN1N2_ri_2, dcN1N2_ri_3);


		CVec3 dcN1N2_rl_1 = thrust::make_tuple<double>(
			(B1 * thrust::get<0>(C2_rl) + C2 * thrust::get<0>(B1_rl)) - (B2 * thrust::get<0>(C1_rl) + C1 * thrust::get<0>(B2_rl)),
			//(B1*C2_rl(1) + C2*B1_rl(1))-(B2*C1_rl(1) + C1*B2_rl(1)), 
			-(A1 * thrust::get<0>(C2_rl) + C2 * thrust::get<0>(A1_rl)) + (A2 * thrust::get<0>(C1_rl) + C1 * thrust::get<0>(A2_rl)) , //CHANGE(9/19): misplaced parenthesis
			//-(A1*C2_rl(1)+C2*A1_rl(1))+(A2*C1_rl(1)+C1*A2_rl(1)), 
			(A1 * thrust::get<0>(B2_rl) + B2 * thrust::get<0>(A1_rl)) - (A1 * thrust::get<0>(B1_rl) + B1 * thrust::get<0>(A1_rl)));
			//(A1*B2_rl(1)+B2*A1_rl(1))-(A1*B1_rl(1)+B1*A1_rl(1));...
			
		CVec3 dcN1N2_rl_2 = thrust::make_tuple<double>(
			(B1 * thrust::get<1>(C2_rl) + C2 * thrust::get<1>(B1_rl)) - (B2 * thrust::get<1>(C1_rl) + C1 * thrust::get<1>(B2_rl)),
			//(B1*C2_rl(2) + C2*B1_rl(2))-(B2*C1_rl(2) + C1*B2_rl(2)), 
			-(A1 * thrust::get<1>(C2_rl) + C2 * thrust::get<1>(A1_rl)) + (A2 * thrust::get<1>(C1_rl) + C1 * thrust::get<1>(A2_rl)) , //CHANGE(9/19): misplaced parenthesis
			//-(A1*C2_rl(2)+C2*A1_rl(2))+(A2*C1_rl(2)+C1*A2_rl(2)), 
			(A1 * thrust::get<1>(B2_rl) + B2 * thrust::get<1>(A1_rl)) - (A1 * thrust::get<1>(B1_rl) + B1 * thrust::get<1>(A1_rl)));
			//(A1*B2_rl(2)+B2*A1_rl(2))-(A1*B1_rl(2)+B1*A1_rl(2));...
			
		CVec3 dcN1N2_rl_3 = thrust::make_tuple<double>(
			(B1 * thrust::get<2>(C2_rl) + C2 * thrust::get<2>(B1_rl)) - (B2 * thrust::get<2>(C1_rl) + C1 * thrust::get<2>(B2_rl)),
			//(B1*C2_rl(3) + C2*B1_rl(3))-(B2*C1_rl(3) + C1*B2_rl(3)), 
			-(A1 * thrust::get<2>(C2_rl) + C2 * thrust::get<2>(A1_rl)) + (A2 * thrust::get<2>(C1_rl) + C1 * thrust::get<2>(A2_rl)) , //CHANGE(9/19): misplaced parenthesis		
			//-(A1*C2_rl(3)+C2*A1_rl(3))+(A2*C1_rl(3)+C1*A2_rl(3)), 
			(A1 * thrust::get<2>(B2_rl) + B2 * thrust::get<2>(A1_rl)) - (A1 * thrust::get<2>(B1_rl) + B1 * thrust::get<2>(A1_rl)));			
			//(A1*B2_rl(3)+B2*A1_rl(3))-(A1*B1_rl(3)+B1*A1_rl(3))];

		Mat_3x3	dcN1N2_rl = thrust::make_tuple<CVec3>(dcN1N2_rl_1, dcN1N2_rl_2, dcN1N2_rl_3);

				double N1N2_nN1nN2 = CVec3_dot(N1,N2) / (nN1*nN2*nN1*nN2);
				CVec3 COSE_rj = CVec3_plus(CVec3_scalermult( 1/(nN1*nN2), DN1N2_rj), CVec3_scalermult( -1.0*N1N2_nN1nN2, PnN1nN2_rj));
					//CVec3_scalermult(nN1 * nN2, DN1N2_rj) , CVec3_scalermult( N1N2_nN1nN2 , PnN1nN2_rj ) );
					//CHANGE(9/14): rewritting the computation to match the matlab version.
	
					//COSE_rjx = ((nN1*nN2)*DN1N2_rj(1) - dot(N1,N2)*PnN1nN2_rj(1))/(nN1*nN2)^2; 
					//COSE_rjy = ((nN1*nN2)*DN1N2_rj(2) - dot(N1,N2)*PnN1nN2_rj(2))/(nN1*nN2)^2;
					//COSE_rjz = ((nN1*nN2)*DN1N2_rj(3) - dot(N1,N2)*PnN1nN2_rj(3))/(nN1*nN2)^2;
	
				CVec3 COSE_rk = CVec3_plus(CVec3_scalermult( 1/(nN1*nN2), DN1N2_rk), CVec3_scalermult( -1.0*N1N2_nN1nN2, PnN1nN2_rk));
					//CVec3_scalermult(nN1 * nN2, DN1N2_rk) , CVec3_scalermult( N1N2_nN1nN2 , PnN1nN2_rk ) );
					//CHANGE(9/14): rewritting the computation to match the matlab version.
	
					//COSE_rkx = ((nN1*nN2)*DN1N2_rk(1) - dot(N1,N2)*PnN1nN2_rk(1))/(nN1*nN2)^2;
					//COSE_rky = ((nN1*nN2)*DN1N2_rk(2) - dot(N1,N2)*PnN1nN2_rk(2))/(nN1*nN2)^2;
					//COSE_rkz = ((nN1*nN2)*DN1N2_rk(3) - dot(N1,N2)*PnN1nN2_rk(3))/(nN1*nN2)^2;
					
				CVec3 COSE_ri = CVec3_plus(CVec3_scalermult( 1/(nN1*nN2), DN1N2_ri), CVec3_scalermult( -1.0*N1N2_nN1nN2, PnN1nN2_ri));
					//CVec3_scalermult(nN1 * nN2, DN1N2_ri) , CVec3_scalermult( N1N2_nN1nN2 , PnN1nN2_ri ) );
					//CHANGE(9/14): rewritting the computation to match the matlab version.
					
					//COSE_rix = ((nN1*nN2)*DN1N2_ri(1) - dot(N1,N2)*PnN1nN2_ri(1))/(nN1*nN2)^2;
					//COSE_riy = ((nN1*nN2)*DN1N2_ri(2) - dot(N1,N2)*PnN1nN2_ri(2))/(nN1*nN2)^2;
					//COSE_riz = ((nN1*nN2)*DN1N2_ri(3) - dot(N1,N2)*PnN1nN2_ri(3))/(nN1*nN2)^2;
					
				CVec3 COSE_rl = CVec3_plus(CVec3_scalermult( 1/(nN1*nN2), DN1N2_rl), CVec3_scalermult( -1.0*N1N2_nN1nN2, PnN1nN2_rl));
					//CVec3_scalermult(nN1 * nN2, DN1N2_rl) , CVec3_scalermult( N1N2_nN1nN2 , PnN1nN2_rl ) );
					//CHANGE(9/14): rewritting the computation to match the matlab version
					
					//COSE_rlx = ((nN1*nN2)*DN1N2_rl(1) - dot(N1,N2)*PnN1nN2_rl(1))/(nN1*nN2)^2;
					//COSE_rly = ((nN1*nN2)*DN1N2_rl(2) - dot(N1,N2)*PnN1nN2_rl(2))/(nN1*nN2)^2;
					//COSE_rlz = ((nN1*nN2)*DN1N2_rl(3) - dot(N1,N2)*PnN1nN2_rl(3))/(nN1*nN2)^2;

			double SINE_rjx = 1/(nN1*nN2) * (CVec3_dot( CVec3_cross(N1,N2) ,thrust::get<0>(dUD_rj))
								+ CVec3_dot( thrust::get<0>(dcN1N2_rj), unitDir ))
								- (1/(nN1*nN2*nN1*nN2)) * CVec3_dot( CVec3_cross(N1, N2), unitDir) * thrust::get<0>(PnN1nN2_rj);
								//(nN1*nN2*(dot(cross(N1,N2), dUD_rj(1,:)) + dot(dcN1N2_rj(1,:), UD)) - dot(cross(N1,N2),UD)*PnN1nN2_rj(1))/(nN1*nN2)^2;
			
			double SINE_rjy = 1/(nN1*nN2) * (CVec3_dot( CVec3_cross(N1,N2) ,thrust::get<1>(dUD_rj))
								+ CVec3_dot( thrust::get<1>(dcN1N2_rj), unitDir ))
								- (1/(nN1*nN2*nN1*nN2)) * CVec3_dot( CVec3_cross(N1, N2), unitDir) * thrust::get<1>(PnN1nN2_rj);			
								//SINE_rjy = (nN1*nN2*(dot(cross(N1,N2), dUD_rj(2,:)) + dot(dcN1N2_rj(2,:), UD)) - dot(cross(N1,N2),UD)*PnN1nN2_rj(2))/(nN1*nN2)^2;
			
			double SINE_rjz = 1.0/(nN1*nN2) * (CVec3_dot( CVec3_cross(N1,N2) ,thrust::get<2>(dUD_rj))
								+ CVec3_dot( thrust::get<2>(dcN1N2_rj), unitDir ))
								- (1.0/(nN1*nN2*nN1*nN2)) * CVec3_dot( CVec3_cross(N1, N2), unitDir) * thrust::get<2>(PnN1nN2_rj);
								//1/(nN1*nN2)*(dot(cross(N1,N2), dUD_rj(3,:)) + dot(dcN1N2_rj(3,:), UD))...
    							//- dot(cross(N1,N2),UD)*PnN1nN2_rj(3) / (nN1*nN2)^2;
			double SINE_rkx = 1.0/(nN1*nN2) * (CVec3_dot( CVec3_cross(N1,N2) ,thrust::get<0>(dUD_rk))
								+ CVec3_dot( thrust::get<0>(dcN1N2_rk), unitDir ))
								- (1.0/(nN1*nN2*nN1*nN2)) * CVec3_dot( CVec3_cross(N1, N2), unitDir) * thrust::get<0>(PnN1nN2_rk);
								//(nN1*nN2*(dot(cross(N1,N2), dUD_rk(1,:)) + dot(dcN1N2_rk(1,:), UD)) - dot(cross(N1,N2),UD)*PnN1nN2_rk(1))/(nN1*nN2)^2;

			double SINE_rky = 1.0/(nN1*nN2) * (CVec3_dot( CVec3_cross(N1,N2) ,thrust::get<1>(dUD_rk))
								+ CVec3_dot( thrust::get<1>(dcN1N2_rk), unitDir ))
								- (1/(nN1*nN2*nN1*nN2)) * CVec3_dot( CVec3_cross(N1, N2), unitDir) * thrust::get<1>(PnN1nN2_rk);
								//(nN1*nN2*(dot(cross(N1,N2), dUD_rk(2,:)) + dot(dcN1N2_rk(2,:), UD)) - dot(cross(N1,N2),UD)*PnN1nN2_rk(2))/(nN1*nN2)^2;
							
			double SINE_rkz = 1.0/(nN1*nN2) * (CVec3_dot( CVec3_cross(N1,N2) ,thrust::get<2>(dUD_rk))
								+ CVec3_dot( thrust::get<2>(dcN1N2_rk), unitDir ))
								- (1/(nN1*nN2*nN1*nN2)) * CVec3_dot( CVec3_cross(N1, N2), unitDir) * thrust::get<2>(PnN1nN2_rk);
								//(nN1*nN2*(dot(cross(N1,N2), dUD_rk(3,:)) + dot(dcN1N2_rk(3,:), UD)) - dot(cross(N1,N2),UD)*PnN1nN2_rk(3))/(nN1*nN2)^2;

			double SINE_rix = 1.0/(nN1*nN2) * (CVec3_dot( CVec3_cross(N1,N2) ,thrust::get<0>(dUD_ri))
								+ CVec3_dot( thrust::get<0>(dcN1N2_ri), unitDir ))
								- (1/(nN1*nN2*nN1*nN2)) * CVec3_dot( CVec3_cross(N1, N2), unitDir) * thrust::get<0>(PnN1nN2_ri);
								//(nN1*nN2*(dot(cross(N1,N2), dUD_ri(1,:)) + dot(dcN1N2_ri(1,:), UD)) - dot(cross(N1,N2),UD)*PnN1nN2_ri(1))/(nN1*nN2)^2;

			double SINE_riy = 1.0/(nN1*nN2) * (CVec3_dot( CVec3_cross(N1,N2) ,thrust::get<1>(dUD_ri))
								+ CVec3_dot( thrust::get<1>(dcN1N2_ri), unitDir ))
								- (1/(nN1*nN2*nN1*nN2)) * CVec3_dot( CVec3_cross(N1, N2), unitDir) * thrust::get<1>(PnN1nN2_ri);
								//(nN1*nN2*(dot(cross(N1,N2), dUD_ri(2,:)) + dot(dcN1N2_ri(2,:), UD)) - dot(cross(N1,N2),UD)*PnN1nN2_ri(2))/(nN1*nN2)^2;

			double SINE_riz = 1.0/(nN1*nN2) * (CVec3_dot( CVec3_cross(N1,N2) ,thrust::get<2>(dUD_ri))
								+ CVec3_dot( thrust::get<2>(dcN1N2_ri), unitDir ))
								- (1/(nN1*nN2*nN1*nN2)) * CVec3_dot( CVec3_cross(N1, N2), unitDir) * thrust::get<2>(PnN1nN2_ri);
								//(nN1*nN2*(dot(cross(N1,N2), dUD_ri(3,:)) + dot(dcN1N2_ri(3,:), UD)) - dot(cross(N1,N2),UD)*PnN1nN2_ri(3))/(nN1*nN2)^2;

			double SINE_rlx = 1.0/(nN1*nN2) * (CVec3_dot( CVec3_cross(N1,N2) ,thrust::get<0>(dUD_rl))
								+ CVec3_dot( thrust::get<0>(dcN1N2_rl), unitDir ))
								- (1/(nN1*nN2*nN1*nN2)) * CVec3_dot( CVec3_cross(N1, N2), unitDir) * thrust::get<0>(PnN1nN2_rl);
								//(nN1*nN2*(dot(cross(N1,N2), dUD_rl(1,:)) + dot(dcN1N2_rl(1,:), UD)) - dot(cross(N1,N2),UD)*PnN1nN2_rl(1))/(nN1*nN2)^2;

			double SINE_rly = 1.0/(nN1*nN2) * (CVec3_dot( CVec3_cross(N1,N2) ,thrust::get<1>(dUD_rl))
								+ CVec3_dot( thrust::get<1>(dcN1N2_rl), unitDir ))
								- (1/(nN1*nN2*nN1*nN2)) * CVec3_dot( CVec3_cross(N1, N2), unitDir) * thrust::get<1>(PnN1nN2_rl);
								//(nN1*nN2*(dot(cross(N1,N2), dUD_rl(2,:)) + dot(dcN1N2_rl(2,:), UD)) - dot(cross(N1,N2),UD)*PnN1nN2_rl(2))/(nN1*nN2)^2;
			
			double SINE_rlz = 1.0/(nN1*nN2) * (CVec3_dot( CVec3_cross(N1,N2) ,thrust::get<2>(dUD_rl))
								+ CVec3_dot( thrust::get<2>(dcN1N2_rl), unitDir ))
								- (1/(nN1*nN2*nN1*nN2)) * CVec3_dot( CVec3_cross(N1, N2), unitDir) * thrust::get<2>(PnN1nN2_rl);

			//Apply Forces
			//using counting iterator we place in the 4 spots of the x-y-z vectors for i,j,k,l
			//you want to reset the force since we are writing to a temporary vector
			//must also write to keys for reduce step. 
			if ((nN1 != 0.0) && (nN2 !=0.0)) {
				idKey[place] = id_i;
				//issue must be in cos_ri
				forceXAddr[place] = spring_constant * (cos(angle_0) * (thrust::get<0>(COSE_ri)) + spring_constant * sin(angle_0) * SINE_rix);
				forceYAddr[place] = spring_constant * (cos(angle_0) * (thrust::get<1>(COSE_ri)) + spring_constant * sin(angle_0) * SINE_riy);
				forceZAddr[place] = spring_constant * (cos(angle_0) * (thrust::get<2>(COSE_ri)) + spring_constant * sin(angle_0) * SINE_riz);

				idKey[place+1] = id_j;
				forceXAddr[place+1] = spring_constant * (cos(angle_0) * (thrust::get<0>(COSE_rj)) + spring_constant * sin(angle_0) * SINE_rjx);
				forceYAddr[place+1] = spring_constant * (cos(angle_0) * (thrust::get<1>(COSE_rj)) + spring_constant * sin(angle_0) * SINE_rjy);
				forceZAddr[place+1] = spring_constant * (cos(angle_0) * (thrust::get<2>(COSE_rj)) + spring_constant * sin(angle_0) * SINE_rjz);

				idKey[place+2] = id_k;
				forceXAddr[place+2] = spring_constant * (cos(angle_0) * (thrust::get<0>(COSE_rk)) + spring_constant * sin(angle_0) * SINE_rkx);
				forceYAddr[place+2] = spring_constant * (cos(angle_0) * (thrust::get<1>(COSE_rk)) + spring_constant * sin(angle_0) * SINE_rky);
				forceZAddr[place+2] = spring_constant * (cos(angle_0) * (thrust::get<2>(COSE_rk)) + spring_constant * sin(angle_0) * SINE_rkz);

				idKey[place+3] = id_l;
				forceXAddr[place+3] = spring_constant * (cos(angle_0) * (thrust::get<0>(COSE_rl)) + spring_constant * sin(angle_0) * SINE_rlx);
				forceYAddr[place+3] = spring_constant * (cos(angle_0) * (thrust::get<1>(COSE_rl)) + spring_constant * sin(angle_0) * SINE_rly);
				forceZAddr[place+3] = spring_constant * (cos(angle_0) * (thrust::get<2>(COSE_rl)) + spring_constant * sin(angle_0) * SINE_rlz);
			}

			//due to symmetry, we the angle does not matter if it is the convex or concave one. 
			double cosAngle = CVec3_dot(N1, N2) / (nN1*nN2);
			
			if (cosAngle > 1.0) {
				cosAngle = 1.0;
			}
			else if (cosAngle < -1.0){
				cosAngle = -1.0;
			}

			double theta_current = acos( cosAngle );
			
			energy = spring_constant * (1 - cos(0.0 - theta_current) );
			
		}
		return energy;
	}
};





#endif /* BENDINGTRIANGLES_H_*/