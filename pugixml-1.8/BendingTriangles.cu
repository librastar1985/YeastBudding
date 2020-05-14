
#include "System.h"
#include "SystemStructures.h"
#include "BendingTriangles.h"

void ComputeCosTriangleSprings(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    BendingTriangleInfoVecs& bendingTriangleInfoVecs) {
/*
if (generalParams.iteration % 1 == 0) {
	
	std::cout<<" iteration: "<< generalParams.iteration << std::endl;
        unsigned id_l, id_j;
        
        unsigned T1 = coordInfoVecs.edges2Triangles_1[0];
        unsigned T2 = coordInfoVecs.edges2Triangles_2[0];

		//these id's are accurate
		unsigned id_k = coordInfoVecs.edges2Nodes_1[0];
        unsigned id_i = coordInfoVecs.edges2Nodes_2[0];
        
		if (T1 != T2) {
			//we need to compute rl and rj from the two involved triangles. 

			
			//one of these is j, but it cannot be equal to id_i or id_k
			unsigned n1T1 = coordInfoVecs.triangles2Nodes_1[T1];
			unsigned n2T1 = coordInfoVecs.triangles2Nodes_2[T1];
			unsigned n3T1 = coordInfoVecs.triangles2Nodes_3[T1];
			if ((n1T1 != id_i) && (n1T1 != id_k)) {
				id_j = n1T1;
			}
			else if ((n2T1 != id_i) && (n2T1 != id_k)) {
				id_j = n2T1;
			}
			else if ((n3T1 != id_i) && (n3T1 != id_k)) {
				id_j = n3T1;
			}

			//one of these is l, find it
			unsigned n1T2 = coordInfoVecs.triangles2Nodes_1[T2];
			unsigned n2T2 = coordInfoVecs.triangles2Nodes_2[T2];
			unsigned n3T2 = coordInfoVecs.triangles2Nodes_3[T2];
			if ((n1T2 != id_i) && (n1T2 != id_k)) {
				id_l = n1T2;
			}
			else if ((n2T2 != id_i) && (n2T2 != id_k)) {
				id_l = n2T2;
			}
			else if ((n3T2 != id_i) && (n3T2 != id_k)) {
				id_l = n3T2;
			}

			std::cout<<"i: "<< id_i << std::endl;
			std::cout<<"k: "<< id_k << std::endl;
			std::cout<<"l: "<< id_l << std::endl;
			std::cout<<"j: "<< id_j << std::endl;
			CVec3 ri = thrust::make_tuple<double>(coordInfoVecs.nodeLocX[id_i], coordInfoVecs.nodeLocY[id_i], coordInfoVecs.nodeLocZ[id_i]);
			CVec3 rj = thrust::make_tuple<double>(coordInfoVecs.nodeLocX[id_j], coordInfoVecs.nodeLocY[id_j], coordInfoVecs.nodeLocZ[id_j]);
			CVec3 rk = thrust::make_tuple<double>(coordInfoVecs.nodeLocX[id_k], coordInfoVecs.nodeLocY[id_k], coordInfoVecs.nodeLocZ[id_k]);
			CVec3 rl = thrust::make_tuple<double>(coordInfoVecs.nodeLocX[id_l], coordInfoVecs.nodeLocY[id_l], coordInfoVecs.nodeLocZ[id_l]);
			std::cout<<"ri: "<< thrust::get<0>(ri)<<" "<< thrust::get<1>(ri)<<" "<< thrust::get<2>(ri)<<std::endl;
            std::cout<<"rj: "<< thrust::get<0>(rj)<<" "<< thrust::get<1>(rj)<<" "<< thrust::get<2>(rj)<<std::endl;
            std::cout<<"rk: "<< thrust::get<0>(rk)<<" "<< thrust::get<1>(rk)<<" "<< thrust::get<2>(rk)<<std::endl;
            std::cout<<"rl: "<< thrust::get<0>(rl)<<" "<< thrust::get<1>(rl)<<" "<< thrust::get<2>(rl)<<std::endl;
            
			CVec3 rjk = CVec3_plus(rk, CVec3_scalermult(-1.0, rj) );
			CVec3 rji = CVec3_plus(ri, CVec3_scalermult(-1.0, rj) );
			CVec3 rli = CVec3_plus(ri, CVec3_scalermult(-1.0, rl) );		
			CVec3 rlk = CVec3_plus(rk, CVec3_scalermult(-1.0, rl) );		
			CVec3 rki = CVec3_plus(ri, CVec3_scalermult(-1.0, rk) );

            std::cout<<"rjk: "<< thrust::get<0>(rjk)<<" "<< thrust::get<1>(rjk)<<" "<< thrust::get<2>(rjk)<<std::endl;
            std::cout<<"rji: "<< thrust::get<0>(rji)<<" "<< thrust::get<1>(rji)<<" "<< thrust::get<2>(rji)<<std::endl;
            std::cout<<"rli: "<< thrust::get<0>(rli)<<" "<< thrust::get<1>(rli)<<" "<< thrust::get<2>(rli)<<std::endl;
            std::cout<<"rlk: "<< thrust::get<0>(rlk)<<" "<< thrust::get<1>(rlk)<<" "<< thrust::get<2>(rlk)<<std::endl;
            std::cout<<"rki: "<< thrust::get<0>(rki)<<" "<< thrust::get<1>(rki)<<" "<< thrust::get<2>(rki)<<std::endl;
			double nrki = sqrt(CVec3_dot(rki, rki));
			//nrki = sqrt(sum(rki.^2));

			CVec3 unitDir =  thrust::make_tuple<double>(thrust::get<0>(rki)/nrki,
														thrust::get<1>(rki)/nrki,
														thrust::get<2>(rki)/nrki);

			//UD is the unit direction we use to check if the cross product is pointing
			//in the right direction.
			double inv_nrki_sq = 1.0/ (CVec3_dot(rki, rki)); //CHANGE(9/13): removing sqrt, this corresponds to norm^2
			CVec3 zero_vec = thrust::make_tuple<double>(0.0, 0.0, 0.0);
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
			std::cout<< "dUD_rj_1"<< thrust::get<0>(thrust::get<0>(dUD_rj))<< " "<< thrust::get<1>(thrust::get<0>(dUD_rj))<< " "<< thrust::get<2>(thrust::get<0>(dUD_rj)) <<std::endl;
			std::cout<< "dUD_rj_2"<< thrust::get<0>(thrust::get<1>(dUD_rj))<< " "<< thrust::get<1>(thrust::get<1>(dUD_rj))<< " "<< thrust::get<2>(thrust::get<1>(dUD_rj)) <<std::endl;
			std::cout<< "dUD_rj_3"<< thrust::get<0>(thrust::get<2>(dUD_rj))<< " "<< thrust::get<1>(thrust::get<2>(dUD_rj))<< " "<< thrust::get<2>(thrust::get<2>(dUD_rj)) <<std::endl;
						
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
				
					std::cout<< "dUD_rk_1"<< thrust::get<0>(thrust::get<0>(dUD_rk))<< " "<< thrust::get<1>(thrust::get<0>(dUD_rk))<< " "<< thrust::get<2>(thrust::get<0>(dUD_rk)) <<std::endl;
					std::cout<< "dUD_rk_2"<< thrust::get<0>(thrust::get<1>(dUD_rk))<< " "<< thrust::get<1>(thrust::get<1>(dUD_rk))<< " "<< thrust::get<2>(thrust::get<1>(dUD_rk)) <<std::endl;
					std::cout<< "dUD_rk_3"<< thrust::get<0>(thrust::get<2>(dUD_rk))<< " "<< thrust::get<1>(thrust::get<2>(dUD_rk))<< " "<< thrust::get<2>(thrust::get<2>(dUD_rk)) <<std::endl;
					
			
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
			
			std::cout<< "dUD_ri_1"<< thrust::get<0>(thrust::get<0>(dUD_ri))<< " "<< thrust::get<1>(thrust::get<0>(dUD_ri))<< " "<< thrust::get<2>(thrust::get<0>(dUD_ri)) <<std::endl;
			std::cout<< "dUD_ri_2"<< thrust::get<0>(thrust::get<1>(dUD_ri))<< " "<< thrust::get<1>(thrust::get<1>(dUD_ri))<< " "<< thrust::get<2>(thrust::get<1>(dUD_ri)) <<std::endl;
			std::cout<< "dUD_ri_3"<< thrust::get<0>(thrust::get<2>(dUD_ri))<< " "<< thrust::get<1>(thrust::get<2>(dUD_ri))<< " "<< thrust::get<2>(thrust::get<2>(dUD_ri)) <<std::endl;
					
			
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
						std::cout<< "dUD_rl_1"<< thrust::get<0>(thrust::get<0>(dUD_rl))<< " "<< thrust::get<1>(thrust::get<0>(dUD_rl))<< " "<< thrust::get<2>(thrust::get<0>(dUD_rl)) <<std::endl;
						std::cout<< "dUD_rl_2"<< thrust::get<0>(thrust::get<1>(dUD_rl))<< " "<< thrust::get<1>(thrust::get<1>(dUD_rl))<< " "<< thrust::get<2>(thrust::get<1>(dUD_rl)) <<std::endl;
						std::cout<< "dUD_rl_3"<< thrust::get<0>(thrust::get<2>(dUD_rl))<< " "<< thrust::get<1>(thrust::get<2>(dUD_rl))<< " "<< thrust::get<2>(thrust::get<2>(dUD_rl)) <<std::endl;
					

			
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

			
			std::cout<<"N1: "<< thrust::get<0>(N1)<<" "<< thrust::get<1>(N1)<<" "<< thrust::get<2>(N1)<<std::endl;
            std::cout<<"N2: "<< thrust::get<0>(N2)<<" "<< thrust::get<1>(N2)<<" "<< thrust::get<2>(N2)<<std::endl;
            

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
			//C2_rl = [-rlk(2)+rli(2) , -rli(1)+rlk(1) , 0];

            std::cout<<"A1_rk: "<< thrust::get<0>(A1_rk)<<" "<< thrust::get<1>(A1_rk)<<" "<< thrust::get<2>(A1_rk)<<std::endl;
            std::cout<<"A1_rj: "<< thrust::get<0>(A1_rj)<<" "<< thrust::get<1>(A1_rj)<<" "<< thrust::get<2>(A1_rj)<<std::endl;
            std::cout<<"A1_ri: "<< thrust::get<0>(A1_ri)<<" "<< thrust::get<1>(A1_ri)<<" "<< thrust::get<2>(A1_ri)<<std::endl;
            std::cout<<"A1_rl: "<< thrust::get<0>(A1_rl)<<" "<< thrust::get<1>(A1_rl)<<" "<< thrust::get<2>(A1_rl)<<std::endl;
            std::cout<<"B1_rk: "<< thrust::get<0>(B1_rk)<<" "<< thrust::get<1>(B1_rk)<<" "<< thrust::get<2>(B1_rk)<<std::endl;
            std::cout<<"B1_rj: "<< thrust::get<0>(B1_rj)<<" "<< thrust::get<1>(B1_rj)<<" "<< thrust::get<2>(B1_rj)<<std::endl;
            std::cout<<"B1_ri: "<< thrust::get<0>(B1_ri)<<" "<< thrust::get<1>(B1_ri)<<" "<< thrust::get<2>(B1_ri)<<std::endl;
            std::cout<<"B1_rl: "<< thrust::get<0>(B1_rl)<<" "<< thrust::get<1>(B1_rl)<<" "<< thrust::get<2>(B1_rl)<<std::endl;
            std::cout<<"C1_rk: "<< thrust::get<0>(C1_rk)<<" "<< thrust::get<1>(C1_rk)<<" "<< thrust::get<2>(C1_rk)<<std::endl;
            std::cout<<"C1_rj: "<< thrust::get<0>(C1_rj)<<" "<< thrust::get<1>(C1_rj)<<" "<< thrust::get<2>(C1_rj)<<std::endl;
            std::cout<<"C1_ri: "<< thrust::get<0>(C1_ri)<<" "<< thrust::get<1>(C1_ri)<<" "<< thrust::get<2>(C1_ri)<<std::endl;
            std::cout<<"C1_rl: "<< thrust::get<0>(C1_rl)<<" "<< thrust::get<1>(C1_rl)<<" "<< thrust::get<2>(C1_rl)<<std::endl;
            std::cout<<"A2_rk: "<< thrust::get<0>(A2_rk)<<" "<< thrust::get<1>(A2_rk)<<" "<< thrust::get<2>(A2_rk)<<std::endl;
            std::cout<<"A2_rj: "<< thrust::get<0>(A2_rj)<<" "<< thrust::get<1>(A2_rj)<<" "<< thrust::get<2>(A2_rj)<<std::endl;
            std::cout<<"A2_ri: "<< thrust::get<0>(A2_ri)<<" "<< thrust::get<1>(A2_ri)<<" "<< thrust::get<2>(A2_ri)<<std::endl;
            std::cout<<"A2_rl: "<< thrust::get<0>(A2_rl)<<" "<< thrust::get<1>(A2_rl)<<" "<< thrust::get<2>(A2_rl)<<std::endl;
            std::cout<<"B2_rk: "<< thrust::get<0>(B2_rk)<<" "<< thrust::get<1>(B2_rk)<<" "<< thrust::get<2>(B2_rk)<<std::endl;
            std::cout<<"B2_rj: "<< thrust::get<0>(B2_rj)<<" "<< thrust::get<1>(B2_rj)<<" "<< thrust::get<2>(B2_rj)<<std::endl;
            std::cout<<"B2_ri: "<< thrust::get<0>(B2_ri)<<" "<< thrust::get<1>(B2_ri)<<" "<< thrust::get<2>(B2_ri)<<std::endl;
            std::cout<<"B2_rl: "<< thrust::get<0>(B2_rl)<<" "<< thrust::get<1>(B2_rl)<<" "<< thrust::get<2>(B2_rl)<<std::endl;
            std::cout<<"C2_rk: "<< thrust::get<0>(C2_rk)<<" "<< thrust::get<1>(C2_rk)<<" "<< thrust::get<2>(C2_rk)<<std::endl;
            std::cout<<"C2_rj: "<< thrust::get<0>(C2_rj)<<" "<< thrust::get<1>(C2_rj)<<" "<< thrust::get<2>(C2_rj)<<std::endl;
            std::cout<<"C2_ri: "<< thrust::get<0>(C2_ri)<<" "<< thrust::get<1>(C2_ri)<<" "<< thrust::get<2>(C2_ri)<<std::endl;
            std::cout<<"C2_rl: "<< thrust::get<0>(C2_rl)<<" "<< thrust::get<1>(C2_rl)<<" "<< thrust::get<2>(C2_rl)<<std::endl;
             

			//Derivative of the dot product of normal vectors

			
			CVec3 At =	CVec3_plus(
					CVec3_scalermult( A1 , A2_rj ),
					CVec3_scalermult( A2 , A1_rj ));
			CVec3 Bt =	CVec3_plus(
					CVec3_scalermult( B1 , B2_rj ),
					CVec3_scalermult( B2 , B1_rj ));
			CVec3 Ct =	CVec3_plus(
					CVec3_scalermult( C1 , C2_rj ),
					CVec3_scalermult( C2 , C1_rj )) ;
			
			CVec3 ABt = CVec3_plus(At,Bt);
			
			CVec3 DN1N2_rj = CVec3_plus(ABt,Ct);
					//A1*A2_rj + A2*A1_rj + B1*B2_rj + B2*B1_rj + C1*C2_rj + C2*C1_rj;

			//DN1N2 := dot product dot(N1,N2), and "_rj" represents the partial 
			//derivative with respect to 1st, 2nd, and 3rd component of rj;
			//i.e. rj(1), rj(2), rj(3) being the x,y,z component of rj.

			CVec3 DN1N2_rk = 
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
                std::cout<<"DN1N2_rj:    "<< thrust::get<0>(DN1N2_rj)<<" "<< thrust::get<1>(DN1N2_rj)<<" "<< thrust::get<2>(DN1N2_rj)<<std::endl;
                std::cout<<"DN1N2_rk:    "<< thrust::get<0>(DN1N2_rk)<<" "<< thrust::get<1>(DN1N2_rk)<<" "<< thrust::get<2>(DN1N2_rk)<<std::endl;
                std::cout<<"DN1N2_ri:    "<< thrust::get<0>(DN1N2_ri)<<" "<< thrust::get<1>(DN1N2_ri)<<" "<< thrust::get<2>(DN1N2_ri)<<std::endl;
                std::cout<<"DN1N2_rl:    "<< thrust::get<0>(DN1N2_rl)<<" "<< thrust::get<1>(DN1N2_rl)<<" "<< thrust::get<2>(DN1N2_rl)<<std::endl;
                
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
				
                std::cout<<"PnN1nN2_rj: "<< thrust::get<0>(PnN1nN2_rj)<<" "<< thrust::get<1>(PnN1nN2_rj)<<" "<< thrust::get<2>(PnN1nN2_rj)<<std::endl;
                std::cout<<"PnN1nN2_rk: "<< thrust::get<0>(PnN1nN2_rk)<<" "<< thrust::get<1>(PnN1nN2_rk)<<" "<< thrust::get<2>(PnN1nN2_rk)<<std::endl;
                std::cout<<"PnN1nN2_ri: "<< thrust::get<0>(PnN1nN2_ri)<<" "<< thrust::get<1>(PnN1nN2_ri)<<" "<< thrust::get<2>(PnN1nN2_ri)<<std::endl;
                std::cout<<"PnN1nN2_rl: "<< thrust::get<0>(PnN1nN2_rl)<<" "<< thrust::get<1>(PnN1nN2_rl)<<" "<< thrust::get<2>(PnN1nN2_rl)<<std::endl;
				
				std::cout<<"dcN1N2_rj_1: "<< thrust::get<0>(dcN1N2_rj_1)<<" "<< thrust::get<1>(dcN1N2_rj_1)<<" "<< thrust::get<2>(dcN1N2_rj_1)<<std::endl;
                std::cout<<"dcN1N2_rj_2: "<< thrust::get<0>(dcN1N2_rj_2)<<" "<< thrust::get<1>(dcN1N2_rj_2)<<" "<< thrust::get<2>(dcN1N2_rj_2)<<std::endl;
                std::cout<<"dcN1N2_rj_3: "<< thrust::get<0>(dcN1N2_rj_3)<<" "<< thrust::get<1>(dcN1N2_rj_3)<<" "<< thrust::get<2>(dcN1N2_rj_3)<<std::endl;
                
	
				std::cout<<"dcN1N2_rk_1: "<< thrust::get<0>(dcN1N2_rk_1)<<" "<< thrust::get<1>(dcN1N2_rk_1)<<" "<< thrust::get<2>(dcN1N2_rk_1)<<std::endl;
                std::cout<<"dcN1N2_rk_2: "<< thrust::get<0>(dcN1N2_rk_2)<<" "<< thrust::get<1>(dcN1N2_rk_2)<<" "<< thrust::get<2>(dcN1N2_rk_2)<<std::endl;
                std::cout<<"dcN1N2_rk_3: "<< thrust::get<0>(dcN1N2_rk_3)<<" "<< thrust::get<1>(dcN1N2_rk_3)<<" "<< thrust::get<2>(dcN1N2_rk_3)<<std::endl;
                
				std::cout<<"dcN1N2_ri_1: "<< thrust::get<0>(dcN1N2_ri_1)<<" "<< thrust::get<1>(dcN1N2_ri_1)<<" "<< thrust::get<2>(dcN1N2_ri_1)<<std::endl;
                std::cout<<"dcN1N2_ri_2: "<< thrust::get<0>(dcN1N2_ri_2)<<" "<< thrust::get<1>(dcN1N2_ri_2)<<" "<< thrust::get<2>(dcN1N2_ri_2)<<std::endl;
                std::cout<<"dcN1N2_ri_3: "<< thrust::get<0>(dcN1N2_ri_3)<<" "<< thrust::get<1>(dcN1N2_ri_3)<<" "<< thrust::get<2>(dcN1N2_ri_3)<<std::endl;
    
				std::cout<<"dcN1N2_rl_1: "<< thrust::get<0>(dcN1N2_rl_1)<<" "<< thrust::get<1>(dcN1N2_rl_1)<<" "<< thrust::get<2>(dcN1N2_rl_1)<<std::endl;
                std::cout<<"dcN1N2_rl_2: "<< thrust::get<0>(dcN1N2_rl_2)<<" "<< thrust::get<1>(dcN1N2_rl_2)<<" "<< thrust::get<2>(dcN1N2_rl_2)<<std::endl;
                std::cout<<"dcN1N2_rl_3: "<< thrust::get<0>(dcN1N2_rl_3)<<" "<< thrust::get<1>(dcN1N2_rl_3)<<" "<< thrust::get<2>(dcN1N2_rl_3)<<std::endl;
				
				
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
					
					std::cout<<"N1N2_nN1nN2:" <<N1N2_nN1nN2<<std::endl;
					
					std::cout<<"nN2: " <<nN2<<std::endl; 
					
					std::cout<<"nN1: " <<nN1<<std::endl;
					std::cout<<"abc1: " <<A1<<"" <<B1<< " " <<C1<<std::endl;
					std::cout<<"abc2: " <<A2<<"" <<B2<< " " <<C2<<std::endl;
					
					std::cout<<"COSE_rj: "<< thrust::get<0>(COSE_rj)<<" "<< thrust::get<1>(COSE_rj)<<" "<< thrust::get<2>(COSE_rj)<<std::endl;
					std::cout<<"COSE_rk: "<< thrust::get<0>(COSE_rk)<<" "<< thrust::get<1>(COSE_rk)<<" "<< thrust::get<2>(COSE_rk)<<std::endl;
					std::cout<<"COSE_ri: "<< thrust::get<0>(COSE_ri)<<" "<< thrust::get<1>(COSE_ri)<<" "<< thrust::get<2>(COSE_ri)<<std::endl;
					std::cout<<"COSE_rl: "<< thrust::get<0>(COSE_rl)<<" "<< thrust::get<1>(COSE_rl)<<" "<< thrust::get<2>(COSE_rl)<<std::endl;
			
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
								//(nN1*nN2*(dot(cross(N1,N2), dUD_rl(3,:)) + dot(dcN1N2_rl(3,:), UD)) - dot(cross(N1,N2),UD)*PnN1nN2_rl(3))/(nN1*nN2)^2;

			std::cout<<"Sine_rj: "<< (SINE_rjx)<<" "<< (SINE_rjy)<<" "<< (SINE_rjz)<<std::endl;
			std::cout<<"Sine_rk: "<< (SINE_rkx)<<" "<< (SINE_rky)<<" "<< (SINE_rkz)<<std::endl;
			std::cout<<"Sine_ri: "<< (SINE_rix)<<" "<< (SINE_riy)<<" "<< (SINE_riz)<<std::endl;
			std::cout<<"Sine_rl: "<< (SINE_rlx)<<" "<< (SINE_rly)<<" "<< (SINE_rlz)<<std::endl;
			
			double angle_0 = 1.5707963267/2.0;
			double spring_constant = bendingTriangleInfoVecs.spring_constant;
			double place_0_x = spring_constant * (cos(angle_0) * (thrust::get<0>(COSE_ri)) + spring_constant * sin(angle_0) * SINE_rix);
			double place_0_y = spring_constant * (cos(angle_0) * (thrust::get<1>(COSE_ri)) + spring_constant * sin(angle_0) * SINE_riy);
			double place_0_z = spring_constant * (cos(angle_0) * (thrust::get<2>(COSE_ri)) + spring_constant * sin(angle_0) * SINE_riz);

			double place_1_x = spring_constant * (cos(angle_0) * (thrust::get<0>(COSE_rj)) + spring_constant * sin(angle_0) * SINE_rjx);
			double place_1_y = spring_constant * (cos(angle_0) * (thrust::get<1>(COSE_rj)) + spring_constant * sin(angle_0) * SINE_rjy);
			double place_1_z = spring_constant * (cos(angle_0) * (thrust::get<2>(COSE_rj)) + spring_constant * sin(angle_0) * SINE_rjz);

			double place_2_x = spring_constant * (cos(angle_0) * (thrust::get<0>(COSE_rk)) + spring_constant * sin(angle_0) * SINE_rkx);
			double place_2_y = spring_constant * (cos(angle_0) * (thrust::get<1>(COSE_rk)) + spring_constant * sin(angle_0) * SINE_rky);
			double place_2_z = spring_constant * (cos(angle_0) * (thrust::get<2>(COSE_rk)) + spring_constant * sin(angle_0) * SINE_rkz);

			double place_3_x = spring_constant * (cos(angle_0) * (thrust::get<0>(COSE_rl)) + spring_constant * sin(angle_0) * SINE_rlx);
			double place_3_y = spring_constant * (cos(angle_0) * (thrust::get<1>(COSE_rl)) + spring_constant * sin(angle_0) * SINE_rly);
			double place_3_z = spring_constant * (cos(angle_0) * (thrust::get<2>(COSE_rl)) + spring_constant * sin(angle_0) * SINE_rlz);
			std::cout<<"Res i : " <<id_i << " " << (place_0_x)<<" "<< (place_0_y)<<" "<< (place_0_z)<<std::endl;
			std::cout<<"Res j : " <<id_j << " " << (place_1_x)<<" "<< (place_1_y)<<" "<< (place_1_z)<<std::endl;
			std::cout<<"Res k : " <<id_k << " " << (place_2_x)<<" "<< (place_2_y)<<" "<< (place_2_z)<<std::endl;
			std::cout<<"Res l : " <<id_l << " " << (place_3_x)<<" "<< (place_3_y)<<" "<< (place_3_z)<<std::endl;
			
												
        }
    
	}*/
    
    thrust::counting_iterator<unsigned> elemId(0); 

	//bendingTriangleInfoVecs.initial_angle = 1.5707963267/2.0;
	thrust::fill(bendingTriangleInfoVecs.tempNodeForceXReduced.begin(),bendingTriangleInfoVecs.tempNodeForceXReduced.end(),0.0);
	thrust::fill(bendingTriangleInfoVecs.tempNodeForceYReduced.begin(),bendingTriangleInfoVecs.tempNodeForceYReduced.end(),0.0);
	thrust::fill(bendingTriangleInfoVecs.tempNodeForceZReduced.begin(),bendingTriangleInfoVecs.tempNodeForceZReduced.end(),0.0);
	thrust::fill(bendingTriangleInfoVecs.tempNodeForceXUnreduced.begin(),bendingTriangleInfoVecs.tempNodeForceXUnreduced.end(),0.0);
	thrust::fill(bendingTriangleInfoVecs.tempNodeForceYUnreduced.begin(),bendingTriangleInfoVecs.tempNodeForceYUnreduced.end(),0.0);
	thrust::fill(bendingTriangleInfoVecs.tempNodeForceZUnreduced.begin(),bendingTriangleInfoVecs.tempNodeForceZUnreduced.end(),0.0);

    //apply force to temporary vectors.
    bendingTriangleInfoVecs.bending_triangle_energy= 
    thrust::transform_reduce(
        thrust::make_zip_iterator(
            thrust::make_tuple(
				elemId,
                coordInfoVecs.edges2Triangles_1.begin(),
                coordInfoVecs.edges2Triangles_2.begin(),
                coordInfoVecs.edges2Nodes_1.begin(),
                coordInfoVecs.edges2Nodes_2.begin())),
        thrust::make_zip_iterator(
            thrust::make_tuple(
				elemId, 
                coordInfoVecs.edges2Triangles_1.begin(),
                coordInfoVecs.edges2Triangles_2.begin(),
                coordInfoVecs.edges2Nodes_1.begin(),
                coordInfoVecs.edges2Nodes_2.begin())) + coordInfoVecs.num_edges,
        CosBendingFunctor(
            bendingTriangleInfoVecs.spring_constant,
            bendingTriangleInfoVecs.initial_angle,        
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data()),
            
            thrust::raw_pointer_cast(bendingTriangleInfoVecs.tempNodeIdUnreduced.data()),
            thrust::raw_pointer_cast(bendingTriangleInfoVecs.tempNodeForceXUnreduced.data()),
            thrust::raw_pointer_cast(bendingTriangleInfoVecs.tempNodeForceYUnreduced.data()),
            thrust::raw_pointer_cast(bendingTriangleInfoVecs.tempNodeForceZUnreduced.data()),
            thrust::raw_pointer_cast(coordInfoVecs.triangles2Nodes_1.data()), 
            thrust::raw_pointer_cast(coordInfoVecs.triangles2Nodes_2.data()), 
            thrust::raw_pointer_cast(coordInfoVecs.triangles2Nodes_3.data())),
		0.0, thrust::plus<double>() );
	
/*	for (unsigned i = 0; i < bendingTriangleInfoVecs.tempNodeIdUnreduced.size(); i++) {

		std::cout<<"id: "<< bendingTriangleInfoVecs.tempNodeIdUnreduced[i]<<std::endl;
		std::cout<< "unreduced F_x: "<< bendingTriangleInfoVecs.tempNodeForceXUnreduced[i]<<std::endl;
		std::cout<< "unreduced F_y: "<< bendingTriangleInfoVecs.tempNodeForceYUnreduced[i]<<std::endl;
		std::cout<< "unreduced F_z: "<< bendingTriangleInfoVecs.tempNodeForceZUnreduced[i]<<std::endl;
	}*/
    //now we have un reduced forces. Sort by id and reduce. 
    //key, then value. Each vector returns sorted		
    thrust::sort_by_key(bendingTriangleInfoVecs.tempNodeIdUnreduced.begin(), bendingTriangleInfoVecs.tempNodeIdUnreduced.end(),
        thrust::make_zip_iterator(
            thrust::make_tuple(
                bendingTriangleInfoVecs.tempNodeForceXUnreduced.begin(),
                bendingTriangleInfoVecs.tempNodeForceYUnreduced.begin(),
                bendingTriangleInfoVecs.tempNodeForceZUnreduced.begin())), thrust::less<unsigned>());
    
    unsigned endKey = thrust::get<0>(
        thrust::reduce_by_key(
            bendingTriangleInfoVecs.tempNodeIdUnreduced.begin(), 
            bendingTriangleInfoVecs.tempNodeIdUnreduced.end(),
        thrust::make_zip_iterator(
            thrust::make_tuple(
                bendingTriangleInfoVecs.tempNodeForceXUnreduced.begin(),
                bendingTriangleInfoVecs.tempNodeForceYUnreduced.begin(),
                bendingTriangleInfoVecs.tempNodeForceZUnreduced.begin())),
			
			bendingTriangleInfoVecs.tempNodeIdReduced.begin(),
        thrust::make_zip_iterator(
            thrust::make_tuple(
                bendingTriangleInfoVecs.tempNodeForceXReduced.begin(),
                bendingTriangleInfoVecs.tempNodeForceYReduced.begin(),
                bendingTriangleInfoVecs.tempNodeForceZReduced.begin())),
		thrust::equal_to<unsigned>(), CVec3Add())) - bendingTriangleInfoVecs.tempNodeIdReduced.begin();//binary_pred, binary_op 
		
    /*	for (unsigned i = 0; i < bendingTriangleInfoVecs.tempNodeIdReduced.size(); i++) {

			std::cout<<"id: "<< bendingTriangleInfoVecs.tempNodeIdReduced[i]<<std::endl;
			std::cout<< "reduced F_x: "<< bendingTriangleInfoVecs.tempNodeForceXReduced[i]<<std::endl;
			std::cout<< "reduced F_y: "<< bendingTriangleInfoVecs.tempNodeForceYReduced[i]<<std::endl;
			std::cout<< "reduced F_z: "<< bendingTriangleInfoVecs.tempNodeForceZReduced[i]<<std::endl;
		}*/
     //apply reduced force to all nodes. 
    thrust::for_each(
        thrust::make_zip_iterator(//1st begin
            thrust::make_tuple(
                bendingTriangleInfoVecs.tempNodeIdReduced.begin(),
                bendingTriangleInfoVecs.tempNodeForceXReduced.begin(),
                bendingTriangleInfoVecs.tempNodeForceYReduced.begin(),
                bendingTriangleInfoVecs.tempNodeForceZReduced.begin())),
        thrust::make_zip_iterator(//1st end
            thrust::make_tuple(
                bendingTriangleInfoVecs.tempNodeIdReduced.begin(),
                bendingTriangleInfoVecs.tempNodeForceXReduced.begin(),
                bendingTriangleInfoVecs.tempNodeForceYReduced.begin(),
                bendingTriangleInfoVecs.tempNodeForceZReduced.begin())) + endKey,
        AddForceFunctor (
            thrust::raw_pointer_cast(coordInfoVecs.nodeForceX.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeForceY.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeForceZ.data())));

};
