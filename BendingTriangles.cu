
#include "System.h"
#include "SystemStructures.h"
#include "BendingTriangles.h"

void ComputeCosTriangleSprings(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    BendingTriangleInfoVecs& bendingTriangleInfoVecs) {
	
    
    thrust::counting_iterator<int> elemId(0); 

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
            generalParams.SCALE_TYPE,
            generalParams.scaling_pow,
            generalParams.gausssigma,
            generalParams.hilleqnconst,
            generalParams.hilleqnpow,
            thrust::raw_pointer_cast(coordInfoVecs.scaling_per_edge.data()),
            bendingTriangleInfoVecs.spring_constant,
            bendingTriangleInfoVecs.spring_constant_weak,
            thrust::raw_pointer_cast(generalParams.edges_in_upperhem.data()),
            thrust::raw_pointer_cast(generalParams.boundaries_in_upperhem.data()),
            bendingTriangleInfoVecs.initial_angle, 
            //bendingTriangleInfoVecs.initial_angle_bud,       
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
	
/*	for (int i = 0; i < bendingTriangleInfoVecs.tempNodeIdUnreduced.size(); i++) {

		std::cout<<"id: "<< bendingTriangleInfoVecs.tempNodeIdUnreduced[i]<<std::endl;
		std::cout<< "unreduced F_x: "<< bendingTriangleInfoVecs.tempNodeForceXUnreduced[i]<<std::endl;
		std::cout<< "unreduced F_y: "<< bendingTriangleInfoVecs.tempNodeForceYUnreduced[i]<<std::endl;
		std::cout<< "unreduced F_z: "<< bendingTriangleInfoVecs.tempNodeForceZUnreduced[i]<<std::endl;
	}*/
    //now we have un reduced forces. Sort by id and reduce. 
    //key, then value. Each vector returns sorted		
    thrust::sort_by_key(bendingTriangleInfoVecs.tempNodeIdUnreduced.begin(), bendingTriangleInfoVecs.tempNodeIdUnreduced.begin() + (bendingTriangleInfoVecs.factor*coordInfoVecs.num_edges),
        thrust::make_zip_iterator(
            thrust::make_tuple(
                bendingTriangleInfoVecs.tempNodeForceXUnreduced.begin(),
                bendingTriangleInfoVecs.tempNodeForceYUnreduced.begin(),
                bendingTriangleInfoVecs.tempNodeForceZUnreduced.begin())), thrust::less<int>());
    
    int endKey = thrust::get<0>(
        thrust::reduce_by_key(
            bendingTriangleInfoVecs.tempNodeIdUnreduced.begin(), 
            bendingTriangleInfoVecs.tempNodeIdUnreduced.begin() + (bendingTriangleInfoVecs.factor*coordInfoVecs.num_edges),
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
		thrust::equal_to<int>(), CVec3Add())) - bendingTriangleInfoVecs.tempNodeIdReduced.begin();//binary_pred, binary_op 
		//std::cout<<"endKey_bend = "<<endKey<<std::endl;
    /*	for (int i = 0; i < bendingTriangleInfoVecs.tempNodeIdReduced.size(); i++) {

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

        //    std::cout<<"Force from bend on node 36 = "<<coordInfoVecs.nodeForceX[35]<<" "<<coordInfoVecs.nodeForceY[35]<<" "<<coordInfoVecs.nodeForceZ[35]<<std::endl;
        //std::cout<<"BEND"<<std::endl;
};
