
#include "System.h"
#include "SystemStructures.h"
#include "LinearSprings.h"

void ComputeLinearSprings(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    LJInfoVecs& ljInfoVecs) {

        thrust::fill(linearSpringInfoVecs.tempNodeForceXReduced.begin(),linearSpringInfoVecs.tempNodeForceXReduced.end(),0.0);
        thrust::fill(linearSpringInfoVecs.tempNodeForceYReduced.begin(),linearSpringInfoVecs.tempNodeForceYReduced.end(),0.0);
        thrust::fill(linearSpringInfoVecs.tempNodeForceZReduced.begin(),linearSpringInfoVecs.tempNodeForceZReduced.end(),0.0);
        thrust::fill(linearSpringInfoVecs.tempNodeForceXUnreduced.begin(),linearSpringInfoVecs.tempNodeForceXUnreduced.end(),0.0);
        thrust::fill(linearSpringInfoVecs.tempNodeForceYUnreduced.begin(),linearSpringInfoVecs.tempNodeForceYUnreduced.end(),0.0);
        thrust::fill(linearSpringInfoVecs.tempNodeForceZUnreduced.begin(),linearSpringInfoVecs.tempNodeForceZUnreduced.end(),0.0);
    
    
        thrust::counting_iterator<unsigned> edgeIdBegin(0);
        thrust::counting_iterator<unsigned> edgeIdEnd(coordInfoVecs.num_edges);

    //std::cout<<"pre linear spring: " <<coordInfoVecs.nodeForceX.size()<<std::endl;
   /* unsigned id =ljInfoVecs.node_id_close[0];
	std::cout<<"partPos: " << coordInfoVecs.nodeLocX[id]<< " "<< coordInfoVecs.nodeLocY[id] << " "<< coordInfoVecs.nodeLocZ[id] << std::endl;
	std::cout<<"partForce: " << coordInfoVecs.nodeForceX[id]<< " "<< coordInfoVecs.nodeForceY[id] << " "<< coordInfoVecs.nodeForceZ[id] << std::endl;
      */  

    linearSpringInfoVecs.linear_spring_energy=
    thrust::transform_reduce(  
        thrust::make_zip_iterator(
            thrust::make_tuple(
                edgeIdBegin,
                coordInfoVecs.edges2Nodes_1.begin(),
                coordInfoVecs.edges2Nodes_2.begin(),
                linearSpringInfoVecs.edge_initial_length.begin())),
        thrust::make_zip_iterator( 
            thrust::make_tuple(
                edgeIdEnd,
                coordInfoVecs.edges2Nodes_1.end(),
                coordInfoVecs.edges2Nodes_2.end(), 
                linearSpringInfoVecs.edge_initial_length.end())),
        LinearSpringFunctor(
            linearSpringInfoVecs.spring_constant, 
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data()),

            thrust::raw_pointer_cast(linearSpringInfoVecs.tempNodeIdUnreduced.data()),
            thrust::raw_pointer_cast(linearSpringInfoVecs.tempNodeForceXUnreduced.data()),
            thrust::raw_pointer_cast(linearSpringInfoVecs.tempNodeForceYUnreduced.data()),
            thrust::raw_pointer_cast(linearSpringInfoVecs.tempNodeForceZUnreduced.data()) ),
        0.0, thrust::plus<double>() ); 
      

    //now we have un reduced forces. Sort by id and reduce. 
    //key, then value. Each vector returns sorted		
    thrust::sort_by_key(linearSpringInfoVecs.tempNodeIdUnreduced.begin(), linearSpringInfoVecs.tempNodeIdUnreduced.end(),
        thrust::make_zip_iterator(
            thrust::make_tuple(
                linearSpringInfoVecs.tempNodeForceXUnreduced.begin(),
                linearSpringInfoVecs.tempNodeForceYUnreduced.begin(),
                linearSpringInfoVecs.tempNodeForceZUnreduced.begin())), thrust::less<unsigned>());

   /* std::cout<<"mid1 linear spring: " <<coordInfoVecs.nodeForceX.size()<<std::endl;
    std::cout<<"partPos: " << coordInfoVecs.nodeLocX[id]<< " "<< coordInfoVecs.nodeLocY[id] << " "<< coordInfoVecs.nodeLocZ[id] << std::endl;
	std::cout<<"partForce: " << coordInfoVecs.nodeForceX[id]<< " "<< coordInfoVecs.nodeForceY[id] << " "<< coordInfoVecs.nodeForceZ[id] << std::endl;
	*/
    unsigned endKey = thrust::get<0>(
        thrust::reduce_by_key(
            linearSpringInfoVecs.tempNodeIdUnreduced.begin(), 
            linearSpringInfoVecs.tempNodeIdUnreduced.end(),
        thrust::make_zip_iterator(
            thrust::make_tuple(
                linearSpringInfoVecs.tempNodeForceXUnreduced.begin(),
                linearSpringInfoVecs.tempNodeForceYUnreduced.begin(),
                linearSpringInfoVecs.tempNodeForceZUnreduced.begin())),
            linearSpringInfoVecs.tempNodeIdReduced.begin(),
        thrust::make_zip_iterator(
            thrust::make_tuple(
                linearSpringInfoVecs.tempNodeForceXReduced.begin(),
                linearSpringInfoVecs.tempNodeForceYReduced.begin(),
                linearSpringInfoVecs.tempNodeForceZReduced.begin())),
        thrust::equal_to<unsigned>(), CVec3Add())) - linearSpringInfoVecs.tempNodeIdReduced.begin();//binary_pred, binary_op 
    
/*
    std::cout<<"mid2 linear spring: " <<coordInfoVecs.nodeForceX.size()<<std::endl;
    std::cout<<"partPos: " << coordInfoVecs.nodeLocX[id]<< " "<< coordInfoVecs.nodeLocY[id] << " "<< coordInfoVecs.nodeLocZ[id] << std::endl;
	std::cout<<"partForce: " << coordInfoVecs.nodeForceX[id]<< " "<< coordInfoVecs.nodeForceY[id] << " "<< coordInfoVecs.nodeForceZ[id] << std::endl;
	*/
    //apply reduced force to all nodes. 
    thrust::for_each(
        thrust::make_zip_iterator(//1st begin
            thrust::make_tuple(
                linearSpringInfoVecs.tempNodeIdReduced.begin(),
                linearSpringInfoVecs.tempNodeForceXReduced.begin(),
                linearSpringInfoVecs.tempNodeForceYReduced.begin(),
                linearSpringInfoVecs.tempNodeForceZReduced.begin())),
        thrust::make_zip_iterator(//1st end
            thrust::make_tuple(
                linearSpringInfoVecs.tempNodeIdReduced.begin(),
                linearSpringInfoVecs.tempNodeForceXReduced.begin(),
                linearSpringInfoVecs.tempNodeForceYReduced.begin(),
                linearSpringInfoVecs.tempNodeForceZReduced.begin())) + endKey,
        AddForceFunctor (
            thrust::raw_pointer_cast(coordInfoVecs.nodeForceX.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeForceY.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeForceZ.data())));
    /*
    std::cout<<"post linear spring: " <<coordInfoVecs.nodeForceX.size()<<std::endl;
    std::cout<<"partPos: " << coordInfoVecs.nodeLocX[id]<< " "<< coordInfoVecs.nodeLocY[id] << " "<< coordInfoVecs.nodeLocZ[id] << std::endl;
	std::cout<<"partForce: " << coordInfoVecs.nodeForceX[id]<< " "<< coordInfoVecs.nodeForceY[id] << " "<< coordInfoVecs.nodeForceZ[id] << std::endl;
	*/

};