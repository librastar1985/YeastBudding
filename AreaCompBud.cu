
#include "System.h"
#include "SystemStructures.h"
#include "AreaCompBud.h"

void ComputeAreaBud(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs//,
    //LinearSpringInfoVecs& linearSpringInfoVecs,
    //LJInfoVecs& ljInfoVecs
    ) {  
    
        thrust::counting_iterator<int> triangleIdBegin(0);
       // thrust::counting_iterator<int> triangleIdEnd(coordInfoVecs.num_triangles);//coordInfoVecs.triangles2Nodes_1.size());

    generalParams.current_bud_area=
    thrust::transform_reduce(  
        thrust::make_zip_iterator(
            thrust::make_tuple(
                triangleIdBegin,
                coordInfoVecs.triangles2Nodes_1.begin(),
                coordInfoVecs.triangles2Nodes_2.begin(),
                coordInfoVecs.triangles2Nodes_3.begin())),
        thrust::make_zip_iterator( 
            thrust::make_tuple(
                triangleIdBegin,
                coordInfoVecs.triangles2Nodes_1.begin(),
                coordInfoVecs.triangles2Nodes_2.begin(), 
                coordInfoVecs.triangles2Nodes_3.begin())) + coordInfoVecs.num_triangles,
        AreaBudCompFunctor(
            //linearSpringInfoVecs.spring_constant, 
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data()),
            thrust::raw_pointer_cast(generalParams.triangles_in_upperhem.data()) ),
        0.0, thrust::plus<double>() ); 
        //This sum is the part without the absolute value and factor of (1/6) in the formula
        

};