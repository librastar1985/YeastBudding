
#include "System.h"
#include "SystemStructures.h"
#include "VolumeComp.h"

void ComputeVolume(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    LJInfoVecs& ljInfoVecs) {  
    
        thrust::counting_iterator<int> triangleIdBegin(0);
       // thrust::counting_iterator<int> triangleIdEnd(coordInfoVecs.num_triangles);//coordInfoVecs.triangles2Nodes_1.size());

    generalParams.current_total_volume=
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
        VolumeCompFunctor(
            linearSpringInfoVecs.spring_constant, 
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data()) ),
        0.0, thrust::plus<double>() ); 
        //This sum is the part without the absolute value and factor of (1/6) in the formula
        
        generalParams.true_current_total_volume = sqrt(generalParams.current_total_volume*generalParams.current_total_volume)/6.0;
        generalParams.volume_energy = generalParams.volume_spring_constant*(generalParams.true_current_total_volume - generalParams.eq_total_volume)*
                                        (generalParams.true_current_total_volume - generalParams.eq_total_volume)/
                                        (2.0*generalParams.Rmin*generalParams.Rmin*generalParams.Rmin*generalParams.eq_total_volume);
        //std::cout<<"Current volume = "<<generalParams.current_total_volume<<std::endl;

};