#include "System.h"
#include "MemRepulsionEnergy.h"



//WARNING: function must not reset coordInfoVecs.nodeForceX etc. 
double ComputeMemRepulsionEnergy(
    CoordInfoVecs& coordInfoVecs,
	LinearSpringInfoVecs& linearSpringInfoVecs, 
	CapsidInfoVecs& capsidInfoVecs,
    GeneralParams& generalParams,
    AuxVecs& auxVecs) {    
    
    thrust::counting_iterator<int> nodeIdBegin(0);
	thrust::counting_iterator<int> nodeIdEnd(generalParams.maxNodeCount);


    linearSpringInfoVecs.memrepulsion_energy=
    thrust::transform_reduce(  
        thrust::make_zip_iterator(
            thrust::make_tuple(
                nodeIdBegin,
                coordInfoVecs.nodeLocX.begin(),
                coordInfoVecs.nodeLocY.begin(),
                coordInfoVecs.nodeLocZ.begin() )),
        thrust::make_zip_iterator( 
            thrust::make_tuple(
                nodeIdBegin,
                coordInfoVecs.nodeLocX.begin(),
                coordInfoVecs.nodeLocY.begin(),
                coordInfoVecs.nodeLocZ.begin() )) + generalParams.maxNodeCount,
        MemRepulsionEnergyFunctor(
            generalParams.Rmin,
            generalParams.abs_Rmin,
            linearSpringInfoVecs.spring_constant_rep1,
            linearSpringInfoVecs.spring_constant_rep2,
            generalParams.maxNodeCount,

            thrust::raw_pointer_cast(coordInfoVecs.nndata1.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nndata2.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nndata3.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nndata4.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nndata5.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nndata6.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nndata7.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nndata8.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nndata9.data()),
            //thrust::raw_pointer_cast(coordInfoVecs.nndata10.data()),
            //thrust::raw_pointer_cast(coordInfoVecs.nndata11.data()),
            //thrust::raw_pointer_cast(coordInfoVecs.nndata12.data()),


            thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data()),          

            thrust::raw_pointer_cast(auxVecs.id_value_expanded.data()),
            thrust::raw_pointer_cast(auxVecs.keyBegin.data()),
            thrust::raw_pointer_cast(auxVecs.keyEnd.data()) ),
        0.0, thrust::plus<double>() ); 

return linearSpringInfoVecs.memrepulsion_energy;

   
};

