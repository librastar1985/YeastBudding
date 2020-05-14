#include "System.h"
#include "MemRepulsionSprings_local.h"



//WARNING: function must not reset coordInfoVecs.nodeForceX etc. 
void ComputeMemRepulsionSprings_local(
    CoordInfoVecs& coordInfoVecs,
	LinearSpringInfoVecs& linearSpringInfoVecs, 
	CapsidInfoVecs& capsidInfoVecs,
    GeneralParams& generalParams,
    AuxVecs& auxVecs) {    
    
    /*thrust::fill(coordInfoVecs.tempNodeForceXReduced.begin(),coordInfoVecs.tempNodeForceXReduced.end(),0.0);
    thrust::fill(coordInfoVecs.tempNodeForceYReduced.begin(),coordInfoVecs.tempNodeForceYReduced.end(),0.0);
    thrust::fill(coordInfoVecs.tempNodeForceZReduced.begin(),coordInfoVecs.tempNodeForceZReduced.end(),0.0);
    thrust::fill(coordInfoVecs.tempNodeForceXUnreduced.begin(),coordInfoVecs.tempNodeForceXUnreduced.end(),0.0);
    thrust::fill(coordInfoVecs.tempNodeForceYUnreduced.begin(),coordInfoVecs.tempNodeForceYUnreduced.end(),0.0);
    thrust::fill(coordInfoVecs.tempNodeForceZUnreduced.begin(),coordInfoVecs.tempNodeForceZUnreduced.end(),0.0);*/

    //CVec4 init(0.0, 0.0, 0.0, 0.0); 
    thrust::counting_iterator<int> begin(0);

    thrust::transform(  
        thrust::make_zip_iterator(
            thrust::make_tuple(
                begin,
                auxVecs.id_bucket.begin(),
                coordInfoVecs.nodeForceX.begin(),
                coordInfoVecs.nodeForceY.begin(),
                coordInfoVecs.nodeForceZ.begin())),
        
        thrust::make_zip_iterator(
            thrust::make_tuple(
                begin,
                auxVecs.id_bucket.begin(),
                coordInfoVecs.nodeForceX.begin(),
                coordInfoVecs.nodeForceY.begin(),
                coordInfoVecs.nodeForceZ.begin())) + generalParams.maxNodeCount,

        thrust::make_zip_iterator(
            thrust::make_tuple(
                coordInfoVecs.nodeForceX.begin(),
                coordInfoVecs.nodeForceY.begin(),
                coordInfoVecs.nodeForceZ.begin())),

        MemRepulsionSpringFunctor_local(
                generalParams.Rmin,
                generalParams.abs_Rmin,
                //linearSpringInfoVecs.spring_constant_rep1,
                //linearSpringInfoVecs.spring_constant_rep2,
                linearSpringInfoVecs.spring_constant_weak,
                linearSpringInfoVecs.spring_constant,
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
               // thrust::raw_pointer_cast(coordInfoVecs.nndata10.data()),
                //thrust::raw_pointer_cast(coordInfoVecs.nndata11.data()),
                //thrust::raw_pointer_cast(coordInfoVecs.nndata12.data()),

                thrust::raw_pointer_cast(generalParams.nodes_in_upperhem.data()),
                thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
                thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
                thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data()), 

                
                          
                thrust::raw_pointer_cast(auxVecs.id_value_expanded.data()),
                thrust::raw_pointer_cast(auxVecs.keyBegin.data()),
                thrust::raw_pointer_cast(auxVecs.keyEnd.data()))
            );
                     
};

