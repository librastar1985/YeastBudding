#include "System.h"
#include "Nodes2Triangles.h"
#include "SystemStructures.h"



//WARNING: function must not reset coordInfoVecs.nodeForceX etc. 
void ComputeNodes2Triangles(
    CoordInfoVecs& coordInfoVecs,
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
                auxVecs.id_bucket.begin()
            )),
        
        thrust::make_zip_iterator(
            thrust::make_tuple(
                begin,
                auxVecs.id_bucket.begin())) + generalParams.maxNodeCount,

        thrust::make_zip_iterator(
            thrust::make_tuple(
                coordInfoVecs.nodes2Triangles_1.begin(),
                coordInfoVecs.nodes2Triangles_2.begin(),
                coordInfoVecs.nodes2Triangles_3.begin(),
                coordInfoVecs.nodes2Triangles_4.begin(),
                coordInfoVecs.nodes2Triangles_5.begin(),
                coordInfoVecs.nodes2Triangles_6.begin(),
                coordInfoVecs.nodes2Triangles_7.begin(),
                coordInfoVecs.nodes2Triangles_8.begin(),
                coordInfoVecs.nodes2Triangles_9.begin()
               )),

        Nodes2TrianglesFunctor(
            
            coordInfoVecs.num_triangles,

            thrust::raw_pointer_cast(coordInfoVecs.triangles2Nodes_1.data()),
            thrust::raw_pointer_cast(coordInfoVecs.triangles2Nodes_2.data()),
            thrust::raw_pointer_cast(coordInfoVecs.triangles2Nodes_3.data()),              
                          
                thrust::raw_pointer_cast(auxVecs.id_value_expanded.data()),
                thrust::raw_pointer_cast(auxVecs.keyBegin.data()),
                thrust::raw_pointer_cast(auxVecs.keyEnd.data()))
            );
                     
};

