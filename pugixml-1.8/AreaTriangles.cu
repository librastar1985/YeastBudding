#include "System.h"
#include "SystemStructures.h"
#include "AreaTriangles.h"

void ComputeAreaTriangleSprings(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    AreaTriangleInfoVecs& areaTriangleInfoVecs) {

        thrust::counting_iterator<unsigned> elemId(0);
        
        thrust::fill(areaTriangleInfoVecs.tempNodeForceXReduced.begin(),areaTriangleInfoVecs.tempNodeForceXReduced.end(),0.0);
        thrust::fill(areaTriangleInfoVecs.tempNodeForceYReduced.begin(),areaTriangleInfoVecs.tempNodeForceYReduced.end(),0.0);
        thrust::fill(areaTriangleInfoVecs.tempNodeForceZReduced.begin(),areaTriangleInfoVecs.tempNodeForceZReduced.end(),0.0);
        thrust::fill(areaTriangleInfoVecs.tempNodeForceXUnreduced.begin(),areaTriangleInfoVecs.tempNodeForceXUnreduced.end(),0.0);
        thrust::fill(areaTriangleInfoVecs.tempNodeForceYUnreduced.begin(),areaTriangleInfoVecs.tempNodeForceYUnreduced.end(),0.0);
        thrust::fill(areaTriangleInfoVecs.tempNodeForceZUnreduced.begin(),areaTriangleInfoVecs.tempNodeForceZUnreduced.end(),0.0);
    
        
        areaTriangleInfoVecs.area_triangle_energy = thrust::transform_reduce( 
			thrust::make_zip_iterator(
				thrust::make_tuple(
                    elemId,
					coordInfoVecs.triangles2Nodes_1.begin(),
					coordInfoVecs.triangles2Nodes_2.begin(),
					coordInfoVecs.triangles2Nodes_3.begin())),
			thrust::make_zip_iterator(
				thrust::make_tuple(
                    elemId,
					coordInfoVecs.triangles2Nodes_1.begin(),
					coordInfoVecs.triangles2Nodes_2.begin(),
					coordInfoVecs.triangles2Nodes_3.begin())) + coordInfoVecs.num_triangles,
            AreaSpringFunctor( 
                areaTriangleInfoVecs.initial_area,
                areaTriangleInfoVecs.spring_constant,
                thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
                thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
                thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data()),
                thrust::raw_pointer_cast(areaTriangleInfoVecs.tempNodeIdUnreduced.data()),
                thrust::raw_pointer_cast(areaTriangleInfoVecs.tempNodeForceXUnreduced.data()),
                thrust::raw_pointer_cast(areaTriangleInfoVecs.tempNodeForceYUnreduced.data()),
                thrust::raw_pointer_cast(areaTriangleInfoVecs.tempNodeForceZUnreduced.data())),
            0.0, thrust::plus<double>() );
        
        
        //now we have un reduced forces. Sort by id and reduce. 
        //key, then value. Each vector returns sorted		
		thrust::sort_by_key(areaTriangleInfoVecs.tempNodeIdUnreduced.begin(), areaTriangleInfoVecs.tempNodeIdUnreduced.end(),
			thrust::make_zip_iterator(
				thrust::make_tuple(
					areaTriangleInfoVecs.tempNodeForceXUnreduced.begin(),
					areaTriangleInfoVecs.tempNodeForceYUnreduced.begin(),
					areaTriangleInfoVecs.tempNodeForceZUnreduced.begin())), thrust::less<unsigned>());

        unsigned endKey = thrust::get<0>(
            thrust::reduce_by_key(
                areaTriangleInfoVecs.tempNodeIdUnreduced.begin(), 
                areaTriangleInfoVecs.tempNodeIdUnreduced.end(),
            thrust::make_zip_iterator(
                thrust::make_tuple(
                    areaTriangleInfoVecs.tempNodeForceXUnreduced.begin(),
                    areaTriangleInfoVecs.tempNodeForceYUnreduced.begin(),
                    areaTriangleInfoVecs.tempNodeForceZUnreduced.begin())),
                areaTriangleInfoVecs.tempNodeIdReduced.begin(),
            thrust::make_zip_iterator(
                thrust::make_tuple(
                    areaTriangleInfoVecs.tempNodeForceXReduced.begin(),
                    areaTriangleInfoVecs.tempNodeForceYReduced.begin(),
                    areaTriangleInfoVecs.tempNodeForceZReduced.begin())),
            thrust::equal_to<unsigned>(), CVec3Add())) - areaTriangleInfoVecs.tempNodeIdReduced.begin();//binary_pred, binary_op 
        

        //apply reduced force to all nodes. 
        thrust::for_each(
            thrust::make_zip_iterator(//1st begin
                thrust::make_tuple(
                    areaTriangleInfoVecs.tempNodeIdReduced.begin(),
                    areaTriangleInfoVecs.tempNodeForceXReduced.begin(),
                    areaTriangleInfoVecs.tempNodeForceYReduced.begin(),
                    areaTriangleInfoVecs.tempNodeForceZReduced.begin())),
            thrust::make_zip_iterator(//1st end
                thrust::make_tuple(
                    areaTriangleInfoVecs.tempNodeIdReduced.begin(),
                    areaTriangleInfoVecs.tempNodeForceXReduced.begin(),
                    areaTriangleInfoVecs.tempNodeForceYReduced.begin(),
                    areaTriangleInfoVecs.tempNodeForceZReduced.begin())) + endKey,
            AddForceFunctor (
                thrust::raw_pointer_cast(coordInfoVecs.nodeForceX.data()),
                thrust::raw_pointer_cast(coordInfoVecs.nodeForceY.data()),
                thrust::raw_pointer_cast(coordInfoVecs.nodeForceZ.data())));
            
};