
#include "System.h"
#include "LJSprings_LJ.h"



//WARNING: function must not reset coordInfoVecs.nodeForceX etc. 
void ComputeLJSprings_LJ(
    CoordInfoVecs& coordInfoVecs,
    LJInfoVecs& ljInfoVecs,  
    GeneralParams& generalParams) {    
     
    CVec4 init(0.0, 0.0, 0.0, 0.0); 
    thrust::counting_iterator<int> begin(0);
    thrust::counting_iterator<int> end(generalParams.maxNodeCountLJ);
CVec4 temp = thrust::transform_reduce(  
                    begin,end,
                    LJSpringFunctor_LJ(
                        ljInfoVecs.Rcutoff_LJ, 
                        ljInfoVecs.Rmin_LJ,
                        ljInfoVecs.epsilon_LJ_rep1,
                        ljInfoVecs.epsilon_LJ_rep2,
                        ljInfoVecs.LJ_PosX,
                        ljInfoVecs.LJ_PosY,
                        ljInfoVecs.LJ_PosZ,
                        thrust::raw_pointer_cast(ljInfoVecs.LJ_PosX_all.data()),
                        thrust::raw_pointer_cast(ljInfoVecs.LJ_PosY_all.data()),
                        thrust::raw_pointer_cast(ljInfoVecs.LJ_PosZ_all.data()),
                        thrust::raw_pointer_cast(ljInfoVecs.forceX_all.data()),
                        thrust::raw_pointer_cast(ljInfoVecs.forceY_all.data()),
                        thrust::raw_pointer_cast(ljInfoVecs.forceZ_all.data())),
                    init, CVec4Add() );//unary, binary 
     
    //temp stores forces on lj particle and last entry as energy
    ljInfoVecs.forceX = thrust::get<0>(temp);
    ljInfoVecs.forceY = thrust::get<1>(temp);
    ljInfoVecs.forceZ = thrust::get<2>(temp);
    //std::cout<<ljInfoVecs.forceZ<<std::endl;

    ljInfoVecs.lj_energy_LJ = thrust::get<3>(temp);
      
    //std::cout<<"lj force from other lj "<< ljInfoVecs.forceX<< " "<<  ljInfoVecs.forceY << " "<<  ljInfoVecs.forceZ << std::endl;
	//std::cout<<"lj points "<< ljInfoVecs.LJ_PosX<< " "<<  ljInfoVecs.LJ_PosY << " "<<  ljInfoVecs.LJ_PosZ << std::endl;

};

