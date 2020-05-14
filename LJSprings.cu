
#include "System.h"
#include "LJSprings.h"



//WARNING: function must not reset coordInfoVecs.nodeForceX etc. 
void ComputeLJSprings(
    CoordInfoVecs& coordInfoVecs,
    LJInfoVecs& ljInfoVecs,  
    GeneralParams& generalParams) {    
     
    CVec4 init(0.0, 0.0, 0.0, 0.0);  
    thrust::counting_iterator<int> begin(0);
    thrust::counting_iterator<int> end(generalParams.maxNodeCount);
CVec4 temp = thrust::transform_reduce(  
                    begin,end,
                    LJSpringFunctor(
                        ljInfoVecs.Rcutoff_M, 
                        ljInfoVecs.Rmin_M,
                        ljInfoVecs.epsilon_M_att1,
                        ljInfoVecs.epsilon_M_att2,
                        ljInfoVecs.epsilon_M_rep1,
                        ljInfoVecs.epsilon_M_rep2,
                        thrust::raw_pointer_cast(generalParams.nodes_in_upperhem.data()),
                        ljInfoVecs.LJ_PosX,
                        ljInfoVecs.LJ_PosY,
                        ljInfoVecs.LJ_PosZ,
                        thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
                        thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
                        thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data()),
                        thrust::raw_pointer_cast(coordInfoVecs.nodeForceX.data()),
                        thrust::raw_pointer_cast(coordInfoVecs.nodeForceY.data()),
                        thrust::raw_pointer_cast(coordInfoVecs.nodeForceZ.data())),
                    init, CVec4Add() );//unary, binary 
     
    //temp stores forces on lj particle and last entry as energy
    ljInfoVecs.forceX = thrust::get<0>(temp);
    ljInfoVecs.forceY = thrust::get<1>(temp);
    ljInfoVecs.forceZ = thrust::get<2>(temp);

    ljInfoVecs.lj_energy_M = thrust::get<3>(temp);
      
    //std::cout<<"lj force from membrane "<< ljInfoVecs.forceX<< " "<<  ljInfoVecs.forceY << " "<<  ljInfoVecs.forceZ << std::endl;
    //std::cout<<"lj energy from membrane "<<ljInfoVecs.lj_energy_M<<std::endl;

};

void AdvanceLJParticle(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LJInfoVecs& ljInfoVecs) {
    
    
    ljInfoVecs.LJ_PosX = ljInfoVecs.LJ_PosX + generalParams.dt * ljInfoVecs.forceX;
    ljInfoVecs.LJ_PosY = ljInfoVecs.LJ_PosY + generalParams.dt * ljInfoVecs.forceY;
    ljInfoVecs.LJ_PosZ = ljInfoVecs.LJ_PosZ + generalParams.dt * ljInfoVecs.forceZ;
};
