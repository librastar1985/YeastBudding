
#include "System.h"
#include "LJEnergy.h"



//WARNING: function must not reset coordInfoVecs.nodeForceX etc. 
double ComputeLJEnergy(
    CoordInfoVecs& coordInfoVecs,
    LJInfoVecs& ljInfoVecs,  
    GeneralParams& generalParams) {    
     
    thrust::counting_iterator<int> begin(0);
    thrust::counting_iterator<int> end(generalParams.maxNodeCount);
 
    ljInfoVecs.lj_energy = thrust::transform_reduce(  
                    begin,end,
                    LJEnergyFunctor(
                        ljInfoVecs.Rcutoff, 
                        ljInfoVecs.Rmin,
                        ljInfoVecs.epsilon,
                        ljInfoVecs.LJ_PosX,
                        ljInfoVecs.LJ_PosY,
                        ljInfoVecs.LJ_PosZ,
                        thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
                        thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
                        thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data())),
                    0.0, thrust::plus<double>() );//unary, binary 
     
    //temp stores forces on lj particle and last entry as energy
    //ljInfoVecs.forceX = thrust::get<0>(temp);
    //ljInfoVecs.forceY = thrust::get<1>(temp);
    //ljInfoVecs.forceZ = thrust::get<2>(temp);

    return ljInfoVecs.lj_energy;

};
