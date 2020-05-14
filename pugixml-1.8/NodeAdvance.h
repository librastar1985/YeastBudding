#ifndef NODEADVANCE_H_
#define NODEADVANCE_H_

#include "SystemStructures.h"

void AdvancePositions(
	CoordInfoVecs& coordInfoVecs,
	GeneralParams& generalParams,
	DomainParams& domainParams);




//advances  polymer positions and  using polymer timestep
struct SaxpyFunctorPrimary : public thrust::binary_function<UCVec3, CVec3, CVec3> {
	double dt;
	double mass;
	unsigned maxNodeCount;
	double domainLengthX;
	double domainLengthY;
	double domainLengthZ;

	__host__ __device__
		//
		SaxpyFunctorPrimary(
			double& _dt, 
			double& _mass,
			unsigned& _maxNodeCount,
			double& _domainLengthX,
			double& _domainLengthY,
			double& _domainLengthZ) :
		dt(_dt),
		mass(_mass),
		maxNodeCount(_maxNodeCount),
		domainLengthX(_domainLengthX),
		domainLengthY(_domainLengthY),
		domainLengthZ(_domainLengthZ) {}

	__device__
		CVec3 operator()(const UCVec3 &p3, const CVec3 &f3) {

		bool isFixed = thrust::get<0>(p3);//true if fixed, false if movable. 
//not using fixed nodes for now 
		if (true) {
			double xLocNew = thrust::get<1>(p3) + dt/mass * thrust::get<0>(f3);
			double yLocNew = thrust::get<2>(p3) + dt/mass * thrust::get<1>(f3);
			double zLocNew = thrust::get<3>(p3) + dt/mass * thrust::get<2>(f3);


			return thrust::make_tuple(xLocNew, yLocNew, zLocNew);
		}
		else {
			return thrust::make_tuple(thrust::get<1>(p3),thrust::get<2>(p3),thrust::get<3>(p3) );
		}
	}                                 

};

#endif /*NODEADVANCE_H_*/
