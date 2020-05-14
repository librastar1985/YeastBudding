#ifndef SYSTEMSTRUCTURES_H_
#define SYSTEMSTRUCTURES_H_

#include <memory>
#include <cmath>
#include <random>
#include <algorithm>
#include <iterator>
#include <thrust/transform_reduce.h>
#include <thrust/functional.h>

#include <thrust/random.h>
#include <thrust/extrema.h>
#include <thrust/sort.h>
#include <thrust/for_each.h>
#include <thrust/transform.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/tuple.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/remove.h>
#include <thrust/copy.h>
#include <thrust/execution_policy.h>
#include <thrust/pair.h>
#include <thrust/unique.h>
#include <thrust/remove.h>
#include <thrust/binary_search.h>
#include <thrust/reduce.h>
#include <thrust/replace.h>
#include <thrust/gather.h>
#include <thrust/random/uniform_real_distribution.h>
#include <stdint.h>
#include <thrust/sequence.h>


typedef thrust::tuple<int, bool, double> Tubd;
typedef thrust::tuple<int, bool> Tub;
typedef thrust::tuple<int, double> Tud;
typedef thrust::tuple<bool, double> Tbd;


typedef thrust::tuple<int, int, double> Tuud;

typedef thrust::tuple<int, int, int, int, int, int, int> Tuuuuuuu;
typedef thrust::tuple<int, int, int, int, double> Tuuuud;
typedef thrust::tuple<int, int, int,int, int> Tuuuuu;
typedef thrust::tuple<int, int, int,int> Tuuuu;
typedef thrust::tuple<int, int, int, double> Tuuud;
typedef thrust::tuple<int, int, int> Tuuu;
typedef thrust::tuple<int, int> Tuu;


typedef thrust::tuple<int, double, double, double> Tuddd;
typedef thrust::tuple<double, double, double, int> Tdddu;
typedef thrust::tuple<double, double> Tdd;

typedef thrust::tuple<bool, double, double, double, double, double, double> BoolCVec6;
typedef thrust::tuple<int, double, double, double, double, double, double,double, double, double> UCVec9;

typedef thrust::tuple<int, double, double, double, double, double, double> UCVec6;
typedef thrust::tuple<int, int, double, double, double, double> U2CVec4;
typedef thrust::tuple<int, int, double, double, double> U2CVec3;
typedef thrust::tuple<int, double, double, double> UCVec3;
typedef thrust::tuple<bool, double, double, double> BoolCVec3;


typedef thrust::tuple<double, double, double, double, double, double, double, double, double> CVec9;
typedef thrust::tuple<double, double, double, double, double, double, double, double> CVec8;
typedef thrust::tuple<double, double, double, double, double, double, double> CVec7;
typedef thrust::tuple<double, double, double, double, double, double> CVec6;
typedef thrust::tuple<double, double, double, double, double> CVec5;
typedef thrust::tuple<double, double, double, double> CVec4;
typedef thrust::tuple<double, double, double> CVec3;
typedef thrust::tuple<double, double> CVec2;

typedef thrust::tuple<CVec3,CVec3,CVec3> Mat_3x3;

struct GeneralParams;
struct DomainParams;
struct AuxVecs;
struct CoordInfoVecs;
struct CapsidInfoVecs;
struct BendingTriangleInfoVecs;
struct AreaTriangleInfoVecs;
struct LinearSpringInfoVecs;
struct LJInfoVecs;


struct HostSetInfoVecs {

	thrust::host_vector<double> scaling_per_edge;
	thrust::host_vector<int> nodes_in_upperhem;
	thrust::host_vector<int> triangles_in_upperhem;
	thrust::host_vector<int> edges_in_upperhem;
	thrust::host_vector<int> edges_in_upperhem_list;
	thrust::host_vector<int> boundaries_in_upperhem;

	thrust::host_vector<int> nodes2Triangles_1;
	thrust::host_vector<int> nodes2Triangles_2;
	thrust::host_vector<int> nodes2Triangles_3;
	thrust::host_vector<int> nodes2Triangles_4;
	thrust::host_vector<int> nodes2Triangles_5;
	thrust::host_vector<int> nodes2Triangles_6;
	thrust::host_vector<int> nodes2Triangles_7;
	thrust::host_vector<int> nodes2Triangles_8;
	thrust::host_vector<int> nodes2Triangles_9;
	//thrust::host_vector<int> nodes2Triangles_10;
	//thrust::host_vector<int> nodes2Triangles_11;
	//thrust::host_vector<int> nodes2Triangles_12;

	thrust::host_vector<int> nndata1;
	thrust::host_vector<int> nndata2;
	thrust::host_vector<int> nndata3;
	thrust::host_vector<int> nndata4;
	thrust::host_vector<int> nndata5;
	thrust::host_vector<int> nndata6;
	thrust::host_vector<int> nndata7;
	thrust::host_vector<int> nndata8;
	thrust::host_vector<int> nndata9;
	//thrust::host_vector<int> nndata10;
	//thrust::host_vector<int> nndata11;
	//thrust::host_vector<int> nndata12;
	
	thrust::host_vector<double> capsidNodeLocX;
	thrust::host_vector<double> capsidNodeLocY;
	thrust::host_vector<double> capsidNodeLocZ;
	
	thrust::host_vector<bool> isNodeFixed;
	
 	thrust::host_vector<double> nodeLocX;
	thrust::host_vector<double> nodeLocY;
	thrust::host_vector<double> nodeLocZ;
	
	thrust::host_vector<double> nodeForceX;
	thrust::host_vector<double> nodeForceY;
	thrust::host_vector<double> nodeForceZ;
	
	//LOCAL COORDS
	//indices of each triangle
	thrust::host_vector<int> triangles2Nodes_1;
	thrust::host_vector<int> triangles2Nodes_2;
	thrust::host_vector<int> triangles2Nodes_3;

	
	//indices of each edge
	thrust::host_vector<int> edges2Nodes_1;
	thrust::host_vector<int> edges2Nodes_2;

	//indices of 2 triangle on each edge
	thrust::host_vector<int> edges2Triangles_1;
	thrust::host_vector<int> edges2Triangles_2;

	//indices of edges on each triangle.
	thrust::host_vector<int> triangles2Edges_1;
	thrust::host_vector<int> triangles2Edges_2;
	thrust::host_vector<int> triangles2Edges_3;

	thrust::host_vector<double> edge_initial_length;
};

struct AddForceFunctor {
	double* forceXAddr;
	double* forceYAddr;
	double* forceZAddr;

	__host__ __device__
	//
		AddForceFunctor(
				double* _forceXAddr,
				double* _forceYAddr,
				double* _forceZAddr) :
			forceXAddr(_forceXAddr),
			forceYAddr(_forceYAddr),
			forceZAddr(_forceZAddr) {}

	__device__
	void operator() (const Tuddd& u1d3) {
			int idToAssign = thrust::get<0>(u1d3);
			if (!isnan(thrust::get<1>(u1d3)) && !isnan(thrust::get<2>(u1d3)) && !isnan(thrust::get<3>(u1d3))) {

			forceXAddr[idToAssign] += thrust::get<1>(u1d3);
			forceYAddr[idToAssign] += thrust::get<2>(u1d3);
			forceZAddr[idToAssign] += thrust::get<3>(u1d3);
			}

	}

};

struct AddForceFunctorAlt {
	int maxNodeCount;
	double* forceXAddr;
	double* forceYAddr;
	double* forceZAddr;

	__host__ __device__
	//
		AddForceFunctorAlt(
				int& _maxNodeCount,
				double* _forceXAddr,
				double* _forceYAddr,
				double* _forceZAddr) :
			maxNodeCount(_maxNodeCount),
			forceXAddr(_forceXAddr),
			forceYAddr(_forceYAddr),
			forceZAddr(_forceZAddr) {}

	__device__
	void operator() (const Tuddd& u1d3) {
		int idToAssign = thrust::get<0>(u1d3);
		if (idToAssign < maxNodeCount){
			forceXAddr[idToAssign] += thrust::get<1>(u1d3);
			forceYAddr[idToAssign] += thrust::get<2>(u1d3);
			forceZAddr[idToAssign] += thrust::get<3>(u1d3);
		}

	}

};
  
struct UCVec3Add : public thrust::binary_function<UCVec3, UCVec3, UCVec3> {
	__host__ __device__ 
		UCVec3 operator()(const UCVec3 &vec1, const UCVec3 &vec2) {
		return thrust::make_tuple(
            thrust::get<0>(vec1) + thrust::get<0>(vec2),
            thrust::get<1>(vec1) + thrust::get<1>(vec2),
			thrust::get<2>(vec1) + thrust::get<2>(vec2),
			thrust::get<3>(vec1) + thrust::get<3>(vec2));
	}
};
struct CVec3Add : public thrust::binary_function<CVec3, CVec3, CVec3> {
	__host__ __device__ 
		CVec3 operator()(const CVec3 &vec1, const CVec3 &vec2) {
		return thrust::make_tuple(
            thrust::get<0>(vec1) + thrust::get<0>(vec2),
			thrust::get<1>(vec1) + thrust::get<1>(vec2),
			thrust::get<2>(vec1) + thrust::get<2>(vec2));
	}
};

struct CVec4Add : public thrust::binary_function<CVec4, CVec4, CVec4> {
	__host__ __device__ 
		CVec4 operator()(const CVec4 &vec1, const CVec4 &vec2) {
		return thrust::make_tuple(
            thrust::get<0>(vec1) + thrust::get<0>(vec2),
			thrust::get<1>(vec1) + thrust::get<1>(vec2),
			thrust::get<2>(vec1) + thrust::get<2>(vec2),
			thrust::get<3>(vec1) + thrust::get<3>(vec2));
	}
};

struct CVec3NormBinary {
	__host__ __device__

		double operator() (const CVec3& vec1, const CVec3& vec2) {
		//divide force by fiber cross section to get stress
		return fabs(
			((thrust::get<0>(vec1) - thrust::get<0>(vec2))) +
			((thrust::get<1>(vec1) - thrust::get<1>(vec2))) +
			((thrust::get<2>(vec1) - thrust::get<2>(vec2))));
	}
};

struct CVec3NormUnary {
	__host__ __device__
		double operator() (const CVec3& vec) {
		//divide force by fiber cross section to get stress
		return (
			sqrt(thrust::get<0>(vec) * thrust::get<0>(vec) +
			thrust::get<1>(vec) * thrust::get<1>(vec) +
			thrust::get<2>(vec) * thrust::get<2>(vec)));
	}
};


///////////////////////////////////////////////
//random number generators
struct psrnormgen {
    
    double a, b;

    __host__ __device__ 
	psrnormgen(
		double _a, 
		double _b) : 
		a(_a), 
		b(_b) {}
 
    __device__ double operator()(const int n) const
    {
        thrust::default_random_engine rng(n);
        thrust::normal_distribution<double> dist(a, b);
        rng.discard(n);
        return dist(rng);
    }
 
};
struct psrunifgen {
    
    double a, b;

    __host__ __device__ 
	psrunifgen(
		double _a, 
		double _b) : 
		a(_a), 
		b(_b) {}
 
    __device__ double operator()(const int n) const
    {
        thrust::default_random_engine rng(n);
        thrust::uniform_real_distribution<double> dist(a, b);
        rng.discard(n);
        return dist(rng);
    }
 
};


struct TorsionAngleFunctor {
	double* locXAddr;
	double* locYAddr;
	double* locZAddr;
	__host__ __device__
		TorsionAngleFunctor(
			double* _locXAddr,
			double* _locYAddr,
			double* _locZAddr) :
			locXAddr(_locXAddr),
			locYAddr(_locYAddr),
			locZAddr(_locZAddr) {}

	__device__
	double operator() (const Tuuu &u3) {
		int indexLeft = thrust::get<0>(u3);
		int indexCenter = thrust::get<1>(u3);
		int indexRight = thrust::get<2>(u3);

		double distLCX = locXAddr[indexLeft] - locXAddr[indexCenter];
		double distLCY = locYAddr[indexLeft] - locYAddr[indexCenter];
		double distLCZ = locZAddr[indexLeft] - locZAddr[indexCenter];
		double distRCX = locXAddr[indexRight] - locXAddr[indexCenter];
		double distRCY = locYAddr[indexRight] - locYAddr[indexCenter];
		double distRCZ = locZAddr[indexRight] - locZAddr[indexCenter];

		//lengths between left & center, right & center
		double lenLC = sqrt(distLCX*distLCX + distLCY*distLCY + distLCZ*distLCZ);
		double lenRC = sqrt(distRCX*distRCX + distRCY*distRCY + distRCZ*distRCZ);


		//normalized dot product
		double cosTheta = (distLCX*distRCX + distLCY*distRCY + distLCZ*distRCZ) / (lenLC * lenRC);
	
		//rounding errors
		if (cosTheta < -1.0) 
			cosTheta = -1.0;
		else if (cosTheta > 1.0) 
			cosTheta = 1.0;

		//double currentAngle = acos(cosTheta);

		return acos(cosTheta);
	}
};



//used to calculate strain
struct AveStrainFunctor {
  __host__ __device__ 
  double operator() (const Tbd &b1d1) {
	  bool isStrainNode = thrust::get<0>(b1d1);
	  if (isStrainNode) 
		  return thrust::get<1>(b1d1);
	  else 
    	return 0.0; 
  }
};

struct NormFunctor {
	
		__host__ __device__
		double operator() (const CVec3& vec) {
		//divide force by fiber cross section to get stress
		double result = sqrt(thrust::get<0>(vec) * thrust::get<0>(vec) +
			thrust::get<1>(vec) * thrust::get<1>(vec) +
			thrust::get<2>(vec) * thrust::get<2>(vec));
		return result;

	}
};

//returns true if is greater than level
struct IsGreaterThanLevel {
	double limit;

	__host__ __device__ 
		IsGreaterThanLevel(
			double& _limit) : 
			limit(_limit) {}
		
	__device__
	//replaces value with 1 if returns true 
		bool operator() (double zPos) {
			return (zPos > limit);//((1-percentPull) * networkLength));
		}
}; 

//returns true if less than llevel.
struct IsLessThanLevel {
		double limit;

	__host__ __device__ 
		IsLessThanLevel(
			double& _limit) : 
			limit(_limit) {}
		
	__device__
	//replaces value with 1 if returns true 
		bool operator() (double zPos) {
			return (zPos < limit);//((1-percentPull) * networkLength));
		}
};

struct tupleEqual {
	
  __host__ __device__
    bool operator()(Tuu x, Tuu y)
    {
      return ( (x.get<0>() == y.get<0>()) && (x.get<1>() == y.get<1>()) );
    }
};

//return true if not equal to input value.
struct IsEqualToOne {

	__host__ __device__ 
	bool operator() (const int& x) {
		return (x != 1);
	}
};

//return true if not equal to input value.
struct isNotEqualZero {

	__host__ __device__ 
	bool operator() (const int& x) {
		return (x != 0);
	}
};

//return true if equal to input value.
struct isEqualZero {
	__host__ __device__ 
	bool operator() (const int& x) {
		return (x == 0);
	}
};

struct is_greater_than {
	int limit;

	 __host__ __device__
	is_greater_than(int& _limit) : limit(_limit) {}
  __device__
  bool operator()(const int& x) {
	if ( x > limit ) {
    	return true;
	}
	else {
		return false;
	}
  }
};
struct is_less_than {
	int limit;

	 __host__ __device__
	is_less_than(int& _limit) : limit(_limit) {}
  __device__
  bool operator()(const int& x) {
	if ( x < limit ) {
    	return true;
	}
	else {
		return false;
	}
  }
};

struct CVec3InnerProduct {
    __host__ __device__ 
    double operator() (CVec3& v1 ) {
        return (thrust::get<0>(v1) * thrust::get<0>(v1) + 
                thrust::get<1>(v1) * thrust::get<1>(v1) + 
                thrust::get<2>(v1) * thrust::get<2>(v1) );
    }
};

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//FUNCTIONS USED IN AREA TRIANGLES AND BENDING TRIANGLES. 

__host__ __device__
inline CVec3 CVec3_cross(CVec3 v1, CVec3 v2) {
	return thrust::make_tuple(
		thrust::get<1>(v1)*thrust::get<2>(v2) - thrust::get<2>(v1)*thrust::get<1>(v2),
		-(thrust::get<0>(v1)*thrust::get<2>(v2) - thrust::get<2>(v1)*thrust::get<0>(v2)),
		thrust::get<0>(v1)*thrust::get<1>(v2) - thrust::get<1>(v1)*thrust::get<0>(v2));
};

__host__ __device__
inline CVec3 CVec3_mult (CVec3 v1, CVec3 v2) {
	return thrust::make_tuple(thrust::get<0>(v1)*thrust::get<0>(v2),
		thrust::get<1>(v1)*thrust::get<1>(v2),
		thrust::get<2>(v1)*thrust::get<2>(v2));
};

__host__ __device__
inline CVec3 CVec3_scalermult (double c, CVec3 v2) {
	return thrust::make_tuple(c*thrust::get<0>(v2),
		c*thrust::get<1>(v2),
		c*thrust::get<2>(v2));
};

__host__ __device__
inline CVec3 CVec3_div (CVec3 v1, CVec3 v2) {
	return thrust::make_tuple(thrust::get<0>(v1) / thrust::get<0>(v2),
		thrust::get<1>(v1) / thrust::get<1>(v2),
		thrust::get<2>(v1) / thrust::get<2>(v2));
};

__host__ __device__
inline double CVec3_dot(CVec3 v1, CVec3 v2) {
	return (
		thrust::get<0>(v1)*thrust::get<0>(v2) +
		thrust::get<1>(v1)*thrust::get<1>(v2) +
		thrust::get<2>(v1)*thrust::get<2>(v2));
};

__host__ __device__
inline CVec3 CVec3_plus(CVec3 v1, CVec3 v2) {
	return thrust::make_tuple(
		thrust::get<0>(v1) + thrust::get<0>(v2),
		thrust::get<1>(v1) + thrust::get<1>(v2),
		thrust::get<2>(v1) + thrust::get<2>(v2));
};


__host__ __device__
inline CVec3 CVec3_plus(CVec3 v1, double d) {
	return thrust::make_tuple(
		thrust::get<0>(v1) + d,
		thrust::get<1>(v1) + d,
		thrust::get<2>(v1) + d);
};

__host__ __device__
inline CVec3 CVec3_plus(CVec3 v1, CVec3 v2, CVec3 v3) {
	return thrust::make_tuple(
		thrust::get<0>(v1) + thrust::get<0>(v2) + thrust::get<0>(v3),
		thrust::get<1>(v1) + thrust::get<1>(v2) +  thrust::get<1>(v3),
		thrust::get<2>(v1) + thrust::get<2>(v2) +  thrust::get<2>(v3));
};

__host__ __device__
inline CVec3 CVec3_plus(CVec3 v1, CVec3 v2, CVec3 v3, CVec3 v4, CVec3 v5, CVec3 v6) {
	return thrust::make_tuple(
		thrust::get<0>(v1) + thrust::get<0>(v2) + thrust::get<0>(v3) + thrust::get<0>(v4) + thrust::get<0>(v5) + thrust::get<0>(v6),
		thrust::get<1>(v1) + thrust::get<1>(v2) + thrust::get<1>(v3) + thrust::get<1>(v4) + thrust::get<1>(v5) + thrust::get<1>(v6),		
		thrust::get<2>(v1) + thrust::get<2>(v2) + thrust::get<2>(v3) + thrust::get<2>(v4) + thrust::get<2>(v5) + thrust::get<2>(v6));
};

__host__ __device__
inline CVec3 CVec3_subtract(CVec3 v1, CVec3 v2) {
	return thrust::make_tuple(thrust::get<0>(v1) - thrust::get<0>(v2),
		thrust::get<1>(v1) - thrust::get<1>(v2),
		thrust::get<2>(v1) - thrust::get<2>(v2));
};

__host__ __device__
inline CVec3 CVec3_minus(CVec3 v1, CVec3 v2) {
	return thrust::make_tuple(
		thrust::get<0>(v1) - thrust::get<0>(v2),
		thrust::get<1>(v1) - thrust::get<1>(v2),
		thrust::get<2>(v1) - thrust::get<2>(v2));
};

#endif /* SYSTEMSTRUCTURES_H_*/