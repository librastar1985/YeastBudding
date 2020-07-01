#ifndef BUCKETSCHEME_H_
#define BUCKETSCHEME_H_

#include "SystemStructures.h"


void initDimensionBucketScheme(
	CoordInfoVecs& coordInfoVecs,
	CapsidInfoVecs& capsidInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams);
	
void buildBucketScheme(
	CoordInfoVecs& coordInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams,
	CapsidInfoVecs& capsidInfoVecs);
	
void extendBucketScheme(
	CoordInfoVecs& coordInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs);





/*
* Functor to compute neighbor buckets(center bucket included) of a node.
* @param input1 bucket index of node
* @param input2 pick from the sequence, which is also global rank of the node
*
* @return output1 bucket indices of node ( all neighbors and the center bucket) of node
* @return output2 global rank of node

example with 3x3x3 grid. with node at position (1,1,1). Notice there are 27 possibilities.
0: (1,1,1)
1: (0,0,1)
2: (1,0,1)
The words left & right denote x change, top and bottom denote y change and upper & lower denote z change
 */

struct NeighborFunctor : public thrust::unary_function<Tuu, Tuu> {
	int numOfBucketsInXDim;
	int numOfBucketsInYDim;
	int numOfBucketsInZDim;

	__host__ __device__ NeighborFunctor(
		int _numOfBucketsInXDim,
		int _numOfBucketsInYDim,
		int _numOfBucketsInZDim ) :
		numOfBucketsInXDim(_numOfBucketsInXDim),
		numOfBucketsInYDim(_numOfBucketsInYDim),
		numOfBucketsInZDim(_numOfBucketsInZDim) {}

	__device__ int operator()(const Tuu &v) {
		int relativeRank = thrust::get<1>(v) % 27;	//27 = 3^3. Takes global node id for calculation

		//takes global bucket id for calculation
		if ((thrust::get<1>(v) == (27*479+1))) {
			//coment

		}
		int area = numOfBucketsInXDim * numOfBucketsInYDim;
		__attribute__ ((unused)) int volume = area * numOfBucketsInZDim;

		int xPos = thrust::get<0>(v) % numOfBucketsInXDim;	//col
		int xPosLeft = xPos - 1;
		int xPosRight = xPos + 1;
		if (xPos == 0) {
			//wraparound int
			xPosLeft = numOfBucketsInXDim-1;
		}
		if (xPosRight >= numOfBucketsInXDim) {
			xPosRight = 0;
		}


		int zPos = thrust::get<0>(v) / area; //z divide by area
		int zPosUp = zPos + 1;
		int zPosLow = zPos - 1;
		if (zPos == 0 ) {
			//wraparound int
			zPosLow = numOfBucketsInZDim-1;
		}
		if (zPosUp >= numOfBucketsInZDim) {
			zPosUp = 0;
		}

		int yPos = (thrust::get<0>(v) - zPos * area) / numOfBucketsInXDim;	//row
		int yPosTop = yPos + 1;
		int yPosBottom = yPos - 1;

		if (yPos == 0) {
			//wraparound unsigend
			yPosBottom = numOfBucketsInYDim-1;
		}
		if (yPosTop >= numOfBucketsInYDim) {
			yPosTop = 0;
		}

		switch (relativeRank) {
		//middle cases
		case 0:
			return thrust::get<0>(v);
			//break;
		case 1:{
				int topLeft = xPosLeft + yPosTop * numOfBucketsInXDim + zPos * area;
				return (topLeft);
				//break;
		}
		case 2:{
				int top = xPos + yPosTop * numOfBucketsInXDim + zPos * area;
			return (top);
			//break;
		}
		case 3:{
				int topRight = xPosRight + yPosTop * numOfBucketsInXDim + zPos * area;
			return topRight;
			//break;
		}
		case 4:{
				int right = xPosRight + yPos * numOfBucketsInXDim + zPos * area;
			return right;
			//break;
		}
		case 5:{
				int bottomRight = xPosRight + yPosBottom * numOfBucketsInXDim + zPos * area;
			return bottomRight;
			//break;
		}
		case 6:{
				int bottom = xPos + yPosBottom * numOfBucketsInXDim + zPos * area;
			return bottom;
			//break;
		}
		case 7:{
				int bottomLeft = xPosLeft + yPosBottom * numOfBucketsInXDim + zPos * area;
			return bottomLeft;
			//break;
		}
		case 8:{
				int left = xPosLeft + yPos * numOfBucketsInXDim + zPos * area;
			return left;
			//break;
		}
		//lower Z cases
		case 9:{
				int lowerCenter = xPos + yPos * numOfBucketsInXDim +  zPosLow * area;
			return lowerCenter;
			//break;
		}
		case 10:{
				int lowerTopLeft = xPosLeft + yPosTop * numOfBucketsInXDim + zPosLow* area;
			return lowerTopLeft;
			//break;
		}
		case 11:{
				int lowerTop = xPos + yPosTop * numOfBucketsInXDim + zPosLow * area;
			return (lowerTop);
			//break;
		}
		case 12:{
				int lowerTopRight = xPosRight + yPosTop * numOfBucketsInXDim  + zPosLow * area;
			return lowerTopRight;
			//break;
		}
		case 13:{
				int lowerRight = xPosRight + yPos * numOfBucketsInXDim + zPosLow * area;
			return (lowerRight);
			//break;
		}
		case 14:{
				int lowerBottomRight = xPosRight + yPosBottom * numOfBucketsInXDim + zPosLow * area;
			return (lowerBottomRight);
			//break;
		}
		case 15:{
				int lowerBottom = xPos + yPosBottom * numOfBucketsInXDim + zPosLow * area;
			return lowerBottom;
			//break;
		}
		case 16:{
				int lowerBottomLeft = xPosLeft + yPosBottom * numOfBucketsInXDim  + zPosLow * area;
			return lowerBottomLeft;
			//break;
		}
		case 17:{
				int lowerLeft = xPosLeft + yPos * numOfBucketsInXDim + zPosLow * area;
			return lowerLeft;
			//break;
		}
		//upper Z cases
		case 18:{
				int upperCenter = xPos + yPos * numOfBucketsInXDim +  zPosUp * area;
			return (upperCenter);
			//break;
		}
		case 19:{
				int upperTopLeft = xPosLeft + yPosTop * numOfBucketsInXDim + zPosUp * area;
			return (upperTopLeft);
			//break;
		}
		case 20:{
				int upperTop = xPos + yPosTop * numOfBucketsInXDim + zPosUp * area;
			return (upperTop);
			//break;
		}
		case 21:{
				int upperTopRight = xPosRight + yPosTop * numOfBucketsInXDim  + zPosUp * area;
			return (upperTopRight);
			//break;
		}
		case 22:{
				int upperRight = xPosRight + yPos * numOfBucketsInXDim + zPosUp * area;
			return (upperRight);
			//break;
		}
		case 23:{
				int upperBottomRight = xPosRight + yPosBottom * numOfBucketsInXDim + zPosUp * area;
			return (upperBottomRight);
			//break;
		}
		case 24:{
				int upperBottom = xPos + yPosBottom * numOfBucketsInXDim + zPosUp * area;
			return (upperBottom);
			//break;
		}
		case 25:{
				int upperBottomLeft = xPosLeft + yPosBottom * numOfBucketsInXDim  + zPosUp * area;
			return (upperBottomLeft);
			//break;
		}
		case 26:{
				int upperLeft = xPosLeft + yPos * numOfBucketsInXDim + zPosUp * area;
			return (upperLeft);
			//break;
		}
		default:{
			int default_Id=ULONG_MAX;
			return (default_Id);
		}

		}
	}
};


struct BucketIndexer {
	double minX;
	double maxX;
	double minY;
	double maxY;
	double minZ;
	double maxZ;

	int XBucketCount;
	int YBucketCount;
	int ZBucketCount;
	double unitLen;

	__host__ __device__

	BucketIndexer(
		double _minX,
		double _maxX,
		double _minY,
		double _maxY,
		double _minZ,
		double _maxZ,
		int _XBucketCount,
		int _YBucketCount,
		int _ZBucketCount,
		double _unitLen) :

		minX(_minX),
		maxX(_maxX),
		minY(_minY),
		maxY(_maxY),
		minZ(_minZ),
		maxZ(_maxZ),
		XBucketCount(_XBucketCount),
		YBucketCount(_YBucketCount),
		ZBucketCount(_ZBucketCount),
		unitLen(_unitLen) {}

	__device__ 
	Tuu operator()(const Tdddu& v) {

			int id = thrust::get<3>(v);

			int x = static_cast<int>((thrust::get<0>(v) - minX) / unitLen);
			int y = static_cast<int>((thrust::get<1>(v) - minY) / unitLen);
			int z = static_cast<int>((thrust::get<2>(v) - minZ) / unitLen);


			// return the bucket's linear index and node's global index
			//return thrust::make_tuple(z * XSize * YSize + y * XSize + x, thrust::get<4>(v));
			int bucket = z * XBucketCount * YBucketCount + y * XBucketCount + x;
			//try to make it so bucket does not return int32Max
			if (bucket == ULONG_MAX) {
				bucket = 0;
			}
			return thrust::make_tuple(bucket, id);

	}
};

template <typename InputIterator1,
          typename InputIterator2,
          typename OutputIterator>
OutputIterator expand(InputIterator1 first1,
                      InputIterator1 last1,
                      InputIterator2 first2,
                      OutputIterator output) {
  typedef typename thrust::iterator_difference<InputIterator1>::type difference_type;
  
  difference_type input_size  = thrust::distance(first1, last1);
  difference_type output_size = thrust::reduce(first1, last1);

  // scan the counts to obtain output offsets for each input element
  thrust::device_vector<difference_type> output_offsets(input_size, 0);
  thrust::exclusive_scan(first1, last1, output_offsets.begin()); 

  // scatter the nonzero counts into their corresponding output positions
  thrust::device_vector<difference_type> output_indices(output_size, 0);
  thrust::scatter_if
    (thrust::counting_iterator<difference_type>(0),
     thrust::counting_iterator<difference_type>(input_size),
     output_offsets.begin(),
     first1,
     output_indices.begin());

  // compute max-scan over the output indices, filling in the holes
  thrust::inclusive_scan
    (output_indices.begin(),
     output_indices.end(),
     output_indices.begin(),
     thrust::maximum<difference_type>());

  // gather input values according to index array (output = first2[output_indices])
  OutputIterator output_end = output; thrust::advance(output_end, output_size);
  thrust::gather(output_indices.begin(),
                 output_indices.end(),
                 first2,
                 output);

  // return output + output_size
  thrust::advance(output, output_size);
  return output;
};





#endif