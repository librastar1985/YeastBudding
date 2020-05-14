#include "BucketScheme.h"
#include "System.h"


void initDimensionBucketScheme(
	CoordInfoVecs& coordInfoVecs,
	CapsidInfoVecs& capsidInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams) {
	
	double membrane_minX = (*(thrust::min_element(coordInfoVecs.nodeLocX.begin(), coordInfoVecs.nodeLocX.end())));
	double membrane_maxX = (*(thrust::max_element(coordInfoVecs.nodeLocX.begin(), coordInfoVecs.nodeLocX.end())));
	double membrane_minY = (*(thrust::min_element(coordInfoVecs.nodeLocY.begin(), coordInfoVecs.nodeLocY.end())));
	double membrane_maxY = (*(thrust::max_element(coordInfoVecs.nodeLocY.begin(), coordInfoVecs.nodeLocY.end())));
	double membrane_minZ = (*(thrust::min_element(coordInfoVecs.nodeLocZ.begin(), coordInfoVecs.nodeLocZ.end())));
	double membrane_maxZ = (*(thrust::max_element(coordInfoVecs.nodeLocZ.begin(), coordInfoVecs.nodeLocZ.end())));
	
	double capsid_minX = 0.0;
	double capsid_maxX = 0.0;
	double capsid_minY = 0.0;
	double capsid_maxY = 0.0;
	double capsid_minZ = 0.0;
	double capsid_maxZ = 0.0;

	if (capsidInfoVecs.nodeLocX.size() > 0){
		capsid_minX = (*(thrust::min_element(capsidInfoVecs.nodeLocX.begin(), capsidInfoVecs.nodeLocX.end())));
		capsid_maxX = (*(thrust::max_element(capsidInfoVecs.nodeLocX.begin(), capsidInfoVecs.nodeLocX.end())));
		capsid_minY = (*(thrust::min_element(capsidInfoVecs.nodeLocY.begin(), capsidInfoVecs.nodeLocY.end())));
		capsid_maxY = (*(thrust::max_element(capsidInfoVecs.nodeLocY.begin(), capsidInfoVecs.nodeLocY.end())));
		capsid_minZ = (*(thrust::min_element(capsidInfoVecs.nodeLocZ.begin(), capsidInfoVecs.nodeLocZ.end())));
		capsid_maxZ = (*(thrust::max_element(capsidInfoVecs.nodeLocZ.begin(), capsidInfoVecs.nodeLocZ.end())));
	}

	domainParams.minX = min(capsid_minX, membrane_minX);
	domainParams.maxX = max(capsid_maxX, membrane_maxX);
	domainParams.minY = min(capsid_minY, membrane_minY);
	domainParams.maxY = max(capsid_maxY, membrane_maxY);
	domainParams.minZ = min(capsid_minZ, membrane_minZ);
	domainParams.maxZ = max(capsid_maxZ, membrane_maxZ);

	if (generalParams.iteration == 0) {
		domainParams.originMinX = domainParams.minX;
		domainParams.originMaxX = domainParams.maxX;		
		domainParams.originMinY = domainParams.minY;
		domainParams.originMaxY = domainParams.maxY;
		domainParams.originMinZ = domainParams.minZ;
		domainParams.originMaxZ = domainParams.maxZ;
		std::cout<<"BSOD Zminmax: "<<domainParams.minZ <<" "<< domainParams.maxZ	<< std::endl;
		std::cout<<"BSOD Yminmax: "<<domainParams.minY <<" "<< domainParams.maxY	<< std::endl;
 		std::cout<<"BSOD Xminmax: "<<domainParams.minX <<" "<< domainParams.maxX	<< std::endl;
 
	}
	
	int temp_bucket_count = domainParams.XBucketCount * domainParams.YBucketCount * domainParams.ZBucketCount;

 	if (temp_bucket_count != domainParams.totalBucketCount) {
		std::cout<<"grid: "<< domainParams.gridSpacing <<std::endl;
		domainParams.XBucketCount = ceil((domainParams.maxX - domainParams.minX) / domainParams.gridSpacing + 1.0);
		domainParams.YBucketCount = ceil((domainParams.maxY - domainParams.minY) / domainParams.gridSpacing + 1.0);
		domainParams.ZBucketCount = ceil((domainParams.maxZ - domainParams.minZ) / domainParams.gridSpacing + 1.0);
		
		domainParams.totalBucketCount = domainParams.XBucketCount * domainParams.YBucketCount * domainParams.ZBucketCount;
		std::cout<<"total bucket count: "<< domainParams.totalBucketCount<<std::endl;


		std::cout<<"maxX: "<< domainParams.maxX<<  std::endl;
		std::cout<<"minZ: "<< domainParams.minZ<<  std::endl;
		std::cout<<"maxZ: "<< domainParams.maxZ<<  std::endl;
		std::cout << "totalBucketsX" << domainParams.XBucketCount << std::endl;
		std::cout << "totalBucketsY" << domainParams.YBucketCount << std::endl;
		std::cout << "totalBucketsZ" << domainParams.ZBucketCount << std::endl;
		
		auxVecs.keyBegin.resize(domainParams.totalBucketCount);
		auxVecs.keyEnd.resize(domainParams.totalBucketCount);
		
	}
	thrust::fill(auxVecs.keyBegin.begin(),auxVecs.keyBegin.end(),0);
	thrust::fill(auxVecs.keyEnd.begin(),auxVecs.keyEnd.end(),0);

}

void extendBucketScheme(
	CoordInfoVecs& coordInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs) {
		
	//memory is already allocated. 
	int endIndexExpanded = (auxVecs.endIndexid_bucket) * 27;

	//test for removing copies. 
	int valuesCount = auxVecs.id_value.size();
	thrust::fill(auxVecs.id_bucket_expanded.begin(),auxVecs.id_bucket_expanded.end(),0);
	thrust::fill(auxVecs.id_value_expanded.begin(),auxVecs.id_value_expanded.end(),0);

	

	/*
	* beginning of constant iterator
	*/
	thrust::constant_iterator<int> first(27);
	/** 
	* end of constant iterator.
	* the plus sign only indicate movement of position, not value.
	* e.g. movement is 5 and first iterator is initialized as 9
	* result array is [9,9,9,9,9];
	*/
	thrust::constant_iterator<int> last = first + (auxVecs.endIndexid_bucket); // this is NOT numerical addition!

	expand(first, last,
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.id_bucket.begin(),
				auxVecs.id_value.begin())),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.id_bucket_expanded.begin(),
				auxVecs.id_value_expanded.begin())));
		

	thrust::counting_iterator<int> countingBegin(0);

	thrust::transform(
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.id_bucket_expanded.begin(),
				countingBegin)),												  
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.id_bucket_expanded.begin(),
				countingBegin)) + endIndexExpanded,
			auxVecs.id_bucket_expanded.begin(),
				
		NeighborFunctor(
			domainParams.XBucketCount,
			domainParams.YBucketCount,
			domainParams.ZBucketCount));



	//int numberOfOutOfRange = thrust::count_if(auxVecs.id_bucket_expanded.begin(),
	//	auxVecs.id_bucket_expanded.end(), is_greater_than(domainParams.totalBucketCount) );
	//int numberInsideRange = endIndexExpanded - numberOfOutOfRange;

	//int endIndexSearch = endIndexExpanded - numberOfOutOfRange;

	thrust::sort_by_key(auxVecs.id_bucket_expanded.begin(),
		auxVecs.id_bucket_expanded.end(),
		auxVecs.id_value_expanded.begin());

//	std::cout<<"out of range: " << numberOfOutOfRange << std::endl;
	
//	std::cout<<"in of range: " << numberInsideRange << std::endl;

/*	auxVecs.id_bucket_expanded.erase(
			auxVecs.id_bucket_expanded.begin() + numberInsideRange,
			auxVecs.id_bucket_expanded.end());

	auxVecs.id_value_expanded.erase(
			auxVecs.id_value_expanded.begin() + numberInsideRange,
			auxVecs.id_value_expanded.end());


*/
	
	thrust::counting_iterator<int> search_begin(0);

	thrust::lower_bound(auxVecs.id_bucket_expanded.begin(),
		auxVecs.id_bucket_expanded.end(), search_begin,
		search_begin + domainParams.totalBucketCount,
		auxVecs.keyBegin.begin());

	thrust::upper_bound(auxVecs.id_bucket_expanded.begin(),
		auxVecs.id_bucket_expanded.end(),search_begin,
		search_begin + domainParams.totalBucketCount, 
		auxVecs.keyEnd.begin()); 

}	   

 
void buildBucketScheme(
	CoordInfoVecs& coordInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams,
	CapsidInfoVecs& capsidInfoVecs) {


	thrust::counting_iterator<int> indexBucketBegin(0);
	// takes counting iterator and coordinates
	// return tuple of keys and values
	// transform the points to their bucket indices
	thrust::transform(
		thrust::make_zip_iterator(
			thrust::make_tuple(
				coordInfoVecs.nodeLocX.begin(),
				coordInfoVecs.nodeLocY.begin(),
				coordInfoVecs.nodeLocZ.begin(),
				indexBucketBegin)),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				coordInfoVecs.nodeLocX.begin(),
				coordInfoVecs.nodeLocY.begin(),
				coordInfoVecs.nodeLocZ.begin(),
				indexBucketBegin)) + generalParams.maxNodeCount,
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.id_bucket.begin(),
				auxVecs.id_value.begin())),			    
		BucketIndexer(
			domainParams.minX, domainParams.maxX, domainParams.minY,
			domainParams.maxY, domainParams.minZ, domainParams.maxZ,
			domainParams.XBucketCount,
			domainParams.YBucketCount,
			domainParams.ZBucketCount,
			domainParams.gridSpacing));

	//transform capsid points to buckets. 
	thrust::counting_iterator<int> indexCapsidBucketBegin(0);

/*	thrust::transform(
		thrust::make_zip_iterator(
			thrust::make_tuple(
				capsidInfoVecs.nodeLocX.begin(),
				capsidInfoVecs.nodeLocY.begin(),
				capsidInfoVecs.nodeLocZ.begin(),
				indexCapsidBucketBegin)),		
		thrust::make_zip_iterator(
			thrust::make_tuple(
				capsidInfoVecs.nodeLocX.begin(),
				capsidInfoVecs.nodeLocY.begin(),
				capsidInfoVecs.nodeLocZ.begin(),
				indexCapsidBucketBegin)) + capsidInfoVecs.maxNodeCount,			
		thrust::make_zip_iterator(
			thrust::make_tuple(
				capsidInfoVecs.id_bucket.begin(),
				capsidInfoVecs.id_value.begin())),		    
		BucketIndexer(
			domainParams.minX, domainParams.maxX, domainParams.minY,
			domainParams.maxY, domainParams.minZ, domainParams.maxZ,
			domainParams.XBucketCount,
			domainParams.YBucketCount,
			domainParams.ZBucketCount,
			domainParams.gridSpacing))*/

thrust::sort_by_key(auxVecs.id_value.begin(),
		auxVecs.id_value.begin() + generalParams.maxNodeCount,
		auxVecs.id_bucket.begin());


										

};


