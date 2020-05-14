#ifndef EDGESWAP_TEST_H_
#define EDGESWAP_TEST_H_

#include "SystemStructures.h"

class Edgeswap {
    std::vector<bool> boundary_node;
    std::vector<unsigned> nndata;

    public:
    Edgeswap(CoordInfoVecs& coordInfoVecs);

	std::vector<bool> DomainBd (CoordInfoVecs& coordInfoVecs);
	std::vector<unsigned> Number_of_Neighbor(CoordInfoVecs& coordInfoVecs);
    void edge_swap (unsigned iedge, 
        GeneralParams& generalParams,
        CoordInfoVecs& coordInfoVecs,
        LinearSpringInfoVecs& linearSpringInfoVecs,
        BendingTriangleInfoVecs& bendingTriangleInfoVecs);
	
};

#endif