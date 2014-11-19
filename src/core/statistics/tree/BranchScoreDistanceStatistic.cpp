#include "BranchScoreDistanceStatistic.h"
#include "TreeBipartitions.h"
#include <algorithm>

#include <boost/dynamic_bitset.hpp>
#include <math.h>

using namespace RevBayesCore;

//
//  BranchScoreDistanceStatistic.cpp
//  RevBayesCore
//
//  Created by Michael Landis on 16/08/14.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

BranchScoreDistanceStatistic::BranchScoreDistanceStatistic(const TypedDagNode<TimeTree> *t1, const TypedDagNode<TimeTree> *t2) : TypedFunction< double >( new double(0.0) ),
tree1( t1 ), tree2(t2)
{
    // add the tree parameter as a parent
    addParameter( tree1 );
    addParameter( tree2 );
    
    update();
}


BranchScoreDistanceStatistic::~BranchScoreDistanceStatistic( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



BranchScoreDistanceStatistic* BranchScoreDistanceStatistic::clone( void ) const
{
    
    return new BranchScoreDistanceStatistic( *this );
}


void BranchScoreDistanceStatistic::update( void )
{
    
    //const TopologyNode& r = tree->getValue().getRoot();
    TreeBipartitions tbp1(tree1);
    TreeBipartitions tbp2(tree2);
    std::vector<boost::dynamic_bitset<> > bipartitions1 = tbp1.getValue();
    std::vector<boost::dynamic_bitset<> > bipartitions2 = tbp2.getValue();
    const std::vector<double> ages1 = tbp1.getBipartitionAges();
    const std::vector<double> ages2 = tbp2.getBipartitionAges();
    
    bool found;
    *value = 0.0;
    for (size_t i = 0; i< bipartitions1.size(); ++i) {
        found = false;
        for (size_t j = 0; j < bipartitions2.size(); ++j) {
            if (bipartitions1[i] == bipartitions2[j]) {
                found = true;
                *value += (ages1[i]-ages2[j])*(ages1[i]-ages2[j]);
                break;
            }
        }
        if (!found) {
            *value += ages1[i]*ages1[i];
        }
    }
    for (size_t i = 0; i< bipartitions2.size(); ++i) {
        found = false;
        for (size_t j = 0; j < bipartitions1.size(); ++j) {
            if (bipartitions2[i] == bipartitions1[j]) {
                found = true;
                *value += (ages1[i]-ages2[j])*(ages1[i]-ages2[j]);
                break;
            }
        }
        if (!found) {
            *value += ages2[i]*ages2[i];
        }
    }
    *value = pow(*value, 0.5);
    
}



void BranchScoreDistanceStatistic::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == tree1)
    {
        tree1 = static_cast<const TypedDagNode<TimeTree>* >( newP );
    }
    if (oldP == tree2)
    {
        tree2 = static_cast<const TypedDagNode<TimeTree>* >( newP );
    }
    
}
