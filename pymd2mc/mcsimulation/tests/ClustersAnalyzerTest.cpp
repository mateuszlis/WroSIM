/*
 * ClustersAnalyzerTest.cpp
 *
 *  Created on: 14-11-2011
 *      Author: lisu
 */
// std
#include <iostream>

// stl
#include <vector>

// Google test tools
#include "gtest/gtest.h"

// project-local
#include "ClustersAnalyzer.h"

// Most of the tests use EXPECT_TRUE instead of appropriate gtest macros
// because it was basically ported from previous test framework.

TEST( ClusterAnalyzerTest, Construct )
{
    ClustersAnalyzer* analyzer = new ClustersAnalyzer( 10, 10.0 );
    EXPECT_TRUE( analyzer != NULL );
    delete analyzer;
}

TEST( ClusterAnalyzerTest, RegisterAtom )
{
    ClustersAnalyzer* analyzer = new ClustersAnalyzer( 1, 1. );
    vector< Distance > v;
    EXPECT_TRUE( analyzer->neighborPairs.size() == 0 );  //should be empty
    //create set of distances
    for ( int i = 1 ; i < 20 ; ++i )
    {
        v.push_back( Distance( 0.11 * i, i ) );
    }
    analyzer->registerAtom( 0, v, 0 );
    EXPECT_TRUE( analyzer->neighborPairs.size() == 9 ); //EXPECT_TRUE count of pairs
    analyzer->registerAtom( 0, v, 1 );
    EXPECT_TRUE( analyzer->neighborPairs.size() == 9 );
    /// validate two randomly chosen
    EXPECT_TRUE( analyzer->neighborPairs[ ClustersAnalyzer::Pair( 0, 1 )] == ClustersAnalyzer::Pair( 2, 1 ) );
    EXPECT_TRUE( analyzer->neighborPairs[ ClustersAnalyzer::Pair( 0, 2 )] == ClustersAnalyzer::Pair( 2, 1 ) );
    // now add one to fit distThr
    v[9].d = 0.9;
    analyzer->registerAtom( 0, v, 2 );
    EXPECT_TRUE( analyzer->neighborPairs[ ClustersAnalyzer::Pair( 0, 10 )] == ClustersAnalyzer::Pair( 1, 2 ) );
    // and move one away
    v[1].d = 1.1;
    analyzer->registerAtom( 0, v, 3 );
    EXPECT_TRUE( analyzer->neighborPairs[ ClustersAnalyzer::Pair( 0, 2 )] == ClustersAnalyzer::Pair( 3, 2 ) );
    // and move it againg close
    v[1].d = 0.1;
    analyzer->registerAtom( 0, v, 4 );
    EXPECT_TRUE( analyzer->neighborPairs[ ClustersAnalyzer::Pair( 0, 2 )] == ClustersAnalyzer::Pair( 1, 4 ) );
    delete analyzer;
}

TEST( ClusterAnalyzerTest, IsClustered )
{
    ClustersAnalyzer* analyzer = new ClustersAnalyzer( 2, 1. );
    vector< Distance > v;
    for ( int i = 1 ; i < 20 ; ++i )
    {
        v.push_back( Distance( 0.11 * i, i ) );
    }
    analyzer->registerAtom( 0, v, 0 );
    EXPECT_TRUE( ! analyzer->isClustered( 0, 0 ) );
    analyzer->registerAtom( 0, v, 1 );
    EXPECT_TRUE( ! analyzer->isClustered( 0, 1 ) );
    analyzer->registerAtom( 0, v, 2 );
    EXPECT_TRUE( analyzer->isClustered( 0, 2 ) );// now it is
    
    vector< Distance > v2;
    v2.push_back( Distance( 0.11, 0 ) );
    analyzer->registerAtom( 1, v2, 0 ); 
    EXPECT_TRUE( !analyzer->isClustered( 1, 0 ) );
    analyzer->registerAtom( 1, v2, 1 ); 
    EXPECT_TRUE( !analyzer->isClustered( 1, 1 ) );
    analyzer->registerAtom( 1, v2, 2 ); 
    EXPECT_TRUE( analyzer->isClustered( 1, 2 ) );
    v2[0].at2Ind = 1.1;
    analyzer->registerAtom( 1, v2, 3 ); 
    EXPECT_TRUE( !analyzer->isClustered( 1, 3 ) );
    v2[0].at2Ind = 0.1;
    analyzer->registerAtom( 1, v2, 4 ); 
    EXPECT_TRUE( !analyzer->isClustered( 1, 4 ) );
    analyzer->registerAtom( 1, v2, 5 ); 
    EXPECT_TRUE( !analyzer->isClustered( 1, 5 ) );
    analyzer->registerAtom( 1, v2, 6 ); 
    EXPECT_TRUE( analyzer->isClustered( 1, 6 ) );
    delete analyzer;
}
TEST( ClusterAnalyzerTest, IsMixed1 )
{
    ClustersAnalyzer* analyzer = new ClustersAnalyzer( 2, 1. );
    Atom* atoms = new Atom[2];
    atoms[0] = Atom("A", "A", 0,0,0 );
    atoms[1] = Atom("A", "A", 0,0,0 ); // consider two same atoms

    vector< Distance > v2;
    v2.push_back( Distance( 0.11, 0 ) );
    analyzer->registerAtom( 1, v2, 0, atoms ); 
    EXPECT_TRUE( !analyzer->isClustered( 1, 0 ) );
    EXPECT_TRUE( !analyzer->isInMixedCluster( 1, 0 ) );
    analyzer->registerAtom( 1, v2, 1 , atoms); 
    EXPECT_TRUE( !analyzer->isClustered( 1, 1 ) );
    EXPECT_TRUE( !analyzer->isInMixedCluster( 1, 1 ) );
    analyzer->registerAtom( 1, v2, 2, atoms ); 
    EXPECT_TRUE( analyzer->isClustered( 1, 2 ) );
    EXPECT_TRUE( !analyzer->isInMixedCluster( 1, 2 ) );
    v2[0].at2Ind = 1.1;
    analyzer->registerAtom( 1, v2, 3, atoms ); 
    EXPECT_TRUE( !analyzer->isClustered( 1, 3 ) );
    EXPECT_TRUE( !analyzer->isInMixedCluster( 1, 3 ) );
    v2[0].at2Ind = 0.1;
    analyzer->registerAtom( 1, v2, 4, atoms ); 
    EXPECT_TRUE( !analyzer->isClustered( 1, 4 ) );
    EXPECT_TRUE( !analyzer->isInMixedCluster( 1, 4 ) );
    analyzer->registerAtom( 1, v2, 5, atoms ); 
    EXPECT_TRUE( !analyzer->isClustered( 1, 5 ) );
    EXPECT_TRUE( !analyzer->isInMixedCluster( 1, 5 ) );
    analyzer->registerAtom( 1, v2, 6, atoms ); 
    EXPECT_TRUE( analyzer->isClustered( 1, 6 ) );
    EXPECT_TRUE( !analyzer->isInMixedCluster( 1, 6 ) );
    delete[] atoms;
    delete analyzer;
}

TEST( ClusterAnalyzerTest, IsMixed2 )
{
    ClustersAnalyzer* analyzer = new ClustersAnalyzer( 2, 1. );
    Atom* atoms = new Atom[2];
    atoms[0] = Atom("A", "A", 0,0,0 );
    atoms[1] = Atom("B", "A", 0,0,0 ); // consider two same atoms

    vector< Distance > v2;
    v2.push_back( Distance( 0.11, 0 ) );
    analyzer->registerAtom( 1, v2, 0, atoms ); 
    EXPECT_TRUE( !analyzer->isClustered( 1, 0 ) );
    EXPECT_TRUE( !analyzer->isInMixedCluster( 1, 0 ) );
    analyzer->registerAtom( 1, v2, 1 , atoms); 
    EXPECT_TRUE( !analyzer->isClustered( 1, 1 ) );
    EXPECT_TRUE( !analyzer->isInMixedCluster( 1, 1 ) );
    analyzer->registerAtom( 1, v2, 2, atoms ); 
    EXPECT_TRUE( analyzer->isClustered( 1, 2 ) );
    EXPECT_TRUE( analyzer->isInMixedCluster( 1, 2 ) );
    v2[0].at2Ind = 1.1;
    analyzer->registerAtom( 1, v2, 3, atoms ); 
    EXPECT_TRUE( !analyzer->isClustered( 1, 3 ) );
    EXPECT_TRUE( !analyzer->isInMixedCluster( 1, 3 ) );
    v2[0].at2Ind = 0.1;
    analyzer->registerAtom( 1, v2, 4, atoms ); 
    EXPECT_TRUE( !analyzer->isClustered( 1, 4 ) );
    EXPECT_TRUE( !analyzer->isInMixedCluster( 1, 4 ) );
    analyzer->registerAtom( 1, v2, 5, atoms ); 
    EXPECT_TRUE( !analyzer->isClustered( 1, 5 ) );
    EXPECT_TRUE( !analyzer->isInMixedCluster( 1, 5 ) );
    analyzer->registerAtom( 1, v2, 6, atoms ); 
    EXPECT_TRUE( analyzer->isClustered( 1, 6 ) );
    EXPECT_TRUE( analyzer->isInMixedCluster( 1, 6 ) );
    delete[] atoms;
    delete analyzer;
}

