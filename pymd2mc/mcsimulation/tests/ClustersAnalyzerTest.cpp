/*
 * ClustersAnalyzerTest.cpp
 *
 *  Created on: 14-11-2011
 *      Author: lisu
 */

#include <iostream>
using namespace std;

#include <assert.h>
#include <vector>

#include "ClustersAnalyzer.h"

void testClustersAnalyzer();
void testRegisterAtom();
void testIsClustered();
void testIsMixed1();
void testIsMixed2();

int main()
{
    testClustersAnalyzer();
    testRegisterAtom();
    testIsClustered();
    testIsMixed1();
    testIsMixed2();

    return 0;
}

void testClustersAnalyzer()
{
    ClustersAnalyzer* analyzer = new ClustersAnalyzer( 10, 10.0 );
    assert( analyzer != NULL );
    delete analyzer;
}

void testRegisterAtom()
{
    ClustersAnalyzer* analyzer = new ClustersAnalyzer( 1, 1. );
    vector< Distance > v;
    assert( analyzer->neighborPairs.size() == 0 );  //should be empty
    //create set of distances
    for ( int i = 1 ; i < 20 ; ++i )
    {
        v.push_back( Distance( 0.11 * i, i ) );
    }
    analyzer->registerAtom( 0, v, 0 );
    assert( analyzer->neighborPairs.size() == 9 ); //assert count of pairs
    analyzer->registerAtom( 0, v, 1 );
    assert( analyzer->neighborPairs.size() == 9 );
    /// validate two randomly chosen
    assert( analyzer->neighborPairs[ ClustersAnalyzer::Pair( 0, 1 )] == ClustersAnalyzer::Pair( 2, 1 ) );
    assert( analyzer->neighborPairs[ ClustersAnalyzer::Pair( 0, 2 )] == ClustersAnalyzer::Pair( 2, 1 ) );
    // now add one to fit distThr
    v[9].d = 0.9;
    analyzer->registerAtom( 0, v, 2 );
    assert( analyzer->neighborPairs[ ClustersAnalyzer::Pair( 0, 10 )] == ClustersAnalyzer::Pair( 1, 2 ) );
    // and move one away
    v[1].d = 1.1;
    analyzer->registerAtom( 0, v, 3 );
    assert( analyzer->neighborPairs[ ClustersAnalyzer::Pair( 0, 2 )] == ClustersAnalyzer::Pair( 3, 2 ) );
    // and move it againg close
    v[1].d = 0.1;
    analyzer->registerAtom( 0, v, 4 );
    assert( analyzer->neighborPairs[ ClustersAnalyzer::Pair( 0, 2 )] == ClustersAnalyzer::Pair( 1, 4 ) );
    delete analyzer;
}

void testIsClustered()
{
    ClustersAnalyzer* analyzer = new ClustersAnalyzer( 2, 1. );
    vector< Distance > v;
    for ( int i = 1 ; i < 20 ; ++i )
    {
        v.push_back( Distance( 0.11 * i, i ) );
    }
    analyzer->registerAtom( 0, v, 0 );
    assert( ! analyzer->isClustered( 0, 0 ) );
    analyzer->registerAtom( 0, v, 1 );
    assert( ! analyzer->isClustered( 0, 1 ) );
    analyzer->registerAtom( 0, v, 2 );
    assert( analyzer->isClustered( 0, 2 ) );// now it is
    
    vector< Distance > v2;
    v2.push_back( Distance( 0.11, 0 ) );
    analyzer->registerAtom( 1, v2, 0 ); 
    assert( !analyzer->isClustered( 1, 0 ) );
    analyzer->registerAtom( 1, v2, 1 ); 
    assert( !analyzer->isClustered( 1, 1 ) );
    analyzer->registerAtom( 1, v2, 2 ); 
    assert( analyzer->isClustered( 1, 2 ) );
    v2[0].at2Ind = 1.1;
    analyzer->registerAtom( 1, v2, 3 ); 
    assert( !analyzer->isClustered( 1, 3 ) );
    v2[0].at2Ind = 0.1;
    analyzer->registerAtom( 1, v2, 4 ); 
    assert( !analyzer->isClustered( 1, 4 ) );
    analyzer->registerAtom( 1, v2, 5 ); 
    assert( !analyzer->isClustered( 1, 5 ) );
    analyzer->registerAtom( 1, v2, 6 ); 
    assert( analyzer->isClustered( 1, 6 ) );
}
void testIsMixed1()
{
    ClustersAnalyzer* analyzer = new ClustersAnalyzer( 2, 1. );
    Atom* atoms = new Atom[2];
    atoms[0] = Atom("A", "A", 0,0,0 );
    atoms[1] = Atom("A", "A", 0,0,0 ); // consider two same atoms

    vector< Distance > v2;
    v2.push_back( Distance( 0.11, 0 ) );
    analyzer->registerAtom( 1, v2, 0, atoms ); 
    assert( !analyzer->isClustered( 1, 0 ) );
    assert( !analyzer->isInMixedCluster( 1, 0 ) );
    analyzer->registerAtom( 1, v2, 1 , atoms); 
    assert( !analyzer->isClustered( 1, 1 ) );
    assert( !analyzer->isInMixedCluster( 1, 1 ) );
    analyzer->registerAtom( 1, v2, 2, atoms ); 
    assert( analyzer->isClustered( 1, 2 ) );
    assert( !analyzer->isInMixedCluster( 1, 2 ) );
    v2[0].at2Ind = 1.1;
    analyzer->registerAtom( 1, v2, 3, atoms ); 
    assert( !analyzer->isClustered( 1, 3 ) );
    assert( !analyzer->isInMixedCluster( 1, 3 ) );
    v2[0].at2Ind = 0.1;
    analyzer->registerAtom( 1, v2, 4, atoms ); 
    assert( !analyzer->isClustered( 1, 4 ) );
    assert( !analyzer->isInMixedCluster( 1, 4 ) );
    analyzer->registerAtom( 1, v2, 5, atoms ); 
    assert( !analyzer->isClustered( 1, 5 ) );
    assert( !analyzer->isInMixedCluster( 1, 5 ) );
    analyzer->registerAtom( 1, v2, 6, atoms ); 
    assert( analyzer->isClustered( 1, 6 ) );
    assert( !analyzer->isInMixedCluster( 1, 6 ) );
}
void testIsMixed2()
{
    ClustersAnalyzer* analyzer = new ClustersAnalyzer( 2, 1. );
    Atom* atoms = new Atom[2];
    atoms[0] = Atom("A", "A", 0,0,0 );
    atoms[1] = Atom("B", "A", 0,0,0 ); // consider two same atoms

    vector< Distance > v2;
    v2.push_back( Distance( 0.11, 0 ) );
    analyzer->registerAtom( 1, v2, 0, atoms ); 
    assert( !analyzer->isClustered( 1, 0 ) );
    assert( !analyzer->isInMixedCluster( 1, 0 ) );
    analyzer->registerAtom( 1, v2, 1 , atoms); 
    assert( !analyzer->isClustered( 1, 1 ) );
    assert( !analyzer->isInMixedCluster( 1, 1 ) );
    analyzer->registerAtom( 1, v2, 2, atoms ); 
    assert( analyzer->isClustered( 1, 2 ) );
    assert( analyzer->isInMixedCluster( 1, 2 ) );
    v2[0].at2Ind = 1.1;
    analyzer->registerAtom( 1, v2, 3, atoms ); 
    assert( !analyzer->isClustered( 1, 3 ) );
    assert( !analyzer->isInMixedCluster( 1, 3 ) );
    v2[0].at2Ind = 0.1;
    analyzer->registerAtom( 1, v2, 4, atoms ); 
    assert( !analyzer->isClustered( 1, 4 ) );
    assert( !analyzer->isInMixedCluster( 1, 4 ) );
    analyzer->registerAtom( 1, v2, 5, atoms ); 
    assert( !analyzer->isClustered( 1, 5 ) );
    assert( !analyzer->isInMixedCluster( 1, 5 ) );
    analyzer->registerAtom( 1, v2, 6, atoms ); 
    assert( analyzer->isClustered( 1, 6 ) );
    assert( analyzer->isInMixedCluster( 1, 6 ) );
}

