/*
 * ClustersAnalyzerTest.cpp
 *
 *  Created on: 14-11-2011
 *      Author: lisu
 */

#include <iostream>
using namespace std;

#include <assert.h>

#include "ClustersAnalyzer.h"

void testClustersAnalyzer();

int main()
{
    testClustersAnalyzer();
    return 0;
}

void testClustersAnalyzer()
{
    ClustersAnalyzer* analyzer = new ClustersAnalyzer( 10, 10.0 );
    assert( analyzer != NULL );
    delete analyzer;
}

