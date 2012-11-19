// std
#include <iostream>
#include <assert.h>
#include <tr1/memory>

// stl
#include <set>

// Google test tools
#include "gtest/gtest.h"
#include "gmock/gmock.h"

// project-local
#include "ProteinTriangularLattice.h"

// tests
#include "testUtils.h"
using namespace testUtils;


TEST( ProteinTriangularLattice, Construct )
{
    static const int LIPID_B_COUNT( 20 );
    std::tr1::shared_ptr< ProteinTriangularLattice > latt( new ProteinTriangularLattice( 100, 10, LIPID_B_COUNT, 0 ) );
    TriangularLattice* localPtr( latt.get() ); // calcSum function accepts TriangularLattice* - this is stupid
    EXPECT_EQ( calcSum( localPtr ), LIPID_B_COUNT * LIPID_B );
}

TEST( ProteinTriangularLattice, Construct_withProteins_Random )
{
    static const int LIPID_B_COUNT( 10 );
    static const int PROTEIN_COUNT( 3 );
    std::tr1::shared_ptr< ProteinTriangularLattice > latt( new ProteinTriangularLattice( 100, 10, LIPID_B_COUNT, PROTEIN_COUNT ) );
    TriangularLattice* localPtr( latt.get() ); // calcSum function accepts TriangularLattice* - this is stupid
    printLatt( latt->getLattice(), 10, 10 );
    EXPECT_EQ( calcSum( localPtr ),  LIPID_B_COUNT * LIPID_B + PROTEIN_COUNT * PROTEIN_A * ProteinTriangularLattice::mProteinSize );
}

TEST( ProteinTriangularLattice, Construct_withProteins_NONRandom_SingleProtein )
{
    static const int LIPID_B_COUNT( 0 );
    static const int PROTEIN_COUNT( 1 );
    std::tr1::shared_ptr< ProteinTriangularLattice > latt( new ProteinTriangularLattice( 100, 10, LIPID_B_COUNT, PROTEIN_COUNT, false ) );
    TriangularLattice* localPtr( latt.get() ); // calcSum function accepts TriangularLattice* - this is stupid
    EXPECT_EQ( calcSum( localPtr ),  LIPID_B_COUNT * LIPID_B + PROTEIN_COUNT * PROTEIN_A * ProteinTriangularLattice::mProteinSize );
}
TEST( ProteinTriangularLattice, Construct_withProteins_NONRandom )
{
    static const int LIPID_B_COUNT( 10 );
    static const int PROTEIN_COUNT( 3 );
    std::tr1::shared_ptr< ProteinTriangularLattice > latt( new ProteinTriangularLattice( 100, 10, LIPID_B_COUNT, PROTEIN_COUNT, false ) );
    TriangularLattice* localPtr( latt.get() ); // calcSum function accepts TriangularLattice* - this is stupid
    EXPECT_EQ( calcSum( localPtr ),  LIPID_B_COUNT * LIPID_B + PROTEIN_COUNT * PROTEIN_A * ProteinTriangularLattice::mProteinSize );
    printLatt( latt->getLattice(), 10, 10 );
    for ( int i = 0 ; i < PROTEIN_COUNT ; ++i )
    {
        EXPECT_EQ( ( *latt )[ latt->getProteins()[i] ], PROTEIN_A );
    }
}

TEST( ProteinTriangularLattice, moveProteinRight )
{
    static int LIPID_B_COUNT( 0 );
    static const int PROTEIN_COUNT( 1 );
    static const ProteinTriangularLattice::lattIndex proteinPos( 0 );
    static const ProteinTriangularLattice::lattIndex siteOnTheRight( 2 );
    std::tr1::shared_ptr< ProteinTriangularLattice > latt( new ProteinTriangularLattice( 100, 10, LIPID_B_COUNT, PROTEIN_COUNT, false ) );

    // let's edit lattice a little by hand
    LIPID_B_COUNT += 3;
    latt->getLattice()[2]= LIPID_B  ;
    latt->getLattice()[11]= LIPID_B ;
    latt->getLattice()[92]= LIPID_B ;

    printLatt( latt->getLattice(), 10, 10 );
    latt->moveProtein( proteinPos, siteOnTheRight ); 

    printLatt( latt->getLattice(), 10, 10 );

    // verify movement by hand
    EXPECT_EQ( latt->getLattice()[0], PROTEIN_A );
    EXPECT_EQ( latt->getLattice()[1], PROTEIN_A );
    EXPECT_EQ( latt->getLattice()[2], PROTEIN_A );
    EXPECT_EQ( latt->getLattice()[10], PROTEIN_A );
    EXPECT_EQ( latt->getLattice()[11], PROTEIN_A );
    EXPECT_EQ( latt->getLattice()[91], PROTEIN_A );
    EXPECT_EQ( latt->getLattice()[92], PROTEIN_A );
    EXPECT_EQ( latt->getLattice()[82], LIPID_B );
    EXPECT_EQ( latt->getLattice()[20], LIPID_B );

    TriangularLattice* localPtr( latt.get() ); 
    EXPECT_EQ( calcSum( localPtr ),  LIPID_B_COUNT * LIPID_B + PROTEIN_COUNT * PROTEIN_A * ProteinTriangularLattice::mProteinSize );
    ProteinTriangularLattice::lattIndex proteinPosAfterMove( 1 );// calculated by hand
    EXPECT_EQ( latt->getProteins()[0], proteinPosAfterMove ); 
}
