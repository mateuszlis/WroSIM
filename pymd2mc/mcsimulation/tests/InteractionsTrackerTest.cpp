// std
#include <iostream>
#include <memory>
// Google test tools
#include "gtest/gtest.h"
#include "gmock/gmock.h"

// project-local
#define protected public
#define private public
#include "InteractionsTracker.h"
#include "ProteinTriangularLattice.h"
#undef protected
#undef private
#include "testUtils.h"


TEST( InteractionsTrackerTest, Ctor_wo_Proteins )
{
    static const int LIPID_B_COUNT( 20 );
    std::shared_ptr< ProteinTriangularLattice > latt( new ProteinTriangularLattice( 100, 10, LIPID_B_COUNT, 0 ) );
    InteractionsTracker tracker( latt.get(), nullptr, 0 );

}
TEST( InteractionsTrackerTest, Ctor_Proteins )
{
    static const int LIPID_B_COUNT( 10 );
    static const int PROTEIN_COUNT( 3 );
    std::shared_ptr< ProteinTriangularLattice > latt( new ProteinTriangularLattice( 100, 10, LIPID_B_COUNT, PROTEIN_COUNT ) );
    InteractionsTracker tracker( latt.get(), latt->mProteins, PROTEIN_COUNT );
}

TEST( InteractionsTrackerTest, Ctor_Proteins_Interactions )
{
    static const int LIPID_B_COUNT( 10 );
    static const int PROTEIN_COUNT( 3 );
    std::shared_ptr< ProteinTriangularLattice > latt( new ProteinTriangularLattice( 81, 9, LIPID_B_COUNT, PROTEIN_COUNT, false ) );
    testUtils::printLatt( latt->getLattice(), 9, 9 );
    InteractionsTracker tracker( latt.get(), latt->mProteins, PROTEIN_COUNT );
    testUtils::printLatt( tracker.mpProtInt.get(), 9, 9 );

    EXPECT_DOUBLE_EQ( tracker.getInteraction( 9 ), 1.5 );
    EXPECT_DOUBLE_EQ( tracker.getInteraction( 21 ), 1.833333 );
    EXPECT_NO_THROW( tracker.registerProteinMove( 56, 65 ) );
    EXPECT_DOUBLE_EQ( tracker.getInteraction( 9 ), 1.833333 );
    EXPECT_DOUBLE_EQ( tracker.getInteraction( 21 ), 1.5 );
}
