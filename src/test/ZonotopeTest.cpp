/**
 * @file   zonotope_test.h
 * @brief  test zonotope implementation
 * @author Yunjun Bai
 * @date April 2021
 * @version 1.0 
 * 
 */

#include "test.h"
#include "zonotope/Zonotope.h"

#include "gtest/gtest.h"
#include <iostream>
#include <stdlib.h>

using namespace reachSolver;

template <typename Number>
class ZonotopeTest : public ::testing::Test {
  protected:
    virtual void SetUp(){}

    virtual void TearDown(){}
};

TYPED_TEST(ZonotopeTest, NoParamConstructor){
    Zonotope<TypeParam> z1;
    //EXPECT_TRUE(z1.empty());
    EXPECT_EQ(z1.order(), TypeParam(0));
    EXPECT_EQ(z1.dimension(), (unsigned)0);
    EXPECT_EQ(z1.center().rows(), 0);
    EXPECT_EQ(z1.center().cols(), 1);
    EXPECT_EQ(z1.generators().rows(), 0);
    EXPECT_EQ(z1.generators().cols(), 0);

    // Zonotope<TypeParam> emp = Zonotope<TypeParam>::Empty();
    // EXPECT_TRUE(emp.empty());
}

TYPED_TEST(ZonotopeTest, DimConstructor){
    unsigned int dim = 4;
    Zonotope<TypeParam> z1(dim);
    EXPECT_EQ(z1.dimension(), (unsigned)4);
}

TYPED_TEST(ZonotopeTest, CenGenConstructor){
    Eigen::Matrix<TypeParam, 2, 1> center;
    Eigen::Matrix<TypeParam, 2, 3> gen;

    center << 1, 2;
    gen << 2, 3, 6, 4, 1, 5;

    Zonotope<TypeParam> z1(center, gen);
    EXPECT_EQ(z1.dimension(), (unsigned)2);
    EXPECT_EQ(z1.center(), center);
    EXPECT_EQ(z1.generators(), gen);
}

TYPED_TEST(ZonotopeTest, Display){
    Eigen::Matrix<TypeParam, 2, 1> center;
    Eigen::Matrix<TypeParam, 2, 3> gen;

    center << 1, 2;
    gen << 2, 3, 6, 4, 1, 5;

    Zonotope<TypeParam> z1(center, gen);
    z1.Display();
} 

TYPED_TEST(ZonotopeTest, Plus){
    Eigen::Matrix<TypeParam, 2, 1> center1, center2, cen_sum;
    Eigen::Matrix<TypeParam, 2, 2> gen1, gen2;
    Eigen::Matrix<TypeParam, 2, 4> gen_sum;

    center1 << 1, 2;
    center2 << 2, 2;
    gen1 << 2, 3, 6, 4;
    gen2 << 3, 3, 4, 4;

    Zonotope<TypeParam> z1(center1, gen1);
    Zonotope<TypeParam> z2(center2, gen2);
    z1 = z1.Plus(z2);

    cen_sum << 3, 4;
    gen_sum << 2, 3, 6, 4, 3, 3, 4, 4;

    EXPECT_EQ(z1.center(), cen_sum);
    EXPECT_EQ(z1.generators(), gen_sum);
}
