/**
 * @file   zonotope_test.h
 * @brief  test zonotope implementation
 * @author Yunjun Bai
 * @date April 2021
 * @version 1.0 
 * 
 */

#include "test.h"
#include "../zonotope/Zonotope.h"

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
	Eigen::Matrix<TypeParam, 2, 3> gen;
	Eigen::Matrix<TypeParam, 2, 1> center;

	gen << 2, 3, 6, 4, 1, 5;
	center << 1, 2;

	Zonotope<TypeParam> z2(center, gen);
	EXPECT_EQ(z2.dimension(), (unsigned)2);
	EXPECT_EQ(z2.center(), center);
	EXPECT_EQ(z2.generators(), gen);
}