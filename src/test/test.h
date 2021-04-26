/**
 * @file   test.h
 * @brief  test headers
 * @author Yunjun Bai
 * @date April 2021
 * @version 1.0 
 * 
 */

#pragma once
#include "gtest/gtest.h"

typedef ::testing::Types<double> allTypes;

TYPED_TEST_CASE(ZonotopeTest, allTypes);
