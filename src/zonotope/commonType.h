/**
 * @file   commonType.h
 * @brief 
 * @author 
 * @date
 * @version 1.0
 * 
 */

#pragma once
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/StdVector>
#include <set>

namespace reachSolver{

/// Wrapper for an Eigen::Matrix type with only one column.
template <typename Number>
using Vector_t = Eigen::Matrix<Number, Eigen::Dynamic, 1>;

/// Wrapper for an Eigen::Matrix type.
template <typename Number>
using Matrix_t = Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>;

template <typename Number>
using Vector_set_t = std::set<Vector_t<Number>>;

struct IntervalMatrix
{
    Matrix_t<double> inf;
    Matrix_t<double> sup;
};

struct Interval
{
    double interval[2];
};

struct linerror_p_type
{
    Vector_t<double> x;
    double u;
};

struct linerror_type
{
    Vector_t<double> f0;
    linerror_p_type p;
};
        
enum class REACHABILITY_RESULT {
    SAFE,
    UNKNOWN
};



}  // namespace reachSolver