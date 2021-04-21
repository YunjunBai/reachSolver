/**
 * @file   commonType.h
 * @brief 
 * @author 
 * @date
 * @version 1.0
 * 
 */

#pragma once

namespace reachSolver{

/// Wrapper for an Eigen::Matrix type with only one column.
template <typename Number>
using vector_t = Eigen::Matrix<Number, Eigen::Dynamic, 1>;

/// Wrapper for an Eigen::Matrix type.
template <typename Number>
using matrix_t = Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>;

template <typename Number>
using vector_set_t = std::set<vector_t<Number>>;

enum class REACHABILITY_RESULT {
	SAFE,
	UNKNOWN
};



}  // namespace reachSolver