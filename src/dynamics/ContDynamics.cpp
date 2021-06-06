/**
 * @file   ContDynamics.cpp
 * @brief  ContDynamics class
 * @author Yunjun Bai
 * @date April 2021
 * @version 1.0 
 * 
 * Reference:
 * CORA ../contDynamics/contDynamics/all
 */

#include "ContDynamics.h"

namespace reachSolver{

template class ContDynamics<double>;


template <typename Number>
const std::string ContDynamics<Number>::name() const{
    return name_;
}

template <typename Number>
const size_t ContDynamics<Number>::num_states() const{
    return num_states_;
}


template <typename Number>
const size_t ContDynamics<Number>::num_inputs() const{
    return num_inputs_;
}


} //namespace reachSolver