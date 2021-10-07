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
ContDynamics<Number>::ContDynamics(std::string name){
    name_ = name;
}

template <typename Number>
ContDynamics<Number>::ContDynamics(std::string name, size_t dim, size_t num_inputs, size_t num_outputs){
    name_ = name;
    dim_ = dim;
    num_inputs_ = num_inputs;
    num_outputs_ = num_outputs;
}

template <typename Number>
ContDynamics<Number>::~ContDynamics(){}

template <typename Number>
const std::string ContDynamics<Number>::name() const{
    return name_;
}

template <typename Number>
const size_t ContDynamics<Number>::num_states() const{
    return dim_;
}

template <typename Number>
const size_t ContDynamics<Number>::dim() const{
    return dim_;
}

template <typename Number>
const size_t ContDynamics<Number>::num_inputs() const{
    return num_inputs_;
}


} //namespace reachSolver