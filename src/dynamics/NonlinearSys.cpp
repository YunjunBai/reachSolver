/**
 * @file   NonlinearSys.cpp
 * @brief  NonlinearSys class
 * @author Yunjun Bai
 * @date April 2021
 * @version 1.0 
 * 
 * Reference:
 * CORA ../contDynamics/nonlinearSys/all
 */

#include "NonlinearSys.h"

namespace reachSolver{

template class NonlinearSys<double>;

/*****************************************************************************
*                                                                           *
*                           Constructors and Destructors                    *
*                                                                           *
*****************************************************************************/
template <typename Number>
NonlinearSys<Number>::NonlinearSys(function_type fun_handle):
    mFile_(fun_handle),
    name_("nonlinearSys"){
    std::vector<size_t> tmp_obj = numberOfInputs(fun_handle, 2);
    num_states_ = tmp_obj[0];
    num_inputs_ = tmp_obj[1];   
    assert(fun_handle!=NULL);
}

template <typename Number>
NonlinearSys<Number>::NonlinearSys(std::string name, function_type fun_handle):
    mFile_(fun_handle),
    name_(name){
    std::vector<size_t> tmp_obj = numberOfInputs(fun_handle, 2);
    num_states_ = tmp_obj[0];
    num_inputs_ = tmp_obj[1];   
    assert(fun_handle!=NULL);
}

template <typename Number>
NonlinearSys<Number>::NonlinearSys(function_type fun_handle, size_t num_states, size_t num_inputs):
    mFile_(fun_handle),
    name_("nonlinearSys"),
    num_states_(num_states),
    num_inputs_(num_inputs){ 
    assert(fun_handle!=NULL);
   
}

template <typename Number>
NonlinearSys<Number>::NonlinearSys(std::string name, function_type fun_handle, size_t num_states, size_t num_inputs):
    mFile_(fun_handle),
    name_(name),
    num_states_(num_states),
    num_inputs_(num_inputs){
    assert(fun_handle!=NULL);
   
}

/*****************************************************************************
*                                                                           *
*                       Public Functions on Properties                      *
*                                                                           *
*****************************************************************************/

template <typename Number>
const std::string NonlinearSys<Number>::name() const{
    return name_;
}

template <typename Number>
const size_t NonlinearSys<Number>::num_states() const{
    return num_states_;
}


template <typename Number>
const size_t NonlinearSys<Number>::num_inputs() const{
    return num_inputs_;
}


template <typename Number>
const NonlinearSys<Number>::function_type NonlinearSys<Number>::mFile() const{
    return mFile_;
}


template <typename Number>
const NonlinearSys<Number>::function_type NonlinearSys<Number>::jacobian() const{
    return jacobian_;
}


template <typename Number>
const NonlinearSys<Number>::function_type NonlinearSys<Number>::hessian() const{
    return hessian_;
}


template <typename Number>
const NonlinearSys<Number>::function_type NonlinearSys<Number>::thirdOrderTensor() const{
    return thirdOrderTensor_;
}


template <typename Number>
const NonlinearSys<Number>::function_type NonlinearSys<Number>::tensors() const{
    return tensors_;
}


/*****************************************************************************
*                                                                           *
*                       Compute Reachable Set                               *
*                                                                           *
*****************************************************************************/

template <typename Number>
ReachableSet<Number> NonlinearSys<Number>::createReachSetObject(TimeInt& time_int, TimePoint& time_point){

}

template <typename Number>
int NonlinearSys<Number>::reach(ReachOptions& options, ReachSpecification& spec, ReachableSet<Number> & R){
    
}

}

