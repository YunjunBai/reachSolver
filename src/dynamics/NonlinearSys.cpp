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
    // options preprocessing
    options = checkOptionsReach(options, 0);
    // compute symbolic derivatives
    derivatives(options);

    // obtain factors for initial state and input solution time step
    int r = options.time_step();
    int *factor = new int[options.taylor_terms()+1];
    for(size_t i = 1; i<=options.taylor_terms()+1; i++){
        factor[i] = std::pow(r,i)/std::tgamma(i);
    }
    options.set_factor(factor);
    options.set_t(options.tStart());
    
    // time period
    int tmp_t_count = (options.tFinal()-options.tStart())/options.timeStep();
    int *tmp_vec = new int[tmp_t_count];
    for (size_t i = 0; i <= tmp_t_count; i++)
    {
        tmp_vec[i] = options.tStart()+options.time_step()*i;
    }    
    Vector_t<Number> tVec(tmp_vec);

    // initialize cell-arrays that store the reachable set
    TimeInt time_int = TimeInt(tmp_t_count-1);
    TimePoint time_point = TimePoint(tmp_t_count-1);

    // initialize reachable set computations
    ReachableSet Rnext = ReachableSet();
    try
    {
        Rnext = initReach(options.R0(), options);
    }catch(SetExplosionException e1 )
    {
        std::cout << e.what() << std::endl;
        return -1;
    }catch(std::exception& e)
    {
        std::cout << "other Exception caught" << std::endl;
        return -1;
    }

    // loop over all reachability steps
    for (size_t i = 2; i < tVec.size(); i++)
    {
        time_int.set_set_rs(i-1, Rnext.time_interval());
        time_int.set_time(i-1, );
        time_point.set_set_rs(i-1, Rnext.time_point());
        time_point.set_time(i-1, tVec[i]);
    }
    
} // NonlinearSys<Number>::reach

} //namespace reachSolver

