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
    int res = 1;
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
    std::vector<Number> tVec = std::vector<Number>(tmp_t_count+1);
    for (size_t i = 0; i <= tmp_t_count; i++)
    {
        tVec[i] = options.tStart()+options.time_step()*i;
    }    

    // initialize cell-arrays that store the reachable set
    std::vector<TimeInt> time_int = std::vector<TimeInt>(tmp_t_count);
    std::vector<TimePoint> time_point = std::vector<TimePoint>(tmp_t_count);

    // initialize reachable set computations
    ReachableSet<Number> Rnext = ReachableSet();
    try
    {
        Rnext = initReach(options.R0(), options);
    }catch(SetExplosionException e1 )
    {
        std::cout << e.what() << std::endl;
        return 0;
    }

    // loop over all reachability steps
    Number t;
    for (size_t i = 1; i < tVec.size(); i++)
    {
        time_int.set_set_rs(i-1, Rnext.time_interval());
        time_int.set_time(i-1, tVec[i-1], tVec[i]);
        time_point.set_set_rs(i-1, Rnext.time_point());
        time_point.set_time(i-1, tVec[i]);

        // check specificaiton
        // if(!spec.empty()){
        //     if(!spec.check(Rnext.ti)){
        //         R = createReachSetObject(time_int, time_point);
        //         return 0;
        //     }
        // }

        // increment time and set counter
        t = tVec[i];
        options.set_t(t);
        if(options.verbose() == 1){
            std::cout << t << std::endl;
        }
        
        //  if a trajectory should be tracked
        // if (options.uTransVec != NULL)
        //     options.uTrans = options.uTransVec(:,i);
        // end

        // compute next reachable set
        try
        {
            Rnext = post(Rnext, options);
        }catch(SetExplosionException e1 )
        {
            R = createReachSetObject(time_int, time_point);
            std::cout << e.what() << std::endl;
            return 0;
        }
    }
    // check specificaiton
    // if(!spec.empty()){
    //     if(!spec.check(Rnext.ti)){
    //         res = 0;
    //     }
    // }

    time_int.set_set_rs(tmp_t_count-1, Rnext.time_interval());
    time_int.set_time(tmp_t_count-1, tVec[tmp_t_count-1], tVec[tmp_t_count]);
    time_point.set_set_rs(tmp_t_count-1, Rnext.time_point());
    time_point.set_time(tmp_t_count-1, tVec[tmp_t_count]);

    R = createReachSetObject(timeInt,timePoint);

    return res;
    
} // NonlinearSys<Number>::reach

template <typename Number>
ReachOptions NonlinearSys<Number>::checkOptionsReach(ReachSpecification& options, int hyb){
    
}

template <typename Number>
void NonlinearSys<Number>::derivatives(ReachSpecification& options){
    
}

template <typename Number>
ReachableSet<Number> NonlinearSys<Number>::initReach(Zonotope<Number>& Rinit, ReachSpecification& options){
    
    ReachableSet<Number> Rnext;
    // compute reachable set using the options.alg = 'linRem' algorithm
    if(strcmp(options.alg(), "linRem") == 0){
        Rnext = initReach_linRem(Rinit, options);
    }

    // loop over all parallel sets
    int setCounter = 1;
    ReachableSet<Number> Rtp = ReachableSet<Number>();
    ReachableSet<Number> Rti = ReachableSet<Number>();
    ReachableSet<Number> R0 = ReachableSet<Number>();
    // for split to loop?
    for (size_t i = 0; i < Rinit.length(); i++)
    {
        // compute reachable set of abstraction
        int dimForSplit = linReach(options, Rinit, Rti, Rtp);

        // check if initial set has to be split
        if(dimForSplit == 0){

        }

    }
    

}

template <typename Number>
ReachableSet<Number> NonlinearSys<Number>::post(ReachableSet<Number>& R, ReachSpecification& options){
    // In contrast to the linear system: the nonlinear system has to be constantly initialized due to the linearization procedure
    ReachableSet<Number> Rnext;
    Rnext = initReach(R.time_point(), options);

    // reduce zonotopes
    for (size_t i = 0; i < Rnext.time_point().length(); i++)
    {
        
    }

    // delete redundant reachable sets
    Rnext = Rnext.deleteRedundantSets(R, options);
    return Rnext;
    
}

template <typename Number>
int NonlinearSys<Number>::linReach(ReachSpecification& options, ReachableSet<Number>& Rstart, ReachableSet<Number>& Rti, ReachableSet<Number>& Rtp){
    
}


} //namespace reachSolver

