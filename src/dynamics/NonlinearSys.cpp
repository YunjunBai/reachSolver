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
ReachableSet<Number> NonlinearSys<Number>::createReachSetObject(TimeInt<Number>& time_int, TimePoint<Number>& time_point){
    
}

template <typename Number>
int NonlinearSys<Number>::reach(ReachOptions<Number>& options, ReachSpecification& spec, ReachableSet<Number> & R){
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
    TimeInt<Number> time_int = TimeInt<Number>(tmp_t_count);
    TimePoint<Number> time_point = TimePoint<Number>(tmp_t_count);

    // initialize reachable set computations
    ReachableSetElement<Number> Rinit_element = ReachableSetElement<Number>(options.R0(), Vector_t<Number>(options.maxError));
    std::vector<ReachableSetElement<Number>> Rinit = {Rinit_element};
    ReachableSet<Number> Rnext = ReachableSet();
    try
    {
        Rnext = initReach(Rinit, options);
    }catch(SetExplosionException e1 )
    {
        std::cout << e.what() << std::endl;
        return 0;
    }

    // loop over all reachability steps
    Number t;
    for (size_t i = 1; i < tVec.size(); i++)
    {
        time_int.set_rs(i-1, Rnext.time_interval());
        time_int.set_time(i-1, tVec[i-1], tVec[i]);
        time_point.set_rs(i-1, Rnext.time_point());
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

    time_int.set_rs(tmp_t_count-1, Rnext.time_interval());
    time_int.set_time(tmp_t_count-1, tVec[tmp_t_count-1], tVec[tmp_t_count]);
    time_point.set_rs(tmp_t_count-1, Rnext.time_point());
    time_point.set_time(tmp_t_count-1, tVec[tmp_t_count]);

    R = createReachSetObject(timeInt,timePoint);

    return res;
    
} // NonlinearSys<Number>::reach

template <typename Number>
ReachOptions<Number> NonlinearSys<Number>::checkOptionsReach(ReachOptions<Number>& options, int hyb){
    
}

template <typename Number>
void NonlinearSys<Number>::derivatives(ReachOptions<Number>& options){
    
}

template <typename Number>
ReachableSet<Number> NonlinearSys<Number>::initReach_linRem(std::vector<ReachableSetElement<Number>>& Rinit, ReachOptions<Number>& options){

}

template <typename Number>
std::vector<Zonotope<Number>> NonlinearSys<Number>::split(Zonotope<Number> input, int number){

}

template <typename Number>
ReachableSet<Number> NonlinearSys<Number>::initReach(std::vector<ReachableSetElement<Number>>& Rinit, ReachOptions<Number>& options){
    
    ReachableSet<Number> Rnext = ReachableSet<Number>();
    // compute reachable set using the options.alg = 'linRem' algorithm
    if(strcmp(options.alg(), "linRem") == 0){
        Rnext = initReach_linRem(Rinit, options);
    }

    // loop over all parallel sets
    int setCounter = 1;
    std::vector<ReachableSetElement<Number>> Rtp = std::vector<ReachableSetElement<Number>>();
    std::vector<ReachableSetElement<Number>> Rti = std::vector<ReachableSetElement<Number>>();
    std::vector<ReachableSetElement<Number>> R0 = std::vector<ReachableSetElement<Number>>();
    Zonotope<Number> Rtemp_tp = Zonotope<Number>();
    Zonotope<Number> Rtemp_ti = Zonotope<Number>();
    for (size_t i = 0; i < Rinit.size(); i++)
    {
        // compute reachable set of abstraction
        int dimForSplit = linReach(options, Rinit[i], Rtemp_ti, Rtemp_tp);

        // check if initial set has to be split
        if(dimForSplit == 0){
            Rtp[setCounter].set_rs(Rtemp_tp);
            Rtp[setCounter].set_prev(i);
            Rti[setCounter].set_rs(Rtemp_ti);
            R0[setCounter] = Rinit{i};
            setCounter = setCounter + 1;
        }else{
            std::cout << "split! ...number of parallel sets: " << Rinit.size()+1 << std::endl;

            // split the initial set
            std::vector<Zonotope<Number>> Rtmp = std::vector<Zonotope<Number>>(2);
            Rtmp = split(Rinit[i].rs(),dimForSplit);
            std::vector<ReachableSetElement<Number>> Rsplit = std::vector<ReachableSetElement<Number>>(2);
            Rsplit[0].set_rs(Rtmp[0]);
            Rsplit[1].set_rs(Rtmp[1]);

            // reset the linearization error
            Rsplit[0].set_error(Vector_t<Number>(options.maxError));
            Rsplit[1].set_error(Vector_t<Number>(options.maxError));

            // compute the reachable set for the splitted sets
            ReachableSet<Number> Rres = ReachableSet<Number>();
            Rres = initReach(Rsplit,options);

            for (size_t j = 0; i < Rres.time_point().size(); j++)
            {
                Rtp[setCounter] = Rres.time_point()[j];
                Rtp[setCounter].set_parent(i);
                Rti[setCounter] = Rres.time_interval()[j];
                R0[setCounter] = Rres.R0()[j];
                setCounter = setCounter + 1;
            } // for j
        } // else
    } // for i
    
    Rnext.set_time_point(Rtp);
    Rnext.set_time_interval(Rti);
    Rnext.set_R0(R0);
    return Rnext;
} // initReach

template <typename Number>
ReachableSet<Number> NonlinearSys<Number>::post(ReachableSet<Number>& R, ReachOptions<Number>& options){
    // In contrast to the linear system: the nonlinear system has to be constantly initialized due to the linearization procedure
    ReachableSet<Number> Rnext;
    Rnext = initReach(R.time_point(), options);

    // reduce zonotopes
    for (size_t i = 0; i < Rnext.time_point().size(); i++)
    {
        if(!Rnext.time_point()[i].rs().empty()){
            Rnext.time_point()[i].rs().reduce(options.zonotopeOrder);
            Rnext.time_interval()[i].rs().reduce(options.zonotopeOrder);
        }
    }

    // delete redundant reachable sets
    Rnext = Rnext.deleteRedundantSets(R, options);
    return Rnext;
}

template <typename Number>
LinearSys<Number> NonlinearSys<Number>::linearize(ReachOptions<Number>& options, Zonotope<Number>& R, ReachOptions<Number>& linOptions){

}

template <typename Number>
double NonlinearSys<Number>::linReach_linRem(ReachableSet<Number>& R, Zonotope<Number>& Rinit, Zonotope<Number>& Rdelta, ReachOptions<Number>& options, ReachableSetElement<Number>& Rti, ReachableSetElement<Number>& Rtp){

}

template <typename Number>
int NonlinearSys<Number>::linReach(ReachOptions<Number>& options, ReachableSetElement<Number>& Rstart, ReachableSetElement<Number>& Rti, ReachableSetElement<Number>& Rtp){
    // extract initial set and abstraction error
    Zonotope<Number> Rinit = Rstart.rs();
    Vector_t<Number> abstrerr = Rstart.error();
    Zonotope<Number>& Rti_internal = Zonotope<Number>();
    Zonotope<Number>& Rtp_internal = Zonotope<Number>();

    // linearize the nonlinear system
    ReachOptions<Number> linOptions = ReachOptions<Number>();
    LinearSys<Number> linSys = linearize(options, Rinit, linOptions); 

    // translate Rinit by linearization point
    Zonotope<Number> Rdelta = Rinit.Plus(-linerror_.p.x);

    // compute reachable set of the linearized system
    LinearReachableSet<Number> R = linSys.initReach(Rdelta, linOptions);
    
    // compute reachable set of the abstracted system including the abstraction error using the selected algorithm
    double perfInd;
    if(strcmp(options.alg(), "linRem") == 0){
        perfInd = linReach_linRem(R, Rinit, Rdelta, options, Rti_internal, Rtp_internal);
    }else{
        // loop until the actual abstraction error is smaller than the estimated linearization error
        Rtp_internal = R.time_point();
        Rti_internal = R.time_interval();
        double perfIndCurr = DBL_MAX;
        Zonotope<Number> VerrorDyn;
        Zonotope<Number> VerrorStat;
        perfInd = 0;
        while (perfIndCurr > 1 && perfInd <= 1)
        {
            //  estimate the abstraction error 
            Vector_t<Number> appliedError = abstrerr * 1.1;
            Zonotope<Number> Verror = Zonotope<Number>(appliedError*0, Eigen::DiagonalMatrix<double, appliedError.size()>(appliedError));
            Zonotope<Number> RallError = linSys.error_solution(options, Verror, NULL);

            // compute the abstraction error using the conservative linearization approach described in [1]
            Vector_t<Number> trueError;
            if (strcmp(options.alg(),'lin') == 0){
                // compute overall reachable set including linearization error
                Zonotope<Number> Rmax = Rti_internal + RallError;
                // compute linearization error
                VerrorDyn = Zonotope<Number>();
                trueError = abstrerr_lin(options, Rmax, VerrorDyn);
                VerrorStat = Zonotope<Number>();
            // }else{
            //     // compute overall reachable set including linearization error
            //     Zonotope<Number> Rmax = Rdelta + RallError;
            //     // compute abstraction error
            }
            // compare linearization error with the maximum allowed error
            perfIndCurr = (trueError/appliedError).maxCoeff();    
            perfInd = (trueError/options.maxError).maxCoeff();
            abstrerr = trueError;

            // if any(abstrerr > 1e+100)
            //     throw(SetExplosionException());

        } // while

        // translate reachable sets by linearization point
        Rti_internal = Rti_internal + linError_.p.x;
        Rtp_internal = Rtp_internal + linError_.p.x;

        // compute the reachable set due to the linearization error
        Zonotope<Number> Rerror = linSys.errorSolution(options, VerrorDyn, VerrorStat);
        
        // add the abstraction error to the reachable sets
        Rti_internal = Rti + Rerror;
        Rtp_internal = Rtp_internal + Rerror;
    } // else
    
    // determine the best dimension to split the set in order to reduce the linearization error
    // dimForSplit = [];
    // if (perfInd > 1){
    //     dimForSplit = select(obj,options,Rstart);
    // }
    int dimForSplit = 0;

    // store the linearization error
    Rtp.set_rs(Rtp_internal);
    Rtp.set_error(abstrerr);
    Rti.set_rs(Rti_internal);

    return dimForSplit;
} // linReach


} //namespace reachSolver

