/**
 * @file   LinearSys.cpp
 * @brief  LinearSys class
 * @author Yunjun Bai
 * @date April 2021
 * @version 1.0 
 * 
 * Reference:
 * CORA ../contDynamics/linearSys/all
 */

#include "LinearSys.h"

namespace reachSolver{

/*****************************************************************************
*                                                                           *
*                           Constructors and Destructors                    *
*                                                                           *
*****************************************************************************/
template <typename Number>
LinearSys<Number>::LinearSys(){
}

/*****************************************************************************
*                                                                           *
*                       Compute Reachable Set                               *
*                                                                           *
*****************************************************************************/

template <typename Number>
LinearReachableSet<Number> LinearSys<Number>::initReach(Zonotope<Number>& Rinit, ReachOptions<Number>& options){
    if(options.usekrylovError() == 0){
       //compute in Krylov space
        return initReach_Krylov(Rinit,options);
    }else{
        return initReach_Euclidean(Rinit,options);
    }
}

template <typename Number>
LinearReachableSet<Number> LinearSys<Number>::initReach_Krylov(Zonotope<Number>& Rinit, ReachOptions<Number>& options){
    
    return new LinearReachableSet<Number>();
}

template <typename Number>
LinearReachableSet<Number> LinearSys<Number>::initReach_Euclidean(Zonotope<Number>& Rinit, ReachOptions<Number>& options){
    // compute exponential matrix
    exponential(options);
    // compute time interval error (tie)
    tie(options);
    // compute reachable set due to input
    inputSolution(options);
    // change the time step
    taylor_.timeStep = options.time_step();

    // compute reachable set of first time interval
    Matrix_t<Number> eAt = (A_*options.time_step()).exp();
    // save data to object structure
    taylor_.eAt = eAt;
    

    return new LinearReachableSet<Number>();
}

template <typename Number>
Zonotope<Number> LinearSys<Number>::error_solution(ReachOptions<Number>& options, Zonotope<Number> Vdyn, Zonotope<Number> Vstat){

}

}