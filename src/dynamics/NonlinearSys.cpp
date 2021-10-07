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
// template <typename Number>
// NonlinearSys<Number>::NonlinearSys(function_type fun_handle)
//     :ContDynamics(new String("nonlinearSys"), 10, 1, 1)
//     :mFile_(fun_handle){
//     std::vector<size_t> tmp_obj = numberOfInputs(fun_handle, 2);
//     num_states_ = tmp_obj[0];
//     num_inputs_ = tmp_obj[1];   
//     assert(fun_handle!=NULL);
// }

// template <typename Number>
// NonlinearSys<Number>::NonlinearSys(std::string name, function_type fun_handle):
//     mFile_(fun_handle),
//     name_(name){
//     std::vector<size_t> tmp_obj = numberOfInputs(fun_handle, 2);
//     num_states_ = tmp_obj[0];
//     num_inputs_ = tmp_obj[1];   
//     assert(fun_handle!=NULL);
// }

template <typename Number>
NonlinearSys<Number>::NonlinearSys(function_type fun_handle, size_t num_states, size_t num_inputs)
    :ContDynamics<Number>(std::string("nonlinearSys"), num_states, num_inputs, 1),
    mFile_(fun_handle){ 
    assert(fun_handle!=NULL);
   
}

template <typename Number>
NonlinearSys<Number>::NonlinearSys(std::string name, function_type fun_handle, size_t num_states, size_t num_inputs)
    :ContDynamics<Number>(name, num_states, num_inputs, 1),
    mFile_(fun_handle){
    assert(fun_handle!=NULL);
   
}

template <typename Number>
NonlinearSys<Number>::~NonlinearSys(){}

/*****************************************************************************
*                                                                           *
*                       Public Functions on Properties                      *
*                                                                           *
*****************************************************************************/

// template <typename Number>
// const typename NonlinearSys<Number>::function_type NonlinearSys<Number>::mFile() const{
//     return mFile_;
// }


// template <typename Number>
// const typename NonlinearSys<Number>::function_type NonlinearSys<Number>::jacobian() const{
//     return jacobian_;
// }


// template <typename Number>
// const typename NonlinearSys<Number>::function_type NonlinearSys<Number>::hessian() const{
//     return hessian_;
// }


// template <typename Number>
// const typename NonlinearSys<Number>::function_type NonlinearSys<Number>::thirdOrderTensor() const{
//     return thirdOrderTensor_;
// }


// template <typename Number>
// const typename NonlinearSys<Number>::function_type NonlinearSys<Number>::tensors() const{
//     return tensors_;
// }


/*****************************************************************************
*                                                                           *
*                       Compute Reachable Set                               *
*                                                                           *
*****************************************************************************/

template <typename Number>
ReachableSet<Number> NonlinearSys<Number>::createReachSetObject(TimeInt<Number>& time_int, TimePoint<Number>& time_point){
    return ReachableSet<Number>();
}

template <typename Number>
ReachOptions<Number> NonlinearSys<Number>::checkOptionsReach(ReachOptions<Number>& options, int hyb){
    return options;
}

// template <typename Number>
// void NonlinearSys<Number>::derivatives(ReachOptions<Number>& options){
    
// }


// template <typename Number>
// int NonlinearSys<Number>::reach(ReachOptions<Number>& options, ReachSpecification& spec, ReachableSet<Number> & R){
 
template <typename Number>
int NonlinearSys<Number>::reach(ReachOptions<Number>& options, ReachableSet<Number> & R){
    int res = 1;
    // options preprocessing
    options = checkOptionsReach(options, 0);
    // compute symbolic derivatives
    // derivatives(options);

    // obtain factors for initial state and input solution time step
    int r = options.time_step();
    Number *factor = new double[options.taylor_terms()+1];
    for(size_t i = 1; i<=options.taylor_terms()+1; i++){
        factor[i] = std::pow(r,i)/std::tgamma(i);
    }
    options.set_factor(factor);
    options.set_t(options.tStart());
    
    // time period
    int tmp_t_count = (options.tFinal()-options.tStart())/options.time_step();
    std::vector<Number> tVec = std::vector<Number>(tmp_t_count+1);
    for (size_t i = 0; i <= tmp_t_count; i++)
    {
        tVec[i] = options.tStart()+options.time_step()*i;
    }    

    // initialize cell-arrays that store the reachable set
    TimeInt<Number> time_int = TimeInt<Number>(tmp_t_count);
    TimePoint<Number> time_point = TimePoint<Number>(tmp_t_count);

    // initialize reachable set computations
    ReachableSetElement<Number> Rinit_element = ReachableSetElement<Number>(options.R0(), Vector_t<Number>(options.max_error()));
    std::vector<ReachableSetElement<Number>> Rinit = {Rinit_element};
    ReachableSet<Number> Rnext = ReachableSet<Number>();
    try
    {
        Rnext = initReach(Rinit, options);
    }catch(SetExplosionException e1 )
    {
        std::cout << e1.what() << std::endl;
        return 0;
    }

    // loop over all reachability steps
    Number t;
    for (size_t i = 1; i < tVec.size(); i++)
    {
        time_int.set_rs(i-1, Rnext.time_interval());
        struct Interval tmp_interval;
        tmp_interval.interval[0] = tVec[i-1];
        tmp_interval.interval[1] = tVec[i];
        time_int.set_time(i-1, tmp_interval);
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
            std::cout << e1.what() << std::endl;
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
    struct Interval tmp_interval;
    tmp_interval.interval[0] = tVec[tmp_t_count-1];
    tmp_interval.interval[1] = tVec[tmp_t_count];
    time_int.set_time(tmp_t_count-1, tmp_interval);
    time_point.set_rs(tmp_t_count-1, Rnext.time_point());
    time_point.set_time(tmp_t_count-1, tVec[tmp_t_count]);

    R = createReachSetObject(time_int,time_point);

    return res;
    
} // NonlinearSys<Number>::reach

template <typename Number>
ReachableSet<Number> NonlinearSys<Number>::initReach_linRem(std::vector<ReachableSetElement<Number>>& Rinit, ReachOptions<Number>& options){
    return ReachableSet<Number> ();
}

// template <typename Number>
// std::vector<Zonotope<Number>> splitOneDim(Zonotope<Number> Zbundle, Matrix_t<double> leftLimit, Matrix_t<double> rightLimit, int dim){
//     // split limits for a given dimension
//     Matrix_t<double> leftLimitMod = leftLimit;
//     leftLimitMod[dim] = 0.5*(leftLimit[dim]+rightLimit[dim]);
//     Matrix_t<double>  rightLimitMod = rightLimit;
//     rightLimitMod[dim] = 0.5*(leftLimit[dim]+rightLimit[dim]);

//     // construct zonotopes which are the left and right boxes
//     IntervalMatrix IM_Zleft, IM_Zright;
//     IM_Zleft.inf = leftLimit; IM_Zleft.sup = rightLimitMod;
//     IM_Zright.inf = leftLimitMod; IM_Zright.sup = rightLimit;
//     Zonotope<Number> Zleft = Zonotope<Number>(IM_Zleft);
//     Zonotope<Number> Zright = Zonotope<Number>(IM_Zright);

//     // generate splitted zonotope bundles
//     std::vector<Zonotope<Number>> Zsplit = new std::vector<Zonotope<Number>>(2);
//     Zsplit[0] = Zbundle & Zleft;
//     Zsplit[1] = Zbundle & Zright; 
//     return Zsplit;
// }

// template <typename Number>
// std::vector<Zonotope<Number>> NonlinearSys<Number>::split(Zonotope<Number> Zbundle, int number){
//     std::vector<Zonotope<Number>> Zsplit;
//     // split given dimension
//     if(number > 0){
//         int N = number;
//         // obtain enclosing interval hull
//         IntervalMatrix IH = Zbundle.interval();
//         // obtain limits
//         Matrix_t<double> leftLimit = IH.inf;
//         Matrix_t<double> rightLimit = IH.sup;
//         // split one dimension
//         Zsplit = splitOneDim(Zbundle,leftLimit,rightLimit,N);
//     }else{
//         Zsplit = new std::vector<Zonotope<Number>>();
//     }
//     return Zsplit;
// }

template <typename Number>
ReachableSet<Number> NonlinearSys<Number>::initReach(std::vector<ReachableSetElement<Number>> Rinit, ReachOptions<Number>& options){
    
    ReachableSet<Number> Rnext = ReachableSet<Number>();
    // compute reachable set using the options.alg = 'linRem' algorithm
    if(options.alg() == std::string("linRem")){
        Rnext = initReach_linRem(Rinit, options);
    }

    // loop over all parallel sets
    int setCounter = 1;
    std::vector<ReachableSetElement<Number>> Rtp = std::vector<ReachableSetElement<Number>>();
    std::vector<ReachableSetElement<Number>> Rti = std::vector<ReachableSetElement<Number>>();
    std::vector<ReachableSetElement<Number>> R0 = std::vector<ReachableSetElement<Number>>();
    ReachableSetElement<Number> Rtemp_tp = ReachableSetElement<Number>();
    ReachableSetElement<Number> Rtemp_ti = ReachableSetElement<Number>();
    for (size_t i = 0; i < Rinit.size(); i++)
    {
        // compute reachable set of abstraction
        int dimForSplit = linReach(options, Rinit[i], Rtemp_ti, Rtemp_tp);

        // check if initial set has to be split
        if(dimForSplit == 0){
            Rtp[setCounter] = Rtemp_tp;
            Rtp[setCounter].set_prev(i);
            Rti[setCounter] = Rtemp_ti;
            R0[setCounter] = Rinit[i];
            setCounter = setCounter + 1;
        }else{
            // std::cout << "split! ...number of parallel sets: " << Rinit.size()+1 << std::endl;

            // // split the initial set
            // std::vector<Zonotope<Number>> Rtmp = std::vector<Zonotope<Number>>(2);
            // Rtmp = split(Rinit[i].rs(),dimForSplit);
            // std::vector<ReachableSetElement<Number>> Rsplit = std::vector<ReachableSetElement<Number>>(2);
            // Rsplit[0].set_rs(Rtmp[0]);
            // Rsplit[1].set_rs(Rtmp[1]);

            // // reset the linearization error
            // Rsplit[0].set_error(Vector_t<Number>(options.maxError));
            // Rsplit[1].set_error(Vector_t<Number>(options.maxError));

            // // compute the reachable set for the splitted sets
            // ReachableSet<Number> Rres = ReachableSet<Number>();
            // Rres = initReach(Rsplit,options);

            // for (size_t j = 0; i < Rres.time_point().size(); j++)
            // {
            //     Rtp[setCounter] = Rres.time_point()[j];
            //     Rtp[setCounter].set_parent(i);
            //     Rti[setCounter] = Rres.time_interval()[j];
            //     R0[setCounter] = Rres.R0()[j];
            //     setCounter = setCounter + 1;
            // } // for j
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
        if(!Rnext.time_point()[i].rs().Empty()){
            Rnext.time_point()[i].rs().Reduce(options.zonotope_order());
            Rnext.time_interval()[i].rs().Reduce(options.zonotope_order());
        }
    }

    // delete redundant reachable sets
    Rnext.deleteRedundantSets(R, options);
    return Rnext;
}

template <typename Number>
LinearSys<Number> NonlinearSys<Number>::linearize(ReachOptions<Number>& options, Zonotope<Number>& R, ReachOptions<Number>& linOptions){
    //linearization point p.u of the input is the center of the input set
    linerror_p_type p;
    p.u = options.uTrans_nonlin();

    //obtain linearization point
    // if(options.haslinearizationPoint() == true){
    //     p.x = options.linearizationPoint();
    // }else{
        // linearization point p.x of the state is the center of the last reachable set R translated by 0.5*delta_t*f0
        Vector_t<Number> tmp_vector1 = Vector_t<Number>(1);
        tmp_vector1 << p.u;
        Vector_t<Number> f0prev = mFile_(R.center(), tmp_vector1);
        try
        {
            p.x = R.center() + f0prev*0.5*options.time_step();
        }catch(const std::exception& e1 )
        {
            std::cout << "time step not yet created" << std::endl;
            p.x = R.center();
        }
    // }

    // substitute p into the system equation to obtain the constant input
    Vector_t<Number> tmp_vector2 = Vector_t<Number>(1);
    tmp_vector2 << p.u;
    Vector_t<Number> f0 = mFile_(p.x, tmp_vector2);
    
    // substitute p into the Jacobian with respect to x and u to obtain the system matrix A and the input matrix B
    Matrix_t<Number> A = Matrix_t<Number>();
    Matrix_t<Number> B = Matrix_t<Number>();
    jacobian(this->mFile_ ,p.x, tmp_vector2, A, B);
    Matrix_t<Number> A_lin = A;
    Matrix_t<Number> B_lin = B;
    linOptions = options;
    // if (strcmp(options.alg(),'linRem') == 0){
    //     // in order to compute dA,dB, we use the reachability set computed for one step in initReach
    //     [dA,dB] = lin_error2dAB(options.Ronestep,options.U,obj.hessian,p);
    //     A = matZonotope(A,{dA});
    //     B = matZonotope(B,{dB});
    //     linSys = linParamSys(A,1,'constParam');
    //     linOptions.compTimePoint = 1;
    // }else{
        //set up linearized system
    Matrix_t<Number> tmp_B = Matrix_t<Number>();
    tmp_B(0, 0) = 1;
    LinearSys<Number> linSys = LinearSys<Number>(std::string("linSys"), A, tmp_B); 
        //B=1 as input matrix encountered in uncertain inputs
    // }

    // set up options for linearized system
    linOptions.set_U((options.U()+options.uTrans_nonlin()-p.u)*B);
    Vector_t<Number> Ucenter = linOptions.U().center();
    linOptions.set_U(linOptions.U() - Ucenter);

    linOptions.set_uTrans_lin(Zonotope<Number>((f0 + Ucenter),Eigen::MatrixXd::Zero(f0.rows(),1)));
    //linOptions.uTrans = zonotope(f0 + Ucenter);
    linOptions.set_originContained(0);

    //save constant input
    linerror_.f0=f0;

    //save linearization point
    linerror_.p=p;

    return linSys;
}

template <typename Number>
double NonlinearSys<Number>::linReach_linRem(LinearReachableSet<Number>& R, Zonotope<Number>& Rinit, Zonotope<Number>& Rdelta, ReachOptions<Number>& options, Zonotope<Number>& Rti, Zonotope<Number>& Rtp){
    return 0;
}

template <typename Number>
int NonlinearSys<Number>::linReach(ReachOptions<Number>& options, ReachableSetElement<Number>& Rstart, ReachableSetElement<Number>& Rti, ReachableSetElement<Number>& Rtp){
    // extract initial set and abstraction error
    Zonotope<Number> Rinit = Rstart.rs();
    Vector_t<Number> abstrerr = Rstart.error();
    Zonotope<Number> Rti_internal = Zonotope<Number>();
    Zonotope<Number> Rtp_internal = Zonotope<Number>();

    // linearize the nonlinear system
    ReachOptions<Number> linOptions = ReachOptions<Number>();
    LinearSys<Number> linSys = linearize(options, Rinit, linOptions); 

    // translate Rinit by linearization point
    Zonotope<Number> Rdelta = Rinit - linerror_.p.x;

    // compute reachable set of the linearized system
    LinearReachableSet<Number> R = linSys.initReach(Rdelta, linOptions);
    
    // compute reachable set of the abstracted system including the abstraction error using the selected algorithm
    double perfInd;
    if(options.alg() == std::string("linRem")){
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
            Zonotope<Number> Verror = Zonotope<Number>(appliedError*0, appliedError.asDiagonal());
            Zonotope<Number> RallError = linSys.error_solution(options, Verror, Zonotope<Number>());

            // compute the abstraction error using the conservative linearization approach described in [1]
            Vector_t<Number> trueError;
            if (options.alg() == std::string("lin")){
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
            perfIndCurr = (trueError.array()/appliedError.array()).maxCoeff(); 
            Matrix_t<Number> tmp_maxerror = Matrix_t<Number>();
            tmp_maxerror(0, 0) = options.max_error();
            perfInd = (trueError.array()/tmp_maxerror.array()).maxCoeff();
            abstrerr = trueError;

            // if any(abstrerr > 1e+100)
            //     throw(SetExplosionException());

        } // while

        // translate reachable sets by linearization point
        Rti_internal = Rti_internal + linerror_.p.x;
        Rtp_internal = Rtp_internal + linerror_.p.x;

        // compute the reachable set due to the linearization error
        Zonotope<Number> Rerror = linSys.error_solution(options, VerrorDyn, VerrorStat);
        
        // add the abstraction error to the reachable sets
        Rti_internal = Rti_internal + Rerror;
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

template <typename Number>
Vector_t<Number> NonlinearSys<Number>::abstrerr_lin(ReachOptions<Number>& options, Zonotope<Number> R, Zonotope<Number> VerrorDyn){
    Vector_t<Number> trueError;
    // compute interval of reachable set
    IntervalMatrix IHx = R.interval();
    // compute intervals of total reachable set
    IntervalMatrix totalInt_x;
    totalInt_x.inf = IHx.inf+linerror_.p.x;
    totalInt_x.sup = IHx.sup+linerror_.p.x;
    
    // compute intervals of input
    IntervalMatrix IHu = options.U().interval();
    // translate intervals by linearization point
    IntervalMatrix totalInt_u;
    Matrix_t<Number> tmp_matrix = Matrix_t<Number>();
    tmp_matrix(0, 0) = linerror_.p.u;
    totalInt_u.inf = IHu.inf+tmp_matrix;
    totalInt_u.sup = IHu.sup+tmp_matrix;

    if(options.tensor_order() == 2){
        // obtain maximum absolute values within IHx, IHu
        Matrix_t<Number> dx = IHx.inf.cwiseAbs().cwiseMax(IHx.sup.cwiseAbs());
        Matrix_t<Number> du = IHu.inf.cwiseAbs().cwiseMax(IHu.sup.cwiseAbs());

        // evaluate the hessian matrix with the selected range-bounding technique
        std::vector<IntervalMatrix> H;
        // if(options.lagrangeRem().method()!='interval'){

        // }else{
        //     if(name_ == "nonlinParamSys"){

        //     }else{
                H = hessian(this->mFile_f_, totalInt_x, totalInt_u);
        //     }
        // }

        // calculate the Lagrange remainder (second-order error)
        Vector_t<Number> errorLagr = Eigen::MatrixXd::Zero(H.size(),1);
        Vector_t<Number> dz = Matrix_t<Number>();
        dz << dx,
              du;
        for (size_t i = 0; i < H.size(); i++)
        {
            Matrix_t<Number>  H__inf, H__sup, H_max;
            H__inf = H[i].inf.cwiseAbs();
            H__sup = H[i].sup.cwiseAbs();
            H_max = H__inf.cwiseMax(H__sup);
            errorLagr[i] = 0.5*dz.adjoint()*H_max*dz;
        }
        VerrorDyn = Zonotope<Number>(0*errorLagr, errorLagr.asDiagonal());
        trueError = errorLagr;
        return trueError;
    }else if(options.tensor_order() == 3){
        // TODO not complete

        // // obtain intervals and combined interval z
        // std::vector<IntervalMatrix> dz = new std::vector<IntervalMatrix>(IHx, IHu);
        // // reduce zonotope
        // Zonotope<Number> Rred = R:
        // Rred.Reduce(options.error_order());
        // // combined zonotope (states + input)
        // Z = cartProd(Rred,options.U);

        // // calculate hessian matrix
        // if (name_ == "nonlinParamSys"){
        //     H = obj.hessian(obj.linError.p.x,obj.linError.p.u,options.paramInt);
        // }else{
        //     H = obj.hessian(obj.linError.p.x,obj.linError.p.u);
        // }
        
        return trueError;
    }else{
        std::cout << "No abstraction error computation for chosen tensor order!" << std::endl;
        return trueError;
    }
}

} //namespace reachSolver

