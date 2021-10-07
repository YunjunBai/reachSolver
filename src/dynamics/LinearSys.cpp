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

template class LinearSys<double>;

/*****************************************************************************
*                                                                           *
*                           Constructors and Destructors                    *
*                                                                           *
*****************************************************************************/
template <typename Number>
LinearSys<Number>::LinearSys(Matrix_t<Number>& A, Matrix_t<Number>& B)
    :ContDynamics<Number>(std::string("linearSys"), A.rows(), B.cols(), 1),
    A_(A), B_(B), c_(0, 0), C_(Eigen::MatrixXd::Identity(1,1)), D_(0, 0), k_(0, 0){ 
   
}

template <typename Number>
LinearSys<Number>::LinearSys(Matrix_t<Number>& A, Matrix_t<Number>& B, Matrix_t<Number>& c)
    :ContDynamics<Number>(std::string("linearSys"), A.rows(), B.cols(), 1),
    A_(A), B_(B), c_(c), C_(Eigen::MatrixXd::Identity(1,1)), D_(0, 0), k_(0, 0){ 
   
}

template <typename Number>
LinearSys<Number>::LinearSys(Matrix_t<Number>& A, Matrix_t<Number>& B, Matrix_t<Number>& c, Matrix_t<Number>& C)
    :ContDynamics<Number>(std::string("linearSys"), A.rows(), B.cols(), C.rows()),
    A_(A), B_(B), c_(c), C_(C), D_(0, 0), k_(0, 0){ 
   
}

template <typename Number>
LinearSys<Number>::LinearSys(Matrix_t<Number>& A, Matrix_t<Number>& B, Matrix_t<Number>& c, Matrix_t<Number>& C, Matrix_t<Number>& D)
    :ContDynamics<Number>(std::string("linearSys"), A.rows(), B.cols(), C.rows()),
    A_(A), B_(B), c_(c), C_(C), D_(D), k_(0, 0){ 
   
}

template <typename Number>
LinearSys<Number>::LinearSys(Matrix_t<Number>& A, Matrix_t<Number>& B, Matrix_t<Number>& c, Matrix_t<Number>& C, Matrix_t<Number>& D, Matrix_t<Number>& k)
    :ContDynamics<Number>(std::string("linearSys"), A.rows(), B.cols(), C.rows()),
    A_(A), B_(B), c_(c), C_(C), D_(D), k_(k){ 
   
}

template <typename Number>
LinearSys<Number>::LinearSys(std::string name, Matrix_t<Number>& A, Matrix_t<Number>& B)
    :ContDynamics<Number>(name, A.rows(), B.cols(), 1),
    A_(A), B_(B), c_(0, 0), C_(Eigen::MatrixXd::Identity(1,1)), D_(0, 0), k_(0, 0){ 
   
}

template <typename Number>
LinearSys<Number>::LinearSys(std::string name, Matrix_t<Number>& A, Matrix_t<Number>& B, Matrix_t<Number>& c)
    :ContDynamics<Number>(name, A.rows(), B.cols(), 1),
    A_(A), B_(B), c_(c), C_(Eigen::MatrixXd::Identity(1,1)), D_(0, 0), k_(0, 0){ 
   
}

template <typename Number>
LinearSys<Number>::LinearSys(std::string name, Matrix_t<Number>& A, Matrix_t<Number>& B, Matrix_t<Number>& c, Matrix_t<Number>& C)
    :ContDynamics<Number>(name, A.rows(), B.cols(), C.rows()),
    A_(A), B_(B), c_(c), C_(C), D_(0, 0), k_(0, 0){ 
   
}

template <typename Number>
LinearSys<Number>::LinearSys(std::string name, Matrix_t<Number>& A, Matrix_t<Number>& B, Matrix_t<Number>& c, Matrix_t<Number>& C, Matrix_t<Number>& D)
    :ContDynamics<Number>(name, A.rows(), B.cols(), C.rows()),
    A_(A), B_(B), c_(c), C_(C), D_(D), k_(0, 0){ 
   
}

template <typename Number>
LinearSys<Number>::LinearSys(std::string name, Matrix_t<Number>& A, Matrix_t<Number>& B, Matrix_t<Number>& c, Matrix_t<Number>& C, Matrix_t<Number>& D, Matrix_t<Number>& k)
    :ContDynamics<Number>(name, A.rows(), B.cols(), C.rows()),
    A_(A), B_(B), c_(c), C_(C), D_(D), k_(k){ 
   
}

template <typename Number>
LinearSys<Number>::~LinearSys(){}

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
    
    return LinearReachableSet<Number>();
}

template <typename Number>
void LinearSys<Number>::exponential(ReachOptions<Number>& options){
    // load data from object/options structure
    Matrix_t<Number> A = A_;
    Matrix_t<Number> A_abs = A_.cwiseAbs();
    size_t taylorTerms = options.taylor_terms();
    size_t dim_local = ContDynamics<Number>::dim();
    Number *factors = options.factor();
    
    // initialize
    std::vector<Matrix_t<Number>> Apower = std::vector<Matrix_t<Number>>(taylorTerms+1);
    Apower[0] = A;
    std::vector<Matrix_t<Number>> Apower_abs = std::vector<Matrix_t<Number>>(taylorTerms+1);
    Apower_abs[0] = A_abs;
    Matrix_t<Number> M = Eigen::MatrixXd::Identity(dim_local,dim_local);

    // compute powers for each term and sum of these
    for (size_t i = 0; i < taylorTerms; i++)
    {
        Apower[i+1] = Apower[i] * A;
        Apower_abs[i+1] = Apower_abs[i] * A_abs;
        M = M + Apower_abs[i] * factors[i];
    }
    // determine error due to finite Taylor series
    Matrix_t<Number> tmp = A_abs*options.time_step();
    Matrix_t<Number> tmp_exp = tmp.array().exp();
    Matrix_t<Number> W = tmp_exp - M;
    // compute absolute value of W for numerical stability
    W =  W.cwiseAbs();
    IntervalMatrix E;
    E.inf = -W;
    E.sup = W;

    // write to object structure
    taylor_.powers = Apower;
    taylor_.error = E;
}

template <typename Number>
void LinearSys<Number>::tie(ReachOptions<Number>& options){
    // load data from object/options structure
    std::vector<Matrix_t<Number>> Apower = taylor_.powers;
    size_t taylorTerms = options.taylor_terms();
    Number* rbyfac = options.factor();
    size_t dim_local = ContDynamics<Number>::dim();

    // initialize Asum
    Matrix_t<Number> Asum_pos = Eigen::MatrixXd::Zero(dim_local, dim_local);
    Matrix_t<Number> Asum_neg = Eigen::MatrixXd::Zero(dim_local, dim_local);
    
    for (size_t i = 2; i <= taylorTerms; i++)
    {
        // compute factor
        double exp1 = -i/(i-1);
        double exp2 = -1/(i-1);
        double factor = (std::pow(i,exp1)-std::pow(i,exp2))*rbyfac[i-1];
        
        // init Apos, Aneg
        Matrix_t<Number> Apos = Eigen::MatrixXd::Zero(dim_local, dim_local);
        Matrix_t<Number> Aneg = Eigen::MatrixXd::Zero(dim_local, dim_local);

        // obtain positive and negative parts
        for (size_t j = 0; j < dim_local; j++)
        {
            for (size_t k = 0; k < dim_local; k++)
            {
                if(Apower[i-1](j, k)>0){
                    Apos(j, k) = Apower[i-1](j, k);
                }else{
                    Aneg(j, k) = Apower[i-1](j, k);
                }
            }
        }
        // compute powers; factor is always negative
        Asum_pos = Asum_pos+factor*Aneg;
        Asum_neg = Asum_neg+factor*Apos;       
    }
    // instantiate interval matrix
    IntervalMatrix Asum;
    Asum.inf = Asum_neg;
    Asum.sup = Asum_pos;

    // write to object structure
    taylor_.F.inf = Asum.inf+taylor_.error.inf;
    taylor_.F.sup = Asum.sup+taylor_.error.sup;
}

template <typename Number>
void LinearSys<Number>::inputSolution(ReachOptions<Number>& options){
    // set of possible inputs
    Zonotope<Number> V = options.U()*B_;

    // compute vTrans 
    Zonotope<Number> vTrans = options.uTrans_lin()*B_;

    // consider constant input
    if (c_ != Matrix_t<Number>(0, 0)){
        vTrans = vTrans + c_;
    }

    Matrix_t<Number> A = A_;
    std::vector<Matrix_t<Number>> Apower = taylor_.powers;
    IntervalMatrix E = taylor_.error;
    size_t taylorTerms = options.taylor_terms();
    double r = options.time_step();
    size_t dim = A.size();
    Number *factors = options.factor();

    // init Vsum
    Zonotope<Number> Vsum = V*r;
    Matrix_t<Number> Asum = r*Eigen::MatrixXd::Identity(dim, dim);
    // compute higher order terms
    for (size_t i = 0; i < taylorTerms; i++)
    {
        Vsum = Vsum+V*(Apower[i]*factors[i+1]);
        Asum = Asum+Apower[i]*factors[i+1];
    }

    // compute overall solution
    Zonotope<Number> inputSolV = Vsum+V*E*r;

    // compute solution due to constant input
    IntervalMatrix eAtInt;
    eAtInt.inf = Asum+E.inf*r;
    eAtInt.sup = Asum+E.sup*r;
    Zonotope<Number> inputSolVtrans = vTrans*eAtInt;

    // compute additional uncertainty if origin is not contained in input set
    Zonotope<Number> inputCorr;
    if(options.originContained() != 0){
        Matrix_t<Number> inputCorr = Eigen::MatrixXd::Zero(dim,1);
    }else{
        // compute inputF
        inputTie(options);
        IntervalMatrix inputF = taylor_.inputF;
        inputCorr = vTrans*inputF;
    }

    // write to object structure
    taylor_.V = V;
    Matrix_t<Number> SolV_Z = Matrix_t<Number>(inputSolV.center().size(), Eigen::Dynamic);
    SolV_Z << inputSolV.center(), inputSolV.generators();
    if(SolV_Z.any()){
        taylor_.RV = inputSolV;
    }else{
        taylor_.RV = Zonotope<Number>(dim);
    }
    Matrix_t<Number> SolVtrans_Z = Matrix_t<Number>(inputSolVtrans.center().size(), Eigen::Dynamic);
    SolVtrans_Z << inputSolVtrans.center(), inputSolV.generators();
    if(SolVtrans_Z.any()){
        taylor_.Rtrans = inputSolVtrans;
    }else{
        taylor_.Rtrans = Zonotope<Number>(dim);
    }
    taylor_.inputCorr = inputCorr;
    taylor_.eAtInt = eAtInt;
}

template <typename Number>
void LinearSys<Number>::inputTie(ReachOptions<Number>& options){
    // load data from object structure
    std::vector<Matrix_t<Number>> Apower = taylor_.powers;
    IntervalMatrix E = taylor_.error;
    size_t taylorTerms = options.taylor_terms();
    double r = options.time_step();
    size_t dim = A_.size();

    // initialize Asum
    Matrix_t<Number> Asum_pos = Eigen::MatrixXd::Zero(dim, dim);
    Matrix_t<Number> Asum_neg = Eigen::MatrixXd::Zero(dim, dim);
    
    for (size_t i = 2; i <= taylorTerms+1; i++)
    {
        // compute factor
        double exp1 = -i/(i-1);
        double exp2 = -1/(i-1);
        double factor = (std::pow(i,exp1)-std::pow(i,exp2))*options.factor()[i-1];
        
        // init Apos, Aneg
        Matrix_t<Number> Apos = Eigen::MatrixXd::Zero(dim, dim);
        Matrix_t<Number> Aneg = Eigen::MatrixXd::Zero(dim, dim);

        // obtain positive and negative parts
        for (size_t j = 0; j < dim; j++)
        {
            for (size_t k = 0; k < dim; k++)
            {
                if(Apower[i-2](j, k)>0){
                    Apos(j, k) = Apower[i-2](j, k);
                }else{
                    Aneg(j, k) = Apower[i-2](j, k);
                }
            }
        }
        // compute powers; factor is always negative
        Asum_pos = Asum_pos+factor*Aneg;
        Asum_neg = Asum_neg+factor*Apos;       
    }

    // instantiate interval matrix
    IntervalMatrix Asum;
    Asum.inf = Asum_neg;
    Asum.sup = Asum_pos;

    // compute error due to finite Taylor series according to internal document "Input Error Bounds in Reachability Analysis"
    IntervalMatrix Einput;
    Einput.inf = Asum_neg*r;
    Einput.sup = Asum_pos*r;

    // write to object structure
    taylor_.inputF.inf = Asum.inf+Einput.inf;
    taylor_.inputF.sup = Asum.sup+Einput.sup;
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
    Matrix_t<Number> eAt = (A_*options.time_step()).array().exp();
    // save data to object structure
    taylor_.eAt = eAt;

    IntervalMatrix F = taylor_.F;
    Zonotope<Number> RV = taylor_.RV;
    Zonotope<Number> inputCorr = taylor_.inputCorr;
    Zonotope<Number> Rtrans = taylor_.Rtrans;

    // first time step homogeneous solution
    Zonotope<Number> Rhom_tp = Rinit*eAt + Rtrans;
    Zonotope<Number> Rhom = Rinit.enclose(Rhom_tp)+Rinit*F+inputCorr;

    // reduce zonotope
    Rhom.Reduce(options.zonotope_order());
    Rhom_tp.Reduce(options.zonotope_order());
    RV.Reduce(options.zonotope_order());

    // save homogeneous and particulate solution
    options.set_Rhom(Rhom);
    options.set_Rhom_tp(Rhom_tp);
    options.set_Raux(RV);
    // if (strcmp(options.linAlg,'wrapping-free') == 0){
    //     options.Rpar=interval(RV);
    // }else{
        options.set_Rpar(RV);
    // }
    options.set_Rtrans(taylor_.Rtrans);

    Zonotope<Number> Rtotal;
    Zonotope<Number> Rtotal_tp;
    // total solution
    // if isa(Rinit,'mptPolytope'){
    //     // convert zonotopes to polytopes
    //     Radd=mptPolytope(RV);
    //     Rtotal=Rhom+Radd;
    //     Rtotal_tp=Rhom_tp+Radd;
    // }else{
        // original computation
        Rtotal = Rhom+RV;
        Rtotal_tp = Rhom_tp+RV;
    // }

    // write results to reachable set struct Rfirst
    LinearReachableSet<Number> Rfirst = LinearReachableSet<Number>();
    Rfirst.set_time_point(Rtotal_tp);
    Rfirst.set_time_interval(Rtotal);

    return Rfirst;
}

template <typename Number>
Zonotope<Number> LinearSys<Number>::error_solution(ReachOptions<Number>& options, Zonotope<Number> Vdyn, Zonotope<Number> Vstat){
    Zonotope<Number> errorStat;
    if(Vstat.Empty()){
        errorStat = Zonotope<Number>();
    }else{
        errorStat = Vstat*taylor_.eAtInt;
    }
    // load data from object/options structure
    std::vector<Matrix_t<Number>> Apower = taylor_.powers;
    IntervalMatrix E = taylor_.error;
    size_t taylorTerms = options.taylor_terms();
    double r = options.time_step();
    Number *factors = options.factor();

    // initialize Asum
    Zonotope<Number> Asum=Vdyn*r;
    for (size_t i = 0; i < taylorTerms; i++)
    {
        // compute powers
        Zonotope<Number> ApowerV = Vdyn*(factors[i+1]*Apower[i]);
        // compute sums
        Asum = Asum+ApowerV;
    }
    // get error due to finite Taylor series
    Zonotope<Number> F;
    // if(isa(Vdyn,'zonoBundle')){
    //     F=E*Vdyn.Z{1}*r;
    // }else{
        F=Vdyn*E*r;
    // }
    // Compute error solution (dyn. + stat.)
    Zonotope<Number> Rerror = Asum+F+errorStat;
    return Rerror;    
}

}