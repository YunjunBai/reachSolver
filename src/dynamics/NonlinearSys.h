/**
 * @file   NonlinearSys.h
 * @brief  NonlinearSys class representation
 * @author Yunjun Bai
 * @date April 2021
 * @version 1.0 
 * 
 * Reference:
 * CORA ../contDynamics/nonlinearSys/all
 */

#pragma once

#include "ContDynamics.h"
#include "LinearSys.h"
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
// #include "derivatives.h"
namespace reachSolver{


/**
 * @class      NonlinearSys 
 * @brief      Class for NonlinearSys.
 * @ingroup    structure
 * @{
 */
template <typename Number>
class NonlinearSys : private ContDynamics<Number>{
  private:
    typedef Vector_t<Number> (*function_type)(Vector_t<Number> param1, Vector_t<Number> param2);
     
    Vector_t<Number> (*mFile_)(Vector_t<Number> vector1, Vector_t<Number> vector2);
    // autodiff::dual2nd (*mFile_f[6])(autodiff::ArrayXdual2nd& vector1, autodiff::ArrayXdual2nd& vector2);
    Number (*mFile_f_[6])(Vector_t<Number> vector1, Vector_t<Number> vector2);
    void jacobian_(Vector_t<Number> vector1, Vector_t<Number> vector2, Matrix_t<Number> matrix1, Matrix_t<Number> matrix2);
    std::vector<IntervalMatrix> hessian_(IntervalMatrix im1, IntervalMatrix im2);
    function_type thirdOrderTensor_;
    function_type tensors_;

    linerror_type linerror_;

    // used in constructor
    std::vector<size_t> numberOfInputs(function_type fun_handle, size_t inpArgs);

    // used in reach: create an object of class reachSet that stores the reachable set
    ReachableSet<Number> createReachSetObject(TimeInt<Number>& time_int, TimePoint<Number>& time_point);

    // used in initReach
    ReachableSet<Number> initReach_linRem(std::vector<ReachableSetElement<Number>>& Rinit, ReachOptions<Number>& options);

    // Splits a zonotope bundle into two zonotope bundles. This is done for one or every generator resulting in n possible splits where n is the system dimension; it is also possible to use a splitting hyperplane
    std::vector<Zonotope<Number>> split(Zonotope<Number> input, int number);

    // used in linReach
    LinearSys<Number> linearize(ReachOptions<Number>& options, Zonotope<Number>& R, ReachOptions<Number>& linOptions);

    double linReach_linRem(LinearReachableSet<Number>& R, Zonotope<Number>& Rinit, Zonotope<Number>& Rdelta, ReachOptions<Number>& options, Zonotope<Number>& Rti, Zonotope<Number>& Rtp);

  public:
    /****************************************************************************
    *                                                                           *
    *                           Constructors and Destructors                    *
    *                                                                           *
    *****************************************************************************/

    /**
     * @brief Constructor with no params
     */
    //NonlinearSys();

    /**
     * @brief Constructor with function handle to the dynamic equation
     * @param fun_handle function handle to the dynamic equation
     */
    // explicit NonlinearSys(function_type fun_handle);

    /**
     * @brief Constructor with two params
     * @param name name of dynamics
     * @param fun_handle function handle to the dynamic equation
     */
    // NonlinearSys(std::string name, function_type fun_handle);

    /**
     * @brief Constructor with three params
     * @param fun_handle function handle to the dynamic equation
     * @param num_states number of states
     * @param num_inputs number of inputs
     */
    NonlinearSys(function_type fun_handle, size_t num_states, size_t num_inputs);

    /**
     * @brief Constructor with four params
     * @param name name of dynamics
     * @param fun_handle function handle to the dynamic equation
     * @param num_states number of states
     * @param num_inputs number of inputs
     */
    NonlinearSys(std::string name, function_type fun_handle, size_t num_states, size_t num_inputs);

    /**
     * @brief Copy Constructor - constructs a NonlinearSys from an existing one.
     * @param other Another NonlinearSys, from which a new NonlinearSys is constructed
     */
    NonlinearSys(const NonlinearSys& other) = default;

    virtual ~NonlinearSys();

    /*****************************************************************************
    *                                                                           *
    *                       Public Functions on Properties                      *
    *                                                                           *
    *****************************************************************************/

    /**
     * @brief Get the mFile handle
     * @return mFile handle
     */
    // const function_type mFile() const;

    /**
     * @brief Get the jacobian
     * @return jacobian handle
     */
    // const function_type jacobian() const;

    /**
     * @brief Get the hessian
     * @return hessian handle
     */
    // const function_type hessian() const;

    /**
     * @brief Get the thirdOrderTensor
     * @return thirdOrderTensor handle
     */
    // const function_type thirdOrderTensor() const;

    /**
     * @brief Get the tensors
     * @return tensors handle
     */
    // const function_type tensors() const;

    /**
     * @brief Get the linerror
     * @return linerror
     */
    const linerror_type linerror() const;

    /**
     * @brief Set the linerror
     * @param linerror
     */
    void set_linerror(linerror_type linerror) const;

    /*****************************************************************************
    *                                                                           *
    *                       Compute Reachable Set                               *
    *                                                                           *
    *****************************************************************************/

    /**
     * @brief computes the reachable continuous set for the entire time horizon of a continuous system
     * @param options options for the computation of reachable sets
     * @param spec object of class specification
     * @param R object of class reachSet storing the computed reachable set
     * @return 1 if specifications are satisfied, 0 if not
     */
    // int reach(ReachOptions<Number>& options, ReachSpecification& spec, ReachableSet<Number> & R);
    int reach(ReachOptions<Number>& options, ReachableSet<Number> & R);

    /**
     * @brief checks if all necessary options are there and have valid values
     * @param options options for nonlinearSys
     * @param hyb called from hybrid Automaton (0/1)
     * @return options for nonlinearSys
     */
    ReachOptions<Number> checkOptionsReach(ReachOptions<Number>& options, int hyb);

    /**
     * @brief computes multivariate derivatives (jacobians, hessians, etc.) of nonlinear systems in a symbolic way; the result is stored in m-files and passed by a handle
     * @param options options struct
     */
    void derivatives(ReachOptions<Number>& options);

    /**
     * @brief computes the reachable continuous set for the first time step
     * @param Rinit initial reachable set
     * @param options struct containing the algorithm settings
     * @return first reachable set
     */
    ReachableSet<Number> initReach(std::vector<ReachableSetElement<Number>> Rinit, ReachOptions<Number>& options);

    /**
     * @brief computes the reachable continuous set for one time step of a nonlinear system by overapproximative linearization
     * @param R reachable set of the previous time step
     * @param options options for the computation of the reachable set
     * @return reachable set of the next time step
     */
    ReachableSet<Number> post(ReachableSet<Number>& R, ReachOptions<Number>& options);

    /**
     * @brief computes the reachable set after linearization
     * @param options struct with algorithm settings
     * @param Rstart initial reachable set
     * @param Rti reachable set for time interval
     * @param Rtp reachable set for time point
     * @return dimForSplit - dimension that is split to reduce the lin. error
     */
    int linReach(ReachOptions<Number>& options, ReachableSetElement<Number>& Rstart, ReachableSetElement<Number>& Rti, ReachableSetElement<Number>& Rtp);

    /**
     * @brief computes the solution due to the linearization error
     * @param options - options struct
     * @param R - reachable set (time-interval solution from linearized system + estimated set of abstraction errors)
     * @param Verrordyn - abstraction error (zonotope)
     * @return trueError - abstraction error (interval)
     */
    Vector_t<Number> abstrerr_lin(ReachOptions<Number>& options, Zonotope<Number> R, Zonotope<Number> VerrorDyn);

};

/** @} */
}