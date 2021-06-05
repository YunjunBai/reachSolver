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

#include "BasicObject.h"
#include <vector> 
#include <string>
#include <cmath>
#include <exception>
#include "ReachOptions.h"
#include "ReachableSet.h"
#include "ReachSpecification.h"
namespace reachSolver{

struct SetExplosionException : public std::exception
{
    const char * what () const throw ()
    {
        return "Set Explosion Exception!";
    }
};

/**
 * @class      NonlinearSys 
 * @brief      Class for NonlinearSys.
 * @ingroup    structure
 * @{
 */
template <typename Number>
class NonlinearSys{
  private:
    typedef Number (*function_type)(Number param1, Number param2);

    std::string name_;
    size_t num_states_;
    size_t num_inputs_;
    //size_t num_outputs_;

    function_type mFile_;
    function_type jacobian_;
    function_type hessian_;
    function_type thirdOrderTensor_;
    function_type tensors_;

    ReachableSet<Number> createReachSetObject(TimeInt& time_int, TimePoint& time_point);

    std::vector<size_t> numberOfInputs(function_type fun_handle, size_t inpArgs);

    ReachableSet<Number> initReach_linRem(Zonotope<Number>& Rinit, ReachSpecification& options);

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
    explicit NonlinearSys(function_type fun_handle);

    /**
     * @brief Constructor with two params
     * @param name name of dynamics
     * @param fun_handle function handle to the dynamic equation
     */
    NonlinearSys(std::string name, function_type fun_handle);

    /**
     * @brief Constructor with three params
     * @param fun_handle number of inputs
     * @param num_states number of states
     * @param num_inputs function handle to the dynamic equation
     */
    NonlinearSys(function_type fun_handle, size_t num_states, size_t num_inputs);

    /**
     * @brief Constructor with four params
     * @param name name of dynamics
     * @param fun_handle function handle to the dynamic equation
     * @param num_states number of states
     * @param num_inputs function handle to the dynamic equation
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
    const std::string name() const;

    const size_t num_states() const;

    const size_t num_inputs() const;

    /**
     * @brief Get the mFile handle
     * @return mFile handle
     */
    const function_type mFile() const;

    /**
     * @brief Get the jacobian
     * @return jacobian handle
     */
    const function_type jacobian() const;

    /**
     * @brief Get the hessian
     * @return hessian handle
     */
    const function_type hessian() const;

    /**
     * @brief Get the thirdOrderTensor
     * @return thirdOrderTensor handle
     */
    const function_type thirdOrderTensor() const;

    /**
     * @brief Get the tensors
     * @return tensors handle
     */
    const function_type tensors() const;

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
    int reach(ReachOptions& options, ReachSpecification& spec, ReachableSet<Number> & R);

    /**
     * @brief checks if all necessary options are there and have valid values
     * @param options options for nonlinearSys
     * @param hyb called from hybrid Automaton (0/1)
     * @return options for nonlinearSys
     */
    ReachOptions checkOptionsReach(ReachSpecification& options, int hyb);

    /**
     * @brief computes multivariate derivatives (jacobians, hessians, etc.) of nonlinear systems in a symbolic way; the result is stored in m-files and passed by a handle
     * @param options options struct
     */
    void derivatives(ReachSpecification& options);

    /**
     * @brief computes the reachable continuous set for the first time step
     * @param Rinit initial reachable set
     * @param options struct containing the algorithm settings
     * @return first reachable set
     */
    ReachableSet<Number> initReach(Zonotope<Number>& Rinit, ReachSpecification& options);

    /**
     * @brief computes the reachable continuous set for one time step of a nonlinear system by overapproximative linearization
     * @param R reachable set of the previous time step
     * @param options options for the computation of the reachable set
     * @return reachable set of the next time step
     */
    ReachableSet<Number> post(ReachableSet<Number>& R, ReachSpecification& options);

    /**
     * @brief computes the reachable set after linearization
     * @param options struct with algorithm settings
     * @param Rstart initial reachable set
     * @param Rti reachable set for time interval
     * @param Rtp reachable set for time point
     * @return dimForSplit - dimension that is split to reduce the lin. error
     */
    int linReach(ReachSpecification& options, ReachableSet<Number>& Rstart, ReachableSet<Number>& Rti, ReachableSet<Number>& Rtp);

};

/** @} */
}