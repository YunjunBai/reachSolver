/**
 * @file   ContDynamics.h
 * @brief  ContDynamics class representation
 * @author Yunjun Bai
 * @date April 2021
 * @version 1.0 
 * 
 * Reference:
 * CORA ../contDynamics/ContDynamics/all
 */

#pragma once

#include "zonotope/BasicObject.h"
#include <vector> 
#include <string>
#include <cmath>
#include <exception>
#include <float.h>
#include "ReachOptions.h"
#include "ReachableSet.h"
// #include "ReachSpecification.h"
#include "TimeRelated.h"

namespace reachSolver{

struct SetExplosionException : public std::exception
{
    const char * what () const throw ()
    {
        return "Set Explosion Exception!";
    }
};

/**
 * @class      ContDynamics 
 * @brief      Class for ContDynamics.
 * @ingroup    structure
 * @{
 */
template <typename Number>
class ContDynamics{
  private:
    //typedef Number (*function_type)(Number param1, Number param2);

    std::string name_;
    size_t dim_;
    size_t num_inputs_;
    size_t num_outputs_;

    // used in constructor
    //std::vector<size_t> numberOfInputs(function_type fun_handle, size_t inpArgs);

    // used in initReach
    //ReachableSet<Number> initReach_linRem(std::vector<ReachableSetElement<Number>>& Rinit, ReachOptions<Number>& options);

    //std::vector<Zonotope<Number>> split(Zonotope<Number> input, int number);

  public:
    /****************************************************************************
    *                                                                           *
    *                           Constructors and Destructors                    *
    *                                                                           *
    *****************************************************************************/

    /**
     * @brief Constructor with no params
     */
    // ContDynamics();

    /**
     * @brief Constructor with function handle to the dynamic equation
     * @param name name of dynamics
     */
    explicit ContDynamics(std::string name);

    /**
     * @brief Constructor with four params
     * @param name name of dynamics
     * @param dim number of states
     * @param num_inputs number of inputs
     * @param num_outputs number of outputs
     */
    ContDynamics(std::string name, size_t dim, size_t num_inputs, size_t num_outputs);

    virtual ~ContDynamics();

    /*****************************************************************************
    *                                                                           *
    *                       Public Functions on Properties                      *
    *                                                                           *
    *****************************************************************************/
    const std::string name() const;

    const size_t num_states() const;
    const size_t dim() const;

    const size_t num_inputs() const;

    /*****************************************************************************
    *                                                                           *
    *                       Compute Reachable Set                               *
    *                                                                           *
    *****************************************************************************/

};

/** @} */
}