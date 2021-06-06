/**
 * @file   LinearSys.h
 * @brief  LinearSys class representation
 * @author Yunjun Bai
 * @date April 2021
 * @version 1.0 
 * 
 * Reference:
 * CORA ../contDynamics/linearSys/all
 */

#pragma once

#include "ContDynamics.h"

namespace reachSolver{

/**
 * @class      LinearSys 
 * @brief      Class for LinearSys.
 * @ingroup    structure
 * @{
 */
template <typename Number>
class LinearSys : private ContDynamics{
  private:

    Matrix_t<Number> A; //system matrix
    Matrix_t<Number> B; //input matrix
    Matrix_t<Number> c; //constant input
    Matrix_t<Number> C; //output matrix
    Matrix_t<Number> D; //throughput matrix
    Matrix_t<Number> k; //output offset
    Matrix_t<Number> taylor; 
    Matrix_t<Number> krylov; 

  public:
    /****************************************************************************
    *                                                                           *
    *                           Constructors and Destructors                    *
    *                                                                           *
    *****************************************************************************/

    /**
     * @brief Constructor with no params
     */
    LinearSys();

    virtual ~LinearSys();

    /*****************************************************************************
    *                                                                           *
    *                       Public Functions on Properties                      *
    *                                                                           *
    *****************************************************************************/


    /*****************************************************************************
    *                                                                           *
    *                       Compute Reachable Set                               *
    *                                                                           *
    *****************************************************************************/

    /**
     * @brief computes the reachable continuous set for the first time step
     * @param Rinit initial reachable set
     * @param options struct containing the algorithm settings
     * @return first reachable set
     */
    LinearReachableSet<Number> initReach(Zonotope<Number>& Rinit, ReachOptions<Number>& options);
};

/** @} */
}