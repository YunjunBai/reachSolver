/**
 * @file   ReachableOptions.h
 * @brief  ReachableOptions class representation
 * @author Yunjun Bai
 * @date April 2021
 * @version 1.0 
 * 
 */

#pragma once

#include "BasicObject.h"
#include "Zonotope.h"
namespace reachSolver{

/**
 * @class      ReachableOptions 
 * @brief      Class for ReachableOptions
 * @ingroup    structure
 * @{
 */

typedef double NumberType;
class ReachOptions{
  private:
    int t_;
    int tStart_;
    int tFinal_;
    Zonotope<NumberType> R0_;
    Zonotope<NumberType> U_;

    int time_step_;
    int taylor_terms_;
    int zonotope_order_;
    int intermediate_order_;
    int error_order_;
    int alg_;
    int tensor_order_;

    NumberType *factor_;

  public:
    /*****************************************************************************
    *                                                                           *
    *                           Constructors and Destructors                    *
    *                                                                           *
    *****************************************************************************/

    /**
     * @brief Constructor with no params
     */
    ReachOptions();

    /*****************************************************************************
    *                                                                           *
    *                       Public Functions on Properties                      *
    *                                                                           *
    *****************************************************************************/


    /**
     * @brief Get the t
     * @return the t
     */
    const int t() const;

    /**
     * @brief Replaces the current t with the parameter
     * @param t
     */
    void set_t(int t);

    /**
     * @brief Get the tStart
     * @return the tStart
     */
    const int tStart() const;

    /**
     * @brief Replaces the current tStart with the parameter
     * @param tStart
     */
    void set_tStart(int tStart);

    /**
     * @brief Get the tFinal
     * @return the tFinal
     */
    const int tFinal() const;

    /**
     * @brief Replaces the current tFinal with the parameter
     * @param tFinal
     */
    void set_tFinal(int tFinal);
    
    /**
     * @brief Get the R0
     * @return the R0
     */
    const Zonotope<NumberType> R0() const;

    /**
     * @brief Replaces the R0 with the parameter
     * @param R0
     */
    void set_R0(Zonotope<NumberType> R0);
    
    /**
     * @brief Get the U
     * @return the U
     */
    const Zonotope<NumberType> U() const;

    /**
     * @brief Replaces the U with the parameter 
     * @param U
     */
    void set_U(Zonotope<NumberType> U);
    
    /**
     * @brief Get the time_step
     * @return the time_step
     */
    const int time_step() const;

    /**
     * @brief Replaces the time_step with the parameter
     * @param time_step
     */
    void set_time_step(int time_step);
    
    /**
     * @brief Get the taylor_terms
     * @return the taylor_terms
     */
    const int taylor_terms() const;

    /**
     * @brief Replaces the taylor_terms with the parameter
     * @param taylor_terms
     */
    void set_taylor_terms(int taylor_terms);
    
    /**
     * @brief Get the zonotope_order
     * @return the zonotope_order
     */
    const int zonotope_order() const;

    /**
     * @brief Replaces the zonotope_order with the parameter
     * @param zonotope_order
     */
    void set_zonotope_order(int zonotope_order);
    
    /**
     * @brief Get the intermediate_order
     * @return the intermediate_order
     */
    const int intermediate_order() const;

    /**
     * @brief Replaces the intermediate_order with the parameter
     * @param intermediate_order
     */
    void set_intermediate_order(int intermediate_order);
    
    /**
     * @brief Get the error_order
     * @return the error_order
     */
    const int error_order() const;

    /**
     * @brief Replaces the error_order with the parameter
     * @param error_order
     */
    void set_error_order(int error_order);
    
    /**
     * @brief Get the alg
     * @return the alg
     */
    const int alg() const;

    /**
     * @brief Replaces the alg with the parameter
     * @param alg
     */
    void set_alg(int alg);

    /**
     * @brief Get the tensor_order
     * @return the tensor_order
     */
    const int tensor_order() const;

    /**
     * @brief Replaces the tensor_order with the parameter
     * @param tensor_order
     */
    void set_tensor_order(int tensor_order);

    /**
   * @brief Get the factor
   * @return the factor
   */
    const NumberType* factor() const;

    /**
   * @brief Replaces the factor with the parameter
   * @param factor
   */
    void set_factor(NumberType* factor);
    
};
/** @} */
}