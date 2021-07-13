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

template <typename Number>
class ReachOptions{
  private:
    double t_;
    double tStart_;
    double tFinal_;
    Zonotope<Number> R0_;
    Zonotope<Number> U_;
    Vector_t<Number> uTrans_;
    Vector_t<Number> Rtrans_;

    double time_step_;
    int taylor_terms_;
    int zonotope_order_;
    int intermediate_order_;
    int error_order_;
    int tensor_order_;

    Number *factor_;

    int verbose_;
    std::string alg_;
    int max_error_;
    int originContained_;

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
    const double t() const;

    /**
     * @brief Replaces the current t with the parameter
     * @param t
     */
    void set_t(double t);

    /**
     * @brief Get the tStart
     * @return the tStart
     */
    const double tStart() const;

    /**
     * @brief Replaces the current tStart with the parameter
     * @param tStart
     */
    void set_tStart(double tStart);

    /**
     * @brief Get the tFinal
     * @return the tFinal
     */
    const double tFinal() const;

    /**
     * @brief Replaces the current tFinal with the parameter
     * @param tFinal
     */
    void set_tFinal(double tFinal);
    
    /**
     * @brief Get the R0
     * @return the R0
     */
    const Zonotope<Number> R0() const;

    /**
     * @brief Replaces the R0 with the parameter
     * @param R0
     */
    void set_R0(Zonotope<Number> R0);
    
    /**
     * @brief Get the U
     * @return the U
     */
    const Zonotope<Number> U() const;

    /**
     * @brief Replaces the U with the parameter 
     * @param U
     */
    void set_U(Zonotope<Number> U);
    
    /**
     * @brief Get the time_step
     * @return the time_step
     */
    const double time_step() const;

    /**
     * @brief Replaces the time_step with the parameter
     * @param time_step
     */
    void set_time_step(double time_step);
    
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
    const Number* factor() const;

    /**
   * @brief Replaces the factor with the parameter
   * @param factor
   */
    void set_factor(Number* factor);

    /**
     * @brief Get the verbose
     * @return the verbose
     */
    const int verbose() const;

    /**
     * @brief Replaces the verbose with the parameter
     * @param verbose
     */
    void set_verbose(int verbose);

    /**
     * @brief Get the alg
     * @return the alg
     */
    const std::string alg() const;

    /**
     * @brief Replaces the alg with the parameter
     * @param alg
     */
    void set_alg(std::string alg);

    /**
     * @brief Get the max_error
     * @return the max_error
     */
    const int max_error() const;

    /**
     * @brief Replaces the max_error with the parameter
     * @param max_error
     */
    void set_max_error(int max_error);
  
    /** 
     * @brief Get the uTrans
     * @return the uTrans
     */
    const Vector_t<Number> uTrans() const;

    /**
     * @brief Replaces the max_error with the parameter
     * @param max_error
     */
    void set_uTrans(Vector_t<Number> uTrans);

    /**
     * @brief Get the originContained
     * @return the originContained
     */
    const int originContained() const;

    /**
     * @brief Replaces the originContained with the parameter
     * @param originContained
     */
    void set_originContained(int originContained);

};
/** @} */
}