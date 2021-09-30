/**
 * @file   ReachOptions.h
 * @brief  ReachOptions class representation
 * @author Yunjun Bai
 * @date April 2021
 * @version 1.0 
 * 
 */

#pragma once

#include "zonotope/BasicObject.h"
#include "zonotope/Zonotope.h"
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
    Zonotope<Number> uTrans_lin_;
    Number uTrans_nonlin_;
    Zonotope<Number> Rtrans_;
    Zonotope<Number> Rhom_;
    Zonotope<Number> Rhom_tp_;
    Zonotope<Number> Raux_;
    Zonotope<Number> Rpar_;

    double time_step_;
    size_t taylor_terms_;
    size_t zonotope_order_;
    size_t intermediate_order_;
    size_t error_order_;
    size_t tensor_order_;

    Number *factor_;

    int verbose_;
    std::string alg_;
    int max_error_;
    int originContained_;

    int reductionInterval_;
    int usekrylovError_;

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
    const size_t taylor_terms() const;

    /**
     * @brief Replaces the taylor_terms with the parameter
     * @param taylor_terms
     */
    void set_taylor_terms(size_t taylor_terms);
    
    /**
     * @brief Get the zonotope_order
     * @return the zonotope_order
     */
    const size_t zonotope_order() const;

    /**
     * @brief Replaces the zonotope_order with the parameter
     * @param zonotope_order
     */
    void set_zonotope_order(size_t zonotope_order);
    
    /**
     * @brief Get the intermediate_order
     * @return the intermediate_order
     */
    const size_t intermediate_order() const;

    /**
     * @brief Replaces the intermediate_order with the parameter
     * @param intermediate_order
     */
    void set_intermediate_order(size_t intermediate_order);
    
    /**
     * @brief Get the error_order
     * @return the error_order
     */
    const size_t error_order() const;

    /**
     * @brief Replaces the error_order with the parameter
     * @param error_order
     */
    void set_error_order(size_t error_order);
    
    /**
     * @brief Get the tensor_order
     * @return the tensor_order
     */
    const size_t tensor_order() const;

    /**
     * @brief Replaces the tensor_order with the parameter
     * @param tensor_order
     */
    void set_tensor_order(size_t tensor_order);

    /**
     * @brief Get the factor
     * @return the factor
     */
    Number* factor() const;

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
    const Zonotope<Number> uTrans_lin() const;

    /**
     * @brief Replaces the uTrans with the parameter
     * @param uTrans
     */
    void set_uTrans_lin(Zonotope<Number> uTrans_lin);
  
    /** 
     * @brief Get the uTrans
     * @return the uTrans
     */
    const Number uTrans_nonlin() const;

    /**
     * @brief Replaces the uTrans with the parameter
     * @param uTrans
     */
    void set_uTrans_nonlin(Number uTrans_nonlin);

    /** 
     * @brief Get the Rhom
     * @return the Rhom
     */
    const Zonotope<Number> Rhom() const;

    /**
     * @brief Replaces the Rhom with the parameter
     * @param Rhom
     */
    void set_Rhom(Zonotope<Number> Rhom);
    
    /** 
     * @brief Get the Rhom_tp
     * @return the Rhom_tp
     */
    const Zonotope<Number> Rhom_tp() const;

    /**
     * @brief Replaces the Rhom_tp with the parameter
     * @param Rhom_tp
     */
    void set_Rhom_tp(Zonotope<Number> Rhom_tp);
    
    /** 
     * @brief Get the Raux
     * @return the Raux
     */
    const Zonotope<Number> Raux() const;

    /**
     * @brief Replaces the Raux with the parameter
     * @param Raux
     */
    void set_Raux(Zonotope<Number> Raux);

    /** 
     * @brief Get the Rpar
     * @return the Rpar
     */
    const Zonotope<Number> Rpar() const;

    /**
     * @brief Replaces the Rpar with the parameter
     * @param Rpar
     */
    void set_Rpar(Zonotope<Number> Rpar);

    /** 
     * @brief Get the Rtrans
     * @return the Rtrans
     */
    const Zonotope<Number> Rtrans() const;

    /**
     * @brief Replaces the Rtrans with the parameter
     * @param Rtrans
     */
    void set_Rtrans(Zonotope<Number> Rtrans);
    
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

    /**
     * @brief Get the reductionInterval
     * @return the reductionInterval
     */
    const int reductionInterval() const;

    /**
     * @brief Replaces the reductionInterval with the parameter
     * @param reductionInterval
     */
    void set_reductionInterval(int reductionInterval);

    /**
     * @brief Get the usekrylovError
     * @return the usekrylovError
     */
    const int usekrylovError() const;

    /**
     * @brief Replaces the usekrylovError with the parameter
     * @param usekrylovError
     */
    void set_usekrylovError(int usekrylovError);
    
};
/** @} */
}