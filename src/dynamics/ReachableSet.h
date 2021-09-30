/**
 * @file   ReachableSet.h
 * @brief  ReachableSet class representation
 * @author Yunjun Bai
 * @date April 2021
 * @version 1.0 
 * 
 */

#pragma once

#include "zonotope/BasicObject.h"
//#include "TimeRelated.h"
#include "zonotope/Zonotope.h"
namespace reachSolver{

/**
 * @class      ReachableSetElement 
 * @brief      Class for ReachableSetElement.
 * @ingroup    structure
 * @{
 */
template <typename Number>
class ReachableSetElement{
  private:
    Zonotope<Number> set_;
    Vector_t<Number> error_; 
    int prev_;
    int parent_;

  public:
    /****************************************************************************
    *                                                                           *
    *                           Constructors and Destructors                    *
    *                                                                           *
    *****************************************************************************/
    /**
     * @brief Constructor with no params
     */
    ReachableSetElement();

    /**
     * @brief Constructor with two params
     * @param set Zonotope<Number>
     * @param error error vector
     */
    ReachableSetElement(Zonotope<Number> set, Vector_t<Number> error);

    /*****************************************************************************
    *                                                                           *
    *                       Public Functions on Properties                      *
    *                                                                           *
    *****************************************************************************/
    /**
     * @brief Get the set
     * @return the set
     */
    Zonotope<Number> rs() const;

    /**
     * @brief Replaces the current set with the parameter
     * @param set
     */
    void set_rs(Zonotope<Number>& set);

    /**
     * @brief Get the error
     * @return the error
     */
    Vector_t<Number> error() const;

    /**
     * @brief Replaces the current error with the parameter
     * @param time
     */
    void set_error(Vector_t<Number>& error);

    /**
     * @brief Get the prev
     * @return the prev
     */
    const int prev() const;

    /**
     * @brief Replaces the current prev with the parameter
     * @param prev
     */
    void set_prev(int prev);

    /**
     * @brief Get the parent
     * @return the parent
     */
    const int parent() const;

    /**
     * @brief Replaces the current parent with the parameter
     * @param parent
     */
    void set_parent(int parent);
};
/** @} */

/**
 * @class      ReachableSet 
 * @brief      Class for ReachableSet.
 * @ingroup    structure
 * @{
 */
template <typename Number>
class ReachableSet{
  private:
    std::vector<ReachableSetElement<Number>> time_point_;
    std::vector<ReachableSetElement<Number>> time_interval_;
    std::vector<ReachableSetElement<Number>> R0_;
    int parent_rs_; //index of parent reachable set
    int loc_; //index of the location

    int internalCount_;

  public:
    /****************************************************************************
    *                                                                           *
    *                           Constructors and Destructors                    *
    *                                                                           *
    *****************************************************************************/
    /**
     * @brief Constructor with no params
     */
    ReachableSet();
    
    /**
     * @brief Constructor with TimePoint
     * @param time_point struct with fields .set and .time storing the time point reachable set
     */
    explicit ReachableSet(std::vector<ReachableSetElement<Number>>& time_point);

    /**
     * @brief Constructor with two params
     * @param time_point struct with fields .set and .time storing the time point reachable set
     * @param parent index of the parent reachable set
     */
    ReachableSet(std::vector<ReachableSetElement<Number>>& time_point, int parent);

    /**
     * @brief Constructor with two params
     * @param time_point struct with fields .set and .time storing the time point reachable set
     * @param time_interval struct with fields .set, .time, and .algebraic (nonlinDAsys) storing the time interval reachable set
     */
    ReachableSet(std::vector<ReachableSetElement<Number>>& time_point, std::vector<ReachableSetElement<Number>>& time_interval);

    /**
     * @brief Constructor with three params
     * @param time_point struct with fields .set and .time storing the time point reachable set
     * @param parent index of the parent reachable set
     * @param loc index of the location
     */
    ReachableSet(std::vector<ReachableSetElement<Number>>& time_point, int parent, int loc);

    /**
     * @brief Constructor with three params
     * @param time_point struct with fields .set and .time storing the time point reachable set
     * @param time_interval struct with fields .set, .time, and .algebraic (nonlinDAsys) storing the time interval reachable set
     * @param parent index of the parent reachable set
     */
    ReachableSet(std::vector<ReachableSetElement<Number>>& time_point, std::vector<ReachableSetElement<Number>>& time_interval, int parent);

    /**
     * @brief Constructor with four params
     * @param time_point struct with fields .set and .time storing the time point reachable set
     * @param time_interval struct with fields .set, .time, and .algebraic (nonlinDAsys) storing the time interval reachable set
     * @param parent index of the parent reachable set
     * @param loc index of the location
     */
    ReachableSet(std::vector<ReachableSetElement<Number>>& time_point, std::vector<ReachableSetElement<Number>>& time_interval, int parent, int loc);

    /**
     * @brief Copy Constructor - constructs a ReachableSet from an existing one.
     * @param other Another ReachableSet, from which a new ReachableSet is constructed
     */
    ReachableSet(const ReachableSet<Number>& other) = default;

    virtual ~ReachableSet();

    /*****************************************************************************
    *                                                                           *
    *                       Public Functions on Properties                      *
    *                                                                           *
    *****************************************************************************/

    /**
     * @brief Get the time_point
     * @return time_point
     */
    const std::vector<ReachableSetElement<Number>> time_point() const;


    /**
     * @brief Replaces the current time_point with the parameter
     * @param time_point
     */
    void set_time_point(std::vector<ReachableSetElement<Number>>& time_point);

    /**
     * @brief Get the time_interval
     * @return time_interval
     */
    const std::vector<ReachableSetElement<Number>> time_interval() const;

    /**
     * @brief Replaces the current time_interval with the parameter
     * @param time_interval
     */
    void set_time_interval(std::vector<ReachableSetElement<Number>>& time_interval);

    /**
     * @brief Get the time_interval
     * @return time_interval
     */
    const std::vector<ReachableSetElement<Number>> R0() const;

    /**
     * @brief Replaces the current R0 with the parameter
     * @param R0
     */
    void set_R0(std::vector<ReachableSetElement<Number>>& R0);

    /**
     * @brief Get the parent_rs
     * @return parent_rs
     */
    const int parent_rs() const;

    /**
     * @brief Replaces the current parent_rs with the parameter
     * @param parent_rs
     */
    void set_parent_rs(int parent_rs);

    /**
     * @brief Get the loc
     * @return loc
     */
    const int loc() const;

    /**
     * @brief Replaces the current loc with the parameter
     * @param loc
     */
    void set_loc(int loc);

    /**
     * @brief Get the internalCount
     * @return internalCount
     */
    const int internalCount() const;

    /**
     * @brief Replaces the current internalCount with the parameter
     * @param internalCount
     */
    void set_internalCount(int internalCount);

    /*****************************************************************************
    *                                                                           *
    *                       Other Functions                                     *
    *                                                                           *
    *****************************************************************************/
    
    /**
     * @brief delete reachable sets that are already covered by other sets
     * @param Rold reachable sets of previous time steps
     * @param options options for the computation of the reachable set
     * @return R - reachable sets
     */
    // ReachableSet<Number> deleteRedundantSets(ReachableSet<Number> Rold, ReachSpecification& options);
    ReachableSet<Number> deleteRedundantSets(ReachableSet<Number> Rold, ReachOptions<Number>& options);
};
/** @} */

/**
 * @class      ReachableSet 
 * @brief      Class for ReachableSet.
 * @ingroup    structure
 * @{
 */
template <typename Number>
class LinearReachableSet{
  private:
    Zonotope<Number> time_point_;
    Zonotope<Number> time_interval_;

  public:
    /****************************************************************************
    *                                                                           *
    *                           Constructors and Destructors                    *
    *                                                                           *
    *****************************************************************************/
    /**
     * @brief Constructor with no params
     */
    LinearReachableSet();

    /**
     * @brief Constructor with two params
     * @param time_point struct with fields .set and .time storing the time point reachable set
     * @param time_interval struct with fields .set, .time, and .algebraic (nonlinDAsys) storing the time interval reachable set
     */
    LinearReachableSet(Zonotope<Number> & time_point, Zonotope<Number> & time_interval);

    /**
     * @brief Copy Constructor - constructs a ReachableSet from an existing one.
     * @param other Another ReachableSet, from which a new ReachableSet is constructed
     */
    LinearReachableSet(const LinearReachableSet<Number>& other) = default;

    virtual ~LinearReachableSet();

    /*****************************************************************************
    *                                                                           *
    *                       Public Functions on Properties                      *
    *                                                                           *
    *****************************************************************************/

    /**
     * @brief Get the time_point
     * @return time_point
     */
    const Zonotope<Number>  time_point() const;


    /**
     * @brief Replaces the current time_point with the parameter
     * @param time_point
     */
    void set_time_point(Zonotope<Number> & time_point);

    /**
     * @brief Get the time_interval
     * @return time_interval
     */
    const Zonotope<Number> time_interval() const;

    /**
     * @brief Replaces the current time_interval with the parameter
     * @param time_interval
     */
    void set_time_interval(Zonotope<Number> & time_interval);


    /*****************************************************************************
    *                                                                           *
    *                       Other Functions                                     *
    *                                                                           *
    *****************************************************************************/

};
/** @} */
}
