/**
 * @file   ReachableSet.h
 * @brief  ReachableSet class representation
 * @author Yunjun Bai
 * @date April 2021
 * @version 1.0 
 * 
 */

#pragma once

#include "BasicObject.h"
#include "TimeRelated.h"
namespace reachSolver{

/**
 * @class      ReachableSet 
 * @brief      Class for ReachableSet.
 * @ingroup    structure
 * @{
 */
template <typename Number>
class ReachableSet{
  private:
    TimePoint time_point_;
    TimeInt time_interval_;
    int parent_rs_; //index of parent reachable set
    int loc_; //index of the location

  public:
    /****************************************************************************
    *                                                                           *
    *                           Constructors and Destructors                    *
    *                                                                           *
    *****************************************************************************/
    /**
     * @brief Constructor with TimePoint
     * @param time_point struct with fields .set and .time storing the time point reachable set
     */
    explicit ReachableSet(TimePoint time_point);

    /**
     * @brief Constructor with two params
     * @param time_point struct with fields .set and .time storing the time point reachable set
     * @param parent index of the parent reachable set
     */
    ReachableSet(TimePoint time_point, int parent);

    /**
     * @brief Constructor with two params
     * @param time_point struct with fields .set and .time storing the time point reachable set
     * @param time_interval struct with fields .set, .time, and .algebraic (nonlinDAsys) storing the time interval reachable set
     */
    ReachableSet(TimePoint time_point, TimeInt time_interval);

    /**
     * @brief Constructor with three params
     * @param time_point struct with fields .set and .time storing the time point reachable set
     * @param parent index of the parent reachable set
     * @param loc index of the location
     */
    ReachableSet(TimePoint time_point, int parent, int loc);

    /**
     * @brief Constructor with three params
     * @param time_point struct with fields .set and .time storing the time point reachable set
     * @param time_interval struct with fields .set, .time, and .algebraic (nonlinDAsys) storing the time interval reachable set
     * @param parent index of the parent reachable set
     */
    ReachableSet(TimePoint time_point, TimeInt time_interval, int parent);

    /**
     * @brief Constructor with four params
     * @param time_point struct with fields .set and .time storing the time point reachable set
     * @param time_interval struct with fields .set, .time, and .algebraic (nonlinDAsys) storing the time interval reachable set
     * @param parent index of the parent reachable set
     * @param loc index of the location
     */
    ReachableSet(TimePoint time_point, TimeInt time_interval, int parent, int loc);

    /**
     * @brief Copy Constructor - constructs a ReachableSet from an existing one.
     * @param other Another ReachableSet, from which a new ReachableSet is constructed
     */
    ReachableSet(const ReachableSet& other) = default;

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
    const TimePoint time_point() const;

    /**
     * @brief Get the time_interval
     * @return time_interval
     */
    const TimeInt time_interval() const;

    /**
     * @brief Get the parent_rs
     * @return parent_rs
     */
    const int parent_rs() const;

    /**
     * @brief Get the loc
     * @return loc
     */
    const int loc() const;
};
/** @} */
}