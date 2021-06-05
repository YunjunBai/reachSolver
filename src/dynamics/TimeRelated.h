/**
 * @file   TimeRelated.h
 * @brief  TimeRelated class representation
 * @author Yunjun Bai
 * @date April 2021
 * @version 1.0 
 * 
 */

#pragma once

#include "BasicObject.h"
#include "Zonotope.h"
namespace reachSolver{

typedef double NumberType;
/**
 * @class      TimePoint 
 * @brief      Class for TimePoint.
 * @ingroup    structure
 * @{
 */
class TimePoint{
  private:
    std::vector<Zonotope<NumberType>> set_;
    std::vector<double> time_;
  public:
    
    /*****************************************************************************
    *                                                                           *
    *                           Constructors and Destructors                    *
    *                                                                           *
    *****************************************************************************/

    /**
     * @brief Constructor with no params
     */
    TimePoint();

    /**
     * @brief Constructor with size
     */
    TimePoint(int size);

    /*****************************************************************************
    *                                                                           *
    *                       Public Functions on Properties                      *
    *                                                                           *
    *****************************************************************************/
    /**
     * @brief Get the set(i)
     * @return the set
     */
    const Zonotope<NumberType> set(int i) const;

    /**
     * @brief Replaces the current set with the parameter
     * @param set
     */
    void set_set_rs(int i, Zonotope<NumberType>& set);

    /**
     * @brief Get the time
     * @return the time
     */
    const double time(int i) const;

    /**
     * @brief Replaces the current time with the parameter
     * @param time
     */
    void set_time(int i, double time);
};
/** @} */

/**
 * @class      TimeInt 
 * @brief      Class for TimeInt.
 * @ingroup    structure
 * @{
 */
class TimeInt{
  private:
    struct Interval
    {
        double interval[2];
    };
    std::vector<Zonotope<NumberType>> set_;
    std::vector<Interval> time_;
  public:

    /*****************************************************************************
    *                                                                           *
    *                           Constructors and Destructors                    *
    *                                                                           *
    *****************************************************************************/

    /**
     * @brief Constructor with no params
     */
    TimeInt();

    /**
     * @brief Constructor with size
     */
    TimeInt(int size);

    /*****************************************************************************
    *                                                                           *
    *                       Public Functions on Properties                      *
    *                                                                           *
    *****************************************************************************/
    /**
     * @brief Get the set(i)
     * @return the set
     */
    const Zonotope<NumberType> set(int i) const;

    /**
     * @brief Replaces the current set with the parameter
     * @param set
     */
    void set_set_rs(int i, Zonotope<NumberType>& set);

    /**
     * @brief Get the time
     * @return the time
     */
    const Interval time(int i) const;

    /**
     * @brief Replaces the current time with the parameter
     * @param time
     */
    void set_time(int i, double time_start, double time_end);
};
/** @} */

}