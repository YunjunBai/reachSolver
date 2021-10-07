/**
 * @file   TimeRelated.cpp
 * @brief  TimeRelated class
 * @author Yunjun Bai
 * @date April 2021
 * @version 1.0 
 * 
 */

#include "TimeRelated.h"

namespace reachSolver{

template class TimePoint<double>;
template class TimeInt<double>;

/*****************************************************************************
*                                                                           *
*                           TimePoint Class                                 *
*                                                                           *
*****************************************************************************/

template <typename Number>
TimePoint<Number>::TimePoint()
    :set_(std::vector<std::vector<ReachableSetElement<Number>>>()),
    time_(std::vector<double>()){

}

template <typename Number>
TimePoint<Number>::TimePoint(int size)
    :set_(std::vector<std::vector<ReachableSetElement<Number>>>(size)),
    time_(std::vector<double>(size)){
        
}

template <typename Number>
const std::vector<ReachableSetElement<Number>> TimePoint<Number>::rs(int i) const{
    return set_[i];
}

template <typename Number>
void TimePoint<Number>::set_rs(int i, std::vector<ReachableSetElement<Number>> set){
    set_[i] = set;
}

template <typename Number>
const double TimePoint<Number>::time(int i) const{
    return time_[i];
}

template <typename Number>
void TimePoint<Number>::set_time(int i, double time){
    time_[i] = time;
}


/*****************************************************************************
*                                                                           *
*                           TimeInt Class                                 *
*                                                                           *
*****************************************************************************/

template <typename Number>
TimeInt<Number>::TimeInt()
    :set_(std::vector<std::vector<ReachableSetElement<Number>>>()),
    time_(std::vector<Interval>()){

}

template <typename Number>
TimeInt<Number>::TimeInt(int size)
    :set_(std::vector<std::vector<ReachableSetElement<Number>>>(size)),
    time_(std::vector<Interval>(size)){
        
}

template <typename Number>
const std::vector<ReachableSetElement<Number>> TimeInt<Number>::rs(int i) const{
    return set_[i];
}

template <typename Number>
void TimeInt<Number>::set_rs(int i, std::vector<ReachableSetElement<Number>> set){
    set_[i] = set;
}

template <typename Number>
const Interval TimeInt<Number>::time(int i) const{
    return time_[i];
}

template <typename Number>
void TimeInt<Number>::set_time(int i, Interval time){
    time_[i] = time;
}


}