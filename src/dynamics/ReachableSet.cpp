/**
 * @file   ReachableSet.cpp
 * @brief  ReachableSet class
 * @author Yunjun Bai
 * @date April 2021
 * @version 1.0 
 * 
 */

#include "ReachableSet.h"

namespace reachSolver{

template class ReachableSetElement<double>;
template class ReachableSet<double>;
template class LinearReachableSet<double>;

/*************************Class  ReachableSetElement************************/
template <typename Number>
ReachableSetElement<Number>::ReachableSetElement()
:set_(Zonotope<Number>()),error_(Vector_t<Number>()),prev_(0),parent_(0){
}

template <typename Number>
ReachableSetElement<Number>::ReachableSetElement(Zonotope<Number> set, Vector_t<Number> error)
:set_(set),error_(error),prev_(0),parent_(0){
}

template <typename Number>
Zonotope<Number> ReachableSetElement<Number>::rs() const{
    return set_;
}

template <typename Number>
void ReachableSetElement<Number>::set_rs(Zonotope<Number>& set){
    set_ = set;
}

template <typename Number>
Vector_t<Number> ReachableSetElement<Number>::error() const{
    return error_;
}

template <typename Number>
void ReachableSetElement<Number>::set_error(Vector_t<Number>& error){
    error_ = error;
}

template <typename Number>
const int ReachableSetElement<Number>::prev() const{
    return prev_;
}

template <typename Number>
void ReachableSetElement<Number>::set_prev(int prev){
    prev_ = prev;
}

template <typename Number>
const int ReachableSetElement<Number>::parent() const{
    return parent_;
}

template <typename Number>
void ReachableSetElement<Number>::set_parent(int parent){
    parent_ = parent;
}

/*************************Class  ReachableSet************************/
template <typename Number>
ReachableSet<Number>::ReachableSet()
    :time_point_(std::vector<ReachableSetElement<Number>>()),
    time_interval_(std::vector<ReachableSetElement<Number>>()),
    R0_(std::vector<ReachableSetElement<Number>>()),
    parent_rs_(0), loc_(0), internalCount_(0){
}

template <typename Number>
ReachableSet<Number>::ReachableSet(std::vector<ReachableSetElement<Number>>& time_point)
    :time_point_(time_point),
    time_interval_(std::vector<ReachableSetElement<Number>>()),
    R0_(std::vector<ReachableSetElement<Number>>()),
    parent_rs_(0),loc_(0), internalCount_(0){
}

template <typename Number>
ReachableSet<Number>::ReachableSet(std::vector<ReachableSetElement<Number>>& time_point, int parent)
    :time_point_(time_point),
    time_interval_(std::vector<ReachableSetElement<Number>>()),
    R0_(std::vector<ReachableSetElement<Number>>()),
    parent_rs_(parent), loc_(0), internalCount_(0){
}

template <typename Number>
ReachableSet<Number>::ReachableSet(std::vector<ReachableSetElement<Number>>& time_point, std::vector<ReachableSetElement<Number>>& time_interval)
    :time_point_(time_point),
    time_interval_(time_interval),
    R0_(std::vector<ReachableSetElement<Number>>()),
    parent_rs_(0),loc_(0), internalCount_(0){
}

template <typename Number>
ReachableSet<Number>::ReachableSet(std::vector<ReachableSetElement<Number>>& time_point, int parent, int loc)
    :time_point_(time_point),
    time_interval_(std::vector<ReachableSetElement<Number>>()),
    R0_(std::vector<ReachableSetElement<Number>>()),
    parent_rs_(parent),loc_(loc), internalCount_(0){
}

template <typename Number>
ReachableSet<Number>::ReachableSet(std::vector<ReachableSetElement<Number>>& time_point, std::vector<ReachableSetElement<Number>>& time_interval, int parent)
    :time_point_(time_point),
    time_interval_(time_interval),
    R0_(std::vector<ReachableSetElement<Number>>()),
    parent_rs_(parent),loc_(0), internalCount_(0){
}

template <typename Number>
ReachableSet<Number>::ReachableSet(std::vector<ReachableSetElement<Number>>& time_point, std::vector<ReachableSetElement<Number>>& time_interval, int parent, int loc)
    :time_point_(time_point),
    time_interval_(time_interval),
    R0_(std::vector<ReachableSetElement<Number>>()),
    parent_rs_(parent),loc_(loc), internalCount_(0){
}

template <typename Number>
ReachableSet<Number>::~ReachableSet(){}

template <typename Number>
const std::vector<ReachableSetElement<Number>> ReachableSet<Number>::time_point() const{
    return time_point_;
}

template <typename Number>
void ReachableSet<Number>::set_time_point(std::vector<ReachableSetElement<Number>>& time_point){
    time_point_ = time_point;
}

template <typename Number>
const std::vector<ReachableSetElement<Number>> ReachableSet<Number>::time_interval() const{
    return time_interval_;
}

template <typename Number>
void ReachableSet<Number>::set_time_interval(std::vector<ReachableSetElement<Number>>& time_interval){
    time_interval_ = time_interval;
}

template <typename Number>
const std::vector<ReachableSetElement<Number>> ReachableSet<Number>::R0() const{
    return R0_;
}

template <typename Number>
void ReachableSet<Number>::set_R0(std::vector<ReachableSetElement<Number>>& R0){
    R0_ = R0;
}

template <typename Number>
const int ReachableSet<Number>::parent_rs() const{
    return parent_rs_;
}

template <typename Number>
void ReachableSet<Number>::set_parent_rs(int parent_rs){
    parent_rs_ = parent_rs;
}

template <typename Number>
const int ReachableSet<Number>::loc() const{
    return loc_;
}

template <typename Number>
void ReachableSet<Number>::set_loc(int loc){
    loc_ = loc;
}

template <typename Number>
const int ReachableSet<Number>::internalCount() const{
    return internalCount_;
}

template <typename Number>
void ReachableSet<Number>::set_internalCount(int internalCount){
    internalCount_ = internalCount;
}

template <typename Number>
void ReachableSet<Number>::deleteRedundantSets(ReachableSet<Number> Rold, ReachOptions<Number>& options){
    // set reduction method
    // redMethod='pca';
    if(Rold.internalCount() == 0){
        this->internalCount_ = 3;
    }else{
        this->internalCount_ = Rold.internalCount()+1;
    }

    if(this->internalCount_ == options.reductionInterval()){
        // reset internal count
        this->internalCount_ = 1;
        // overapproximate reachable set of time point(!) by parallelpipeds and
        // save them as polytopes
        // R.P=[];
        // for i=1:length(R.tp)
        //     R.tp{i}.set=reduce(R.tp{i}.set,redMethod,1);
        //     %generate polytope
        //     R.P{i}=polytope(R.tp{i}.set,options);
    }else if(this->internalCount_ == 2){
        // just ignore for using polytope
    }

}

/*************************Class  LinearReachableSet************************/

template <typename Number>
LinearReachableSet<Number>::LinearReachableSet()
:time_point_(Zonotope<Number>()),time_interval_(Zonotope<Number>()){
}

template <typename Number>
LinearReachableSet<Number>::LinearReachableSet(Zonotope<Number> & time_point, Zonotope<Number> & time_interval)
:time_point_(time_point),time_interval_(time_interval){
}

template <typename Number>
LinearReachableSet<Number>::~LinearReachableSet(){}

template <typename Number>
const Zonotope<Number> LinearReachableSet<Number>::time_point() const{
    return time_point_;
}

template <typename Number>
void LinearReachableSet<Number>::set_time_point(Zonotope<Number>& time_point){
    time_point_ = time_point;
}

template <typename Number>
const Zonotope<Number> LinearReachableSet<Number>::time_interval() const{
    return time_interval_;
}

template <typename Number>
void LinearReachableSet<Number>::set_time_interval(Zonotope<Number>& time_interval){
    time_interval_ = time_interval;
}

} // namespace