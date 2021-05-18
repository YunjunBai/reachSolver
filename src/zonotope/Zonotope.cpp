/**
 * @file   Zonotope.cpp
 * @brief  Zonotope implementation
 * @author Yunjun
 * @date April 2021
 * @version 1.0
 * 
 * Reference:
 *   CORA ../contSet/zonotope/all.m 
 */

#include "Zonotope.h"

namespace reachSolver{

/*****************************************************************************
*                                                                           *
*                           Constructors and Destructors                    *
*                                                                           *
*****************************************************************************/
/**
 * @brief Constructor with no params
 */
template <typename Number, typename Converter, typename Setting>
Zonotope<Number,Converter,Setting>::Zonotope():
dimension_(0),center_(Vector_t<Number>::Matrix()),generators_( Matrix_t<Number>::Matrix()){ }

/**
 * @brief Constructor with dimension
 * @param dimension Dimensionality of Zonotope
 */
template <typename Number, typename Converter, typename Setting>
Zonotope<Number,Converter,Setting>::Zonotope(size_t dimension):
    dimension_(dimension),
    center_(Vector_t<Number>::Zero(dimension)), 
    generators_(Matrix_t<Number>::Zero(dimension, 1)) {
    assert(dimension!=0);
}

/**
 * @brief Constructor with center and generators.
 * @param center A  vector
 * @param generators A  matrix
 */
template <typename Number, typename Converter, typename Setting>
Zonotope<Number,Converter,Setting>::Zonotope(const Vector_t<Number>& center, const Matrix_t<Number>& generators):
    dimension_(center.rows()),
    center_(center),
    generators_(generators){
    assert(center.rows()==generators.rows());
}

template <typename Number, typename Converter, typename Setting>
Zonotope<Number,Converter,Setting>::~Zonotope(){}

/*****************************************************************************
*                                                                           *
*                            Public Functions on Set Properties             *
*                                                                           *
*****************************************************************************/

/**
 * @brief Get the dimension of Zonotope
 * @return the dimension
 */
template <typename Number, typename Converter, typename Setting>
size_t Zonotope<Number,Converter,Setting>::dimension() const{
    return dimension_;
}

/**
 * @brief Get the current center
 * @return center a nx1 matrix
 */
template <typename Number, typename Converter, typename Setting>
const Vector_t<Number>& Zonotope<Number,Converter,Setting>::center() const{
    return center_;
}

/**
 * @brief Replaces the current center with the parameter center
 * @param center a nx1 matrix
 */
template <typename Number, typename Converter, typename Setting>
void Zonotope<Number,Converter,Setting>::set_center(const Vector_t<Number>& center){
    if (dimension_ == 0){
        dimension_ = center.rows();
        generators_ = matrix_t<Number>::Zero(dimension_, 1);
    }
    assert((std::size_t)center.rows() == dimension_ && "Center has to have same dimensionality as zonotope." );
    center_ = center;
    // uniteEqualVectors();
    // removeEmptyGenerators();
}

/**
 * @brief Get the current generators
 * @return center a nxm matrix
 */
template <typename Number, typename Converter, typename Setting>
const Matrix_t<Number>& Zonotope<Number,Converter,Setting>::generators() const{
    return generators_;
}

/**
 * @brief Replaces the current matrix of generators with the parameter generators
 * @param generators a nxm matrix
 */
template <typename Number, typename Converter, typename Setting>
void Zonotope<Number,Converter,Setting>::set_generators(const Matrix_t<Number>& generators){
    if (dimension_ == 0) {
        dimension_ = generators.rows();
        center_ = vector_t<Number>::Zero(dimension_);
    }
    assert((std::size_t)generators.rows() == dimension_ && "Generators have to have same dimensionality as zonotope" );
    generators_ = generators;
    // uniteEqualVectors();
    // removeEmptyGenerators();
}

/**
 * @brief Add generators to Zonotope. Simply performs setGenerators if generators was previously not initialized.
 * @param generators a nxm matrix
 * @return true if able to add generators
 */
template <typename Number, typename Converter, typename Setting>
bool Zonotope<Number,Converter,Setting>::AddGenerators(const Matrix_t<Number>& generators){

}

/**
 * @brief Get the order
 * @return zonotope order
 */
template <typename Number, typename Converter, typename Setting>
Number Zonotope<Number,Converter,Setting>::order() const{

}

/**
 * @brief Number of generators
 * @return number of generators
 */
template <typename Number, typename Converter, typename Setting>
size_t Zonotope<Number,Converter,Setting>::numGenerators() const{

}

/**
 * @brief Removes zero generators in generator matrix
 */
template <typename Number, typename Converter, typename Setting>
void Zonotope<Number,Converter,Setting>::DeleteZeroGenerators(){

}

/**
 * @brief Changes the dimension of a Zonotope. if new_dim > old dim, new rows are initialized with null
 * @param new_dim The new dimension of the Zonotope
 * @return True, if change in dimension was successful
 */
template <typename Number, typename Converter, typename Setting>
bool Zonotope<Number,Converter,Setting>::ChangeDimension(size_t new_dim){

}

/**
 * @brief Reduces the order of a zonotope
 * @param limitOrder order of reduced zonotope
 */
template <typename Number, typename Converter, typename Setting>
void Zonotope<Number,Converter,Setting>::Reduce(unsigned limitOrder){

}

/**
 * @brief Clears the generators and center of the Zonotope and sets dimensionality to zero
 */
template <typename Number, typename Converter, typename Setting>
void Zonotope<Number,Converter,Setting>::Clear(){

}

/**
 * @brief display the zonotope
 */
template <typename Number, typename Converter, typename Setting>
void Zonotope<Number,Converter,Setting>::Display() const{
    
}

} // namespace reachSolver

