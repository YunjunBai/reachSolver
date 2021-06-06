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

template class Zonotope<double>;
/*****************************************************************************
*                                                                           *
*                           Constructors and Destructors                    *
*                                                                           *
*****************************************************************************/
/**
 * @brief Constructor with no params
 */
template <typename Number>
Zonotope<Number>::Zonotope():
dimension_(0),center_(0, 1),generators_(0, 0){ }

/**
 * @brief Constructor with dimension
 * @param dimension Dimensionality of Zonotope
 */
template <typename Number>
Zonotope<Number>::Zonotope(size_t dimension):
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
template <typename Number>
Zonotope<Number>::Zonotope(const Vector_t<Number>& center, const Matrix_t<Number>& generators):
    dimension_(center.rows()),
    center_(center),
    generators_(generators){
    assert(center.rows()==generators.rows());
}

template <typename Number>
Zonotope<Number>::~Zonotope(){}

/*****************************************************************************
*                                                                           *
*                            Public Functions on Set Properties             *
*                                                                           *
*****************************************************************************/

/**
 * @brief Get the dimension of Zonotope
 * @return the dimension
 */
template <typename Number>
size_t Zonotope<Number>::dimension() const{
    return dimension_;
}

/**
 * @brief Get the current center
 * @return center a nx1 matrix
 */
template <typename Number>
const Vector_t<Number>& Zonotope<Number>::center() const{
    return center_;
}

/**
 * @brief Replaces the current center with the parameter center
 * @param center a nx1 matrix
 */
template <typename Number>
void Zonotope<Number>::set_center(const Vector_t<Number>& center){
    if (dimension_ == 0){
        dimension_ = center.rows();
        generators_ = Matrix_t<Number>::Zero(dimension_, 1);
    }
    assert((std::size_t)center.rows() == dimension_ && "Center has to have same dimensionality as zonotope." );
    center_ = center;
    // uniteEqualVectors();
    DeleteZeroGenerators();
}

/**
 * @brief Get the current generators
 * @return center a nxm matrix
 */
template <typename Number>
const Matrix_t<Number>& Zonotope<Number>::generators() const{
    return generators_;
}

/**
 * @brief Replaces the current matrix of generators with the parameter generators
 * @param generators a nxm matrix
 */
template <typename Number>
void Zonotope<Number>::set_generators(const Matrix_t<Number>& generators){
    if (dimension_ == 0){
        dimension_ = generators.rows();
        center_ = Vector_t<Number>::Zero(dimension_);
    }
    assert((std::size_t)generators.rows() == dimension_ && "Generators have to have same dimensionality as zonotope" );
    generators_ = generators;
    // uniteEqualVectors();
    DeleteZeroGenerators();
}

/**
 * @brief Add generators to Zonotope. Simply performs setGenerators if generators was previously not initialized.
 * @param generators a nxm matrix
 * @return true if able to add generators
 */
template <typename Number>
bool Zonotope<Number>::AddGenerators(const Matrix_t<Number>& generators){
    if (dimension_ == 0){
        dimension_ = generators.rows();
        center_ = Vector_t<Number>::Zero(dimension_);
        generators_ = generators;
        return true;
    }
    assert((std::size_t)generators.rows() == dimension_ && "Added generators have to have same dimensionality as zonotope" );
    Matrix_t<Number> tmp = generators_;
    generators_.resize( tmp.rows(), generators.cols() + tmp.cols() );
    generators_ << tmp, generators;

    // uniteEqualVectors();
    DeleteZeroGenerators();
    return true;
}

/**
 * @brief Get the order
 * @return zonotope order
 */
template <typename Number>
Number Zonotope<Number>::order() const{
    // object empty.
    if (generators_.rows() == 0){
        return Number( 0 );
    }
    return Number(generators_.cols()) / Number(generators_.rows());
}

/**
 * @brief Number of generators
 * @return number of generators
 */
template <typename Number>
size_t Zonotope<Number>::numGenerators() const{
    return generators_.cols();
}

template <typename Number>
void Zonotope<Number>::RemoveGenerator(unsigned int colomn){
    Eigen::Index numRows = generators_.rows();
    Eigen::Index numCols = generators_.cols() - 1;

    if (colomn < numCols) {
        generators_.block(0, colomn, numRows, numCols - colomn) =
                generators_.block(0, colomn + 1, numRows, numCols - colomn);
    }

    generators_.conservativeResize(numRows, numCols);
}

/**
 * @brief Removes zero generators in generator matrix
 */
template <typename Number>
void Zonotope<Number>::DeleteZeroGenerators(){
    Vector_t<Number> zero_vector = Vector_t<Number>::Zero(dimension_);

    std::vector<unsigned> zeroIndex;
    for (unsigned i = 0; i < generators_.cols(); i++) {
        if (generators_.col(i) == zero_vector) {
            zeroIndex.push_back( i );
        }
    }

    for (std::vector<unsigned>::reverse_iterator r_it = zeroIndex.rbegin();r_it != zeroIndex.rend(); ++r_it) {
        RemoveGenerator(*r_it);
    }
}

/**
 * @brief Changes the dimension of a Zonotope. if new_dim > old dim, new rows are initialized with null
 * @param new_dim The new dimension of the Zonotope
 * @return True, if change in dimension was successful
 */
template <typename Number>
bool Zonotope<Number>::ChangeDimension(size_t new_dim){
    assert( new_dim != 0 && "Cannot change dimensionality of zonotope to zero" );
    if (new_dim == dimension_){
        return false;
    }else{
        center_.conservativeResize(new_dim, Eigen::NoChange);
        generators_.conservativeResize(new_dim, Eigen::NoChange);

        // If new dim > old dim, initialize all new rows to zero
        for (unsigned i = dimension_; i < new_dim; i++){
            center_.row(i).setZero();
            generators_.row(i).setZero();
        }

        dimension_ = new_dim;
        return true;
    }
}

/**
 * @brief Reduces the order of a zonotope
 * @param limitOrder order of reduced zonotope
 */
template <typename Number>
void Zonotope<Number>::Reduce(unsigned limitOrder){
    
}

/**
 * @brief Clears the generators and center of the Zonotope and sets dimensionality to zero
 */
template <typename Number>
void Zonotope<Number>::Clear(){
    generators_.resize(0, 0);
    center_.resize(0, 1);
    dimension_ = 0;
}

/**
 * @brief display the zonotope
 */
template <typename Number>
void Zonotope<Number>::Display() const{
    std::cout << "This Zonotope's dimension is " << dimension_ << std::endl;
    std::cout << "This Zonotope's center is \n" << center_ << std::endl;
    std::cout << "This Zonotope's generators are \n" << generators_ << std::endl;
}

/*****************************************************************************
*                                                                           *
*                          Basic Set Operations                               *
*                                                                           *
*****************************************************************************/

/**
 * @brief implement the linear maps of a set, i.e., "*" operator
 * @param matrix a matrix
 * @return a  zonotope = matrix * this zonotope
 */
template <typename Number>
Zonotope<Number> Zonotope<Number>::Times(const Matrix_t<Number>& matrix) const{
    assert(matrix.cols() == center_.rows());
    assert(matrix.cols() == center_.rows());
    Zonotope<Number> result;
    result.set_center(matrix * center_ );
    result.set_generators(matrix * generators_);
    return result;
}

/**
 * @brief Get the minkowski addition of two zonotope,i.e., "+" operator
 * @param another_zonotope 
 * @return a  zonotope = zonotope1 + zonotope2
 */
template <typename Number>
Zonotope<Number> Zonotope<Number>::Plus(const Zonotope& another_zonotope) const{
    assert(dimension_ == another_zonotope.dimension());
    Zonotope<Number> sum;
    sum.set_center(this->center_ + another_zonotope.center());
    Matrix_t<Number> sum_generators;
    sum_generators.resize(dimension_, generators_.cols() + another_zonotope.generators().cols());
    sum_generators << generators_, another_zonotope.generators();
    sum.set_generators(sum_generators);
    //sum.uniteEqualVectors();
    sum.DeleteZeroGenerators();
    //sum.reduce();
    return sum;
}

/**
 * @brief Get the enclosure for the convex hull of two zonotope
 * @param another_zonotope 
 * @return a  zonotope enclosing the convex hull
 */
template <typename Number>
Zonotope<Number> Zonotope<Number>::Convexhull(const Zonotope& another_zonotope) const{
    
}


} // namespace reachSolver

