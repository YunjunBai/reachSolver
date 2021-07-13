/**
 * @file   Zonotope.h
 * @brief  Zonotope class representation
 * @author Yunjun Bai
 * @date April 2021
 * @version 1.0 
 * 
 * Reference:
 *   CORA ../contSet/zonotope/zonotope.m 
 */
#pragma once

#include "BasicObject.h"

namespace reachSolver {

/**
 * @class      Zonotope 
 * @brief      Class for Zonotopes.
 * @tparam     Number     The used number type.
 * @ingroup    structure
 * @{
 */
template <typename Number>
class Zonotope : private BasicObject {
  private:
    size_t dimension_;
    Vector_t<Number> center_;
    Matrix_t<Number> generators_;
    
    /**
     * @brief Remove a generator
     * @param colomn column of generator
     */
    void RemoveGenerator(unsigned int colomn);

  public:
    /*****************************************************************************
    *                                                                           *
    *                           Constructors and Destructors                    *
    *                                                                           *
    *****************************************************************************/

    /**
     * @brief Constructor with no params
     */
    Zonotope();

    /**
     * @brief Constructor with dimension
     * @param dimension Dimensionality of Zonotope
     */
    explicit Zonotope(size_t dimension);

    /**
     * @brief Constructor with center and generators.
     * @param center A vector 
     * @param generators A  matrix
     */
    Zonotope(const Vector_t<Number>& center, const Matrix_t<Number>& generators);

    /**
     * @brief Copy Constructor - constructs a zonotope from an existing one.
     * @param other Another Zonotope, from which a new zonotope is constructed
     */
    Zonotope(const Zonotope& other) = default;

    virtual ~Zonotope();

    /*****************************************************************************
    *                                                                           *
    *                       Public Functions on Properties                      *
    *                                                                           *
    *****************************************************************************/

    /**
     * @brief Get the dimension of Zonotope
     * @return the dimension
     */
    size_t dimension() const;

    /**
     * @brief Get the current center
     * @return center a nx1 matrix
     */
    const Vector_t<Number>& center() const;

    /**
     * @brief Replaces the current center with the parameter center
     * @param center a nx1 matrix
     */
    void set_center(const Vector_t<Number>& center);

    /**
     * @brief Get the current generators
     * @return center a nxm matrix
     */
    const Matrix_t<Number>& generators() const;

    /**
     * @brief Replaces the current matrix of generators with the parameter generators
     * @param generators a nxm matrix
     */
    void set_generators(const Matrix_t<Number>& generators);

    /**
     * @brief Add generators to Zonotope. Simply performs setGenerators if generators was previously not initialized.
     * @param generators a nxm matrix
     * @return true if able to add generators
     */
    bool AddGenerators(const Matrix_t<Number>& generators);

    /**
     * @brief Get the order
     * @return zonotope order
     */
    Number order() const;

    /**
     * @brief Number of generators
     * @return number of generators
     */
    size_t numGenerators() const;

    /**
     * @brief Removes zero generators in generator matrix
     */
    void DeleteZeroGenerators();

    /**
     * @brief Changes the dimension of a Zonotope. if new_dim > old dim, new rows are initialized with null
     * @param new_dim The new dimension of the Zonotope
     * @return True, if change in dimension was successful
     */
    bool ChangeDimension(size_t new_dim);

    /**
     * @brief Reduces the order of a zonotope
     * @param limitOrder order of reduced zonotope
     */
    void Reduce(unsigned limitOrder);
    
    /**
     * @brief Whether the Zonotope is empty
     * @return 0 if empty, 1 if not empty
     */
    int Empty();

    /**
     * @brief Clears the generators and center of the Zonotope and sets dimensionality to zero
     */
    void Clear();

    /**
     * @brief display the zonotope
     */
    void Display() const;

    /**
     * @brief Judge whether two zonotope is equal
     * @param another_zonotope 
     * @return True, if equal
     */
    bool operator==(const Zonotope<Number>& another_zonotope) const {
        if ( this->dimension_ != another_zonotope.dimension() ) {
            return false;
        }
        if ( this->center_ != another_zonotope.center() ) {
            return false;
        }
        if ( this->generators_ != another_zonotope.generators() ) {
            return false;
        }
        return true;
    }

    /*****************************************************************************
    *                                                                           *
    *                          Basic Set Operations                             *
    *                                                                           *
    *****************************************************************************/

    /**
     * @brief implement the linear maps of a set, i.e., "*" operator
     * @param matrix a matrix
     * @return a  zonotope = matrix * this zonotope
     */
    Zonotope Times(const Matrix_t<Number>& matrix) const;

    /**
     * @brief implement the linear maps of a set, i.e., "*" operator
     * @param matrix a matrix
     * @return a  zonotope = matrix * this zonotope
     */
    Zonotope operator*(const Matrix_t<Number>& matrix) const;


    /**
     * @brief implement the times of zonotope and intervalmatrix
     * @param int_matrix an IntervalMatrix
     * @return a  zonotope = int_matrix * this zonotope
     */
    Zonotope Times(const IntervalMatrix int_matrix) const;

    /**
     * @brief implement the times of zonotope and intervalmatrix
     * @param int_matrix an IntervalMatrix
     * @return a  zonotope = int_matrix * this zonotope
     */
    Zonotope operator*(const IntervalMatrix int_matrix) const;

    /**
     * @brief Get the minkowski addition of two zonotope,i.e., "+" operator
     * @param another_zonotope 
     * @return a  zonotope = zonotope1 + zonotope2
     */
    Zonotope Plus(const Zonotope<Number>& another_zonotope) const;

    /**
     * @brief Get the addition of a zonotope and a vector
     * @param vector 
     * @return the sum
     */
    Zonotope Plus(const Vector_t<Number>& vector) const;

    /**
     * @brief Get the minkowski addition of two zonotope,i.e., "+" operator
     * @param another_zonotope 
     * @return the sum
     */
    Zonotope operator+(const Zonotope<Number>& another_zonotope) const;

    /**
     * @brief Get the addition of a zonotope and a vector
     * @param another_zonotope 
     * @return the sum
     */
    Zonotope operator+(const Vector_t<Number>& vector) const;

    /**
     * @brief Get the enclosure for the convex hull of two zonotope
     * @param another_zonotope 
     * @return a zonotope enclosing the convex hull
     */
    Zonotope Convexhull(const Zonotope& another_zonotope) const;

    /**
     * @brief Generates a zonotope that encloses a zonotopes and its linear transformation
     * @param another_zonotope second zonotope object, satisfying Z2 = (M * Z1) + Zplus
     * @return zonotope, that encloses Z1 and Z2
     */
    Zonotope enclose(const Zonotope& another_zonotope) const;

    
};

/** @} */

}  // namespace reachSolver
