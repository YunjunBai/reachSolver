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
 * @tparam     Converter  The converter.
 * @ingroup    structure
 * @{
 */
template <typename Number, typename Converter, typename Setting>
class Zonotope : private BasicObject {
  private:
    size_t dimension_;
    vector_t<Number> center_;
    matrix_t<Number> generators_;
    
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
	 * @param center A (nx1) vector
	 * @param generators A (nxm) vector
	 */
	Zonotope(const vector_t<Number>& center, const matrix_t<Number>& generators);

	/**
	 * @brief Copy Constructor - constructs a zonotope from an existing one.
	 * @param other Another Zonotope, from which a new zonotope is constructed
	 */
	Zonotope(const Zonotope& other) = default;

	/**
	 * @brief Copy Constructor - constructs a 2D-zonotope of from an existing ND one.
	 * @param other : Another Zonotope, from which a new zonotope is constructed
	 * @param d1 : 1st dimension (0 <= d1 < other.dimension)
	 * @param d2 : 2nd dimension (0 <= d2 < other.dimension) d1!=d2
	 */
	Zonotope(const Zonotope& other, unsigned d1, unsigned d2);

	virtual ~Zonotope();

    /*****************************************************************************
	*                                                                           *
	*                           Public Functions                                *
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
	const vector_t<Number>& center() const;

    /**
	 * @brief Replaces the current center with the parameter center
	 * @param center a nx1 matrix
	 */
	void set_center(const vector_t<Number>& center);

    /**
	 * @brief Get the current generators
	 * @return center a nxm matrix
	 */
	const matrix_t<Number>& generators() const;

    /**
	 * @brief Replaces the current matrix of generators with the parameter generators
	 * @param generators a nxm matrix
	 */
	void set_generators(const matrix_t<Number>& generators);

    /**
	 * @brief Add generators to Zonotope. Simply performs setGenerators if generators was previously not initialized.
	 * @param generators a nxm matrix
	 * @return true if able to add generators
	 */
	bool AddGenerators(const matrix_t<Number>& generators);

    /**
	 * @brief Get the order
	 * @return zonotope order
	 */
	Number order() const;

    /**
	 * @brief Number of generators
	 * @return number of generators
	 */
	size_t size() const;

    /**
	 * @brief Removes empty (null) columns in generator matrix
	 */
	void RemoveEmptyGenerators();

	/**
     * @brief It's important to do it, so we can reduce the necessary amount of calls of corners!
     */
	void UniteEqualVectors();

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
	 * @brief Clears the generators and center of the Zonotope and sets dimensionality to zero
	 */
	void Clear();

    /**
	 * @brief Print the zonotope
	 */
	void Print() const;

    /**
	 * @brief Judge whether two zonotope is equal
	 * @param another_zonotope 
	 * @return True, if equal
	 */
    bool operator==(const ZonotopeT<Number, Converter, Setting>& another_zonotope) const {
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
	*                           Algorithm Functions                             *
	*                                                                           *
	*****************************************************************************/
    /**
	 * @brief Get the enclosure for the convex hull of two zonotope
	 * @param another_zonotope 
	 * @return a  zonotope enclosing the convex hull
	 */
    Zonotope Convexhull( const Zonotope& another_zonotope ) const;
};

/** @} */

}  // namespace reachSolver
