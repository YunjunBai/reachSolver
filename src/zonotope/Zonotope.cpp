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
     Zonotope<Number,Converter,Setting>::Zonotope()
        :dimension_(0),center_((0,1)),generators_(0,0;0,0) { }

    /**
	 * @brief Constructor with dimension
	 * @param dimension Dimensionality of Zonotope
	 */
    template <typename Number, typename Converter, typename Setting>
	Zonotope<Number,Converter,Setting>:: Zonotope(size_t dimension):
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
        dimension_(center.row()),
        center_(center),
        generators_(generators){
        assert(center.row()==generators.row());
    }

    template <typename Number, typename Converter, typename Setting>
	Zonotope<Number,Converter,Setting>:: ~Zonotope(){}


    
} // namespace reachSolver

