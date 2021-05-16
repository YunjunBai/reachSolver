/**
 * @file   BasicObject.h
 * @brief  Basic class of these set representations
 * @author Yunjun Bai
 * @date April 2021
 * @version 1.0
 * 
 * @defgroup structure
 * Reference:
 *   CORA ../contSet/contSet/contSet.m 
 */

#pragma once
#include <cstddef>
#include <iostream>
#include "commonType.h"

namespace reachSolver {

/**
* @class   BasicObject
* @brief 
* @ingroup structure
* @{
*
* detailed description for the class.
*/
class BasicObject{

  protected:

    BasicObject() = default;
    BasicObject(const BasicObject& in) = default;
    BasicObject(BasicObject& in) = default;

    BasicObject& operator=(const BasicObject& in) = default;
    BasicObject& operator=(BasicObject&& in) = default;
};
/** @} */

}  // namespace reachSolver