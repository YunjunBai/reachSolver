/**
 * @file   BasicObject.h
 * @brief 
 * @author 
 * @date
 * @version 1.0
 * 
 * @defgroup structure
 * Reference:
 *   CORA ../contSet/contSet/contSet.m 
 */

#pragma once
#include <cstddef>
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