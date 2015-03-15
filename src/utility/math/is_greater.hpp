/*
 * ### autorights BEGIN ###
 * siegel / siegelcl is an application that computes Siegel sets and cusp
 * lists for the action of SU(n,1;O) on HCn. Its development was supported
 * by EPSRC and University College London.
 * 
 * Copyright (C) 2009  Brian Tyler
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ### autorights END ###
 */

/**
 *  @file     is_greater.hpp
 *  @brief    This is a header implementation file for is_greater
 *  @note     Safe floating point greater than comparison
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-08
 */


#ifndef _SG_IS_GREATER_H
#define _SG_IS_GREATER_H 1

// Global includes
#include <functional>

#include <boost/call_traits.hpp>

// Local includes
#include "utility/math/is_equal.hpp"


namespace sg
{
namespace utility
{
namespace math
{
/**
 *  @class    is_greater
 *  @brief    A greater than functor for floating point types
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-03
 *  @param    Type The floating point type to compare.
 */
template <class _Type> struct is_greater
  : public std::binary_function<typename boost::call_traits<_Type>::param_type,
                                typename boost::call_traits<_Type>::param_type,
                                bool>
{
  //! Object type
  typedef is_greater<_Type> self_type;
  //! The value type
  typedef _Type value_type;
  
  
  /**
   *  @brief  Determine if one floating point value is strictly greater than
   *          another, accounting for machine precision.
   *  @param  lhs The lefthand floating point number
   *  @param  rhs The righthand floating point number
   *  @return \c True if lhs > rhs + precision::zero.
   */
  bool operator( )
      ( typename self_type::first_argument_type __lhs,
        typename self_type::second_argument_type __rhs ) const
  {
    typedef is_equal<value_type> _ie;
    
    return bool( __lhs > __rhs && !_ie( )( __lhs, __rhs ) );
  }
};

}// math
}// utility
}// sg
#endif
