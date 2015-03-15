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
 *  @file     square.hpp
 *  @brief    This is a header implementation file for square
 *  @note     Template squareing functor
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-08
 */


#ifndef _SG_SQUARE_H
#define _SG_SQUARE_H 1

// Global includes
#include <functional>

#include <boost/call_traits.hpp>


namespace sg
{
namespace utility
{
namespace math
{
/**
 *  @struct   square
 *  @brief    A squareing functor
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-08
 *  @param    Type A numerical type point type
 */
template <class _Type> struct square
  : public std::unary_function
      < typename boost::call_traits<_Type>::param_type, _Type >
{
  //! Object type
  typedef square<_Type> self_type;
  //! The value type
  typedef _Type value_type;
  
  
  /**
   *  @brief  Squares an input value.
   *  @param  value The value to square
   *  @return \c value * \c value
   */
  value_type operator( ) ( typename self_type::argument_type __value ) const {
    return value_type( __value * __value );
  }
};

}// math
}// utility
}// sg
#endif
