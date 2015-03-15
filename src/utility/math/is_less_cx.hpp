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
 *  @file     is_less_cx.hpp
 *  @brief    This is a header implementation file for is_less_cx
 *  @note     Include this file to implement a strict weak ordering on complex
 *            numbers: This is solely for the purpose of being able to use
 *            them in std::set containers.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-17
 */


#ifndef _SG_IS_LESS_CX_H
#define _SG_IS_LESS_CX_H 1

// Global includes
#include <functional>
#include <complex>

#include <boost/call_traits.hpp>

// Local includes
#include "utility/math/is_less.hpp"


namespace sg
{
namespace utility
{
namespace math
{
/**
   *  @class    is_less_cx
   *  @brief    A partial ordering functor for complex types.
   *  @author   Brian Tyler
   *  @version  1.0
   *  @date     2008-03-17
   *  @param    Type The floating point type taken by the complex number.
 */
template <class _Float> struct is_less_cx
  : public std::binary_function<const std::complex<_Float>&,
                                const std::complex<_Float>&,
                                bool>
{
  //! Object type
  typedef is_less_cx<_Float> self_type;
  //! The floating point type
  typedef _Float float_type;
  //! The complex type
  typedef std::complex<float_type> complex_type;
  
  
  /**
   *  @brief  A partial ordering for complex types.
   *  @param  lhs The lefthand complex number
   *  @param  rhs The righthand complex number
   *  @return \c true if \c lhs and \c rhs are partially ordered lhs < rhs.
   *  @note   The ordering says that points which are very close (according to
   *          the precision setting) are considered to be the same.
   */
  bool operator( )
      ( const complex_type& __lhs, const complex_type& __rhs ) const
  {
    typedef is_less<float_type> _il;
    typedef is_equal<float_type> _ie;
    
    return bool(    _il()( __lhs.real(), __rhs.real() )
                 || (    _ie()( __lhs.real(), __rhs.real() )
                      && _il()( __lhs.imag(), __rhs.imag() )
                    )
               );
  }
};

}// math
}// utility
}// sg
#endif
