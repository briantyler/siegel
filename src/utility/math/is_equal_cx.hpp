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
 *  @file     is_equal_cx.hpp
 *  @brief    This is a header implementation file for is_equal_cx.
 *  @note     Include this file to test whether two complex values are close
 *            together.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-03
 */


#ifndef _SG_IS_EQUAL_CX_H
#define _SG_IS_EQUAL_CX_H 1

// Global includes
#include <complex>
#include <functional>

// Local includes
#include "utility/math/is_equal.hpp"


namespace sg
{
namespace utility
{
namespace math
{
/**
 *  @class    is_equal_cx
 *  @brief    A function object for determining if two numercal values are
 *            close.
 *  @author   Brian Tyler
 *  @version  1.1
 *  @date     2008-03-03
 *  @param    Float The floating point type taken by the complex number.
 */
template <class _Float> class is_equal_cx
  :std::binary_function< const std::complex<_Float>&,
                         const std::complex<_Float>&,
                         bool >
{
  public:
    //! The object type
    typedef is_equal_cx<_Float> self_type;
    //! The floating point type
    typedef _Float float_type;
    //! The complex type
    typedef std::complex<float_type> complex_type;
    
    
    /**
     *  @brief  Determine if two complex values are very close
     *          (ie computationally equivalent).
     *  @param  lhs The lefthand variable
     *  @param  rhs The righthand variable
     *  @return \c True if the real and imaginary parts of \c lhs and \c rhs
     *          differ by less than precision zero.
     */
    bool operator( )
        ( const complex_type& __lhs, const complex_type& __rhs ) const
    {
      return bool(    is_equal<_Float>()( __lhs.real(), __rhs.real() )
                   && is_equal<_Float>()( __lhs.imag(), __rhs.imag() )
                 );
    }
};

}// functors
}// utility
}// sg
#endif
