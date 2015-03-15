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
 *  @file     point_inner_product.hpp
 *  @brief    This is a header implementation file for point_inner_product.
 *  @note     Include this file to calculate the hermitian inner product.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-05-11
 */

#ifndef _SG_POINT_INNER_PRODUCT_H
#define _SG_POINT_INNER_PRODUCT_H 1

// Local includes
#include "utility/math/hermitian_inner_product.hpp"


namespace sg
{
namespace utility
{
namespace math
{
  /**
   *  @brief  Computes the Hermitian inner product of two points
   *  @param  _Point The point class
   *  @param  lhs a point
   *  @param  rhs a point
   *  @return The hermitian inner product of two hyperbolic points
   */
template <class _Point>
    typename _Point::zeta_type
    point_inner_product( const _Point& __lhs, const _Point& __rhs)
{
  typedef typename _Point::zeta_type complex_type;
  
  return complex_type(   __lhs.dependent()() + std::conj( __rhs.dependent()() )
                       + hermitian_inner_product( __lhs.zeta_ref_begin(),
                                                  __lhs.zeta_ref_end(),
                                                  __rhs.zeta_ref_begin()
                                                )
                     );
}

}// algorithms
}// geometry
}// sg
#endif
