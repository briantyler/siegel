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
 *  @file     conjugate_multiply.hpp
 *  @brief    This is a header implementation file for conjugate_multiply.
 *  @note     Include this file to perform conjugate / hermitian multiplication
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-11
 */

#ifndef _SG_CONJUGATE_MULTIPLY_H
#define _SG_CONJUGATE_MULTIPLY_H 1

// Global includes
#include<functional>


namespace sg
{
namespace utility
{
namespace math
{
/**
 *  @struct   conjugate_multiply
 *  @brief    A binary function object that performs conjugate multiplication
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-11
 *  @param    _Cx A complex type.
 *
 *  Suppose \f$ z, y \in \mathbb{C} \f$. then this function object returns
 *  \f$ z\overline{y} \f$.
 */
template <class _Cx>
    struct conjugate_multiply
  : public std::binary_function<const _Cx&, const _Cx&, _Cx>
{
  //! Object type
  typedef conjugate_multiply<_Cx> self_type;
  //! Complex Type
  typedef _Cx complex_type;
  
  
  /**
   *  @brief  Performs Hermitian multiplication
   *  @param  lhs A complex number
   *  @param  rhs A complex number
   *  @return The hermitian product: \f$ lhs \times \overline{rhs} \f$
   */
  complex_type operator( )
      ( const complex_type& __lhs, const complex_type& __rhs ) const {
    return complex_type( __lhs * std::conj(__rhs) );
  }
};

}// math
}// utility
}// sg
#endif
