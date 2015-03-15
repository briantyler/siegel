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
 *  @file     hermitian_quadratic_product.hpp
 *  @brief    This is a header implementation file for
              \c hermitian_quadratic_product.
 *  @note     Include this file to compute the Hermitian quadratic product of a
 *            sequence.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-11
 */

#ifndef _SG_HERMITIAN_QUADRATIC_PRODUCT_H
#define _SG_HERMITIAN_QUADRATIC_PRODUCT_H 1

// Global includes
#include <algorithm>
#include <boost/iterator/iterator_traits.hpp>

// Local includes
#include "utility/math/norm_sum.hpp"


namespace sg
{
namespace utility
{
namespace math
{
/**
 *  @brief  Computes Hermitian quadratic product
 *  @param  first1 An input iterator.
 *  @param  last1 An input iterator.
 *  @return The hermitian quadratic product: \f$ \sigma |z_i|^2 \f$.
 */
template <typename _InputIterator1>
    typename boost::iterator_value<_InputIterator1>::type::value_type
    hermitian_quadratic_product
    (_InputIterator1 __first1,_InputIterator1 __last1)
{
  typedef typename boost::iterator_value<_InputIterator1>::type complex_type;
  typedef typename complex_type::value_type value_type;
  typedef norm_sum<complex_type> _ns;
  
  return value_type( std::for_each( __first1, __last1, _ns() ) );
}

}// math
}// utility
}// sg
#endif
