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
 *  @file     hermitian_inner_product.hpp
 *  @brief    This is a header implementation file for \c hermitian_inner_product.
 *  @note     Include this file to calculate the hermitian inner product.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-11
 */

#ifndef _SG_HERMITIAN_INNER_PRODUCT_H
#define _SG_HERMITIAN_INNER_PRODUCT_H 1

// Global includes
#include <numeric>
#include <functional>
#include <boost/iterator/iterator_traits.hpp>

// Local includes
#include "utility/math/conjugate_multiply.hpp"


namespace sg
{
namespace utility
{
namespace math
{
/**
 *  @brief  Computes Hermitian inner product
 *  @param  first1     An input iterator.
 *  @param  last1      An input iterator.
 *  @param  first2     An input iterator.
 *  @return The hermitian inner product:
 *          \f$ \sigma lhs_i \times \overline{rhs}_i \f$
 *  @note   There is a significant speed reduction for non-inbuilt types, this
 *          is due to the fact that std::complex has been optimised for
 *          \c float and \c double. There is not much that can be done about
 *          this except for rewriting an optimised complex class, which is not
 *          something I really want to waste my time on at the moment
 *          (or probably ever). Same goes for all complex calculations.
 *
 *  Applies the operator to the corresponding elements in the two
 *  input ranges.
 *
 */
template <typename _InputIterator1, typename _InputIterator2>
    typename boost::iterator_value<_InputIterator1>::type
    hermitian_inner_product
    ( _InputIterator1 __first1, _InputIterator1 __last1,
      _InputIterator2 __first2 )
{
  typedef typename boost::iterator_value<_InputIterator1>::type complex_type;
  
  typedef conjugate_multiply<complex_type> _cm;
  typedef std::plus<complex_type> _pl;
  
  return complex_type( std::inner_product( __first1, __last1, __first2,
                                           complex_type(0.0,0.0), _pl(), _cm()
                                         )
                     );
}

}// math
}// utility
}// sg
#endif
