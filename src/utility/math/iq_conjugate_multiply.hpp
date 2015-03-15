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
 *  @file     iq_conjugate_multiply.hpp
 *  @brief    This is a header implementation file for \c iq_conjugate_multiply.
 *  @note     Include this file to perform conjugate / hermitian multiplication on
 *            iq_numbers
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-17
 */

#ifndef _SG_IQ_CONJUGATE_MULTIPLY_H
#define _SG_IQ_CONJUGATE_MULTIPLY_H 1

// Global includes
#include <functional>


namespace sg
{
namespace utility
{
namespace math
{
/**
 *  @struct   iq_conjugate_multiply
 *  @brief    A binary function object that performs conjugate multiplication
 *            on iq_numbers.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-17
 *  @param    _IqNumber An iq_number type
 */
template <class _IqNumber>
    struct iq_conjugate_multiply
  : public std::binary_function<const _IqNumber&, const _IqNumber&, _IqNumber>
{
  //! Object type
  typedef iq_conjugate_multiply<_IqNumber> self_type;
  //! IqNumber type
  typedef _IqNumber iq_number_type;
  
  
  /**
   *  @brief  Performs Hermitian multiplication
   *  @param  lhs An \c iq_number
   *  @param  rhs An \c iq_number
   *  @return The hermitian product: \f$ lhs \times \overline{rhs} \f$
   */
  iq_number_type operator( )
      ( const iq_number_type& __lhs, const iq_number_type& __rhs ) const {
    return iq_number_type( __lhs * __rhs.conj() );
  }
};

}// math
}// utility
}// sg
#endif
