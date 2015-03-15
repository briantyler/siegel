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
 *  @file     is_even.hpp
 *  @brief    This is a header implementation file for is_even
 *  @note     Include this file to check if an integer is even or odd
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-09
 */


#ifndef _SG_IS_EVEN_H
#define _SG_IS_EVEN_H 1

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
 *  @class   is_even
 *  @brief   Determines if a number is even
 *  @author  Brian Tyler
 *  @version 1.0
 *  @date    2008-04-09
 *  @param   Integer The integer type to check the parity of.
 *  @note    This might seem slightly stupid, but since -1 % 2 = -1 there is
 *           sense in it
 */
template <class _Integer> struct is_even
  : public std::unary_function<
      typename boost::call_traits<_Integer>::param_type,
      bool>
{
  //! Object type
  typedef is_even<_Integer> self_type;
  //! Value type
  typedef _Integer value_type;
  
  /**
   *  @brief  Return the parity of an integer.
   *  @param  value THe integer to get the parity of.
   *  @return True if \c value is even.
   */
  bool operator()( typename self_type::argument_type __value ) const {
    return bool( (__value % 2) == 0 );
  }
};


/**
 *  @class   is_odd
 *  @brief   Determines if a number is odd
 *  @author  Brian Tyler
 *  @version 1.0
 *  @date    2008-04-09
 *  @param   Integer The integer type to check the parity of.
 */
template <class _Integer> struct is_odd
  : public std::unary_function<
      typename boost::call_traits<_Integer>::param_type,
      bool>
{
  //! Object type
  typedef is_even<_Integer> self_type;
  //! Value type
  typedef _Integer value_type;
  
  
  /**
   *  @brief  Return the parity of an integer.
   *  @param  value THe integer to get the parity of.
   *  @return True if \c value is odd.
   */
  bool operator()( typename self_type::argument_type __value ) const {
    return bool( !is_even<value_type>()(__value) );
  }
};

}// math
}// utility
}// sg
#endif
