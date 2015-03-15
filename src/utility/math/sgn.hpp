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
 *  @file     sgn.hpp
 *  @brief    This is a header implementation file for sgn
 *  @note     Include this file to get the sign of a number
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-11
 */


#ifndef _SG_SGN_H
#define _SG_SGN_H 1

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
 *  @struct   sgn
 *  @brief    Computes the sign of a numerical value.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-11
 *  @param    Type The type of the number to get the sign of.
 */
template <class _Type> struct sgn
  : public std::unary_function< typename boost::call_traits<_Type>::param_type,
                                short>
{
  //! Object type
  typedef sgn<_Type> self_type;
  //! Value type
  typedef _Type value_type;
  
  /**
   *  @brief  Get the sign of an input value.
   *  @param  value The \c value to find the sign of.
   *  @return 1 if > 0, -1 if < 0, 0 otherwise.
   *  @note   This functor is not stable around zero; in other words values very
   *          close to zero are not considered to be zero, if this behaviour is
   *          required, then the \c is_zero functor should be used before
   *          employing the \c sgn functor.
   */
  short operator()( typename self_type::argument_type __value ) const {
    return short( __value > 0 ? 1 : ( __value < 0 ? -1 : 0 ) );
  }
};


#ifdef __GMP_PLUSPLUS__
template <> struct sgn<mpz_class>
  : public std::unary_function< const mpz_class&, short>
{
  //! Object type
  typedef sgn<mpz_class> self_type;
  //! Value type
  typedef mpz_class value_type;
  
  short operator()( const mpz_class& __value ) const {
    return short( ::sgn( __value ) );
  }
};


template <> struct sgn<mpf_class>
  : public std::unary_function< const mpf_class&, short>
{
  //! Object type
  typedef sgn<mpf_class> self_type;
  //! Value type
  typedef mpf_class value_type;
  
  short operator()( const mpf_class& __value ) const {
    return short( ::sgn( __value ) );
  }
};


template <> struct sgn<mpq_class>
  : public std::unary_function< const mpq_class&, short>
{
  //! Object type
  typedef sgn<mpq_class> self_type;
  //! Value type
  typedef mpq_class value_type;
  
  short operator()( const mpq_class& __value ) const {
    return short( ::sgn( __value ) );
  }
};
#endif // __GMP_PLUSPLUS__

}// math
}// utility
}// sg
#endif
