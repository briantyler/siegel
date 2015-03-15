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
 *  @file     abs.hpp
 *  @brief    This is a header implementation file for abs
 *  @note     Include this file to get the absolute value of a number
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-17
 */


#ifndef _SG_ABS_H
#define _SG_ABS_H 1

// Global includes
#include <cmath>
#include <functional>

#include <boost/call_traits.hpp>

// Local includes
#include "utility/math/sgn.hpp"


namespace sg
{
namespace utility
{
namespace math
{
/**
 *  @brief  Computes the absolute value of a numerical value.
 *  @param  Type The type to take the absolute value of:
 *          Specialised for \c double and the gmp classes
 *  @param  value  The value to find the absolute value of.
 *  @return The absolute value of \c value.
 */
template <class _Type> struct abs
  : public std::unary_function< typename boost::call_traits<_Type>::param_type,
                                _Type>
{
  //! Object type
  typedef abs<_Type> self_type;
  //! Value type
  typedef _Type value_type;
  
  value_type operator()( typename self_type::argument_type __value ) const {
    return value_type( sgn<value_type>()(__value) == -1 ? -__value : __value );
  }
};


template <> struct abs<double>
  : public std::unary_function< const double, double >
{
  //! Object type
  typedef abs<double> self_type;
  //! Value type
  typedef double value_type;
  
  const value_type operator()( const double __value ) const {
    return value_type( ::fabs( __value ) );
  }
};


#ifdef __GMP_PLUSPLUS__
template <> struct abs<mpz_class>
  : public std::unary_function< const mpz_class&, mpz_class >
{
  //! Object type
  typedef abs<mpz_class> self_type;
  //! Value type
  typedef mpz_class value_type;
  
  value_type operator()( const mpz_class& __value ) const {
    return value_type( ::abs( __value ) );
  }
};


template <> struct abs<mpf_class>
  : public std::unary_function< const mpf_class&, mpf_class >
{
  //! Object type
  typedef abs<mpf_class> self_type;
  //! Value type
  typedef mpf_class value_type;
  
  value_type operator()( const mpf_class& __value ) const {
    return value_type( ::abs( __value ) );
  }
};


template <> struct abs<mpq_class>
  : public std::unary_function< const mpq_class&, mpq_class >
{
  //! Object type
  typedef abs<mpq_class> self_type;
  //! Value type
  typedef mpq_class value_type;
  
  value_type operator()( const mpq_class& __value ) const {
    return value_type( ::abs( __value ) );
  }
};
#endif

}// math
}// utility
}// sg
#endif
