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
 *  @file     is_zero.hpp
 *  @brief    This is a header implementation file for \c is_zero.
 *  @note     Include this file to test whether a numerical value is close to
 *            zero
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-28
 */


#ifndef _SG_IS_ZERO_H
#define _SG_IS_ZERO_H 1

// Global includes
#include <cmath>
#include <functional>

// Local includes
#include "utility/precision.hpp"


namespace sg
{
namespace utility
{
namespace math
{
/**
 *  @class    is_zero
 *  @brief    A function object for determining a numerical value is close to
 *            zero
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-28
 *  @param    Type The type to check for equality with zero.
 *  @note     Since there is no need for a subtraction this is faster than
 *            using \c is_equal(value,0)
 */
template <class _Type> class is_zero
  : public std::unary_function< const _Type&, bool >
{
  public:
    typedef is_zero<_Type> self_type;
    typedef _Type value_type;
    
    /**
     *  @brief  Determines if a value is close to precision zero.
     *  @param  value The valuse to check for equality with precision::zero
     *  @return \c True if \c value is no larger than precision::zero.
     */
    bool operator( )( const _Type& __value ) const {
      
      static const _Type zero( precision::zero() );
      static const _Type mzero( -zero );
      
      return bool( ( mzero <= __value ) && ( __value <= zero ) );
    }
};


template <>
    class is_zero<double>
  :public std::unary_function< const double, bool >
{
  public:
    typedef is_zero<double> self_type;
    typedef double value_type;
    
    const bool operator( )( const double __value ) const {
      return bool( std::fabs( __value ) <= precision::zero() );
    }
};


template <>
    class is_zero<long double>
  :public std::unary_function< const long double, bool >
{
  public:
    typedef is_zero<long double> self_type;
    typedef long double value_type;
    
    bool operator( )( const long double __value ) const {
      static const long double zero( precision::zero() );
      
      return bool( std::fabs( __value ) <= zero );
    }
};


template <>
    class is_zero<int>
  :public std::unary_function< const int, bool >
{
  public:
    typedef is_zero<int> self_type;
    typedef int value_type;
    
    bool operator( )( const int __value ) const {
      return bool( __value == 0 );
    }
};


template <>
    class is_zero<long>
  :public std::unary_function< const long, bool >
{
  public:
    typedef is_zero<long> self_type;
    typedef long value_type;
    
    bool operator( )( const long __value ) const {
      return bool( __value == 0 );
    }
};


#ifdef __GMP_PLUSPLUS__
template <>
    class is_zero<mpz_class>
  : public std::unary_function< const mpz_class&, bool >
{
  public:
    typedef is_zero<mpz_class> self_type;
    typedef mpz_class value_type;
    
    bool operator( )( const mpz_class& __value ) const {
      static const mpz_class zero(0);
      return bool( __value == zero );
    }
};
#endif

}// functors
}// utility
}// sg
#endif
