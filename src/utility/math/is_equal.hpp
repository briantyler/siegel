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
 *  @file     is_equal.hpp
 *  @brief    This is a header implementation file for is_equal.
 *  @note     Include this file to test whether two numerical values are close
 *            together.
 *  @author   Brian Tyler
 *  @version  1.1
 *  @note     Made \c is_equal into a function object
 *  @date     2008-03-03
 */


#ifndef _SG_IS_EQUAL_H
#define _SG_IS_EQUAL_H 1

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
 *  @class    is_equal
 *  @brief    A function object for determining if two numercal values are
 *            close.
 *  @author   Brian Tyler
 *  @version  1.1
 *  @date     2008-03-03
 *  @param    Type The floating point type to check equality on.
 */
template <class _Type> class is_equal
  : public std::binary_function< const _Type&, const _Type&, bool >
{
  public:
    typedef is_equal<_Type> self_type;
    typedef _Type value_type;
    
    /**
     *  @brief  Determine if two floating point values are very close
     *          (ie computationally equivalent).
     *  @param  lhs The lefthand variable
     *  @param  rhs The righthand variable
     *  @return \c True if \c lhs and \c rhs differ by less than precision zero.
     */
    bool operator( )( const _Type& __lhs, const _Type& __rhs ) const {
      static const _Type zero( precision::zero() );
      static const _Type mzero( -zero );
      const _Type tmp( __lhs - __rhs );
      
      return bool( ( mzero <= tmp ) && ( tmp <= zero ) );
    }
};


template <>
    class is_equal<double>
  :public std::binary_function< const double, const double, bool >
{
  public:
    typedef is_equal<double> self_type;
    typedef double value_type;
    
    bool operator( )( const double __lhs, const double __rhs ) const {
      return bool( std::fabs( __lhs - __rhs ) <= precision::zero() );
    }
};


template <>
    class is_equal<long double>
  :public std::binary_function< const long double, const long double,
                                bool >
{
  public:
    typedef is_equal<long double> self_type;
    typedef long double value_type;
    
    bool operator( )
        ( const long double __lhs, const long double __rhs ) const
    {
      static const long double zero( precision::zero() );
      return bool( std::fabs( __lhs - __rhs ) <= zero );
    }
};


template <>
    class is_equal<int>
  : public std::binary_function< const int, const int, bool >
{
  public:
    typedef is_equal<int> self_type;
    typedef int value_type;
    
    bool operator( )( const int __lhs, const int __rhs ) const {
      return bool( __lhs == __rhs );
    }
};


template <>
    class is_equal<long>
  : public std::binary_function< const long, const long, bool >
{
  public:
    typedef is_equal<long> self_type;
    typedef long value_type;
    
    bool operator( )( const long __lhs, const long __rhs ) const {
      return bool( __lhs == __rhs );
    }
};


#ifdef __GMP_PLUSPLUS__
template <>
    class is_equal<mpz_class>
  : public std::binary_function< const mpz_class&, const mpz_class&, bool >
{
  public:
    typedef is_equal<mpz_class> self_type;
    typedef mpz_class value_type;
    
    bool operator( )( const mpz_class& __lhs, const mpz_class& __rhs ) const {
      return bool( __lhs == __rhs );
    }
};
#endif


}// functors
}// utility
}// sg

#endif
