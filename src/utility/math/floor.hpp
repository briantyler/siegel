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
 *  @file     floor.hpp
 *  @brief    This is a header implementation file for floor
 *  @note     Include this file to send floating point numbers to integer floor.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-17
 */


#ifndef _SG_FLOOR_H
#define _SG_FLOOR_H 1

#include <cmath>
#include <functional>

#include "utility/functors/stream_cast.hpp"
#include "utility/math/is_equal.hpp"
#include "utility/math/truncate.hpp"
#include "utility/precision.hpp"


namespace sg
{
namespace utility
{
namespace math
{
/**
 *  @struct  floor
 *  @brief   Take a floating point number to its floor and convert to an
 *           integer type.
 *  @author  Brian Tyler
 *  @version 1.0
 *  @date    2008-03-17
 *  @param   To The integer type to convert to.
 *  @param   From The floating point type to convert from.
 *  @note    Specialised for native and gmp types.
 */
template <class _To, class _From> struct floor
  : public std::unary_function<const _From&,_To>
{
  
  //! Object type
  typedef floor<_To,_From> self_type;
  //! Type to convert to
  typedef _To to_type;
  //!Type to convert from
  typedef _From from_type;
  
  /**
   *  @brief  Compute the floor of a float and convert to an integral type
   *  @param  value The floating point to find the flooring of
   *  @return The floor of __value cast to an integral type
   */
  _To operator( ) ( const _From& __value ) const {
    // By adding our precision zero we are basically saying that something
    // very close to the floor should stay at the floor
    _To floored = truncate<_To,_From>( )( __value + precision::zero() );
  
    if( __value < 0 ){
      _From tmp;
      functors::stream_cast<_To,_From>( )( floored, tmp );
      if( !math::is_equal<_From>( )( tmp, __value ) )  --floored;
    }
    return floored;
  }
};

template<> struct floor<long,double>
  : public std::unary_function<const double,long>
{
  typedef floor<long,double> self_type;
  typedef long to_type;
  typedef double from_type;
    
  long operator( )( const double __value ) const {
    return long( ::floor( __value + precision::zero() ) );
  }
};


#ifdef __GMP_PLUSPLUS__

template<> struct floor<mpz_class,double>
  : public std::unary_function<const double,mpz_class>
{
  typedef floor<mpz_class,double> self_type;
  typedef mpz_class to_type;
  typedef double from_type;
    
  mpz_class operator( )( const double __value ) const {
    return mpz_class( ::floor( __value + precision::zero() ) );
  }
};

template<> struct floor<long,mpf_class>
  : public std::unary_function<const mpf_class&,long>
{
  typedef floor<long,mpf_class> self_type;
  typedef long to_type;
  typedef mpf_class from_type;
  
  long operator( )( const mpf_class& __value ) const {
    return long( mpf_class( ::floor( __value + precision::zero() ) ).get_si() );
  }
};

template<> struct floor<mpz_class,mpf_class>
  : public std::unary_function<const mpf_class&,mpz_class>
{
  typedef floor<mpz_class,mpf_class> self_type;
  typedef mpz_class to_type;
  typedef mpf_class from_type;
  
  mpz_class operator( )( const mpf_class& __value ) const {
    return mpz_class( ::floor( __value + precision::zero() ) );
  }
};

template<> struct floor<long,mpq_class>
  : public std::unary_function<const mpq_class&,long>
{
  typedef floor<long,mpq_class> self_type;
  typedef long to_type;
  typedef mpq_class from_type;
  
  long operator( )( const mpq_class& __value ) const {
    return long( floor<long,mpf_class>( )( mpf_class( __value ) ) );
  }
};

template<> struct floor<mpz_class,mpq_class>
  : public std::unary_function<const mpq_class&,mpz_class>
{
  typedef floor<mpz_class,mpq_class> self_type;
  typedef mpz_class to_type;
  typedef mpq_class from_type;
  
  mpz_class operator( )( const mpq_class& __value ) const {
    return mpz_class( floor<mpz_class,mpf_class>( )( mpf_class( __value ) ) );
  }
};
#endif //__GMP_PLUSPLUS__


}// math
}// utility
}// sg

#endif
