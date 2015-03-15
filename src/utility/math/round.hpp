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
 *  @file     round.hpp
 *  @brief    This is a header implementation file for round
 *  @note     Include this file to round floating point numbers to integers.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-28
 */


#ifndef _SG_ROUND_H
#define _SG_ROUND_H 1

// Global includes
#include "utility/math/sgn.hpp"
#include "utility/functors/stream_cast.hpp"
#include "utility/precision.hpp"


namespace sg
{
namespace utility
{
namespace math
{
/**
 *  @struct  round
 *  @brief   Round a floating point number and convert to an integer type.
 *  @author  Brian Tyler
 *  @version 1.0
 *  @date    2008-03-17
 *  @param   To The integer type to convert to.
 *  @param   From The floating point type to convert from.
 *  @note    Specialised for native and gmp types.
 */
template <class _To, class _From> struct round
  : public std::unary_function<const _From&, _To>
{
  
  //! Object type
  typedef round<_To,_From> self_type;
  //! Type to convert to
  typedef _To to_type;
  //!Type to convert from
  typedef _From from_type;
  
  /**
   *  @brief  Round a float and convert to an integral type
   *  @param  value The floating point to round
   *  @return The rounded value cast to an integral type
   */
  _To operator( ) ( const _From& __value ) const {
    static const _From phalf(  0.5 + precision::zero() );
    static const _From nhalf( -0.5 + precision::zero() );
  
    _To rounded;
    if( __value < 0 ){
      functors::stream_cast<_From,_To>( _From( __value + phalf ), rounded );
    } else {
      functors::stream_cast<_From,_To>( _From( __value + nhalf ), rounded );
    }
    return rounded;
  }
};

template<> struct round<long,double>
  : public std::unary_function<const double,long>
{
  typedef round<long,double> self_type;
  typedef long to_type;
  typedef double from_type;
    
  long operator( )( const double __value ) const {
    static const double phalf(  0.5 + precision::zero() );
    static const double nhalf( -0.5 + precision::zero() );
    
    return long( __value < 0.0 ? __value + nhalf : __value + phalf );
  }
};

template<> struct round<long,long double>
  : public std::unary_function<const long double,long>
{
  typedef round<long,long double> self_type;
  typedef long to_type;
  typedef long double from_type;
    
  long operator( )( const long double __value ) const {
    static const long double phalf(  0.5 + precision::zero() );
    static const long double nhalf( -0.5 + precision::zero() );
    
    return long( __value < 0.0 ? __value + nhalf : __value + phalf );
  }
};


#ifdef __GMP_PLUSPLUS__

template<> struct round<mpz_class,double>
  : public std::unary_function<const double,mpz_class>
{
  typedef round<mpz_class,double> self_type;
  typedef mpz_class to_type;
  typedef double from_type;
    
  mpz_class operator( )( const double __value ) const {
    static const double phalf(  0.5 + precision::zero() );
    static const double nhalf( -0.5 + precision::zero() );
    
    return mpz_class( __value < 0.0 ? __value + nhalf : __value + phalf );
  }
};

template<> struct round<long,mpf_class>
  : public std::unary_function<const mpf_class&,long>
{
  typedef round<long,mpf_class> self_type;
  typedef long to_type;
  typedef mpf_class from_type;
  
  long operator( )( const mpf_class& __value ) const {
    static const mpf_class phalf(  0.5 + precision::zero() );
    static const mpf_class nhalf( -0.5 + precision::zero() );
  
    long rounded(0);
    switch ( sgn<mpf_class>()( __value ) ){
      case -1:
        rounded = mpf_class( __value + nhalf ).get_si();
        break;
      case 0:
        rounded = 0L;
        break;
      case 1:
        rounded = mpf_class( __value + phalf ).get_si();
        break;
    }
    return rounded;
  }
};

template<> struct round<mpz_class,mpf_class>
  : public std::unary_function<const mpf_class&,mpz_class>
{
  typedef round<mpz_class,mpf_class> self_type;
  typedef mpz_class to_type;
  typedef mpf_class from_type;
    
  mpz_class operator( )( const mpf_class& __value ) const {
    static const mpf_class phalf(  0.5 + precision::zero() );
    static const mpf_class nhalf( -0.5 + precision::zero() );
  
    mpz_class rounded(0);
    switch ( sgn<mpq_class>()( __value ) ){
      case -1:
        rounded = static_cast<mpz_class>( __value + nhalf );
        break;
      case 0:
        rounded = static_cast<mpz_class>(0);
        break;
      case 1:
        rounded = static_cast<mpz_class>( __value + phalf );
        break;
    }
    return rounded;
  }
};

template<> struct round<long,mpq_class>
  : public std::unary_function<const mpq_class&,long>
{
  typedef round<long,mpq_class> self_type;
  typedef long to_type;
  typedef mpq_class from_type;
    
  long operator( )( const mpq_class& __value )  const {
    static const mpq_class phalf(  0.5 + precision::zero() );
    static const mpq_class nhalf( -0.5 + precision::zero() );
  
    long rounded(0);
    switch ( sgn<mpq_class>()( __value ) ){
      case -1:{
        rounded = mpz_class( __value + nhalf ).get_si();
        break;
      }
      case 0:
        rounded = 0L;
        break;
        case 1:{
          rounded = mpz_class( __value + phalf ).get_si();
          break;
        }
    }
    return rounded;
  }
};

template<> struct round<mpz_class,mpq_class>
  : public std::unary_function<const mpq_class&,mpz_class>
{
  typedef round<mpz_class,mpq_class> self_type;
  typedef mpz_class to_type;
  typedef mpq_class from_type;
    
  mpz_class operator( )( const mpq_class& __value ) const {
    static const mpq_class phalf(  0.5 + precision::zero() );
    static const mpq_class nhalf( -0.5 + precision::zero() );
  
    mpz_class rounded(0);
    switch ( sgn<mpq_class>()( __value ) ){
      case -1:{
        rounded = static_cast<mpz_class>( __value + nhalf );
        break;
      }
      case 0:{
        rounded = static_cast<mpz_class>(0);
        break;
      }
      case 1:{
        rounded = static_cast<mpz_class>( __value + phalf );
        break;
      }
    }
    return rounded;
  }
};
#endif // __GMP_PLUSPLUS__

}// math
}// utility
}// sg
#endif
