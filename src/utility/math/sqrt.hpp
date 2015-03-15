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
 *  @file     sqrt.hpp
 *  @brief    This is a header implementation file for sqrt
 *  @note     Include this file to compute real valued integer sqrts.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-03
 */


#ifndef _SG_SQRT_H
#define _SG_SQRT_H 1

// Global includes
#include <cmath>
#include <functional>

#include <boost/mpl/not.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_same.hpp>

// Local includes
#include "utility/functors/stream_cast.hpp"


namespace sg
{
namespace utility
{
namespace math
{
/**
 *  @struct  sqrt
 *  @brief    An integer to floating point square root functor.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-03
 *  @param   _To The floating type to convert to
 *  @param   _From The integer type to convert from.
 */
template <class _To, class _From> struct sqrt
  : public std::unary_function<const _From&, _To>
{
  //! Object type
  typedef sqrt<_To,_From> self_type;
  //! Type to convert to
  typedef _To to_type;
  //! Type to convert from
  typedef _From from_type;
  
  
  /**
   *  @brief  Square roots an integer and convert to a floating point type.
   *  @param  value The input value
   *  @return The square root of \c value converted to a floating point type.
   */
  to_type operator( )( const from_type& __value ) const {
    to_type output;
    double dreal;
    functors::stream_cast<from_type,double>( )( __value, dreal );
    functors::stream_cast<double,to_type>( )( ::sqrt( dreal ), output );
    return output;
  }
};


template <class _From> struct sqrt<long,_From>
{
  // Don't allow conversion to long - sounds like syntax error
  // The assert always evaluates to false, but it is not instantiated until the
  // functor is called.
  BOOST_MPL_ASSERT( ( boost::mpl::not_<boost::is_same<_From, _From> > ) );
};


template <> struct sqrt<double,long>
  : public std::unary_function<const long, double>
{
  typedef sqrt<double,long> self_type;
  typedef double to_type;
  typedef long from_type;
  
  double operator( )( const long __value ) const {
    return double( ::sqrt( static_cast<double>( __value ) ) );
  }
};


#ifdef __GMP_PLUSPLUS__
template <class _From> struct sqrt<mpz_class,_From>
{
  // Don't allow conversion to mpz_class - sounds like syntax error
  BOOST_MPL_ASSERT( ( boost::mpl::not_<boost::is_same<_From, _From> > ) );
};


template <> struct sqrt<double,mpz_class>
  : public std::unary_function<const mpz_class&, double>
{
  typedef sqrt<double,mpz_class> self_type;
  typedef double to_type;
  typedef mpz_class from_type;
  
  double operator( )( const mpz_class& __value ) const {
    return double( ::sqrt( __value.get_d() ) );
  }
};


template <> struct sqrt<mpf_class,long>
  : public std::unary_function<const long, mpf_class>
{
  typedef sqrt<mpf_class,long> self_type;
  typedef mpf_class to_type;
  typedef long from_type;
  
  mpf_class operator( )( const long __value ) const {
    return mpf_class( ::sqrt( static_cast<mpf_class>( __value ) ) );
  }
};


template <> struct sqrt<mpf_class,mpz_class>
  : public std::unary_function<const mpz_class&, mpf_class>
{
  typedef sqrt<mpf_class,mpz_class> self_type;
  typedef mpf_class to_type;
  typedef mpz_class from_type;
  
  mpf_class operator( )( const mpz_class& __value ) const {
    return mpf_class( ::sqrt( static_cast<mpf_class>( __value ) ) );
  }
};
#endif

}// math
}// utility
}// sg
#endif
