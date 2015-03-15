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
 *  @file     is_square.hpp
 *  @brief    An inline header file for the \c is_square class.
 *  @note     Checks whether an integer is a perfect square
 *  @author   Brian Tyler
 *  @version  1.1
 *  @date     2008-04-02
 *  @note     v1.1 made is_square into a function object
 */


#ifndef _SG_IS_SQUARE_H
#define _SG_IS_SQUARE_H 1

// Global includes
#include <functional>
#include <iostream>

#include <boost/call_traits.hpp>

// Local includes
#include "utility/functors/stream_cast.hpp"


namespace sg
{
namespace utility
{
namespace math
{
/**
 *  @class    is_square_candidate
 *  @brief    Checks whether an integer is not a perfect square
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-02
 *  @param    _Integer An integral type; defaults to \c long
 */
  template<class _Integer> struct is_square_candidate
  : public std::unary_function
    <typename boost::call_traits<_Integer>::param_type, bool>
{
  //! The object type
  typedef is_square_candidate<_Integer> self_type;
  //! The integral type
  typedef _Integer integer_type;
  
  /**
   *  @brief  Implements a well known algorithm for testing square candidacy
   *  @param  value The candidate
   *  @return \c false if \c value is NOT a perfect square
   *  @note   A return __value of true does not imply \c value ia a square.
   */
  bool operator( )( typename self_type::argument_type __value ) const {
    // Note: with modern compilers __value & 63UL is a total waste of time
    typedef utility::functors::stream_cast<integer_type, unsigned long> _sc;
    unsigned long tmp;
    _sc( )( __value % (64UL*63UL*65UL*11UL), tmp );
    if( !sq64( tmp % 64UL ) ) return false;
    if( !sq63( tmp % 63UL ) ) return false;
    if( !sq65( tmp % 65UL ) ) return false;
    return bool( sq11( tmp % 11UL ) );
  }

  private:
    static bool sq64( short int __value ) {
      static bool arr[]={
        1,1,0,0,1,0,0,0,0,1, 0,0,0,0,0,0,1,1,0,0, 0,0,0,0,0,1,0,0,0,0,
        0,0,0,1,0,0,1,0,0,0, 0,1,0,0,0,0,0,0,0,1, 0,0,0,0,0,0,0,1,0,0, 0,0,0,0};
      return bool( arr[__value] );
    }
    
    static bool sq63( short int __value ) {
      static bool arr[]={
        1,1,0,0,1,0,0,1,0,1, 0,0,0,0,0,0,1,0,1,0, 0,0,1,0,0,1,0,0,1,0,
        0,0,0,0,0,0,1,1,0,0, 0,0,0,1,0,0,1,0,0,1, 0,0,0,0,0,0,0,0,1,0, 0,0,0};
      return bool( arr[__value] );
    }
    
    static bool sq65( short int __value ) {
      static bool arr[]={
      1,1,0,0,1,0,0,0,0,1, 1,0,0,0,1,0,1,0,0,0, 0,0,0,0,0,1,1,0,0,1,
      1,0,0,0,0,1,1,0,0,1, 1,0,0,0,0,0,0,0,0,1, 0,1,0,0,0,1,1,0,0,0, 0,1,0,0,1};
      return bool( arr[__value] );
    }
    
    static bool sq11( short int __value ) {
      static bool arr[]={1,1,0,1,1,1,0,0,0,1, 0};
      return bool( arr[__value] );
    }
};


/**
 *  @class    is_square
 *  @brief    Checks whether an integer is a perfect square
 *  @author   Brian Tyler
 *  @version  1.1
 *  @date     2008-04-02
 *  @param    _Integer An integral type; defaults to \c long
 */
template<class _Integer> struct is_square
  : public std::binary_function
    <typename boost::call_traits<_Integer>::param_type, _Integer&, bool>
{
  //! The object type
  typedef is_square<_Integer> self_type;
  //! The integral type
  typedef _Integer integer_type;
  
  /**
   *  @brief  Tests an integer to see if it is a perfect square, if it is
   *          compute the root.
   *  @param  value The candidate.
   *  @param  sqrt The integer square root of \c value if the function returns
   *          true, otherwise it is undefined.
   *  @return True if \c value is a perfect square.
   */
  bool operator( )
      ( typename self_type::first_argument_type __value,
        typename self_type::second_argument_type __sqrt ) const
  {
    bool b = false;
    if( is_square_candidate<integer_type>( )( __value ) ) {
       // __value is a good square candidate
      __sqrt = __value;
      integer_type y = ( ( __sqrt + 1 ) / 2 );
      
      while( y < __sqrt ){
        __sqrt = y;
        y = ( ( __sqrt + __value/__sqrt ) /2 );
      }
      
      b = ( __sqrt * __sqrt == __value );
    }
    return b;
  }
  
  /**
   *  @brief  Tests an integer to see if it is a perfect square.
   *  @param  value The candidate.
   *  @return True if \c value is a perfect square.
   */
  bool operator( ) ( typename self_type::first_argument_type __value ) const {
    integer_type t;
    return bool ( operator()( __value, t ) );
  }
};


#ifdef __GMP_PLUSPLUS__
template<> struct is_square<mpz_class>
  : public std::binary_function
    <boost::call_traits<mpz_class>::param_type, mpz_class&, bool>
{
  //! The object type
  typedef is_square<mpz_class> self_type;
  //! The integral type
  typedef mpz_class integer_type;
  
  bool operator( )
      ( self_type::first_argument_type __value,
        self_type::second_argument_type __sqrt ) const
  {
    bool b = false;
    if( mpz_perfect_square_p( __value.get_mpz_t() ) ) {
      mpz_sqrt( __sqrt.get_mpz_t(), __value.get_mpz_t() );
      b = true;
    }
    return b;
  }
};
#endif //__GMP_PLUSPLUS__


}// math
}// utility
}// sg
#endif
