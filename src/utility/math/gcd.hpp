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
 *  @file     gcd.hpp
 *  @brief    An inline header file for the \c gcd class.
 *  @note     Include to compute the extended GCD of a pair of numbers
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-14
 */


#ifndef _SG_GCD_H
#define _SG_GCD_H 1

// Global includes
#include <functional>

// Local includes
#include "utility/math/abs.hpp"
#include "utility/math/sgn.hpp"


namespace sg
{
namespace utility
{
namespace math
{
/**
 *  @struct  gcd_solution
 *  @brief   The solution to a greatest common divisor system.
 *  @author  Brian Tyler
 *  @version 1.0
 *  @date    2008-04-14
 *  @param   _Integer An integral type; defaults to \c long
 */
template <class _Integer = long> struct gcd_solution {
  //! Object type
  typedef gcd_solution<_Integer> self_type;
    //! Integer type
  typedef _Integer integer_type;
  
  //! The coefficient of the first variable a1
  integer_type a;
  //! The coefficient of the second variable b1
  integer_type b;
  //! The gcd of the two variables
  integer_type gcd;
  
  // Note a*a1 + b*b1 = gcd
  
  /**
   *  @brief  Defines the output stream operator
   *  @param  __os An ostream.
   *  @param  __t The object to stream out.
   */
  friend std::ostream& operator<<( std::ostream& __os, const self_type& __t ){
    __os << '[' << __t.a << ',' << __t.b << ',' << __t.gcd << ']';
    return __os;
  }
};


/**
 *  @class    gcd
 *  @brief    Functor implementing Euclid's extended algorithm.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-14
 *  @param    _Integer An integral type; defaults to \c long
 */
template <class _Integer = long> class gcd
  :public std::binary_function<
            const _Integer&, // since we need abs value then it is ok to pass
            const _Integer&, // POD by reference
            gcd_solution<_Integer>
          >
{
  public:
    //! Object type
    typedef gcd<_Integer> self_type;
    //! Integer type
    typedef _Integer integer_type;
    // Types defined by binary_function
    //! First argument type
    typedef typename self_type::first_argument_type first_argument_type;
    //! Second argument type
    typedef typename self_type::second_argument_type second_argument_type;
    //! Solution type
    typedef typename self_type::result_type solution_type;
    //! Result type
    typedef typename self_type::result_type result_type;
    
    
    /**
     *  @brief  Compute the Extended GCD of \c a and \c b
     *  @param  a an integer.
     *  @param  b an integer.
     *  @return A gcd_solution {u,v,g} such that \f$ au +bv = g \f$ where
     *          \f$ g = \gcd(a,b) \f$
     *  @note   Cohen Algorithm 1.3.6
     */
    result_type operator( )
        ( first_argument_type __a, second_argument_type __b ) const
    {
      typedef sgn<integer_type> _sgn_fnc;
      result_type result;
      
      integer_type a( abs<integer_type>()(__a) ), b( abs<integer_type>()(__b) );
      
      result.a = 1;
      result.gcd = a;
      
      if( b == 0 ) {
        result.a = _sgn_fnc()(__a);
        result.b = 0;
        return result;
      }
      
      integer_type v1(0), v3( b ), q, t1, t3;
      
      while ( v3 != 0 ) {
        q = result.gcd / v3;
        t3 = result.gcd % v3;
        t1 = result.a - q*v1;
        result.a = v1;
        result.gcd = v3;
        v1 = t1;
        v3 = t3;
      }
      
      // We now have the gcd and bezout coefficients, it remains to compensate
      // for the sign of a and b
      result.b = ( result.gcd - a * result.a ) / b;
      
      short sa = _sgn_fnc()(__a);
      short sb = _sgn_fnc()(__b);
      
      if( sa == -1 && sb == -1 ) {
        result.b = -result.b;
        result.a = -result.a;
      }
      else if( sa != sb ) {
        if( __b*result.b > __a*result.a ) result.a = -result.a;
        else result.b = -result.b;
      }
      return result;
    }
};

}// math
}// utility
}// sg
#endif
