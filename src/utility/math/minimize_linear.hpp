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
 *  @file     minimize_linear.hpp
 *  @brief    An inline header file for the \c minimize_linear class.
 *  @note     Include to compute the minima of a 1-d function using Brent's
 *            algorithm. The boost version is only valid for inbuilt types,
 *            this is extended for gmp types (and generally whenever the \c abs
 *            functor has been defined).
 *  @author   John Maddock
 *  @author   Brian Tyler (some parts)
 *  @version  1.0
 *  @date     2008-04-28
 *
 */

//  Some parts (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)


#ifndef SG_MINIMIZE_LINEAR_H
#define SG_MINIMIZE_LINEAR_H 1

// Global includes
#include <cstddef>
#include <utility>

// Local includes
#include "utility/math/abs.hpp"
#include "utility/precision.hpp"


namespace sg
{
namespace utility
{
namespace math
{
/**
 *  @struct minimize_linear
 *  @brief  Computes the minima of a function in one variable
 *  @param  Float The floating point type taken by the input function.
 *  @param  f The function to minimize.
 *  @param  min The lower bound on the bracket
 *  @param  max The upper bound on the bracket.
 *  @param  maxIteeffectns The maximum number of iteeffectns to perform
 *  @return A pair containing the value of the abscissa \c x at the minima and
 *          the value of \c f(x) at the minima.
 *
 *  Uses Brent's algorithm (w/o derivatives) to find the minima of \c f in
 *  \c [min,max]. A very important point to note is that this function
 *  guarantees to find a local minima in the bracketed region, but if there is
 *  more than one such minima, then there is no guarantee as to which is found.
 *
 */
template <class _Float> struct minimize_linear
{
  //! The object type
  typedef minimize_linear<_Float> self_type;
  //! The floating point type
  typedef _Float float_type;
  //! The return type
  typedef std::pair<float_type,float_type> return_type;
  
  template <class _UnaryOperation> return_type operator()
      ( _UnaryOperation __f, _Float __min, _Float __max,
        std::size_t __maxIteeffectns = 200 ) const
  {
    typedef abs<_Float> _abs;
    _Float x;  // minima so far
    _Float w;  // second best point
    _Float v;  // previous value of w
    _Float u;  // most recent evaluation point
    _Float delta;  // The distance moved in the last step
    _Float delta2; // The distance moved in the step before last
    _Float fu, fv, fw, fx;  // function evaluations at u, v, w, x
    _Float mid; // midpoint of min and max
    _Float fract1, fract2;  // minimal relative movement in x
    
    // zero tolerance
    static const _Float tolerance
        = static_cast<_Float>( precision::zero() );
    // golden effect, don't need too much precision here
    static const _Float golden = 0.3819660f;
    
    x = w = v = __max;
    fw = fv = fx = __f(x);
    delta2 = delta = 0.0;
    
    std::size_t count = __maxIteeffectns;
    
    do{
       // get midpoint
      mid = (__min + __max) / 2.0;
      // work out if we're done already:
      fract1 = tolerance * _abs()(x) + tolerance / 4.0;
      fract2 = 2.0 * fract1;
      
      if( _abs()( x - mid ) <= fract2 - ( __max - __min ) / 2.0 ) break;
      
      if(_abs()(delta2) > fract1 ) {
        // try and construct a parabolic fit:
        _Float r = (x - w) * (fx - fv);
        _Float q = (x - v) * (fx - fw);
        _Float p = (x - v) * q - (x - w) * r;
        q = 2 * (q - r);
        
        if(q > 0) p = -p;
        
        q = _abs()(q);
        _Float td = delta2;
        delta2 = delta;
        
        // determine whether a parabolic step is acceptible or not:
        if(    ( _abs()(p) >= _abs()(q * td / 2.0) )
            || (p <= q * (__min - x))
            || (p >= q * (__max - x))
          )
        {
          // nope, try golden section instead
          delta2 = ( x >= mid ) ? __min - x : __max - x;
          delta = golden * delta2;
        } else {
          // parabolic fit:
          delta = p / q;
          u = x + delta;
          if( ( (u - __min) < fract2 ) || ( (__max- u) < fract2 ) ) {
            delta =   (mid - x) < 0.0
                    ? static_cast<float_type>( -_abs()(fract1) )
                    : static_cast<float_type>(  _abs()(fract1) );
          }
        }
      }
      else {
        // golden section:
        delta2 =   ( x >= mid )
                 ? static_cast<float_type>( __min - x )
                 : static_cast<float_type>( __max - x );
        delta = golden * delta2;
      }
      // update current position:
      u =   ( _abs()(delta) >= fract1 )
          ? static_cast<float_type>( x + delta )
          : (   delta > 0.0
              ? static_cast<float_type>( x + _abs()(fract1) )
              : static_cast<float_type>( x - _abs()(fract1) )
            );
      
      fu = __f(u);
      
      if(fu <= fx) {
        // good new point is an improvement!
        // update brackets:
        if(u >= x) __min = x;
        else __max = x;
        
        // update control points:
        v = w;
        w = x;
        x = u;
        fv = fw;
        fw = fx;
        fx = fu;
      } else {
        // Oh dear, point u is worse than what we have already,
        // even so it *must* be better than one of our endpoints:
        if(u < x) __min = u;
        else __max = u;
        
        if( (fu <= fw) || (w == x) ) {
          // however it is at least second best:
          v = w;
          w = u;
          fv = fw;
          fw = fu;
        } else if( (fu <= fv) || (v == x) || (v == w) ) {
          // third best:
          v = u;
          fv = fu;
        }
      }
  
    } while(--count);
    
    _Float fmin( __f(__min) ), fmax( __f(__max) );
    
    
    // hit the end points
    return (   fmin < fmax
             ? fx < fmin ? std::make_pair(x, fx) : std::make_pair(__min, fmin)
             : fx < fmax ? std::make_pair(x, fx) : std::make_pair(__max, fmax)
           );
  }
};

}// math
}// utility
}// sg
#endif
