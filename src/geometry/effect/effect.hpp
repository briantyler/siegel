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
 *  @file     effect.hpp
 *  @brief    This is a header implementation file for the \c effect function.
 *  @note     Include this file to compute the effect function at a point.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-29
 */


#ifndef _SG_RATIO_H
#define _SG_RATIO_H 1

// Global includes
#include <cassert>
#include <functional>

// Local includes
#include "utility/math/point_inner_product.hpp"


namespace sg
{
namespace geometry
{
namespace effect
{
/**
 *  @class    effect
 *  @brief    A function object that computes the effect function at a point.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-29
 *  @param    Point The hyperbolic point type.
 */
template <class _Point> struct effect
  : public std::unary_function<const _Point&, typename _Point::float_type>
{
  //! Object type
  typedef effect<_Point> self_type;
  //! Point type
  typedef _Point point_type;
  //! Floating point type
  typedef typename _Point::float_type float_type;
  //! Complex type
  typedef typename _Point::zeta_type complex_type;
  
  /**
   *  @brief  Default constructor.
   *  @note   The functor must be bound to a base point prior to use.
   */
  effect( ) : base_(0) { }
  
  /**
   *  @brief  Base constructor
   *  @param  base The base point which generates the function
   */
  effect( const point_type* __base ) : base_(__base) { }
  
  
  /**
   *  @brief  Bind to a new base point.
   *  @param  base The base point to bind to.
   */
  void bind_base( const point_type* __base ) { base_ = __base; }
  
  
  /**
   *  @brief  Compute the effect function at \c current
   *  @param  current The current point
   *  @return The value of the effect function at current.
   */
  float_type operator() ( const point_type& __current ) const {
    using utility::math::point_inner_product;
    // During debugging ensure the base point is valid
    assert( base_ != 0 );
    
    return float_type( std::norm( point_inner_product<point_type>( __current,
                                                                   *base_
                                                                 )
                                )
                     );
  }
  
  
  private:
    //! The base point which generates the function
    const point_type* base_;
};

}// effect
}// geometry
}// sg
#endif
