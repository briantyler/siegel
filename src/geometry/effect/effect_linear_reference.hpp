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
 *  @file     effect_linear_reference.hpp
 *  @brief    An inline header file for the \c effect_linear class.
 *  @note     Include to compute the linearised effect function
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-23
 */


#ifndef _SG_RATIO_LINEAR_REFERENCE_H
#define _SG_RATIO_LINEAR_REFERENCE_H 1

// Global includes
#include <cassert>
#include <boost/call_traits.hpp>

// Local includes
#include "geometry/effect/effect.hpp"


namespace sg
{
namespace geometry
{
namespace effect
{
/**
 *  @struct   effect_linear_reference
 *  @brief    A unary function object that computes the gradient of the effect
 *            function at a point in space
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-23
 *  @param    Point The hyperbolic point type.
 *  @note     This function is inefficient from the point of minimisation,
 *            since it is called many times with the same gradient and base
 *            points. As such a certain amount of pre-computation can be done
 *            which reduces the amount of calculations performed on each loop.
 *            However I have left it in since it is obvious what it does, so
 *            the notion of linearisation is clear, but it should not be used
 *            in a build enviroment.
 */
template <class _Point> class effect_linear_reference
  : public std::unary_function<
          typename boost::call_traits<typename _Point::float_type>::param_type,
          typename _Point::float_type
    >
{
  public:
    //! Object type
    typedef effect_linear_reference<_Point> self_type;
    //! Point type
    typedef _Point point_type;
    //! The floating point type;
    typedef typename point_type::float_type float_type;
    
  private:
    //! Type of the effect function
    typedef effect<point_type> effect_type;
    
    
  public:
    /**
     *  @brief  Default constructor.
     *  @note   The base, current and gradient points must be bound prior to
     *          use.
     */
    effect_linear_reference()
  : effect_(), current_(0), gradient_(0)
    { }
    
    
    /**
     *  @brief  Bind to a new base point.
     *  @param  base The base point to bind to.
     */
    void bind_base( const point_type* __base ) {
      effect_.bind_base( __base );
    }
    
    /**
     *  @brief  Bind to a new gradient point.
     *  @param  gradient The gradient point to bind to.
     */
    void bind_gradient( const point_type* __gradient ) {
      gradient_ = __gradient;
    }
    
    /**
     *  @brief  Bind to a new current point.
     *  @param  current The current point to bind to.
     */
    void bind_current( const point_type* __current ) {
      current_ = __current;
    }
    
    /**
     *  @brief  Initialize the internal variables.
     *  @note   Call after base, current, or gradient change. This does
     *          nothing, it is just supplied for compatibility.
     */
    void initialize( void ) { }
    
    
    /**
     *  @brief  Compute the value of the function at \c lambda
     *  @param  lambda The argument of the function linearised effect function
     *  @return The value at \c lambda
     */
    float_type operator() ( typename type::argument_type __lambda ) const {
      // Ensure variables valid at debug time.
      assert( current_ != 0 && gradient_ != 0 );
      
      return float_type( effect_( (*current_) + __lambda * (*gradient_) ) );
    }
    
    
  private:
    //! Ratio function type
    effect_type effect_;
    //! Current position
    const point_type* current_;
    //! Gradient at a point
    const point_type* gradient_;
};

}// effect
}// geometry
}// sg
#endif
