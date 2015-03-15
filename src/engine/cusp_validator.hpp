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
 *  @file     cusp_validator.hpp
 *  @brief    An inline header file for the \c cusp_validator class.
 *  @note     The engine for validating cusps.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-05-02
 */


#ifndef _SG_CUSP_VALIDATOR_H
#define _SG_CUSP_VALIDATOR_H 1

// Global includes
#include <cassert>
#include <cstddef>
#include <cmath>

// Local includes
#include "structure/geometric/hyperbolic/cusp.hpp"
#include "geometry/effect/minimize_effect.hpp"
#include "utility/math/is_less_equal.hpp"


namespace sg
{
namespace engine
{
/**
 *  @class    cusp_validator
 *  @brief    An engine for validating cusps to determine if they are effective
 *            on some region.
 *  @param    N The dimension of the space
 *  @param    _Float A floating point type; defaults to \c double
 *  @param    _Integer An integral type; defaults to \c long
 *  @param    _Id The identifier of the field
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-05-02
 */
template < std::size_t N, class _Float = double, class _Integer = long,
           std::size_t _Id = 0> class cusp_validator
{
  public:
    //! The object type
    typedef cusp_validator<N,_Float,_Integer,_Id> self_type;
    //! The float type
    typedef _Float float_type;
    //! The integral type
    typedef _Integer integer_type;
    //! The cusp type
    typedef structure::geometric::hyperbolic::cusp
        <N,float_type,integer_type,_Id> cusp_type;
    //! The point type
    typedef typename cusp_type::point_type point_type;
    //! The field type
    typedef typename cusp_type::field_type field_type;
    
  private:
    //! Less than functor
    typedef utility::math::is_less_equal<float_type> less_eq_fnc;
    //! The type of the effect minimization function
    typedef geometry::effect::minimize_effect<N,float_type> min_effect_fnc;
    
  public:
    //! Type of the bounding space
    typedef typename min_effect_fnc::space_type space_type;
    
    
    /**
     *  @brief  Default constructor
     *  @note   Validator must be bound to a cusp and space prior to use.
     */
    cusp_validator( )
  : cusp_(0), minimizer_( ) { }
    
    
    /**
     *  @brief  Validation operator
     *  @return True if the cusp is capable of affecting the space
     */
    bool operator() ( void ) {
      assert( cusp_ != 0 );
      
      // Note that as long as the final Siegel height is significantly larger
      // than the zero tollerance then this will give all the cusps which
      // contribute to the fundamental domain.
      return bool( less_eq_fnc()( minimizer_(), cusp_->threshold() ) );
    }
    
    
    /**
     *  @brief  Bind to a new cusp.
     *  @param  cusp The cusp to bind to.
     */
    void bind_cusp( const cusp_type& __cusp ) {
      cusp_ = &__cusp;
      minimizer_.bind_base( &( cusp_->point() ) );
    }
    
    
    /**
     *  @brief  Bind to a new space.
     *  @param  space The space to bind to.
     */
    void bind_space( const space_type& __space ) {
      minimizer_.bind_space( __space );
    }
    
    
    /**
     *  @brief  Compute the effect of the cusp on the space.
     *  @return The cusps effect on the space.
     *  @note   The greater the effect, the larger the region in the space that
     *          is repelled by the cusp. Result is only valid after calling the
     *          operator() function.
     */
    float_type effect( void ) const {
      float_type e( cusp_->threshold() - minimizer_.minimum() );
      
      if( e > 0.0 ) e = 2.0 * ::sqrt( e );
      
      return e;
    }
    
    
  private:
    //! Constant reference to the cusp to be minimized
    const cusp_type* cusp_;
    //! The minimizer
    min_effect_fnc minimizer_;
};

}// engine
}// sg
#endif
