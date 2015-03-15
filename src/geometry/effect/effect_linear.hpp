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
 *  @file     effect_linear.hpp
 *  @brief    An inline header file for the \c effect_linear class.
 *  @note     Include to compute the linearised effect function
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-23
 */


#ifndef _SG_RATIO_LINEAR_H
#define _SG_RATIO_LINEAR_H 1

// Global includes
#include <cassert>
#include <boost/call_traits.hpp>

// Local includes
#include "utility/math/hermitian_inner_product.hpp"


namespace sg
{
namespace geometry
{
namespace effect
{
/**
 *  @struct   effect_linear
 *  @brief    A unary function object that computes the gradient of the effect
 *            function at a point in space
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-23
 *  @param    Point The hyperbolic point type.
 */

template <class _Point> class effect_linear
  : public std::unary_function<
          typename boost::call_traits<typename _Point::float_type>::param_type,
          typename _Point::float_type
    >
{
  public:
    //! Object type
    typedef effect_linear<_Point> self_type;
    //! Point type
    typedef _Point point_type;
    //! The floating point type;
    typedef typename point_type::float_type float_type;
    //! The complex type;
    typedef typename point_type::zeta_type complex_type;
    
    
    /**
     *  @brief  Default constructor.
     *  @note   The base, current and gradient points must be bound prior to
     *          use.
     */
    effect_linear( )
  : base_(0), current_(0), gradient_(0) { }

    
    
    /**
     *  @brief  Bind to a new base point.
     *  @param  base The base point to bind to.
     */
    void bind_base( const point_type* __base ) {
      base_ = __base;
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
     *  @note   Call after base, current, or gradient change.
     */
    void initialize( void ) {
      using utility::math::hermitian_inner_product;
      
      // During debugging ensure the base point is valid
      assert( base_ != 0 );
      
      a_ = -(gradient_->dependent().real_ref());
      
      complex_type bg =
          hermitian_inner_product( base_->zeta_ref_begin(),
                                   base_->zeta_ref_end(),
                                   gradient_->zeta_ref_begin()
                                 );
      complex_type bc =
          hermitian_inner_product( base_->zeta_ref_begin(),
                                   base_->zeta_ref_end(),
                                   current_->zeta_ref_begin()
                                 );
      complex_type cg =
          hermitian_inner_product( current_->zeta_ref_begin(),
                                   current_->zeta_ref_end(),
                                   gradient_->zeta_ref_begin()
                                 );
      
      b_ = cg.real() - bg.real();
      c_ = - base_->dependent().real_ref() - current_->dependent().real_ref()
           - bc.real();
      
      s_ = gradient_->r_ref() - bg.imag();
      t_ = current_->r_ref() - base_->r_ref() - bc.imag();
    }
    
    
    /**
     *  @brief  Compute the value of the function at \c lambda
     *  @param  lambda The argument of the function linearised effect function
     *  @return The value at \c lambda
     */
    float_type operator() ( typename self_type::argument_type __lambda ) const {
      // During debugging ensure the base point is valid
      assert( base_ != 0 );
      
      float_type re = __lambda * ( __lambda * a_ + b_ ) + c_;
      float_type im = __lambda * s_ + t_;
    
      return float_type( re*re + im*im );
    }
    
    
    // Access the local variables
    const float_type& a( void ) { return a_; }
    const float_type& b( void ) { return b_; }
    const float_type& c( void ) { return c_; }
    const float_type& s( void ) { return s_; }
    const float_type& t( void ) { return t_; }
    
    
  private:
    //! The base point
    const point_type* base_;
    //! Current position
    const point_type* current_;
    //! Gradient at a point
    const point_type* gradient_;
    
    // Local variables
    float_type a_, b_, c_, s_, t_;
};

}// effect
}// geometry
}// sg
#endif
