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
 *  @file     minimize_effect.hpp
 *  @brief    An inline header file for the \c minimize_effect class.
 *  @note     Include to compute the minima of the multi dimensional effect
 *            function. This is a fairly general multidimensional algorithm
 *            (without derivatives) however it has been specialised to minimize
 *            the effect function and moreover, to perform a large number of
 *            minimizations.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-28
 *
 */

#ifndef SG_MINIMISE_RATIO_H
#define SG_MINIMISE_RATIO_H 1

// Global includes
#include <cassert>
#include <cstddef>
#include <utility>

// Local includes
#include "utility/precision.hpp"
#include "utility/math/abs.hpp"
#include "utility/math/is_zero.hpp"
#include "utility/math/is_equal.hpp"
#include "utility/math/sgn.hpp"
#include "utility/math/minimize_linear.hpp"

#include "geometry/effect/effect.hpp"
#include "geometry/effect/effect_linear.hpp"
#include "structure/geometric/hyperbolic/hyperbolic_space.hpp"
#include "structure/geometric/hyperbolic/hyperbolic_point.hpp"


namespace sg
{
namespace geometry
{
namespace effect
{
/**
 *  @class    minimize_effect
 *  @brief    Minimizes the effect function on a compact Heisenberg region.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-28
 *  @param    N The hyperbolic dimension of the space to minimise over.
 *  @param    Float The floating point type.
 */
template <std::size_t N, class _Float = double> class minimize_effect
{
  public:
    //! Object type
    typedef minimize_effect<N,_Float> self_type;
    //! The type of the point input variable
    typedef _Float float_type;
    //! The type of the bounding space
    typedef structure::geometric::hyperbolic::hyperbolic_space<N,float_type>
        space_type;
    //! The type of the bounding space
    typedef structure::geometric::hyperbolic::hyperbolic_point<N,float_type>
        point_type;
    //! The function type
    typedef effect<point_type> effect_type;
    
    
  private:
    //! The type of the function restricted to one dimension
    typedef effect_linear<point_type> line_type;
    //! The type of the function which minimizes the 1-D restriction
    typedef utility::math::minimize_linear<float_type> line_min_type;
    //! Zero functor
    typedef utility::math::is_zero<float_type> _isz_fnc;
    //! Equality functor
    typedef utility::math::is_equal<float_type> _eq_fnc;
    //! Sign functor
    typedef utility::math::sgn<float_type> _sgn_fnc;
    //! Type of a pair of floats
    typedef std::pair<float_type,float_type> pair_float_type;
    //! Point heisenberg iterator
    typedef typename point_type::iterator pt_iterator;
    //! Constant point heisenberg iterator
    typedef typename point_type::const_iterator pt_const_iterator;
    //! Constant space heisenberg iterator
    typedef typename space_type::const_iterator sp_const_iterator;
    //! Type of a boost array of points
    typedef boost::array<point_type,point_type::heisenberg_size> point_array;
    
  public:
    /**
     *  @brief  Default constructor
     *  @note   The functor must be bound to both a base point and a space
     *          prior to use.
     */
    minimize_effect( )
  : space_(0), base_(0), current_(), gradientArr_(), minimum_(), effect_(),
    line_(), lineMin_()
    { line_.bind_current( &current_ ); }
    
    
    /**
     *  @brief  Space and Base point constructor
     *  @param  space The space which bounds the problem
     *  @param  base The base point which generates the effect function
     */
    minimize_effect( const space_type& __space, const point_type& __base )
  : space_(&__space), base_(&__base), current_(), gradientArr_(), minimum_(),
    effect_(), line_(), lineMin_()
    { bind_base(base_), line_.bind_current( &current_ ); }
    
    
    /**
     *  @brief  Bind to a new base point.
     *  @param  base The base point to bind to.
     */
    void bind_base( const point_type* __base ) {
      base_ = __base;
      effect_.bind_base( __base );
      line_.bind_base( __base );
    }
    
    /**
     *  @brief  Bind to a new space.
     *  @param  space The space to bind to.
     */
    void bind_space( const space_type& __space ) { space_ = &__space; }
    
    
    /**
     *  @brief  Find the minimum of the effect point in the space
     *  @param  maxIteeffectns The maximum number of iteeffectns to go through
     *  @return The minimum value.
     */
    const float_type& operator() ( std::size_t __maxIterations = 200 ) {
      // During debugging ensure the base point and space are valid
      assert( base_ != 0 && space_ != 0 );
      
      // Check if the base point is already in the heisenberg space, if so it
      // is its own minimum.
      if(    space_->contains( *base_ )
          && _eq_fnc()( base_->height(), space_->height()().lower() )
        )
      {
        // Nothing to do since the base point is in the space
        current_ = *base_;
        minimum_ = 0.0;
        
        return minimum_;
      }
      
      // Initialize the functor
      sys_initialize();
      
      // Loop until all iteeffectns are used up, or the minimum is reached
      std::size_t count = __maxIterations;
      
      // Implementation of the conjugate gradient method
      // Note that although we would like to use our knowledge of the actual
      // gradient, due to the restricitons imposed by the boundary, it turns
      // out not be very helpful, and an implementation of a direction set
      // method is quite efficient.
      do {
        // Numerical Recipies states that gradients should be reset after every
        // N + 1 iteeffectns (where N is the dimension of the problem)
        sys_reset_gradients();
        
        float_type min = minimum_;
        for( int i = 0; i != point_type::heisenberg_size; ++i ) {
          sys_conjugate_gradient_loop();
        }
        if( _eq_fnc()( minimum_, min ) ) break;
        
      } while (--count);
      
      return minimum_;
    }
    
    
    /**
     *  @brief  Get a constant reference to the current minimum.
     *  @return A constant reference to the current minimum
     */
    const float_type& minimum( void ) const { return minimum_; }
    
    /**
     *  @brief  Get a constant reference to the current closest point
     *  @return A constant reference to the current closest point
     */
    const point_type& current( void ) const { return current_; }
    
    
  private:
        /**
   *  @brief Initialize the minimizer
         */
    void sys_initialize( void ) {
      // First put the current point somewhere in the space: the closest
      // Euclidean point is a decent guess
      pt_iterator ptIt = current_.begin();
      pt_const_iterator bIt = base_->begin();
      sp_const_iterator spIt = space_->begin();
      for( ; ptIt != current_.heisenberg_end(); ++ptIt, ++spIt ) {
        *ptIt = spIt->closest(*bIt);
      }
      current_.height() = space_->height_ref().lower();
      current_.initialize();
      
      // Then set the starting minimum as the value at current
      minimum_ = effect_( current_ );
    }
    
    
    /**
     *  @brief  Perform a loop of the congjugate gradient method
     */
    void sys_conjugate_gradient_loop( void ) {
      // When reaching the boundary the gradient information can sometimes be
      // misleading as under projection the gradient may point us directly into
      // the wall of a valley
      
      point_type currentTmp = current_;
      
      for( typename point_array::iterator gaIt = gradientArr_.begin();
           gaIt != gradientArr_.end(); ++gaIt )
      {
        // Bind the linear restriction to the new direction vector.
        line_.bind_gradient( gaIt );
        
        // Compute the bounds on lambda (displacement)
        sys_compute_lambda_bounds( gaIt );
        
          // Compute the minimum between these bounds
        line_.initialize();
        lambdaPair_
            = lineMin_( line_, lambdaBounds_.first, lambdaBounds_.second );
        
        // If the displacement is non-zero then move to the new point.
        if( !_isz_fnc()( lambdaPair_.first ) ) {
          current_ += lambdaPair_.first * (*gaIt);
          // Very occasionally current will move slightly ouside the space,
          // when this happens we bring it back by moving to the closest point
          // in the space under a Euclidean metric. Due to the very small
          // margins of error this doesn't invalidate the minimisation.
          sys_correct_current();
        }
        
        // Update the direction set (and repeat).
        *gaIt = (   gaIt != gradientArr_.end() - 1
                  ? *(gaIt + 1) : current_ - currentTmp
                );
      }
      
      // Note that although the linearisation of the effect function is
      // analytically equivalent to the effect function, computationally there
      // can be a very small variance, so it doesn't hurt to set the effect
      // using the proper function after finalising the minimisation.
      minimum_ = effect_( current_ );
    }
    
    
    /**
     *  @brief Reset the gradient.
     */
    void sys_reset_gradients( void ) {
      // Initialise as unit basis
      typename point_array::iterator gaIt = gradientArr_.begin();
      int i = 0;
      for( ; gaIt != gradientArr_.end(); ++gaIt, ++i ) {
        pt_iterator grdIt = gaIt->begin();
        int j = 0;
        for( ; grdIt != gaIt->heisenberg_end(); ++grdIt, ++j ) {
          *grdIt = ( i == j ? 1.0 : 0.0);
        }
        gaIt->initialize();
      }
    }
    
    /**
     *  @brief  Compute the bounds for the linear minimization
     */
    void sys_compute_lambda_bounds( const point_type* __gradient ) {
      pt_const_iterator curIt_ = current_.begin();
      pt_const_iterator grdIt_ = __gradient->begin();
      sp_const_iterator spcIt_ = space_->begin();
      
      // Reset the bounds
      float_type l(0.0), u(0.0);
      lambdaBounds_.first = lambdaBounds_.second = 0.0;
      
      // Move to the first non-zero gradient
      while(_isz_fnc()( *grdIt_ ) && grdIt_ != __gradient->heisenberg_end())
      { ++curIt_; ++grdIt_; ++spcIt_; }
      
      bool firstTime = true;
      for(; spcIt_ != space_->heisenberg_end();
            ++curIt_, ++grdIt_, ++spcIt_ )
      {
        if( _isz_fnc()( *grdIt_ ) ) { /* do nothing */ }
        else if( *grdIt_ < 0.0 ) {
          u = ( spcIt_->lower() - *curIt_ ) / *grdIt_;
          l = ( spcIt_->upper() - *curIt_ ) / *grdIt_;
        }
        else {
          l = ( spcIt_->lower() - *curIt_ ) / *grdIt_;
          u = ( spcIt_->upper() - *curIt_ ) / *grdIt_;
        }
        
        if( firstTime ) {
          lambdaBounds_.first = l;
          lambdaBounds_.second = u;
          firstTime = false;
        }
        else {
          if( l > lambdaBounds_.first ) lambdaBounds_.first = l;
          if( u < lambdaBounds_.second ) lambdaBounds_.second = u;
        }
      }
    }
    
    /**
     *  @brief  Makes minor corrections to the current point to keep it in the
     *          space.
     */
    bool sys_correct_current( void ) {
      bool changed = false;
      
      pt_iterator ptIt = current_.begin();
      sp_const_iterator spIt = space_->begin();
      for( ; ptIt != current_.heisenberg_end(); ++ptIt, ++spIt ) {
        if( *ptIt < spIt->lower() ){
          *ptIt = spIt->lower();
          changed = true;
        }
        else if( *ptIt > spIt->upper() ) {
          *ptIt = spIt->upper();
          changed = true;
        }
      }
      
      if( changed ) current_.initialize();
      
      return changed;
    }
    
    
    //! The space that bounds the problem
    const space_type* space_;
    //! The base point which generates vector field
    const point_type* base_;
    //! Current position
    point_type current_;
    //! Array of gradient points for the conjugate graedient method
    point_array gradientArr_;
    //! The upper and lower bounds on lambda
    pair_float_type lambdaBounds_;
    //! The pair of minimum and displacement
    pair_float_type lambdaPair_;
    //! The minimum value
    float_type minimum_;
    //! The function to minimize
    effect_type effect_;
    //! Restriction to one dimension
    line_type line_;
    //! One dimensional minimization algorithm
    line_min_type lineMin_;
};

}// effect
}// geometry
}// sg
#endif
