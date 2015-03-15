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
 *  @file     lattice_bound.hpp
 *  @brief    An inline header file for the \c lattice_bound class.
 *  @note     \c lattice_bound keeps track of the current bound in a
 *            \c zeta_lattice
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-05-10
 */


#ifndef _SG_LATTICE_BOUND_H
#define _SG_LATTICE_BOUND_H 1

// Global includes
#include <cstddef>
#include <cassert>
#include <cmath>
#include <algorithm>

#include <boost/array.hpp>
#include <boost/call_traits.hpp>

// Local includes
#include "structure/numerical/iq_number.hpp"
#include "utility/math/square.hpp"
#include "utility/math/is_greater_equal.hpp"
#include "utility/functors/stream_cast.hpp"


namespace sg
{
namespace structure
{
namespace geometric
{
namespace lattice
{
namespace detail
{
/**
 *  @class    lattice_bound
 *  @brief    Keeps track of the current bound in a \c zeta_lattice
 *  @author   Brian Tyler
 *  @version  1.1
 *  @date     2008-05-10
 *  @note     \c v1.1 \c 2008-06-05 Modified the engine so that the threshold
 *            and bound does not depend on the norm even when the class number
 *            is non-trivial.
 *  @note     This object is specific to the SU(n,1) problem; whenever a
 *            distance is set (at index i say), all distances below it (ie all
 *            \f$ j < i \f$)  need to be set before getting the lower bounds.
 */
template <std::size_t N, class _Float = double, class _Integer = long,
          std::size_t _Id = 0>
  class lattice_bound
{
  public:
    //! Object type
    typedef lattice_bound<N,_Float,_Integer,_Id> self_type;
    //! Float type
    typedef _Float float_type;
    //! The integer type
    typedef _Integer integer_type;
    //! The imaginary quadratic number type
    typedef numerical::iq_number<float_type,integer_type,_Id> iq_number_type;
    //! The imaginary quadratic field type
    typedef typename iq_number_type::field_type field_type;
    //! The size type
    typedef std::size_t size_type;
    
  private:
    //! Type of the array of distance components and bounds
    typedef boost::array<float_type,N> distance_array_type;
    //! The floating point parameter type
    typedef typename boost::call_traits<float_type>::param_type flt_param_type;
    //! The integer parameter type
    typedef typename boost::call_traits<integer_type>::param_type int_param_type;
    
    
  public:
    /**
     *  @brief  Default constructor.
     */
    lattice_bound( )
  : dilation_(0), height_(0.0) { }
    
    
    /**
     *  @brief  Dilation and height constructor.
     *  @param  dilation The initial dilation factor
     *  @param  height The initial height component
     */
    lattice_bound( int_param_type __dilation, flt_param_type __height )
  : dilation_( __dilation ), height_( __height )
    { }
    
    
    /**
     *  @brief  Initializes the bound.
     *  @note   Should be called after chaning the dilation factor, the height
     *          component or starting a new iteeffectn.
     */
    void initialize( void ) {
      typedef utility::functors::stream_cast<integer_type,float_type> _sc;
      _sc()( dilation_, dilationFlt_ );
      
      sqrtDilationInv_ = 2.0/::sqrt(dilationFlt_);
      changed_ = true;
      lastLoc_ = distance_array_type::static_size;
      std::fill( distanceArr_.begin(), distanceArr_.end(), 0.0 );
    }
    
    
    /**
     *  @brief  Get a constant reference to the current dilation factor.
     *  @return A constant reference to the dilation factor.
     */
    const integer_type& dilation( void ) const { return dilation_; }
    
    /**
     *  @brief  Get a reference to the current dilation factor.
     *  @return A reference to the dilation factor.
     */
    integer_type& dilation( void ) { return dilation_; }
    
    
    /**
     *  @brief  Get a constant reference to the minimum height bound.
     *  @return A constant reference to the minimum height bound.
     */
    const float_type& height( void ) const { return height_; }
    
    /**
     *  @brief  Get a reference to the minimum height bound.
     *  @return A reference to the minimum height bound.
     */
    float_type& height( void ) { return height_; }
    
    
    /**
     *  @brief  Get the current bound at \c loc
     *  @param  loc The location of the bound to return
     *  @return The bound on the lattice region at \c loc
     *  
     *  \f$ \sqrt{\frac{2}{\delta} - h
     *            - \Sum_{i > loc}{|\zeta_i/\beta\delta - z_i|}^{N-1}
     *           } \f$ if the result is totally real, otherwise \f$ 0.0 \f$
     */
    const float_type& operator() ( size_type __loc ) const {
      // Only recompute the bound if the index has changed
      if( changed_ || lastLoc_ != __loc ) {
        if( __loc == N - 1 ) {
          bound_ = sqrtDilationInv_ - height_;
        }
        else {
          bound_ = sqrtDilationInv_ - height_ - distanceArr_[__loc + 1];
        }
        
        if( bound_ < 0.0 ) {
          bound_ = 0.0;
        }
        else {
          bound_ = ::sqrt( bound_ );
        }
        
        // The bound has been computed, so set changed to false
        changed_ = false;
        lastLoc_ = __loc;
      }
      
      return bound_;
    }
    
    
    /**
     *  @brief  Get the lattice component of the r bound.
     *  @return Returns the lattice component of the r bound for the current
     *          lattice point.
     *  
     *  \f$ \sqrt{\detla^2 -
     *            \left(\frac{\delta}{2}
     *              \left( \Sum{|\zeta_i/\beta\delta - z_i|}^{N-1} + height \right)
     *            \right)^2
     *           } \f$ if the result is totally real, otherwise \f$ 0.0 \f$
     */
    float_type r_bound ( void ) const {
      typedef utility::math::square<float_type> _sq;
      // RVO
      float_type bound(dilationFlt_);
      bound -= _sq()( 0.5*dilationFlt_*( distanceArr_.front() + height_ ) );
      
      // Don't want to take the square root of a negative number
      if( bound < 0.0 ) {
        bound = 0.0;
      }
      else {
        bound = ::sqrt( bound );
      }
      
      return bound;
    }
    
    
    /**
     *  @brief  Get the total bound for the current lattice point.
     *  @return A total bound for the current lattice point.
     *  
     *  \f$ \frac{2}{\delta} - h - \Sum{|\zeta_i/\beta\delta - z_i|}^{N-1} \f$
     */
    float_type total_bound ( ) const {
      // Allow negative bounds, this gives a criteria for an invalid lattice
      // point.
      return float_type( sqrtDilationInv_ - height_ - distanceArr_.front() );
    }
    
    
    /**
     *  @brief  Validate a point with respect to the original space
     *  @return True if the current point is not too far from the original
     *          space.
     */
    bool validate ( void ) const {
      typedef utility::math::is_greater_equal<float_type> _ige;
      return _ige()( total_bound(), 0.0 );
    }
    
    
    /**
     *  @brief  Set the distance at \c loc
     *  @param  loc The location of the distance component to update.
     */
    void set_distance( size_type __loc, flt_param_type __distance ) {
      if( __loc == N-1 ) {
        distanceArr_.back() = __distance*__distance;
      }
      else {
        distanceArr_[__loc] = distanceArr_[__loc+1] + __distance*__distance;
      }
      
      // The bound has changed
      changed_ = true;
    }
    
    
  private:
    // dilation_ and height_ could be constant references, but this seems to
    // be a bit unnecessary as it will have no real impact on speed and
    // complicates ownership issues.
    //! The dilation factor of the current cusp
    integer_type dilation_;
    //! The floating point representation of the dilation factor
    float_type dilationFlt_;
    //! The current minimum height of the space
    float_type height_;
    
    //! The bound on the lattice region
    mutable float_type bound_;
    //! Has the bound changed since it was lasted computed?
    mutable bool changed_;
    //! The last location to be computed
    mutable size_type lastLoc_;
    
    //! Square root of the inverse of the dilation factor
    float_type sqrtDilationInv_;
    //! The part attributable to the distance of the point from the space
    distance_array_type distanceArr_;
};

}// detail
}// lattice
}// geometric
}// structure
}// sg
#endif
