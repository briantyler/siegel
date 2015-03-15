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
 *  @file     complex_region.hpp
 *  @brief    An inline header file for the \c complex_region class.
 *  @note     Represents a region in complex space
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-23
 */


#ifndef _SG_COMPLEX_REGION_H
#define _SG_COMPLEX_REGION_H 1

// Global includes
#include <complex>

#include <boost/array.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/operators.hpp>
#include <boost/call_traits.hpp>

// Local includes
#include "utility/functors/stream_cast_tuple.hpp"
#include "structure/geometric/euclidean/real_interval.hpp"


namespace sg
{
namespace structure
{
namespace geometric
{
namespace euclidean
{
/**
 *  @class    complex_region
 *  @brief    A class representing an interval on the real line.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-23
 *  @param    _Float A floating point type; defaults to \c double
 */
template <class _Float = double> class complex_region
  : private boost::equality_comparable< complex_region<_Float> >
{
  public:
    //! Object type
    typedef complex_region<_Float> self_type;
    //! Value type
    typedef _Float float_type;
    //! Interval type
    typedef real_interval<_Float> interval_type;
    //! Real type
    typedef _Float real_type;
    //! Complex type
    typedef std::complex<_Float> complex_type;
   
  private:
    //! The type of the pair of intervals
    typedef boost::tuple<interval_type,interval_type> tuple_type;
    //! Float parameter type
    typedef typename boost::call_traits<float_type>::param_type param_type;
    //! The array of corners
    typedef boost::array<complex_type,4> corner_array_type;
    
  public:
    //! Constant forward corner iterator
    typedef typename corner_array_type::const_iterator const_iterator;
    //! Constant reverse corner iterator
    typedef typename corner_array_type::const_reverse_iterator
        const_reverse_iterator;
  
    /**
     *  @brief  Default constructor
     */
    complex_region( ) : region_( ) { }
    
    /**
     *  @brief    Explicit interval constructor
     *  @param    re the real interval
     *  @param    im the imaginary interval
     */
    // Note __real__ and __imag__ are C keywords, and the compiler doesn't
    // like it if it sees __real or __imag.
    complex_region( const interval_type& __re, const interval_type& __im )
    : region_( __re, __im ) { }
    
    
    /**
     *  @brief  Get a constant reference to the real interval
     */
    const interval_type& real( void ) const { return boost::get<0>( region_ ); }
    
    /**
     *  @brief  Get a reference to the real interval
     */
    interval_type& real( void ) {
      return const_cast<interval_type&>
          ( static_cast<const self_type&>( *this ).real() );
    }
    
    
    /**
     *  @brief  Get a constant reference to the imaginary interval
     */
    const interval_type& imag( void ) const { return boost::get<1>( region_ ); }
    
    /**
     *  @brief  Get a reference to the imaginary interval
     */
    interval_type& imag( void ) {
      return const_cast<interval_type&>
          ( static_cast<const self_type&>( *this ).imag() );
    }
    
    
    /**
     *  @brief  Puts the bound into canonical form
     *  @note   Ensures that \c lower < \c upper in both bounds and builds
     *          the corner points.
     */
    void initialize( void ) {
      initialize_intervals();
      initialize_corners();
    }
    
    /**
     *  @brief  Initializes the sub intervals
     */
    void initialize_intervals( void ) {
      real().initialize( );
      imag().initialize( );
    }
    
    /**
     *  @brief  Initialize (compute) the corner points
     */
    void initialize_corners( void ) {
      corners_[0] = complex_type( real().lower(), imag().lower() );
      corners_[1] = complex_type( real().upper(), imag().lower() );
      corners_[2] = complex_type( real().lower(), imag().upper() );
      corners_[3] = complex_type( real().upper(), imag().upper() );
    }
    
    
    /**
     *  @brief  Set the region by the bottom left and top right corner points.
     */
    void set_bl_tr( const complex_type& __bl, const complex_type& __tr ) {
      real().lower() = __bl.real();
      imag().lower() = __bl.imag();
      
      real().upper() = __tr.real();
      imag().upper() = __tr.imag();
      
      initialize();
    }
    
    
    /**
     *  @brief  Determines if the region contains a point
     *  @param  value The point in complex space
     *  @return True if \c value is contained in the region
     */
    bool contains( const complex_type& __value ) const {
      return bool(    real().contains( __value.real() )
                   && imag().contains( __value.imag() )
                 );
    }
    
    
    /**
     *  @brief  Find the closest point in the region to a given value
     *  @param  value The point in complex space
     *  @return The closest point to \c value in the region
     */
    complex_type closest( const complex_type& __value ) const {
      return complex_type( real().closest( __value.real() ),
                           imag().closest( __value.imag() )
                         );
    }
    
    
    /**
     *  @brief  Compute the distance between a complex point and the nearest
     *          point in the region
     *  @param  value The point in complex space
     *  @return The distance between value and the closest point in the region
     */
    float_type distance( const complex_type& __value ) const {
      return float_type( std::abs( __value - closest(__value) ) );
    }
    
    /**
     *  @brief  Compute the squared distance between a complex point and the
     *          nearest point in the region.
     *  @param  value The point in complex space
     *  @return The squared distance between value and the closest point in
     *          the region
     */
    float_type distance2( const complex_type& __value ) const {
      return float_type( std::norm( __value - closest(__value) ) );
    }
    
    
    /**
     *  @brief  Extend the region by \c value in both dimesnions
     *  @param  value The amount to extend by
     *  @return \c *this after extension
     *  @note   If \c value is negative then the region will be contracted and
     *          the interval bounds may become inverted; invalidating the region.
     */
    self_type& extend( param_type __value ) {
      real().extend( __value );
      imag().extend( __value );
      initialize_corners();
      return *this;
    };
    
    
    /**
     *  @brief Get a constant forward iterator for iterating over the corners
     */
    const_iterator begin( void ) const {
      return const_iterator( corners_.begin() );
    };
    
    /**
     *  @brief Get a one past the end corner iterator
     */
    const_iterator end( void ) const {
      return const_iterator( corners_.end() );
    };
    
    /**
     *  @brief Get a constant reverse iterator for iterating over the corners
     */
    const_reverse_iterator rbegin( void ) const {
      return const_reverse_iterator( corners_.rbegin() );
    };
    
    /**
     *  @brief Get a one before the beginning corner iterator
     */
    const_reverse_iterator rend( void ) const {
      return const_reverse_iterator( corners_.rend() );
    };
    
    
    /**
     *  @brief   Gets the complex number at the bottom lefthand corner of
     *           the region.
     *  @return  The complex number at the bottom left of the region.
     */
    const complex_type& bl( void ) const {
      return corners_[0];
    }
    
    /**
     *  @brief   Gets the complex number at the bottom righthand corner of
     *           the region.
     *  @return  The complex number at the bottom right of the region.
     */
    const complex_type& br( void ) const {
      return corners_[1];
    }
    
    /**
     *  @brief   Gets the complex number at the top lefthand corner of
     *           the region.
     *  @return  The complex number at the top left of the region.
     */
    const complex_type& tl( void ) const {
      return corners_[2];
    }
    
    /**
     *  @brief   Gets the complex number at the bottom righthand corner of
     *           the region.
     *  @return  The complex number at the top right of the region.
     */
    const complex_type& tr( void ) const {
      return corners_[3];
    }
    
    
    /**
     *  @brief  Convert a rectangle to a complex region
     *  @param  rectangle The rectangle to convert from
     */
    template <class _Rectangle>
        void from_rectangle( const _Rectangle& __rectangle )
    {
      real().set( __rectangle.bl().real(), __rectangle.br().real() );
      imag().set( __rectangle.bl().imag(), __rectangle.tl().imag() );
      
      initialize_corners();
    }
    
    
    /**
     *  @brief  Compares two complex_region objects for equality.
     *  @param  rhs A \c complex_region
     */
    bool operator==( const self_type& __rhs ) const {
      return bool( real() == __rhs.real() && imag() == __rhs.imag() );
    }
    
    
    /**
     *  @brief  Defines the output stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::ostream& operator<<( std::ostream& __os, const self_type& __t ) {
      typedef utility::functors::stream_cast_tuple<tuple_type,std::ostream> sc;
      
      sc( )( __t.region_, __os );
      return __os;
    }
    
    /**
     *  @brief  Defines the input stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::istream& operator>>( std::istream& __is, self_type& __t ) {
      typedef utility::functors::stream_cast_tuple<std::istream,tuple_type> sc;
      
      sc( )( __is, __t.region_ );
      return __is;
    }
  
  
  private:
    //! The region
    tuple_type region_;
    //! The corner points
    corner_array_type corners_;
};

}// euclidean
}// geometric
}// structure
}// sg
#endif
