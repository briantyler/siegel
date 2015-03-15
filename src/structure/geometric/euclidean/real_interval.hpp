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
 *  @file     real_interval.hpp
 *  @brief    An inline header file for the \c real_interval class.
 *  @note     Represents an interval on the real line
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-23
 */


#ifndef _SG_REAL_INTERVAL_H
#define _SG_REAL_INTERVAL_H 1

// Global includes
#include <boost/tuple/tuple.hpp>
#include <boost/operators.hpp>
#include <boost/call_traits.hpp>

// Local includes
#include "utility/math/abs.hpp"
#include "utility/math/comparison.hpp"
#include "utility/functors/stream_cast_tuple.hpp"

namespace sg
{
namespace structure
{
namespace geometric
{
namespace euclidean
{
/**
 *  @class    real_interval
 *  @brief    A class representing an interval on the real line.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-23
 *  @param    _Float A floating point type; defaults to \c double
 */
template <class _Float = double> class real_interval
  : private boost::equality_comparable< real_interval<_Float>
    , boost::field_operators< real_interval<_Float>, _Float > >
    //, boost::subtractable< real_interval<_Float>, _Float
    //, boost::multipliable< real_interval<_Float>, _Float
    //, boost::dividable< real_interval<_Float>, _Float
    //> > > > >
{
  public:
    //! Object type
    typedef real_interval<_Float> self_type;
    //! The float type
    typedef _Float float_type;
    
  private:
    //! The tuple type
    typedef boost::tuple<float_type,float_type> tuple_type;
    //! Parameter type
    typedef typename boost::call_traits<float_type>::param_type param_type;
    
    
  public:
    /**
     *  @brief  Default constructor
     */
    real_interval( ) : interval_( 0, 0 ) { }
    
    /**
     *  @brief    Explicit interval constructor
     *  @param    lower the lower bound
     *  @param    upper the upper bound
     *  @warning  It is the responsibility of the programmer to ensure that
     *            \c lower < \c upper if in doubt call \c initialize()
     */
    real_interval( param_type __lower, param_type __upper )
    : interval_( __lower, __upper ) { }
    
    
    /**
     *  @brief  Get a constant reference to the lower bound
     */
    const float_type& lower( void ) const { return boost::get<0>( interval_ ); }
    
    /**
     *  @brief  Get a reference to the lower bound
     */
    float_type& lower( void ) {
      return const_cast<float_type&>(static_cast<const self_type&>(*this).lower());
    }
    
    /**
     *  @brief  Get a constant reference to the upper bound
     */
    const float_type& upper( void ) const { return boost::get<1>( interval_ ); }
    
    /**
     *  @brief  Get a reference to the lower bound
     */
    float_type& upper( void ) {
      return const_cast<float_type&>(static_cast<const self_type&>(*this).upper());
    }
    
    
    /**
     *  @brief    Sets the upper and lower bound
     *  @param    lower the lower bound
     *  @param    upper the upper bound
     *  @note     Ensures that \c lower < \c upper
     */
    void set( param_type __lower, param_type __upper ) {
      if( __lower > __upper ){
        lower() = __upper;
        upper() = __lower;
      } else {
        lower() = __lower;
        upper() = __upper;
      }
    }
    
    
    /**
     *  @brief  Puts the interval into canonical form
     *  @note   Ensures that \c lower < \c upper
     */
    void initialize( void ) {
      if( lower() > upper() ){
        float_type temp = lower();
        lower() = upper();
        upper() = temp;
      }
    }
    
    
    /**
     *  @brief  Get the length of the interval
     *  @return \c upper - \c lower
     */
    float_type length( void ) const {
      return float_type( upper() - lower() );
    };
    
    /**
     *  @brief  Get the midpoint of the interval
     *  @return (upper + lower)/2
     *  @todo   WRITE test script (although it is pretty obvious that it's ok)
     */
    float_type midpoint( void ) const {
      return float_type( (upper() + lower()) * 0.5 );
    };
    
    
    /**
     *  @brief  Determines if the interval [-,-] contains a point
     *  @param  value The point to check
     *  @return True if \c value is contained in the interval
     *  @note   Closed interval
     */
    bool contains( param_type __value ) const {
      typedef utility::math::is_less_equal<float_type> _le;
      
      return bool( _le()(lower(),__value) && _le()(__value,upper()) );
    }
    
    /**
     *  @brief  Determines if the interval [-,-) contains a point
     *  @param  value The point to check
     *  @return True if \c value is contained in the interval
     *  @note   Lower closed, upper open interval
     */
    bool contains_co( param_type __value ) const {
      typedef utility::math::is_less<float_type> _l;
      typedef utility::math::is_less_equal<float_type> _le;
      
      return bool( _le()(lower(),__value) && _l()(__value,upper()) );
    }
    
    /**
     *  @brief  Determines if the interval (-,-] contains a point
     *  @param  value The point to check
     *  @return True if \c value is contained in the interval
     *  @note   Lower open, upper closed interval
     */
    bool contains_oc( param_type __value ) const {
      typedef utility::math::is_less<float_type> _l;
      typedef utility::math::is_less_equal<float_type> _le;
      
      return bool( _l()(lower(),__value) && _le()(__value,upper()) );
    }
    
    /**
     *  @brief  Determines if the interval (-,-) contains a point
     *  @param  value The point to check
     *  @return True if \c value is contained in the interval
     *  @note   Open interval
     */
    bool contains_oo( param_type __value ) const {
      typedef utility::math::is_less<float_type> _l;
      
      return bool( _l()(lower(),__value) && _l()(__value,upper()) );
    }
    
    
    /**
     *  @brief  Find the closest point in the interval to a given value
     *  @param  value The point on the real line
     *  @return The closest point to \c value in the interval
     *  @note   Before this used the is_less functor, there seems no point in
     *          this, but I may have assumed that it does elsewhere (like in the
     *          minimisation class); something to be aware of.
     */
    float_type closest( param_type __value ) const {
      return float_type(   __value < lower()
                         ? lower()
                         : ( upper() < __value ? upper() : __value )
                       );
    }
    
    
    /**
     *  @brief  Find the distance between a point and the closest point in the
     *          interval
     *  @param  value The point on the real line
     *  @return The distance to \c value from the interval
     */
    float_type distance( param_type __value ) const {
      typedef utility::math::abs<float_type> _abs;
      
      return float_type( _abs()( __value - closest(__value) ) );
    }
    
    
    /**
     *  @brief  Extend the interval by \c value at both ends
     *  @param  value The amount to extend by
     *  @return \c *this after extension
     *  @note   If \c value is negative then the interval will be contracted
     *          and the bounds could be reversed.
     */
    self_type& extend( param_type __value ) {
      lower() -= __value;
      upper() += __value;
      return *this;
    };
    
    
    /**
     *  @brief  Increments the interval by a real number.
     *  @param  rhs The real number to increment by.
     *  @return A constant copy of \c *this incremented by \c rhs
     */
    self_type& operator +=( param_type __rhs ) {
      lower() += __rhs;
      upper() += __rhs;
      return *this;
    }
    
    
    /**
     *  @brief  Decrements the interval by a real number.
     *  @param  rhs The real number to decrement by.
     *  @return A constant copy of \c *this decremented by \c rhs
     */
    self_type& operator -=( param_type __rhs ) {
      lower() -= __rhs;
      upper() -= __rhs;
      return *this;
    }
    
    
    /**
     *  @brief  Multiplies the interval by a real number.
     *  @param  rhs The real number to multiply by.
     *  @return A constant copy of \c *this multiplied by \c rhs
     *  @note   This DOESN'T multiply the length of the interval by \c rhs it
     *          multiplies each end point, so if \c rhs < 0 the interval
     *          becomes invalidated and must be re-initialized.
     */
    self_type& operator *=( param_type __rhs ) {
      lower() *= __rhs;
      upper() *= __rhs;
      return *this;
    }
    
    
    /**
     *  @brief  Divides the interval by a real number.
     *  @param  rhs The real number to divide by.
     *  @return A constant copy of \c *this divided by \c rhs
     *  @note   See *= (this function has odd behaviour)
     */
    self_type& operator /=( param_type __rhs ) {
      lower() /= __rhs;
      upper() /= __rhs;
      return *this;
    }
    
    
    /**
     *  @brief  Compares two real_interval objects for equality.
     *  @param  rhs A \c real_interval
     */
    bool operator==( const self_type& __rhs ) const {
      utility::math::is_equal<float_type> _ie;
      return bool(    _ie( lower(), __rhs.lower() )
                   && _ie( upper(), __rhs.upper() )
                 );
    }
    
    
    /**
     *  @brief  Defines the output stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::ostream& operator<<( std::ostream& __os, const self_type& __t ) {
      using utility::functors::stream_cast_tuple;
      typedef stream_cast_tuple<tuple_type,std::ostream> _sc;
      
      _sc( )( __t.interval_, __os );
      return __os;
    }
    
    /**
     *  @brief  Defines the input stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::istream& operator>>( std::istream& __is, self_type& __t ) {
      using utility::functors::stream_cast_tuple;
      typedef stream_cast_tuple<std::istream,tuple_type> _sc;
      
      _sc( )( __is, __t.interval_ );
      return __is;
    }
  
  
  private:
    //! The interval
    tuple_type interval_;
};

}// euclidean
}// geometric
}// structure
}// sg

#endif
