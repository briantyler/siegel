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
 *  @file     interval_data.hpp
 *  @brief    An inline header file for the \c interval_data class.
 *  @note     Represents the division of an interval.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-17
 */


#ifndef _SG_INTERVAL_DATA_H
#define _SG_INTERVAL_DATA_H 1

// Global includes
#include <boost/operators.hpp>
#include <boost/call_traits.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_same.hpp>

// Local includes
#include "utility/math/is_equal.hpp"
#include "utility/math/ceil.hpp"
#include "utility/functors/stream_cast.hpp"
#include "utility/io/string_parser.hpp"


namespace sg
{
namespace structure
{
namespace geometric
{
namespace euclidean
{
namespace detail
{
/**
 *  @class    interval_data
 *  @brief    A class representing the division interval of the real line.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-17
 *  @param    _Float A floating point type; defaults to \c double
 *  @param    _Integer An integer type; defaults to \c long
 */
template <class _Float = double, class _Integer = long> class interval_data
  : private boost::equality_comparable< interval_data<_Float,_Integer> >
{
  public:
    //! Object type
    typedef interval_data<_Float,_Integer> self_type;
    //! Real type
    typedef _Float float_type;
    //! Integer type
    typedef _Integer integer_type;
   
  private:
    //! Integer parameter type
    typedef typename boost::call_traits<integer_type>::param_type
        integer_param_type;
    //! Float parameter type
    typedef typename boost::call_traits<float_type>::param_type
        float_param_type;
    
  public:
    /**
     *  @brief  Sets the stride and resolution.
     *  @param  resolution The resolution.
     *  @param  length The length of the interval.
     */
    void set_data
        ( integer_param_type __resolution, float_param_type __length )
    {
      typedef utility::functors::stream_cast<integer_type,float_type> _sc;
      
      // During debugging ensure that resolution and length are positive.
      assert( __resolution > 0 && __length > 0.0 );
      
      resolution_ = __resolution;
      
      // The inefficiency introduced here by the stream_cast is so
      // insignificant that this comment is a waste of bytes, there is no
      // point being clever and doing some sort of is_convertible meta-magic.
      float_type resolution;
      _sc()( __resolution, resolution );
      stride_ = __length / resolution;
    }
    
    
    /**
     *  @brief  A constant reference to the resolution
     */
    const integer_type& resolution( void ) const { return resolution_; }
    
    /**
     *  @brief  A constant reference to the stride.
     *  @note   The stride is the size of sub-interval such that such that
     *          stride * resolution = length of enire interval.
     */
    const float_type& stride( void ) const { return stride_; }
    
    
    /**
     *  @brief  Get the subinterval of \c interval at \c loc.
     *  @param  Interval The interval type
     *  @param  interval The interval to get the sub interval of.
     *  @param  loc The index of the sub interval.
     *  @return The subinterval of \c interval at \c loc.
     */
    template <class _Interval>
        _Interval subinterval_at
        ( const _Interval& __interval, integer_param_type __loc ) const
    {
      typedef utility::functors::stream_cast<integer_type,float_type> _sc;
      
      // Ensure that the real types match
      BOOST_MPL_ASSERT( (
            boost::is_same<typename _Interval::float_type, float_type>
      ) );
      
      // Again, there is no need to use meta-techniques, it's fast enough
      float_type loc;
      _sc()( __loc, loc );
      
      float_type lower( __interval.lower() + (loc * stride_) );
      return _Interval( lower, lower + stride_ );
    }
    
     
    /**
     *  @brief  Compares two interval_data objects for equality.
     *  @param  rhs A \c interval_data
     */
    bool operator==( const self_type& __rhs ) const {
      typedef utility::math::is_equal<float_type> _ie;
      
      return bool(    resolution_ == __rhs.resolution_
          &&
          _ie( )( stride_, __rhs.stride_ )
                 );
    }
    
    
    /**
     *  @brief  Defines the output stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::ostream& operator<<( std::ostream& __os, const self_type& __t ) {
      typedef utility::functors::stream_cast<integer_type,std::ostream> _sci;
      typedef utility::functors::stream_cast<float_type,std::ostream> _scf;
      
      __os << '[';
      _sci( )( __t.resolution_, __os );
      __os << ',';
      _scf( )( __t.stride_, __os );
      __os << ']';
      
      return __os;
    }
    
    /**
     *  @brief  Defines the input stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::istream& operator>>( std::istream& __is, self_type& __t ) {
      typedef utility::functors::stream_cast<std::string,integer_type> _sci;
      typedef utility::functors::stream_cast<std::string,float_type> _scf;
      typedef utility::io::string_parser _sp;
      
      std::string s;
      __is >> s;
      
      _sp::string_vector_type vs = _sp()(s);
      _sci( )( vs.at(0), __t.resolution_ );
      _scf( )( vs.at(1), __t.stride_ );
      
      return __is;
    }
    
    
    /**
     *  @struct   compute_data
     *  @brief    A binary function object that computes interval data.
     *  @author   Brian Tyler
     *  @version  1.0
     *  @date     2008-03-11
     *  @note     Test this with the Heisenberg slice.
     */
    template <class _Interval> struct compute_data
    {
      typedef compute_data<_Interval> self_type;
      typedef _Interval interval_type;
      typedef interval_data<float_type,integer_type> data_type;
      
      // Ensure that the real types match
      BOOST_MPL_ASSERT( (
            boost::is_same<typename interval_type::float_type, float_type>
      ) );
      
      /**
       *  @brief  Constructor
       *  @note   Initialises the resolution to one
       */
      compute_data( float_param_type __resolution, float_param_type __length )
        : resolution_(1), correction_(__resolution/__length) { }
      
      /**
       *  @brief  Sets the interval data based on interval and starting conditions
       *  @param  interval An interval type
       *  @param  data An interval data type
       */
      void operator( ) ( const interval_type& __interval, data_type& __data ) {
        typedef utility::math::ceil<integer_type,float_type> _ceil;
         
        /* If
         *    resolution = f.p. resolution of a sub-interval a hypercube
         *    length = length of one edge of that hypercube
         * Then we want to get as close as possible to matching the ratio len / res
         * in order to evenly divide the space as best as possible:
         *    [length of this interval] / [resolution of this interval]
         *   ~
         *    len / res
         * Therefore
         *    [resolution of this interval] = (res/len) * [length of this interval]
         * 
         * Choose the ceiling of this value to avoid the possiblity of having a
         * smaller final resolution than the input resolution, this makes it easier
         * to increase the resolution. If this becomes problematic something more
         * sophisticated may be necessary.
         */
        
        float_type length( __interval.length() );
        __data.set_data( _ceil( )( correction_ * length ), length );
        resolution_ *= __data.resolution();
      }
      
      
      /**
       *  @brief  Get the cumulative resolution
       *  @return The cumulative resolution
       */
      operator integer_type ( ) { return resolution_; }
      
      private:
        //! The cumulative resolution
        integer_type resolution_;
        //! The correction factor
        float_type correction_;
    };
    
    
  private:
    //! The resolution of the interval
    integer_type resolution_;
    //! The stride of the interval
    float_type stride_;
    
};

}// detail
}// euclidean
}// geometric
}// structure
}// sg
#endif
