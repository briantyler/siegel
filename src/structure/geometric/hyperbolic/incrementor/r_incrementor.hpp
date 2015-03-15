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
 *  @file     r_incrementor.hpp
 *  @brief    An inline header file for the \c r_incrementor class.
 *  @note     The \c r_incrementor class creates an easy way to increment the
 *            \c r component of a hyperbolic point
 *  @author   Brian Tyler
 *  @version  1.1
 *  @date     2008-03-24
 */


#ifndef _SG_R_INCREMENTOR_H
#define _SG_R_INCREMENTOR_H 1

// Global includes
#include <functional>

#include <boost/operators.hpp>
#include <boost/call_traits.hpp>

// Local includes
#include "utility/math/is_equal.hpp"
#include "utility/functors/stream_cast.hpp"


namespace sg
{
namespace structure
{
namespace geometric
{
namespace hyperbolic
{
namespace incrementor
{
/**
 *  @class    r_incrementor
 *  @brief    Class for efficiently incrementing the \c r component of a
 *            hyperbolic point
 *  @author   Brian Tyler
 *  @version  1.1
 *  @date     2008-03-24
 *  @param    _Point A \c hyperbolic_point class
 */
template <class _Point> class r_incrementor
  : public std::binary_function<_Point&, typename _Point::return_type&, void>
  , private boost::equality_comparable< r_incrementor<_Point> >
{
  public:
    //! Object type
    typedef r_incrementor<_Point> self_type;
    //! Point type
    typedef _Point point_type;
    //! Return type of the point
    typedef typename point_type::return_type return_type;
    //! Parameter type of the return type
    typedef typename boost::call_traits<return_type>::param_type param_type;
    
    
    /**
     *  @brief  Default constructor
     */
    r_incrementor( ) : stride_(0.0) { }
    
    /**
     *  @brief  Stride constructor
     *  @param  stride The stride by which to increment the point.
     */
    r_incrementor( param_type __stride ) : stride_(__stride) { }
    
    
    /**
     *  @brief  Set the stride of the incrementor
     *  @param  stride The stride by which to increment the point.
     */
    void set_stride( param_type __stride ) { stride_ = __stride; }
    
    
    /**
     *  @brief  Get the stride of the incrementor
     *  @return The stride of the incrementor.
     */
    param_type stride( void ) const { return stride_; }
    
    
    /**
     *  @brief  The operator which increments the r component of the point
     *  @param  point The point to increment.
     *  @param  rvalue The \c r component of \c point
     *  @note   Although it might seem pointless to pass the r component of the
     *          point in, since it could just be called in the function, it
     *          makes sense in relation to the zeta incrementor, where the
     *          value could be either the real or imaginary part of any one
     *          of the zeta coordinates of the zeta array. This method turns
     *          out to be surprisingly efficient and is also quite simple.
     */
    void operator( )( point_type& __point, return_type& __rvalue ) const {
      __point.dependent_.imag_common() = (__rvalue += stride_);
    }
    
    
    /**
     *  @brief  Compares two r_incrementor objects for equality.
     *  @param  rhs An \c r_incrementor
     */
    bool operator==( const self_type& __rhs ) const {
      typedef utility::math::is_equal<return_type> _ie;
      return bool(  _ie( )( stride_, __rhs.stride_ ) );
    }
    
    
    /**
     *  @brief  Defines the output stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::ostream& operator<<( std::ostream& __os, const self_type& __t ) {
      typedef utility::functors::stream_cast<return_type,std::ostream> _sc;
      
      _sc( )( __t.stride_, __os );
      return __os;
    }
    
    /**
     *  @brief  Defines the input stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::istream& operator>>( std::istream& __is, self_type& __t ) {
      typedef utility::functors::stream_cast<std::istream,return_type> _sc;
      
      _sc( )( __is, __t.stride_ );
      return __is;
    }
    
    
  private:
    //! The stride to increment by.
    return_type stride_;
};

}// incrementor
}// hyperbolic
}// geometric
}// structure
}// sg
#endif
