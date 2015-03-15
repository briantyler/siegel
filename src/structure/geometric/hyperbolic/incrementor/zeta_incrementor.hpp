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
 *  @file     zeta_incrementor.hpp
 *  @brief    An inline header file for the zeta incrementor classes.
 *  @note     The zeta incrementor classes offer an efficient interface for
 *            iterating the zeta component of hyperbolic points. By directly
 *            accessing the data members of the hyperbolic point, and using
 *            knowledge of the inner hermitian form, points can be incremented
 *            with very little overhead.
 *  @author   Brian Tyler
 *  @version  1.2
 *  @date     2008-03-24
 */


#ifndef _SG_ZETA_INCREMENTOR_H
#define _SG_ZETA_INCREMENTOR_H 1

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
 *  @class    zeta_incrementor
 *  @brief    Class for efficiently incrementing the zeta component of a
 *            hyperbolic point.
 *  @author   Brian Tyler
 *  @version  1.2
 *  @date     2008-03-24
 *  @param    _Point The \c hyperbolic_point type to increment
 */
template <class _Point> class zeta_incrementor
  : public std::binary_function<_Point&, typename _Point::return_type&, void>
  , private boost::equality_comparable< zeta_incrementor<_Point> >
{
  public:
    //! Object type
    typedef zeta_incrementor<_Point> self_type;
    //! Point type
    typedef _Point point_type;
    //! Real type
    typedef typename point_type::return_type return_type;
    //! Parameter type of the return type
    typedef typename boost::call_traits<return_type>::param_type param_type;
    
    /**
     *  @brief Default constructor
     */
    zeta_incrementor( )
    : stride_(0.0), qf_(0.0) { }
    
    
    /**
     *  @brief  Location and Stride constructor
     *  @param  stride The amount by which to increment the component
     */
    zeta_incrementor( param_type __stride )
  : stride_(), qf_() { set_stride(__stride); }
    
    
    /**
     *  @brief  Sets the stride of the incrementor.
     *  @param  stride The amount to increment by on each step.
     */
    void set_stride( param_type __stride ) {
      stride_ = __stride;
      this->compute_qf();
    }
    
    
    /**
     *  @brief  Get the stride of the incrementor
     *  @return The stride of the incrementor.
     */
    param_type stride( void ) const { return stride_; }
    
    
    /**
     *  @brief  Get the amount to increase the dependent variable by, given the
     *          initial value of a particular coordinate.
     *  @param  zvalue The real or imaginary zeta value being incremented.
     *  @return The change in the dependent coordinate.
     */
    return_type dependent( param_type __zvalue ) const {
      return return_type( -( __zvalue * stride_ ) + qf_ );
    }
    
    
    
    /**
     *  @brief  The operator which increments the r component of the point
     *  @param  point The point to increment.
     *  @param  zvalue The real or imaginary zeta value to increment.
     */
    void operator( )( point_type& __point, return_type& __zvalue ) const {
      __point.dependent_.real_ref() += dependent( __zvalue );
      __zvalue += stride_;
    }
      
      
    /**
     *  @brief  Compares two r_incrementor objects for equality.
     *  @param  rhs A \c zeta_incrementor
     */
    bool operator==( const self_type& __rhs ) const {
      typedef utility::math::is_equal<return_type> _ie;
      return bool( _ie( )( stride_, __rhs.stride_ ) );
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
      
      return_type stride;
      _sc( )( __is, stride );
      __t.set_stride( stride );
      
      return __is;
    }
    
    
  private:
    /**
     *  @brief  Computes the quadratic form contribution of \c stride_
     */
    void compute_qf( void ) {
      qf_ = (-0.5)*( stride_ * stride_ );
    }
    
    
    //! The stride to increment by
    return_type stride_;
    //! The quadratic form contribution of \c stride_
    return_type qf_;
};

}// incrementor
}// hyperbolic
}// geometric
}// structure
}// sg
#endif
