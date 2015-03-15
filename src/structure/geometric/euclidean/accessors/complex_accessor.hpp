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
 *  @file     complex_accessor.hpp
 *  @brief    An inline header file for the \c complex_accessor class.
 *  @note     The \c complex_accessor class creates a common interface to
 *            complex numbers for hyperbolic structures.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-26
 */


#ifndef _SG_COMPLEX_ACCESSOR_H
#define _SG_COMPLEX_ACCESSOR_H 1

// Global includes
#include <complex>

#include <boost/operators.hpp>

// Local includes
#include "utility/math/is_equal_cx.hpp"
#include "utility/functors/stream_cast.hpp"


namespace sg
{
namespace structure
{
namespace geometric
{
namespace euclidean
{
namespace accessor
{
/**
 *  @class    complex_accessor
 *  @brief    An accessor class for complex numbers
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-26
 *  @param    _Float - The floating point type of the real and imaginary
 *            components of the complex number.
 */
template <class _Float = double>
    class complex_accessor
  : private boost::equality_comparable< complex_accessor<_Float>
             , boost::equality_comparable< complex_accessor<_Float>,
                                           std::complex<_Float>
                                         >
            >
{
  public:
    //! Object's type
    typedef complex_accessor<_Float> self_type;
    //! Type of a coordinate of zeta
    typedef std::complex<_Float> zeta_type;
    //! Type of the real component of zeta
    typedef _Float real_type;
    //! Type of the imaginary component of zeta
    typedef _Float imag_type;
    //! The common return type of real and imaginary parts
    typedef _Float return_type;
    
    
    /**
     *  @brief  Empty constructor
     */
    complex_accessor( ) : zeta_(0.0,0.0) { }
    
    
    /**
     *  @brief  Zeta constructor
     *  @param  zeta The zeta value to construct from
     */
    complex_accessor( const zeta_type& __zeta ) : zeta_(__zeta) { }
    
    
    /**
     *  @brief   Zeta assignment
     *  @param   zeta The zeta value to assign from
     *  @return  A reference to this.
     */
    self_type& operator=( const zeta_type& __zeta ) {
      zeta_ = __zeta;
      return *this;
    }
    
    
    /**
     *  @brief   Get a constant reference to the real part of zeta
     *  @return  A constant reference to the real part of zeta.
     */
    const real_type& real_ref( void ) const {
      return zeta_.real( );
    }
    
    
    /**
     *  @brief   Get a constant reference to the imaginary part of zeta
     *  @return  A constant reference to the imaginary part of zeta.
     */
    const imag_type& imag_ref( void ) const {
      return zeta_.imag( );
    }
    
    
    /**
     *  @brief   A common interface between real and imaginary parts
     *  @return  An object which is constructed from the real part of zeta
     *           and shares a common type with an imaginary return type.
     */
    const return_type& real_common( void ) const {
      return zeta_.real( );
    }
    
    /**
     *  @brief   A common interface between real and imaginary parts
     *  @return  An object which is constructed from the imag part of zeta
     *           and shares a common type with a real return type.
     */
    const return_type& imag_common( void ) const {
      return zeta_.imag( );
    }
    
    
    /**
     *  @brief  Compares two complex_accessor objects for equality.
     *  @param  rhs A \c complex_accessor
     */
    bool operator==( const self_type& __rhs ) const {
      return bool( *this == __rhs.zeta_ );
    }
    
    
    /**
     *  @brief  Compares two complex_accessor objects for equality.
     *  @param  rhs A \c complex_type
     */
    bool operator==( const zeta_type& __rhs ) const {
      typedef utility::math::is_equal_cx<_Float> _ie;
      return bool( _ie( )( zeta_, __rhs ) );
    }
    
    
    /**
     *  @brief  Defines the output stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::ostream& operator<<( std::ostream& __os, const self_type& __t ) {
      typedef utility::functors::stream_cast<zeta_type,std::ostream> _sc;
      
      _sc( )( __t.zeta_, __os );
      return __os;
    }
    
    /**
     *  @brief  Defines the input stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::istream& operator>>( std::istream& __is, self_type& __t ) {
      typedef utility::functors::stream_cast<std::istream,zeta_type> _sc;
      
      _sc( )( __is, __t.zeta_ );
      return __is;
    }
    
    
  protected:
    //! Zeta component
    zeta_type zeta_;
};

}// accessor
}// euclidean
}// geometric
}// structure
}// sg
#endif
