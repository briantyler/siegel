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
 *  @file     region_accessor.hpp
 *  @brief    An inline header file for the \c region_accessor class.
 *  @note     The \c region_accessor class creates a common interface to
 *            complex regions for heisenberg structures.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-12
 */


#ifndef _SG_REGION_ACCESSOR_H
#define _SG_REGION_ACCESSOR_H 1

// Global includes
#include <boost/operators.hpp>

// Local includes
#include "structure/geometric/euclidean/complex_region.hpp"


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
 *  @class    region_accessor
 *  @brief    An accessor class for real numbers
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-28
 *  @param    _Float A floating point type; defaults to \c double
 */
template <class _Float = double >
    class region_accessor
  : private boost::equality_comparable< region_accessor<_Float>
              , boost::equality_comparable< region_accessor<_Float>,
                                            complex_region<_Float>
                                          >
            >
{
  public:
    //! Object's type
    typedef region_accessor<_Float> self_type;
    //! Type of a coordinate of zeta
    typedef complex_region<_Float> zeta_type;
    //! Type of the real component of zeta
    typedef typename zeta_type::interval_type real_type;
    //! Type of the imaginary component of zeta
    typedef real_type imag_type;
    //! The common return type of real and imaginary parts
    typedef real_type return_type;
    
    
    /**
     *  @brief  Empty constructor
     */
    region_accessor( ) : zeta_( ) { }
    
    
    /**
     *  @brief  Zeta constructor
     *  @param  zeta The zeta value to construct from
     */
    region_accessor( const zeta_type& __zeta ) : zeta_( __zeta ) { }
    
    
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
     *  @note    This needs to be a constant object.
     */
    const return_type& real_common( void ) const {
      return zeta_.real( );
    }
    
    /**
     *  @brief   A common interface between real and imaginary parts
     *  @return  An object which is constructed from the imag part of zeta
     *           and shares a common type with a real return type.
     *  @note    This needs to be a constant object.
     */
    const return_type& imag_common( void ) const {
      return zeta_.imag( );
    }
    
    
    /**
     *  @brief  Compares two region_accessor objects for equality.
     *  @param  rhs A \c region_accessor
     */
    bool operator==( const self_type& __rhs ) const {
      return bool( *this == __rhs.zeta_ );
    }
    
    
    /**
     *  @brief  Compares two region_accessor objects for equality.
     *  @param  rhs A \c complex_type
     */
    bool operator==( const zeta_type& __rhs ) const {
      return bool( zeta_ == __rhs );
    }
    
    
    /**
     *  @brief  Defines the output stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::ostream& operator<<( std::ostream& __os, const self_type& __t ) {
      __os << __t.zeta_;
      return __os;
    }
    
    /**
     *  @brief  Defines the input stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::istream& operator>>( std::istream& __is, self_type& __t ) {
      __is >> __t.zeta_;
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
