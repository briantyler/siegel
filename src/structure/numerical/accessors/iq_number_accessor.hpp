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
 *  @file     iq_number_accessor.hpp
 *  @brief    An inline header file for the \c iq_number_accessor class.
 *  @note     The \c iq_number_accessor class creates a common interface to
 *            complex numbers for heisenberg structures.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-26
 */


#ifndef _SG_IQ_NUMBER_ACCESSOR_H
#define _SG_IQ_NUMBER_ACCESSOR_H 1

// Global includes
#include <boost/operators.hpp>

// Local includes
#include "structure/numerical/iq_number.hpp"


namespace sg
{
namespace structure
{
namespace numerical
{
namespace accessor
{
/**
 *  @class    iq_number_accessor
 *  @brief    An accessor class for complex numbers
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-26
 *  @param    _Float The floating point type defaults to \c double
 *  @param    _Integer The integer type defaults to \c long
 *  @param    _Id The multiple singleton identifier of the iq_field
 */
template <class _Float = double, class _Integer = long, std::size_t _Id = 0>
    class iq_number_accessor
  : private boost::equality_comparable< iq_number_accessor<_Float,_Integer,_Id>,
                boost::equality_comparable<
                    iq_number_accessor<_Float,_Integer,_Id>,
                    iq_number<_Float,_Integer,_Id>
                >
            >
{
  public:
    //! Object type
    typedef iq_number_accessor<_Float,_Integer,_Id> self_type;
    //! Type of a coordinate of zeta
    typedef iq_number<_Float,_Integer,_Id> zeta_type;
    //! Type of the real component of zeta
    typedef typename zeta_type::integer_type real_type;
    //! Type of the imaginary component of zeta
    typedef real_type imag_type;
    //! The common return type of real and imaginary parts
    typedef real_type return_type;
    
    
    /**
     *  @brief  Empty constructor
     */
    iq_number_accessor( ) : zeta_( ) { }
    
    
    /**
     *  @brief  Zeta constructor
     *  @param  zeta The zeta value to construct from
     */
    iq_number_accessor( const zeta_type& __zeta ) : zeta_(__zeta) { }
    
    
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
     *  @brief  Compares two iq_number_accessor objects for equality.
     *  @param  lhs A \c iq_number_accessor
     *  @param  rhs A \c iq_number_accessor
     */
    bool operator==( const self_type& __rhs ) const {
      return bool( *this == __rhs.zeta_ );
    }
    
    
    /**
     *  @brief  Compares two iq_number_accessor objects for equality.
     *  @param  lhs A \c iq_number_accessor
     *  @param  rhs A \c iq_number
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
}// numerical
}// structure
}// sg
#endif
