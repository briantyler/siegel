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
 *  @file     zeta_coordinate.hpp
 *  @brief    An inline header file for the \c zeta_coordinate class.
 *  @note     The \c zeta_coordinate class is something like a traits class,
 *            in that it provides a common interface for all zeta components.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-23
 */


#ifndef _SG_ZETA_COORDINATE_H
#define _SG_ZETA_COORDINATE_H 1

// Global includes
#include <boost/operators.hpp>


namespace sg
{
namespace structure
{
namespace geometric
{
namespace detail
{
/**
 *  @class    zeta_coordinate
 *  @brief    A common interface to the zeta coordinates of different types
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-23
 *  @note     There is no need to overload the \c == operator because
 *            \c zeta_coordinate inherits from _Accessor and this type
 *            type should have it's equality operator overloaded.
 */
template <class _Accessor>
    class zeta_coordinate : public _Accessor
{
  public:
    //! Object's type
    typedef zeta_coordinate<_Accessor> self_type;
    //! The accessor type
    typedef _Accessor accessor_type;
    //! The common return type; a common type that real and imag cast to.
    typedef typename accessor_type::zeta_type zeta_type;
    //! Type of the real component of zeta
    typedef typename accessor_type::real_type real_type;
    //! Type of the imaginary component of zeta
    typedef typename accessor_type::imag_type imag_type;
    //! Type of the common interface between real and imag
    typedef typename accessor_type::return_type return_type;
    
    
    /**
     *  @brief  Empty constructor
     */
    zeta_coordinate( ) : accessor_type( ) { }
    
    
    /**
     *  @brief  Zeta constructor
     *  @param  zeta The zeta value to construct from
     */
    zeta_coordinate( const zeta_type& __zeta ) : accessor_type(__zeta) { }
    
    
    /**
     *  @brief  Zeta assignment
     *  @param  zeta The zeta value to assign from
     *  @return A reference to this.
     */
    self_type& operator=( const zeta_type& __zeta ) {
      accessor_type::operator=( __zeta );
      return *this;
    }
    
    
    /**
     *  @brief   Get a constant reference to zeta
     *  @return  A constant reference to zeta.
     */
    const zeta_type& ref( void ) const {
      return accessor_type::zeta_;
    }
    
    /**
     *  @brief   Get a reference to zeta
     *  @return  A reference to zeta.
     */
    zeta_type& ref( void ) {
      return accessor_type::zeta_;
    }
    
    
    /**
     *  @brief   Get a reference to the common real part of zeta
     *  @return  A reference to the common real part of zeta.
     */
    const return_type& real_common( void ) const {
      return accessor_type::real_common( );
    }
    
    /**
     *  @brief   Get a reference to the common real part of zeta
     *  @return  A reference to the common real part of zeta.
     */
    return_type& real_common( void ) {
      // Cast away constness and use _Accessor
      return const_cast<return_type&>
          ( static_cast<const accessor_type&>( *this ).real_common( ) );
    }
    
    
    /**
     *  @brief   Get a reference to the common imaginary part of zeta
     *  @return  A reference to the common imaginary part of zeta.
     */
    const return_type& imag_common( void ) const {
      return accessor_type::imag_common( );
    }
    
    /**
     *  @brief   Get a reference to the common imaginary part of zeta
     *  @return  A reference to the common imaginary part of zeta.
     */
    return_type& imag_common( void ) {
      // Cast away constness and use _Accessor
      return const_cast<return_type&>
          ( static_cast<const accessor_type&>( *this ).imag_common( )  );
    }
    
    
    /**
     *  @brief   Get a constant reference to the real part of zeta
     *  @return  A constant reference to the real part of zeta.
     */
    const real_type& real_ref( void ) const {
      return accessor_type::real_ref( );
    }
    
    /**
     *  @brief   Get a reference to the real part of zeta
     *  @return  A reference to the real part of zeta.
     */
    real_type& real_ref( void ) {
      // Cast away constness and use _Accessor
      return const_cast<real_type&>
          ( static_cast<const accessor_type&>( *this ).real_ref( ) );
    }
    
    /**
     *  @brief   Get a constant reference to the imaginary part of zeta
     *  @return  A constant reference to the imaginary part of zeta.
     */
    const imag_type& imag_ref( void ) const {
      return accessor_type::imag_ref( );
    }
    
    
    /**
     *  @brief   Get a reference to the imaginary part of zeta
     *  @return  A reference to the imaginary part of zeta.
     */
    imag_type& imag_ref( void ) {
      // Cast away constness and use _Accessor
      return const_cast<imag_type&>
          ( static_cast<const accessor_type&>( *this ).imag_ref( )  );
    }
    
    
    /**
     *  @brief   Returns a constant reference to the underlying zeta
     *           component.
     *  @return  A constant reference to the zeta component.
     */
    const zeta_type& operator ( ) ( void ) const {
      return accessor_type::zeta_;
    }
    
    /**
     *  @brief   Returns a reference to the underlying zeta component.
     *  @return  A reference to the zeta component.
     */
    zeta_type& operator ( ) ( void ) {
      return accessor_type::zeta_;
    }
    
    /**
     *  @brief   Returns a constant reference to the underlying zeta
     *           component.
     *  @return  A constant reference to the zeta component.
     */
    const zeta_type& operator* ( void ) const {
      return accessor_type::zeta_;
    }
    
    /**
     *  @brief   Returns a reference to the underlying zeta component.
     *  @return  A reference to the zeta component.
     */
    zeta_type& operator* ( void ) {
      return accessor_type::zeta_;
    }
    
    
    /**
     *  @brief   Conversion to zeta_type
     *  @return  A copy of zeta
     */
    operator zeta_type ( ) const {
      return zeta_type( accessor_type::zeta_ );
    }
};

}// detail
}// geometric
}// structure
}// sg
#endif
