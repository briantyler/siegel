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
 *  @file     real_coordinate.hpp
 *  @brief    An inline header file for the \c real_coordinate class.
 *  @note     The \c real_coordinate class is a wrapper class; it provides a
 *            common interface for all real coordinates
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-28
 */


#ifndef _SG_REAL_COORDINATE_H
#define _SG_REAL_COORDINATE_H 1

// Global includes
#include <boost/operators.hpp>
#include <boost/call_traits.hpp>


namespace sg
{
namespace structure
{
namespace geometric
{
namespace detail
{
/**
 *  @class    real_coordinate
 *  @brief    A common interface to the real coordinates of different types
 *  @note     The real coordinates are the "r" and "height" coordinates.
 *  @author   Brian Tyler
 *  @version  1.2
 *  @date     2008-02-28
 */

template <class _Accessor> class real_coordinate : public _Accessor
{
  public:
    //! The object type
    typedef real_coordinate<_Accessor> self_type;
    //! Type of the inherited accessor
    typedef _Accessor accessor_type;
    //! The value type of the accessor
    typedef typename accessor_type::value_type value_type;
    //! The common return type
    typedef typename accessor_type::return_type return_type;
    
    
  private:
    //! Value parameter type
    typedef typename boost::call_traits<value_type>::param_type param_type;
    
  public:
    
    
    /**
     *  @brief  Default constructor
     */
    real_coordinate( ) : accessor_type( ) { }
    
    
    /**
     *  @brief  Value constructor
     *  @param  value The value to construct from
     */
    real_coordinate( param_type __value ) : accessor_type( __value ){ }
    
    
    /**
     *  @brief  Assignment operator from value_type
     *  @param  value The value to assign from
     */
    self_type& operator=( param_type __value ) {
      accessor_type::operator=( __value );
      return *this;
    }
    
    
    /**
     *  @brief   A common interface between r and zeta
     *  @return  An object which is constructed from r and shares a common
     *           type with the zeta return type.
     */
    return_type& common( void ) {
      return const_cast<return_type&>( accessor_type::common() );
    }
    
    /**
     *  @brief   A common interface between r and zeta
     *  @return  An object which is constructed from r and shares a common
     *           type with the zeta return type.
     */
    const return_type& common( void ) const {
      return accessor_type::common();
    }
    
    
    /**
     *  @brief Get a constant reference to the raw real value.
     */
    const value_type& ref( void ) const { return accessor_type::value_; }
    
    /**
     *  @brief Get a reference to the raw real value.
     */
    value_type& ref( void ){ return accessor_type::value_; }
    
    
    /**
     *  @brief Get a constant reference to the raw real value.
     */
    const value_type& operator ( ) ( void ) const {
      return accessor_type::value_;
    }
    
    /**
     *  @brief Get a reference to the raw real value.
     */
    value_type& operator ( ) ( void ) { return accessor_type::value_; }
    
    
    /**
     *  @brief Provide conversion operator to value_type
     */
    operator value_type ( ) const {
      return value_type( accessor_type::value_ );
    }
};

}// detail
}// geometric
}// structure
}// sg
#endif
