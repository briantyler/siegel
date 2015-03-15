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
 *  @file     r_incrementor_accessor.hpp
 *  @brief    An inline header file for the \c r_incrementor_accessor class.
 *  @note     The \c r_incrementor_accessor class creates a common interface to
 *            r_incrementors numbers for heisenberg structures.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-24
 */


#ifndef _SG_R_INCREMENTOR_ACCESSOR_H
#define _SG_R_INCREMENTOR_ACCESSOR_H 1

// Global includes
#include <iostream>
#include <string>

#include <boost/operators.hpp>
#include <boost/function.hpp>
#include <boost/ref.hpp>

// Local includes
#include "structure/geometric/hyperbolic/incrementor/r_incrementor.hpp"
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
namespace accessor
{
/**
 *  @class    r_incrementor_accessor
 *  @brief    An accessor class for real numbers
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-24
 *  @note     The common object is constructed from a constant reference to the
 *            incrementor function object. As such altering the common
 *            interface will result in the object no longer being valid.
 *            ( Although it will compile, the decision to allow it to compile
 *            was based on the fact that it makes the common interface easier
 *            to implement ).
 */
template <class _Point>
    class r_incrementor_accessor
  : private boost::equality_comparable< r_incrementor_accessor<_Point>,
              boost::equality_comparable< r_incrementor_accessor<_Point>,
                                          incrementor::r_incrementor<_Point>
                                        >
            >
{
  public:
    //! Object's type
    typedef r_incrementor_accessor<_Point> self_type;
    //! The real type
    typedef incrementor::r_incrementor<_Point> value_type;
    //! The common return type with zeta
    typedef boost::function< typename value_type::result_type
                             ( typename value_type::first_argument_type,
                               typename value_type::second_argument_type
                             )
                           > return_type;
    
    
    /**
     *  @brief  Empty constructor
     */
    r_incrementor_accessor( ) : value_(), common_( boost::cref( value_ ) ) { }
    
    
    /**
     *  @brief  Copy constructor
     *  @param  that The incrementor to construct from
     */
    r_incrementor_accessor( const self_type& __that )
  : value_( __that.value_ ),  common_( boost::cref( value_ ) ) { }
    
    
    /**
     *  @brief  Value constructor
     *  @param  value The real value to construct from
     */
    r_incrementor_accessor( const value_type& __value )
    : value_( __value ),  common_( boost::cref( value_ ) ) { }
    
    
    /**
     *  @brief   Copy assignment
     *  @param   that The incrementor to construct from
     *  @return  A reference to this.
     */
    self_type& operator=( const self_type& __that ) {
      value_ = __that.value_;
      return *this;
    }
    
    /**
     *  @brief   Value assignment
     *  @param   value The real value to assign from
     *  @return  A reference to this.
     */
    self_type& operator=( const value_type& __value ) {
      value_ = __value;
      return *this;
    }
    
    
    /**
     *  @brief   A common interface between r and zeta
     *  @return  An object which is constructed from r and shares a common
     *           type with the zeta return type.
     */
    const return_type& common( void ) const { return common_; }
    
    
    /**
     *  @brief  Compares two r_incrementor_accessor objects for equality.
     *  @param  lhs A \c r_incrementor_accessor
     *  @param  rhs A \c r_incrementor_accessor
     */
    bool operator==( const self_type& __rhs ) const {
      return bool( *this == __rhs.value_ );
    }
    
    
    /**
     *  @brief  Compares a r_incrementor_accessor to a real for equality.
     *  @param  lhs A \c r_incrementor_accessor
     *  @param  rhs A \c value_type
     */
    bool operator==( const value_type& __rhs ) const {
      return bool( value_ == __rhs );
    }
    
    
    /**
     *  @brief  Defines the output stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::ostream& operator<<( std::ostream& __os, const self_type& __t ) {
      typedef utility::functors::stream_cast<value_type,std::ostream> _sc;
      
      _sc( )( __t.value_, __os );
      return __os;
    }
    
    /**
     *  @brief  Defines the input stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::istream& operator>>( std::istream& __is, self_type& __t ) {
      typedef utility::functors::stream_cast<std::istream,value_type> _sc;
      
      _sc( )( __is, __t.value_ );
      return __is;
    }
    
    
  protected:
    //! Incrementor function object
    value_type value_;
    //! Common function object interface
    return_type common_;
};

}// accessor
}// incrementor
}// hyperbolic
}// geometric
}// structure
}// sg
#endif
