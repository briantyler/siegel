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
 *  @file     real_accessor.hpp
 *  @brief    An inline header file for the \c real_accessor class.
 *  @note     The \c real_accessor class creates a common interface to
 *            real numbers for hyperbolic structures.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-28
 */


#ifndef _SG_REAL_ACCESSOR_H
#define _SG_REAL_ACCESSOR_H 1

// Global includes
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
namespace euclidean
{
namespace accessor
{
/**
 *  @class    real_accessor
 *  @brief    An accessor class for real numbers
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-28
 *  @param    _Numerical A real numerical type
 */
template <class _Numerical>
    class real_accessor
  : private boost::equality_comparable< real_accessor<_Numerical>
    , boost::equality_comparable< real_accessor<_Numerical>, _Numerical > >
{
  public:
    //! Object's type
    typedef real_accessor<_Numerical> self_type;
    //! The real type
    typedef _Numerical value_type;
    //! The common return type of with zeta
    typedef value_type return_type;
    
  private:
    //! Value param_type
    typedef typename boost::call_traits<value_type>::param_type param_type;
    
    
  public:
    /**
     *  @brief  Empty constructor
     */
    real_accessor( ) : value_(0) { }
    
    
    /**
     *  @brief  Value constructor
     *  @param  value The real value to construct from
     */
    real_accessor( param_type __value ) : value_(__value) { }
    
    
    /**
     *  @brief   Value assignment
     *  @param   value The real value to assign from
     *  @return  A reference to this.
     */
    self_type& operator=( param_type __value ) {
      value_ = __value;
      return *this;
    }
    
    
    /**
     *  @brief   A common interface between r and zeta
     *  @return  An object which is constructed from r and shares a common
     *           type with the zeta return type.
     */
    const return_type& common( void ) const { return value_; }
    
    
    /**
     *  @brief  Compares two real_accessor objects for equality.
     *  @param  rhs A \c real_accessor
     */
    bool operator==( const self_type& __rhs ) const {
      return bool( *this == __rhs.value_ );
    }
    
    
    /**
     *  @brief  Compares a real_accessor to a real for equality.
     *  @param  rhs A \c value_type
     */
    bool operator==( param_type __rhs ) const {
      typedef utility::math::is_equal<value_type> _ie;
      
      return bool( _ie( )( value_, __rhs ) );
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
    //! _Numerical object
    value_type value_;
};

}// accessor
}// euclidean
}// geometric
}// structure
}// sg
#endif
