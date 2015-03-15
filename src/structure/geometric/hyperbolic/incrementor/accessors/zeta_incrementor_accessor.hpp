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
 *  @file     zeta_incrementor_accessor.hpp
 *  @brief    An inline header file for the \c zeta_incrementor_accessor class.
 *  @note     The \c zeta_incrementor_accessor class creates a common interface
 *            to zeta_incrementors for heisenberg structures.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-24
 */


#ifndef _SG_ZETA_INCREMENTOR_ACCESSOR_H
#define _SG_ZETA_INCREMENTOR_ACCESSOR_H 1

// Global includes
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/operators.hpp>
#include <boost/function.hpp>
#include <boost/ref.hpp>

// Local includes
#include "structure/geometric/hyperbolic/incrementor/zeta_incrementor.hpp"
#include "utility/functors/stream_cast_tuple.hpp"


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
 *  @class    zeta_incrementor_accessor
 *  @brief    An accessor class for complex numbers
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-24
 *  @param    _Point The hyperbolic point type to increment
 */
template <class _Point>
    class zeta_incrementor_accessor
  : private boost::equality_comparable< zeta_incrementor_accessor<_Point>,
              boost::equality_comparable< zeta_incrementor_accessor<_Point>,
                                          incrementor::zeta_incrementor<_Point>
                                        >
              >
{
  public:
    //! Object's type
    typedef zeta_incrementor_accessor<_Point> self_type;
    //!Point type
    typedef _Point point_type;
    //! Type of the real component of zeta
    typedef incrementor::zeta_incrementor<_Point> real_type;
    //! Type of a coordinate of zeta
    typedef boost::tuple<real_type,real_type> zeta_type;
    //! Type of the imaginary component of zeta
    typedef real_type imag_type;
    //! The common return type of real and imaginary parts
    typedef boost::function< typename real_type::result_type
                             ( typename real_type::first_argument_type,
                               typename real_type::second_argument_type
                             )> return_type;
    
  private:
    //! The return type tuple
    typedef boost::tuple<return_type,return_type> return_tuple_type;
    
    
  public:
    /**
     *  @brief  Empty constructor
     */
    zeta_incrementor_accessor( )
    : zeta_(),
      common_( boost::cref( boost::get<0>( zeta_ ) ),
               boost::cref( boost::get<1>( zeta_ ) )
             )
    { }
    
    /**
     *  @brief  Copy constructor
     *  @param  that The accessor to construct from
     */
    zeta_incrementor_accessor( const self_type& __that )
  : zeta_(__that.zeta_),
    common_( boost::cref( boost::get<0>( zeta_ ) ),
             boost::cref( boost::get<1>( zeta_ ) )
           )
    { }
    
    /**
     *  @brief  Zeta constructor
     *  @param  zeta The zeta value to construct from
     */
    zeta_incrementor_accessor( const zeta_type& __zeta )
  : zeta_(__zeta),
    common_( boost::cref( boost::get<0>( zeta_ ) ),
             boost::cref( boost::get<1>( zeta_ ) )
           )
    { }
    
    
    /**
     *  @brief   Zeta assignment
     *  @param   zeta The zeta value to assign from
     *  @return  A reference to this.
     */
    self_type& operator=( const self_type& __that ) {
      zeta_ = __that.zeta_;
      return *this;
    }
    
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
      return boost::get<0>( zeta_ );
    }
    
    
    /**
     *  @brief   Get a constant reference to the imaginary part of zeta
     *  @return  A constant reference to the imaginary part of zeta.
     */
    const imag_type& imag_ref( void ) const {
      return boost::get<1>( zeta_ );
    }
    
    
    /**
     *  @brief   A common interface between real and imaginary parts
     *  @return  An object which is constructed from the real part of zeta
     *           and shares a common type with an imaginary return type.
     */
    const return_type& real_common( void ) const {
      return boost::get<0>( common_ );
    }
    
    /**
     *  @brief   A common interface between real and imaginary parts
     *  @return  An object which is constructed from the imag part of zeta
     *           and shares a common type with a real return type.
     */
    const return_type& imag_common( void ) const {
      return boost::get<1>( common_ );
    }
    
    
    /**
     *  @brief  Compares two zeta_incrementor_accessor objects for equality.
     *  @param  rhs A \c zeta_incrementor_accessor
     */
    bool operator==( const self_type& __rhs ) const {
      return bool( *this == __rhs.zeta_ );
    }
    
    
    /**
     *  @brief  Compares two zeta_incrementor_accessor objects for equality.
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
      typedef utility::functors::stream_cast_tuple<zeta_type,std::ostream> _sc;
      
      _sc( )( __t.zeta_, __os );
      return __os;
    }
    
    /**
     *  @brief  Defines the input stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::istream& operator>>( std::istream& __is, self_type& __t ) {
      typedef utility::functors::stream_cast_tuple<std::istream,zeta_type> _sc;
      
      _sc( )( __is, __t.zeta_ );
      return __is;
    }
    
    
  protected:
    //! Zeta component
    zeta_type zeta_;
    //! The common type tuple
    return_tuple_type common_;
};

}// accessor
}// incrementor
}// hyperbolic
}// geometric
}// structure
}// sg
#endif
