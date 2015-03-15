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
 *  @file     rectangle.hpp
 *  @brief    An inline header file for the \c rectangle class.
 *  @version  1.0
 *  @date     2008-04-04
 */


#ifndef _SG_RECTANGLE_H
#define _SG_RECTANGLE_H 1

// Global includes
#include <cassert>
#include <complex>
#include <algorithm>
#include <functional>

#include <boost/array.hpp>
#include <boost/call_traits.hpp>
#include <boost/algorithm/minmax_element.hpp>
#include <boost/operators.hpp>

// Local includes
#include "utility/functors/stream_cast.hpp"
#include "utility/io/string_parser.hpp"
#include "utility/math/is_equal_cx.hpp"


namespace sg
{
namespace structure
{
namespace geometric
{
namespace euclidean
{
/**
 *  @class    rectangle
 *  @brief    A class modeling a rectagular region in complex space.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-04
 *  @param    _Float A floating point type; defaults to \c double
 */
template <class _Float = double> class rectangle
  : private boost::equality_comparable< rectangle<_Float> >
{
  public:
    //! Object type
    typedef rectangle<_Float> self_type;
    //! Floating point type
    typedef _Float float_type;
    //! The complex type
    typedef std::complex<float_type> complex_type;
    
  private:
    /**
     *  @brief The array type holding the vertices of the rectangle
     *  @note [bl tl br tr] t=top, b=bottom, r=right and l=left
     */
    typedef boost::array<complex_type,4> array_type;
    
  public:
    /**
     *  @brief  Default constructor.
     */
    rectangle() { };
    
    /**
     *  @brief  Corner constructor.
     *  @param  bl The bottom lefthand corner.
     *  @param  tr The top righthand corner.
     */
    rectangle( const complex_type& __bl, const complex_type& __tr ) {
      set( __bl, __tr );
    };
    
    /**
     *  @brief  Set rectangle from its bottom left and top righthand corners.
     *  @param  bl The bottom lefthand corner.
     *  @param  tr The top righthand corner.
     */
    void set( const complex_type& __bl, const complex_type& __tr ) {
      // At debug time ensure the coordinates are in the right order.
      assert( __bl.real() <= __tr.real() && __bl.imag() <= __tr.imag() );
      
      sys_bl() = __bl;
      sys_tr() = __tr;
      
      // Sometimes you just need to have a few dumb lines of code
      sys_tl().real() = bl().real();
      sys_br().real() = tr().real();
      sys_br().imag() = bl().imag();
      sys_tl().imag() = tr().imag();
    }
    
    
    /**
     *  @brief  Multiply the rectangle by a complex number, then make this the
     *          smallest rectangle that contains the transformed shape.
     *  @param  transform The complex number to transform by.
     */
    void transform_contain( const complex_type& __transform ) {
      // Transform the rectangle
      typedef typename array_type::iterator _iter;
      for( _iter it = rectangle_.begin(); it != rectangle_.end(); ++it ) {
        *it *= __transform;
      }
      
      // Get the minimum and maximum real and imaginary components
      typedef std::pair< _iter, _iter > _pi;
      _pi m1 = boost::minmax_element( rectangle_.begin( ),
                                      rectangle_.end( ), prv_compare_re( )
                                    );
      _pi m2 = boost::minmax_element( rectangle_.begin( ),
                                      rectangle_.end( ), prv_compare_im( )
                                    );
      
      // The bottom left point is set from the minimal real and imaginary parts
      // The top right point is set from the maximal real and imaginary parts
      set( complex_type( (m1.first)->real(), (m2.first)->imag() ),
           complex_type( (m1.second)->real(), (m2.second)->imag() )
          );
    }
    
    
  private:
    /**
     *  @struct prv_compare_re
     *  @brief  A helper comparison functor for comparing the real parts
     *          of a complex number.
     */
    struct prv_compare_re{
      bool operator()
          ( const complex_type& __lhs, const complex_type& __rhs ) const
      {
        return bool( __lhs.real() < __rhs.real() );
      }
    };
    
    /**
     *  @struct prv_compare_im
     *  @brief  A helper comparison functor for comparing the imaginary parts
     *          of a complex number.
     */
    struct prv_compare_im{
      bool operator()
          ( const complex_type& __lhs, const complex_type& __rhs ) const
      {
        return bool( __lhs.imag() < __rhs.imag() );
      }
    };
    
    
  public:
    /**
     *  @brief   Get the complex number at the bottom left of the rectangle.
     *  @return  The complex number at the bottom left of the rectangle.
     */
    const complex_type& bl( void ) const {
      return rectangle_[0];
    }
    
    /**
     *  @brief   Get the complex number at the top left of the rectangle.
     *  @return  The complex number at the top left of the rectangle.
     */
    const complex_type& tl( void ) const {
      return rectangle_[1];
    }
    
    /**
     *  @brief   Get the complex number at the bottom right of the rectangle.
     *  @return  The complex number at the bottom right of the rectangle.
     */
    const complex_type& br( void ) const {
      return rectangle_[2];
    }
    
    /**
     *  @brief   Get the complex number at the top right of the rectangle.
     *  @return  The complex number at the top rightof the rectangle.
     */
    const complex_type& tr( void ) const {
      return rectangle_[3];
    }
    
    
    /**
     *  @brief  Convert a complex region to a rectangle
     *  @param  region The region to convert from
     */
    template <class _Region> void from_region( const _Region& __region ) {
      set( __region.bl(), __region.tr() );
    }
    
    
    /**
     *  @brief  Compares two rectangle objects for equality.
     *  @param  rhs A \c rectangle
     */
    bool operator==( const self_type& __rhs ) const {
      typedef utility::math::is_equal_cx<float_type> _ie;
      
      return bool(    _ie( )( bl(), __rhs.bl() )
                   && _ie( )( tr(), __rhs.tr() )
                 );
    }
    
    
    /**
     *  @brief  Defines the output stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::ostream& operator<<( std::ostream& __os, const self_type& __t ){
      typedef utility::functors::stream_cast<complex_type,std::ostream> _sc;
      
      __os << '[';
      _sc()( __t.bl(), __os );
      __os << ',';
      _sc()( __t.tr(), __os );
      __os << ']';
      
      return __os;
    }
    
    /**
     *  @brief  Defines the input stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::istream& operator>>( std::istream& __is, self_type& __t ){
      typedef utility::functors::stream_cast<std::string,complex_type> _sc;
      typedef  utility::io::string_parser _sp;
      
      _sp::string_vector_type vs;
      std::string s;
      __is >> s;
      vs = _sp()( s );
      assert( vs.size() == 2 );
      
      complex_type bl, tr;
      _sc()( vs.at(0), bl );
      _sc()( vs.at(1), tr );
      
      __t.set(bl,tr);
      
      return __is;
    }
    
    
  private:
    // Note that these provide easy access, but they should not be public as
    // it is a bad idea to allow write access to the corners.
    // Do this with const casting the const members in case the underlying data
    // structure changes - makes code more flexible and the overhead is
    // optimised out at compile time when the build flags are on.
    /**
     *  @brief   Get the complex number at the bottom left of the rectangle.
     *  @return  The complex number at the bottom left of the rectangle.
     */
    complex_type& sys_bl( void ) {
      return const_cast<complex_type&>
          ( static_cast<const self_type&>( *this ).bl() );
    }
    
    /**
     *  @brief   Get the complex number at the top left of the rectangle.
     *  @return  The complex number at the top left of the rectangle.
     */
    complex_type& sys_tl( void ) {
      return const_cast<complex_type&>
          ( static_cast<const self_type&>( *this ).tl() );
    }
    
    /**
     *  @brief   Get the complex number at the bottom right of the rectangle.
     *  @return  The complex number at the bottom right of the rectangle.
     */
    complex_type& sys_br( void ) {
      return const_cast<complex_type&>
          ( static_cast<const self_type&>( *this ).br() );
    }
    
    /**
     *  @brief   Get the complex number at the top right of the rectangle.
     *  @return  The complex number at the top rightof the rectangle.
     */
    complex_type& sys_tr( void ) {
      return const_cast<complex_type&>
          ( static_cast<const self_type&>( *this ).tr() );
    }
    
    
    //! The rectangle
    array_type rectangle_;
};

}// euclidean
}// geometric
}// structure
}// sg
#endif
