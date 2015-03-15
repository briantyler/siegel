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
 *  @file     region_data.hpp
 *  @brief    An inline header file for the \c region_data class.
 *  @note     Represents the division of a region in complex space.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-17
 */


#ifndef _SG_REGION_DATA_H
#define _SG_REGION_DATA_H 1

// Global includes
#include <boost/tuple/tuple.hpp>
#include <boost/operators.hpp>

// Local includes
#include "structure/geometric/euclidean/detail/interval_data.hpp"
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
namespace detail
{
/**
 *  @class    region_data
 *  @brief    A class representing the division of a region in complex space.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-17
 *  @param    _Float A floating point type; defaults to \c double
 *  @param    _Integer An integer type; defaults to \c long
 */
template <class _Float = double, class _Integer = long> class region_data
  : private boost::equality_comparable< region_data<_Float,_Integer> >
{
  public:
    //! Object type
    typedef region_data<_Float,_Integer> self_type;
    //! Real type
    typedef _Float float_type;
    //! Integer type
    typedef _Integer integer_type;
    //! Interval data type
    typedef interval_data<float_type,integer_type> interval_data_type;
   
  private:
    //! The tuple type containing the two intervals
    typedef boost::tuple<interval_data_type,interval_data_type> tuple_type;
    
  public:
    /**
     *  @brief  Get a constant reference to the real interval data
     */
    const interval_data_type& real( void ) const {
      return real_;
    }
    
    /**
     *  @brief  Get a reference to the real interval data
     */
    interval_data_type& real( void ) {
      return const_cast<interval_data_type&>
          ( static_cast<const self_type&>(*this).real() );
    }
    
    /**
     *  @brief  Get a constant reference to the imaginary interval data
     */
    const interval_data_type& imag( void ) const {
      return imag_;
    }
    
    /**
     *  @brief  Get a reference to the imaginary interval data
     */
    interval_data_type& imag( void ) {
      return const_cast<interval_data_type&>
          ( static_cast<const self_type&>(*this).imag() );
    }
    
    
    /**
     *  @brief  Defines the output stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::ostream& operator<<( std::ostream& __os, const self_type& __t ) {
      typedef utility::functors::stream_cast<interval_data_type,std::ostream> _sc;
      
      __os << '[';
      _sc( )( __t.real(), __os );
      __os << ',';
      _sc( )( __t.imag(), __os );
      __os << ']';
      
      return __os;
    }
    
    /**
     *  @brief  Defines the input stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::istream& operator>>( std::istream& __is, self_type& __t ) {
      typedef utility::functors::stream_cast<std::string,interval_data_type> _sc;
      typedef utility::io::string_parser _sp;
      
      std::string s;
      __is >> s;
      
      _sp::string_vector_type vs = _sp()(s);
      _sc( )( vs.at(0), __t.real() );
      _sc( )( vs.at(1), __t.imag() );
      
      return __is;
    }
    
    /**
     *  @brief  Compares two region_data objects for equality.
     *  @param  rhs A \c region_data object
     */
    bool operator==( const self_type& __rhs ) const {
      return bool( real_ == __rhs.real_ && imag_ == __rhs.imag_ );
    }
    
  private:
    interval_data_type real_;
    interval_data_type imag_;
};

}// detail
}// euclidean
}// geometric
}// structure
}// sg
#endif
