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
 *  @file     heisenberg_structure.hpp
 *  @brief    An inline header file for the \c heisenberg_structure class.
 *  @note     The \c heisenberg_structure class provides a common interface to
 *            the Heisenberg part of all hyperbolic objects.
 *  @author   Brian Tyler
 *  @version  1.1
 *  @date     2008-02-25
 */


#ifndef _SG_HEISENBERG_STRUCTURE_H
#define _SG_HEISENBERG_STRUCTURE_H 1

// Global includes
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/operators.hpp>

// Local includes
#include "structure/geometric/detail/zeta_array.hpp"
#include "structure/geometric/detail/real_coordinate.hpp"
#include "utility/functors/stream_cast.hpp"
#include "utility/io/string_parser.hpp"

namespace sg
{
namespace structure
{
namespace geometric
{
namespace detail
{
/**
 *  @class    heisenberg_structure
 *  @brief    A common interface to the zeta coordinates of different types
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-25
 *  @param    N The (complex hyperbolic) dimension of the space
 *  @param    ZetaAccessor The accessor for the zeta component
 *  @param    RAccessor The accessor for the r component
 */
template <std::size_t N, class _ZetaAccessor, class _RAccessor>
    class heisenberg_structure
  : public zeta_array<N,_ZetaAccessor>,
    private boost::equality_comparable<
                  heisenberg_structure<N,_ZetaAccessor,_RAccessor>
            >
{
  // Ensure that the return types match
    BOOST_MPL_ASSERT( (
          boost::is_same< typename _RAccessor::return_type
                        , typename _ZetaAccessor::return_type
                        >
                      )
                    );
  public:
    //! Object type
    typedef heisenberg_structure<N,_ZetaAccessor,_RAccessor> self_type;
    //! The Heisenberg structure type
    typedef self_type heisenberg_structure_type;
    
    
  public:
    // Pull the typedefs up
    typedef typename self_type::size_type size_type;
    typedef typename self_type::return_type return_type;
    typedef typename self_type::zeta_array_type zeta_array_type;
    typedef typename self_type::iterator iterator;
    typedef typename self_type::const_iterator const_iterator;
    
    //! The r coordinate type
    typedef real_coordinate<_RAccessor> r_coordinate_type;
    //! The r value type
    typedef typename r_coordinate_type::value_type r_type;
    
    
    /**
     *  @brief  Default constructor.
     */
    heisenberg_structure( ) : zeta_array_type() {
      self_type::add_ptr_at( self_type::heisenberg_size-1, &( r_.common() ) );
    }
    
    
    /**
     *  @brief  Copy constructor.
     *  @param  that The heisenberg structure to copy
     */
    heisenberg_structure( const self_type& __that )
  : zeta_array_type( __that ), r_( __that.r_ )
    {
      self_type::add_ptr_at( self_type::heisenberg_size-1, &( r_.common() ) );
    }
    
    
    /**
     *  @brief  Assignment operator
     *  @param  that The heisenberg structure to assign from
     *  @return \c *this
     */
    heisenberg_structure& operator=( const self_type& __that ) {
      zeta_array_type::operator=( __that );
      r_ = __that.r_;
      return *this;
    }
    
    
    /**
     *  @brief Returns a contant reference to the r coordinate.
     */
    const r_coordinate_type& r( void ) const { return r_; }
    
    /**
     *  @brief Returns a reference to the r coordinate.
     */
    r_coordinate_type& r( void ) { return r_; }
    
    
    /**
     *  @brief Returns a contant reference to the underlying r coordinate data.
     */
    const r_type& r_ref( void ) const { return r_.ref(); }
    
    /**
     *  @brief Returns a reference to the underlying r coordinate data.
     */
    r_type& r_ref( void ) { return r_.ref(); }
    
    
    /**
     *  @brief Returns a contant reference to the common r coordinate data.
     */
    const return_type& r_common( void ) const { return r_.common(); }
    
    /**
     *  @brief Returns a reference to the common r coordinate data.
     */
    return_type& r_common( void ) { return r_.common(); }
    // Now the object looks like a heisenberg space of dimension 2N-1
    
    
    /**
     *  @brief  Provide a sentinal iterator for the end of the common structure.
     *  @return An iterator pointing one past the end of the common structure.
     */
    iterator end( void ) { return heisenberg_end(); }
    
    /**
     *  @brief  Provide a sentinal iterator for the end of the common structure.
     *  @return A const iterator pointing one past the end of the common
     *          structure.
     */
    const_iterator end( void ) const { return heisenberg_end(); }
    
    /**
     *  @brief  Provide a sentinal iterator for the end of the common heisenberg
     *          structure.
     *  @return An iterator pointing one past the end of the common heisenberg
     *          structure.
     */
    iterator heisenberg_end( void ) {
      return iterator(   self_type::pointers_.begin()
                       + self_type::heisenberg_size
                     );
    }
    
    /**
     *  @brief  Provide a constant sentinal iterator for the end of the common
     *          heisenberg structure.
     *  @return A constant iterator pointing one past the end of the common
     *          heisenberg structure.
     */
    const_iterator heisenberg_end( void ) const {
      return const_iterator(   self_type::pointers_.begin()
                             + self_type::heisenberg_size
                           );
    }
    
    
    /**
     *  @brief   Compare two \c heisenberg_structures for equality
     *  @param   rhs A Heisenberg structure
     *  @return  True if they are equal.
     */
    bool operator==( const self_type& __rhs ) const {
      return bool(    (    static_cast<const zeta_array_type&>(*this)
                        ==
                           static_cast<const zeta_array_type&>(__rhs)
                      )
                   && ( r_ == __rhs.r_ ) );
    }
    
    /**
     *  @brief  Defines the output stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::ostream& operator<<( std::ostream& __os, const self_type& __t ){
      using utility::functors::stream_cast;
      typedef stream_cast<r_coordinate_type, std::ostream> _sc;
      
      __os << '[' << static_cast<const zeta_array_type&>(__t) << ',';
      _sc()( __t.r_, __os );
      __os << ']';
      return __os;
    }
    
    /**
     *  @brief  Defines the input stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::istream& operator>>( std::istream& __is, self_type& __t ){
      using utility::functors::stream_cast;
      using utility::io::string_parser;
      using std::string;
      
      typedef stream_cast<string,zeta_array_type> _scz;
      typedef stream_cast<string,r_coordinate_type> _scr;
      
      string s;
      __is >> s;
      
      string_parser::string_vector_type vs = string_parser()(s);
      
      _scz( )( vs.at(0), static_cast<zeta_array_type&>(__t) );
      _scr( )( vs.at(1), __t.r_ );
        
      return __is;
    }
    
    
  private:
    r_coordinate_type r_;
};

}// detail
}// geometric
}// structure
}// sg
#endif
