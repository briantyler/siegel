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
 *  @file     hyperbolic_structure.hpp
 *  @brief    An inline header file for the \c hyperbolic_structure class.
 *  @note     The \c hyperbolic_structure class provides a common interface to
 *            all hyperbolic objects.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-04
 */


#ifndef _SG_HYPERBOLIC_STRUCTURE_H
#define _SG_HYPERBOLIC_STRUCTURE_H 1

// Global includes
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/operators.hpp>

// Local includes
#include "structure/geometric/detail/heisenberg_structure.hpp"
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
 *  @class    hyperbolic_structure
 *  @brief    A common interface to all hyperbolic objects
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-04
 *  @param    N The (complex hyperbolic) dimension of the space
 *  @param    ZetaAccessor The accessor for the zeta component
 *  @param    RAccessor The accessor for the r component
 *  @param    HeightAccessor The accessor for the height component
 *  @param    Dim The real dimension of the object
 */
template < std::size_t N, class _ZetaAccessor,
           class _RAccessor,class _HeightAccessor
         >
    class hyperbolic_structure
  : public heisenberg_structure<N,_ZetaAccessor,_RAccessor>
  , private boost::equality_comparable<
                hyperbolic_structure< N, _ZetaAccessor, _RAccessor,
                                      _HeightAccessor
                                    >
            >
{
  public:
    //! Object type
    typedef hyperbolic_structure <N,_ZetaAccessor,_RAccessor,_HeightAccessor>
        self_type;
    //! The hyperbolic structure type
    typedef self_type hyperbolic_structure_type;
    
    //! The type of the real coordinate
    typedef real_coordinate<_HeightAccessor> height_coordinate_type;
    
    // Ensure that the return types match
    BOOST_MPL_ASSERT( (
          boost::is_same<
                          typename height_coordinate_type::return_type
                        ,
                          typename self_type::return_type
                        >
                      )
                    );
    
  public:
    //! Height type
    typedef typename height_coordinate_type::value_type height_type;
    
    // Pull up some typedefs
    typedef typename self_type::size_type size_type;
    typedef typename self_type::return_type return_type;
    typedef typename self_type::zeta_array_type zeta_array_type;
    typedef typename self_type::heisenberg_structure_type heisenberg_structure_type;
    typedef typename self_type::iterator iterator;
    typedef typename self_type::const_iterator const_iterator;
    
    
    /**
     *  @brief  Default constructor.
     */
    hyperbolic_structure( ) : heisenberg_structure_type() {
      self_type::add_ptr_at( self_type::hyperbolic_size-1, &( height_.common() ) );
    }
    
    
    /**
     *  @brief  Copy constructor.
     *  @param  that The hyperbolic structure to copy
     */
    hyperbolic_structure( const self_type& __that )
  : heisenberg_structure_type( __that ), height_( __that.height_ )
    {
      self_type::add_ptr_at( self_type::hyperbolic_size-1, &( height_.common() ) );
    }
    
    
    /**
     *  @brief  Assignment operator
     *  @param  that The heisenberg structure to assign from
     *  @return \c *this
     */
    hyperbolic_structure& operator=( const self_type& __that ) {
      heisenberg_structure_type::operator=( __that );
      height_ = __that.height_;
      return *this;
    }
    
    
    /**
     *  @brief Returns a contant reference to the height coordinate.
     */
    const height_coordinate_type& height( void ) const { return height_; }
    
    /**
     *  @brief Returns a reference to the height coordinate.
     */
    height_coordinate_type& height( void ) { return height_; }
    
    
    /**
     *  @brief Returns a contant reference to the underlying height
     *         coordinate data.
     */
    const height_type& height_ref( void ) const { return height_.ref(); }
    
    /**
     *  @brief Returns a reference to the underlying height coordinate data.
     */
    height_type& height_ref( void ) { return height_.ref(); }
    
    
    /**
     *  @brief Returns a contant reference to the common height coordinate data.
     */
    const return_type& height_common( void ) const { return height_.common(); }
    
    /**
     *  @brief Returns a reference to the common height coordinate data.
     */
    return_type& height_common( void ) { return height_.common(); }
    
    // Now the object looks like a Hyperbolic object of dimension N
    
    
    /**
     *  @brief  Provide a sentinal iterator for the end of the common structure.
     *  @return An iterator pointing one past the end of the common structure.
     */
    iterator end( void ) { return hyperbolic_end(); }
    
    /**
     *  @brief  Provide a sentinal iterator for the end of the common structure.
     *  @return A const iterator pointing one past the end of the common
     *          structure.
     */
    const_iterator end( void ) const { return hyperbolic_end(); }
    
    /**
     *  @brief  Provide a sentinal iterator for the end of the common
     *          hyperbolic structure.
     *  @return An iterator pointing one past the end of the common
     *          hyperbolic structure.
     */
    iterator hyperbolic_end( void ) {
      return iterator(   self_type::pointers_.begin()
                       + self_type::hyperbolic_size
                     );
    }
    
    /**
     *  @brief  Provide a constant sentinal iterator for the end of the common
     *          hyperbolic structure.
     *  @return A constant iterator pointing one past the end of the common
     *          hyperbolic structure.
     */
    const_iterator hyperbolic_end( void ) const {
      return const_iterator(   self_type::pointers_.begin()
                             + self_type::hyperbolic_size
                           );
    }
    
    
    /**
     *  @brief   Compare two \c hyperbolic_structures for equality
     *  @param   lhs A Hyperbolic structure
     *  @return  True if they are equal.
     */
    bool operator==( const self_type& __rhs ) const  {
      typedef const typename self_type::heisenberg_structure_type& _ht;
      
      return bool(    ( static_cast<_ht>(*this) == static_cast<_ht>(__rhs) )
                   &&
                      ( height_  ==  __rhs.height_ )
                 );
    }
    
    
    /**
     *  @brief  Defines the output stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::ostream& operator<<( std::ostream& __os, const self_type& __t ){
      using utility::functors::stream_cast;
      typedef stream_cast<height_coordinate_type, std::ostream> _sc;
      
      __os << '[' << static_cast<const heisenberg_structure_type&>(__t) << ',';
      _sc()( __t.height_, __os );
      __os << ']';
      return __os;
    }
    
    /**
     *  @brief  Defines the input stream operator
     *  @param  __is An istream.
     *  @param  __t The object to stream out.
     */
    friend std::istream& operator>>( std::istream& __is, self_type& __t ){
      using utility::functors::stream_cast;
      using utility::io::string_parser;
      using std::string;
      
      typedef typename self_type::heisenberg_structure_type _ht;
      
      typedef stream_cast<string,_ht> _scb;
      typedef stream_cast<string,height_coordinate_type> _sch;
      
      string s;
      __is >> s;
      
      string_parser::string_vector_type vs = string_parser( )(s);
      
      _scb( )( vs.at(0), static_cast<_ht&>(__t) );
      _sch( )( vs.at(1), __t.height_ );
        
      return __is;
    }
    
    
  private:
    height_coordinate_type height_;
};

}// detail
}// geometric
}// structure
}// sg
#endif
