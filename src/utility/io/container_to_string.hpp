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
 *  @file     container_to_string.hpp
 *  @brief    An inline header file for \c container_to_string.
 *  @note     Include to convert \b containers to strings.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-06
 */

#ifndef _SG_CONTAINER_TO_STRING_H
#define _SG_CONTAINER_TO_STRING_H 1

// Global includes
#include <string>
#include <functional>
#include <boost/iterator/iterator_traits.hpp>

// Local includes
#include "utility/functors/stream_cast.hpp"


namespace sg
{
namespace utility
{
namespace io
{
/**
 *  @class    container_to_string
 *  @brief    Functor for parsing a container to a string.
 *  @param    _Iterator The type of the input iterator.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-06
 */
  template <class _Iterator> class container_to_string
: public std::binary_function< _Iterator, _Iterator, std::string >
{
  public:
    //! The object type
    typedef container_to_string<_Iterator> self_type;
    //! The value type
    typedef typename boost::iterator_value<_Iterator>::type value_type;
    //! The iterator type
    typedef _Iterator iterator_type;
    //! The string type
    typedef std::string string_type;
    
    /**
     *  @brief  Empty constructor
     *  @note   Default parenthesis and delimeter:
     *          <ul>
     *            <li>\b open "["
     *            <li>\b close "]"
     *            <li>\b delimeter ","
     *          </ul>
     */
    container_to_string( ) : open_('['), close_(']'), delimeter_(',') { }
    
    
    /**
     *  @brief  Parenthesis and delimeter constructor
     *  @param  open The opening parenthesis
     *  @param  close The closing parenthesis
     *  @param  delimeter The delimiter
     */
    container_to_string
        ( const char __open, const char __close, const char __delimeter )
    :open_(__open), close_(__close), delimeter_(__delimeter) { }
    
    
    /**
     *  @brief  Set the opening parenthesis
     *  @param  c The opening parenthesis character.
     */
    void set_open( const char __c ) { open_ = __c; }
    
    /**
     *  @brief  Set the closing parenthesis
     *  @param  c The closing parenthesis character.
     */
    void set_close( const char __c ) { close_ = __c; }
    
    /**
     *  @brief  Set the delimeter
     *  @param  c The delimeter character.
     */
    void set_delimeter( const char __c ) { delimeter_ = __c; }
    
    
    /**
     *  @brief  Converts a container into a string representation.
     *  @param  begin an iterator pointing to the begining of the container.
     *  @param  end an iterator pointing one past the end of the container.
     *  @return A string representation of the container.
     */
    string_type operator( )
         ( iterator_type __begin, iterator_type __end ) const
    {
      typedef functors::stream_cast<value_type, string_type>  _sc;
      
      string_type element, output;
      // Open parenthesis
      output += open_;
      for( iterator_type it = __begin; it != __end; ++it ){
        // Cast the string to an element
        _sc()( *it, element );
        
        // Add the element
        output += element + delimeter_;
      }
      // Remove the final delimeter
      if ( output.size() > 1 ) output.erase( output.end() - 1 );
      
      // Close parenthesis
      output += close_;
      return output;
    }
    
    
  private:
    //!  The opening parenthesis
    char open_;
    
    //! The closing parenthesis
    char close_;
    
    //!The delimeter
    char delimeter_;
};

}// io
}// utility
}// sg
#endif
