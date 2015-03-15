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
 *  @file     string_to_container.hpp
 *  @brief    An inline header file for \c string_to_container.
 *  @note     Include to convert strings to \c stl::containers
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-06
 */

#ifndef _SG_STRING_TO_CONTAINER_H
#define _SG_STRING_TO_CONTAINER_H 1

// Global includes
#include <string>
#include <functional>

// Local includes
#include "utility/functors/stream_cast.hpp"
#include "utility/io/string_parser.hpp"


namespace sg
{
namespace utility
{
namespace io
{
/**
 *  @class    string_to_container
 *  @brief    Functor for parsing a string to an stl container
 *  @param    _Container An STL Container
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-06
 */
template <class _Container> class string_to_container
  :public std::binary_function< const std::string&, _Container&, void >
{
  public:
    //! The object type
    typedef string_to_container<_Container> self_type;
    //! The container type
    typedef _Container container_type;
    //! The value type
    typedef typename container_type::value_type value_type;
    //! The string parser type
    typedef string_parser parser_type;
    //! The string vector type
    typedef typename parser_type::string_vector_type string_vector_type;
    //! The string type
    typedef typename parser_type::string_type string_type;
    
    /**
     *  @brief  Default constructor
     */
    string_to_container( ) : sp_() { }
    
    /**
     *  @brief  \c string_parser constructor
     *  @note   Use this constructor if you require the string to be parsed
     *          using non-standard parenthesis and delimeters.
     */
    string_to_container( const parser_type& __sp ) : sp_(__sp) { }
    
    
    /**
     *  @brief  Access the string_parser object
     *  @return A reference to the string parser which does most of the work.
     */
    parser_type& parser( void ) { return sp_; }
    
    /**
     *  @brief  Access the string_parser object
     *  @return A reference to the string parser which does most of the work.
     */
    const parser_type& parser( void ) const { return sp_; }
    
    
    /**
     *  @brief  Converts a string of homogeneous objects into a container.
     *  @param  input The string to convert to a container.
     *  @return A container formed from the string.
     */
    void operator( )(const string_type& __input, container_type& __container) {
      typedef functors::stream_cast<string_type, value_type> _sc;
      typedef string_parser::string_vector_type::const_iterator c_iter;
      
      string_parser::string_vector_type vs( sp_( __input ) );
      
      // Cast each element into an object then push this into the container
      value_type element;
      __container.clear( );
      
      for( c_iter cit = vs.begin(); cit != vs.end(); ++cit ){
        _sc()( *cit, element );
        __container.push_back( element );
      }
    }
    
    
  private:
    //! The parser which carries out most of the work.
    parser_type sp_;
};

}// io
}// utility
}// sg
#endif
