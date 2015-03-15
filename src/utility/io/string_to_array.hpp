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
 *  @file     string_to_array.hpp
 *  @brief    An inline header file for \c string_to_array.
 *  @note     Include to convert strings to \c boost::arrays
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-11
 */

#ifndef _SG_STRING_TO_ARRAY_H
#define _SG_STRING_TO_ARRAY_H 1

// Global includes
#include <boost/iterator/iterator_traits.hpp>

// Local includes
#include "utility/io/string_parser.hpp"
#include "utility/functors/stream_cast.hpp"


namespace sg
{
namespace utility
{
namespace io
{
/**
 *  @class    string_to_array
 *  @brief    Functor for parsing a string to fixed size container.
 *  @param    Iterator The type of the output iterator.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-11
 */
template <class _Iterator> class string_to_array
{
  public:
    //! The object type
    typedef string_to_array<_Iterator> self_type;
    //! The value type
    typedef typename boost::iterator_value<_Iterator>::type value_type;
    //! The string parser type
    typedef string_parser parser_type;
    //! The string vector type
    typedef typename parser_type::string_vector_type string_vector_type;
    //! The string type
    typedef typename parser_type::string_type string_type;
    
    /**
     *  @brief  Default constructor
     */
    string_to_array( ) : sp_() { }
    
    /**
     *  @brief  \c string_parser constructor
     *  @note   Use this constructor if you require the string to be parsed
     *          using non-standard parenthesis and delimeters.
     */
    string_to_array( const parser_type& __sp ) : sp_(__sp) { }
    
    
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
     *  @brief  Converts a string of homogeneous objects into an array.
     *  @param  input The string to convert to an array.
     *  @return A array formed from the string.
     */
    void operator( )
        ( const string_type& __input, _Iterator __begin, _Iterator __end )
    {
      typedef functors::stream_cast<string_type,value_type> _sc;
      
      string_vector_type vs( sp_( __input ) );
      string_vector_type::const_iterator it = vs.begin();
      for( ; it != vs.end() && __begin != __end; ++it, ++__begin ) {
        _sc()( *it, *__begin );
      }
      
      if( it != vs.end() || __begin != __end ) {
        throw exceptions::bad_vector_input( 1, __input );
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
