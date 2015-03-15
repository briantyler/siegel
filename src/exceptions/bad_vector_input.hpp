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
 *  @file     bad_vector_input.hpp
 *  @brief    An inline header file for bad_vector_input.
 *  @note     Include to throw and handle bad vector inputs.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-03
 */

#ifndef _SG_BAD_VECTOR_INPUT_H
#define _SG_BAD_VECTOR_INPUT_H 1

// Global includes
#include <sstream>
#include <exception>


namespace sg
{
namespace exceptions
{
/**
 *  @class  bad_vector_input
 *  @brief  Indicates that a vector representation of a string is misformed.
 *  @note   Graphically indicates the point in the string that the first
 *          error occurs.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-03
 */

class bad_vector_input : public std::exception
{
  public:
    /**
     *  @brief  Constructs the exception with information about the error.
     *  @param  position The position in the string where the first error
     *          occured.
     *  @param  input The string in which there is an error.
     */
    bad_vector_input( const int __position, const std::string& __input ) throw()
    :position_( __position ), input_( __input ) { }
    
    
    /**
     *  @brief  Virtual Destructor
     */
    virtual ~bad_vector_input( void ) throw() { }
    
    
    /**
     *  @brief  Returns an error message.
     *  @return A C-style error string.
     *  @throw  nothing
     */
    virtual const char* what() const throw(){
      std::stringstream ss;
      ss << "Malformed input vector:" << std::endl;
      ss << error_position_string();
      return ss.str().c_str();
    }
    
    
    /**
     *  @brief  Returns the position at which the error occured.
     *  @throw  nothing
     */
    int position( void ) const throw(){
      return position_;
    }
    
    
  private:
    /**
     *  @brief  Empty constructor.
     *  @note   Private so non-constructable by this method.
     */
    bad_vector_input( void ) throw()
    :position_(0), input_("") { }
    
    
    /**
     *  @brief  Returns a string indicating where the first error occurs.
     *  @throw  nothing
     */
    std::string error_position_string( void ) const throw(){
      /* Error string looks like:
       *
       * [1,2,[3,4,5,[6,7],]]
       * -----------------^
       *
       * Indicating that the comma invalidates the vector
       */
      using std::string;

      string errorString( input_ );
      string errorIndicator( position_-1, '-' );
      errorString += '\n' + errorIndicator + '^';

      return errorString;
    }
    
    
    /**
     *  @var    position_
     *  @note   The position the error occured at.
     *  @var    input_
     *  @note   The input string in which the error occured.
     */
    int position_;
    std::string input_;
    
    
};
  
  
}// exceptions
}// sg

#endif
