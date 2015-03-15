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
 *  @file     string_parser.hpp
 *  @brief    An inline header file for \c string_parser.
 *  @note     Include to parse strings to a vector of strings
 *  @author   Brian Tyler
 *  @version  1.2
 *  @date     2008-02-06
 */

#ifndef _SG_STRING_PARSER_H
#define _SG_STRING_PARSER_H 1

// Global includes
#include <utility>
#include <string>
#include <vector>
#include <set>
#include <stack>
#include <functional>
#include <algorithm>

#include <boost/algorithm/string.hpp>

// Local includes
#include "exceptions/bad_vector_input.hpp"


namespace sg
{
namespace utility
{
namespace io
{
/**
 *  @class    string_parser
 *  @brief    Functor for parsing a string to a vector of strings
 *  @author   Brian Tyler
 *  @version  1.2
 *  @date     2008-02-05
 */
class string_parser
: public std::unary_function< const std::string&,
                              std::vector<std::string> >
{
  public:
    //! The object type
    typedef string_parser self_type;
    //! The string_type
    typedef std::string string_type;
    //! Type of the string vector
    typedef std::vector<string_type> string_vector_type;
    //! Type of the exception thrown by the functor
    typedef exceptions::bad_vector_input exception_type;
    
    /**
     *  @brief  Default constructor
     *  @note   Default parenthesis and delimeter:
     *          <ul>
     *            <li>\b parenthesis "()", "[]", "{}" and "<>"
     *            <li>\b delimeters "," and " "
     *          </ul>
     *  @note  The only character not suitable for a parenthesis is ' '
     *         (white space) because white space is trimmed from the input
     *         string.
     *  
     *  Constructs a new \c string_parser with default parenthesis and
     *  delimeter.
     */
    string_parser( ) { }
    
    
    /**
     *  @brief  Add new parenthesis
     *  @param  left The lefthand parenthesis
     *  @param  right The rightand parenthesis
     *  @return \c true if the parenthesis was added successfully.
     *  @note   Both \c left and \c right must be different and not ' '
     */
    bool add_parenthesis( const char __left, const char __right ) {
      return punctuation_.add_parenthesis(__left,__right);
    }
    
    /**
     *  @brief   Adds a delimeter
     *  @param   c The \c char for the delimter
     */
    bool add_delimiter( const char __delimiter ) {
      return punctuation_.add_delimiter(__delimiter);
    }
    
    
    /**
     *  @brief  Remove parenthesis from the punctuation
     *  @param  left The lefthand parenthesis
     *  @param  right The rightand parenthesis
     *  @return \c true if the parenthesis was removed successfully.
     */
    bool remove_parenthesis( const char __left, const char __right ) {
      return punctuation_.remove_parenthesis(__left,__right);
    }
    
    
    /**
     *  @brief  Remove delimeter from the punctuation
     *  @param  delimeter The delimiter to remove
     *  @return \c true if the delimiter was removed successfully.
     */
    bool remove_delimiter( const char __delimiter ) {
      return punctuation_.remove_delimiter( __delimiter );
    }
    
    
    /**
     *  @brief  Parses a string to a vector of strings.
     *  @param  input The string to parse
     *  @return The \c input string parsed to a vector.
     *  @note   The parenthesis and delimeters can be modified using the
     *          \c add and \c remove functions.
     */
    string_vector_type operator( )( const string_type& __input )
    {
      // Copy the input string and remove white space
      input_ = boost::trim_copy( __input );
      
      // Reset the members.
      depth_ = 0;
      pos_ = 0;
      current_.clear();
      elements_.clear();
      nextRight_ = std::stack<char>();
      
      // Process the input string character by character
      typedef string_type::const_iterator c_iter;
      for( c_iter it = input_.begin(); it != input_.end(); ++it ) {
        sys_process_string( *it );
      }
      
      if( depth_ != 0 ) throw exception_type( pos_, input_ );
      
      return elements_;
    }
    
    
    /**
     *  @brief  Determines the depth of nested vectors in the string.
     *  @param  input The string to check
     *  @return The depth of nested vectors.
     *  @note   The criteria for a valid depth is that the brackets must
     *          balance and there must be the same number of left and right
     *          brackets.
     */
    int depth( const string_type& __input ) const {
      typedef std::string::const_iterator c_iter;
      std::stack<char> nr;
      
      int depth(0), balance(0), position(0);
      
      for( c_iter cit = __input.begin(); cit != __input.end();
           ++cit, ++position ) {
        switch( punctuation_.get_char_type( *cit ) ){
          case Left:
          {
            std::pair<char,bool> p = punctuation_.get_right(*cit);
            if( !p.second ) {
              throw exception_type( position, __input );
            }
            else {
              nr.push( p.first );
            }
            
            if( ++balance > depth ) ++depth;
            break;
          }
          case Right:
            if( *cit != nr.top() ) {
              throw exception_type( pos_, input_ );
            }
            nr.pop();
            
            --balance;
            break;
          default:
            break;
        }
        if( balance < 0 ) {
          throw exception_type( position, __input );
        }
      }
      if( balance != 0 ) {
        throw exception_type( position - 1, __input );
      }
      
      return depth;
    }
    
    
  private:
    /**
     *  @enum   char_type
     *  @brief  A private Enum for recording the type of a character found in
     *          an input string.
     *  @note   The possiblities are
     *          <ul>
     *            <li>\c Left Indicates left parenthesis
     *            <li>\c Right Indicates right parenthesis
     *            <li>\c Delimeter Indicates a delimeter
     *            <li>\c Normal Indicates none of the above
     *          </ul>
     */
    enum char_type {
      Left,
      Right,
      Delimiter,
      Normal
    };
    
    
    /**
     *  @brief   Processes each character of the input string
     *  @param   c The character to process
     */
    void sys_process_string( const char __c ){
      // Increment the position for error recording
      ++pos_;
      
      // Get the character type of the current and previous characters
      char_type ct = punctuation_.get_char_type( __c );
      char_type pt = (pos_==1 ? Normal : punctuation_.get_char_type( prev_ ));
      
      switch( ct ) {
        case Left:
        {
          // Increment the depth
          if( ++depth_ == 1 && pt == Delimiter ) {
            throw exception_type( pos_, input_ );
          }
          
          std::pair<char,bool> p = punctuation_.get_right(__c);
          if( !p.second ) {
            throw exception_type( pos_, input_ );
          }
          else {
            nextRight_.push( p.first );
          }
          
          if( depth_ != 1 ) current_.push_back( __c );
          prev_ = __c;
          return;
        }
        
        case Right:
        {
          // ",]" found: ERROR
          if( depth_ == 0 || pt == Delimiter || __c != nextRight_.top() ) {
            throw exception_type( pos_, input_ );
          }
          
          nextRight_.pop();
          
          if( --depth_ == 0 ) {
            sys_add_element( Right );
            return;
          }
          break;
        }
        
        case Delimiter:
        {
          // ",," found: ERROR
          if( pt == Delimiter ) {
            throw exception_type( pos_, input_ );
          }
          if( depth_ == 1 ) {
            sys_add_element( Delimiter );
            return;
          }
          break;
        }
        
        case Normal: { } // Nothing to do
      }
      
      if( depth_ > 0 ) current_.push_back( __c );
      else throw exception_type( pos_, input_ );
      prev_ = __c;
    }
    
    /**
     *  @brief Add the current string to the elements vector.
     */
    void sys_add_element( const char_type __ct ) {
      boost::trim( current_ );
      if( current_.size() != 0 ) {
        elements_.push_back( current_ );
        current_.clear();
      }
      else if ( __ct != Right ) {
        throw exception_type( pos_, input_ );
      }
    }
    
    /**
     *  @struct parenthesis
     *  @brief  A pair of left and right parentheses.
     */
    struct parenthesis {
      //! The object type
      typedef parenthesis self_type;
      
      //! The lefthand parenthesis
      char left;
      //! The righthand parenthesis
      char right;
      
      /**
       *  @brief  Determins if a character is a left or right parenthesis
       *  @param  c The character to check
       *  @return \c true if c is either \c left or \c right
       */
      bool contains( const char __c ) const {
        return bool( left == __c || right == __c );
      }
      
      /**
       *  @brief  Less than operator
       *  @param  lhs A parenthesis
       *  @param  rhs A parenthesis
       *  @return A strict ordering on \c lhs and \c rhs
       */
      friend bool operator< ( const self_type __lhs, const self_type __rhs )
      {
        return bool(    __lhs.left < __rhs.left
                     || (    __lhs.left == __rhs.left
                          && __lhs.right < __rhs.right )
                   );
      }
    };
    
    
    /**
     *  @struct punctuation
     *  @brief  The puctuation which defines the grammar
     */
    class punctuation {
      public:
            /**
             *  @brief  Default constructor
             *  @note   Default parenthesis and delimeter:
             *          <ul>
             *            <li>\b parenthesis "()", "[]", "{}" and "<>"
             *            <li>\b delimeters "," and " "
             *          </ul>
             */
        punctuation() {
          add_parenthesis('(',')');
          add_parenthesis('[',']');
          add_parenthesis('{','}');
          add_parenthesis('<','>');
          add_delimiter(',');
          add_delimiter(' ');
        }
        
        
        /**
         *  @brief  Add new parenthesis
         *  @param  left The lefthand parenthesis
         *  @param  right The rightand parenthesis
         *  @return \c true if the parenthesis was added successfully.
         *  @note   Both \c left and \c right must be different and not ' '
         */
        bool add_parenthesis( const char __left, const char __right ) {
          bool b =    __left != ' ' && __right != ' ' && __left != __right
                   && !contains(__left) && !contains(__right);
          
          parenthesis p = {__left,__right};
          if( b ) parentheses_.insert( p );
          return b;
        }
        
        /**
         *  @brief  Remove parenthesis from the punctuation
         *  @return \c true if the parenthesis was removed successfully.
         */
        bool remove_parenthesis( const char __left, const char __right ) {
          parenthesis p = {__left,__right};
          return parentheses_.erase( p );
        }
        
        /**
         *  @brief  Add new delimeter
         *  @param  delimeter The delimiter to add
         *  @return \c true if the delimiter was added successfully.
         */
        bool add_delimiter( const char __delimiter ) {
          bool b = !contains(__delimiter);
          if( b ) delimiters_.insert( __delimiter );
          return b;
        }
        
        /**
         *  @brief  Remove delimeter from the punctuation
         *  @param  delimeter The delimiter to remove
         *  @return \c true if the delimiter was removed successfully.
         */
        bool remove_delimiter( const char __delimiter ) {
          return delimiters_.erase( __delimiter );
        }
        
        /**
         *  @brief  Determines if a character is in the punctuation list.
         *  @param  c The character to check
         *  @return True if \c c is a left parenthesis, right parenthesis or a
         *          delimeter.
         */
        bool contains( const char __c ) const {
          typedef parentheses_type::const_iterator c_iter;
          
          bool b = is_delimiter(__c);
          
          if(!b) {
            for( c_iter it = parentheses_.begin(); it != parentheses_.end()
                 ; ++it )
            {
              if( it->contains(__c) ) {
                b = true;
                break;
              }
            }
          }
          return b;
        }
        
        /**
         *  @brief  Determines if a character is a left parenthesis.
         *  @param  left The character to check
         *  @return True if \c left is a left parenthesis
         */
        bool is_left( const char __left ) const {
          typedef parentheses_type::const_iterator c_iter;
          
          bool b = false;
          for( c_iter it = parentheses_.begin(); it != parentheses_.end()
               ; ++it )
          {
            if( it->left == __left ) {
              b = true;
              break;
            }
          }
          return b;
        }
        
        /**
         *  @brief  Determines if a character is a right parenthesis.
         *  @param  right The character to check
         *  @return True if \c right is a right parenthesis
         */
        bool is_right( const char __right ) const {
          typedef parentheses_type::const_iterator c_iter;
          
          bool b = false;
          for( c_iter it = parentheses_.begin(); it != parentheses_.end()
               ; ++it )
          {
            if( it->right == __right ) {
              b = true;
              break;
            }
          }
          return b;
        }
        
        /**
         *  @brief  Determines if a character is a delimiter.
         *  @param  delimiter The character to check
         *  @return True if \c delimiter is a delimiter
         */
        bool is_delimiter( const char __delimiter ) const {
          return bool( delimiters_.count( __delimiter ) );
        }
        
        
        /**
         *  @brief  Determines the character type of a character.
         *  @param  c The character to check
         *  @return The character type of \c c
         */
        char_type get_char_type( const char __c ) const {
          char_type ct = Normal;
          if ( is_delimiter(__c) ) ct = Delimiter;
          else if ( is_left(__c) )  ct = Left;
          else if ( is_right(__c) )  ct = Right;
          return ct;
        }
        
        
        /**
         *  @brief  Get the righthand parenthesis from the left.
         *  @param  left The character to get the righthand pair of.
         *  @return A \c std::pair the second coordinate is a boolean which
         *          determines whether \c left is actually a left parenthesis
         *          and the first is the right parenthesis if it exists.
         */
        std::pair<char,bool> get_right( const char __left ) const {
          typedef parentheses_type::const_iterator c_iter;
          
          bool b = false;
          char c = ' ';
          
          for( c_iter it = parentheses_.begin(); it != parentheses_.end()
               ; ++it )
          {
            if( it->left == __left ) {
              b = true;
              c = it->right;
              break;
            }
          }
          
          return std::make_pair(c,b);
        }
        
      private:
        //! The type of the parentheses
        typedef std::set<parenthesis> parentheses_type;
        //! The type of the delimiters
        typedef std::set<char> char_set_type;
        
        //! Set of parenthesis
        parentheses_type parentheses_;
        //! The set of delimeters
        char_set_type delimiters_;
    };
    
    
    //! The punctuation that defines the grammar.
    punctuation punctuation_;
    
    //! The input string.
    string_type input_;
    //! The depth of the current vector: depth = nested parenthesis.
    int depth_;
    //! Current processing position in the input string.
    int pos_;
    //! The previous character in \c input_
    char prev_;
    //! The next right parenthesis that is expected
    std::stack<char> nextRight_;
    //!  The current element of the output vector
    string_type current_;
    //! The output vector
    string_vector_type elements_;
};

}// io
}// utility
}// sg
#endif
