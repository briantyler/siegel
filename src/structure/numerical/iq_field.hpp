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
 *  @file     iq_field.hpp
 *  @brief    An inline header file for iq_field.
 *  @note     Include this file to work with imaginary quadratic fields.
 *  @author   Brian Tyler
 *  @version  2.1
 *  @date     2008-02-04
 */

#ifndef _SG_IQ_FIELD_H
#define _SG_IQ_FIELD_H 1

// Global includes
#include <cstddef>
#include <iostream>
#include <string>

#include <boost/call_traits.hpp>

// Local includes
#include "utility/math/sqrt.hpp"
#include "utility/math/is_square.hpp"
#include "utility/functors/stream_cast.hpp"
#include "utility/io/string_parser.hpp"


namespace sg
{
namespace structure
{
namespace numerical
{
/**
 *  @class    iq_field
 *  @brief    Represents an imaginary quadratic field.
 *  @param    Float The real type (field of fractions) defaults to \c double
 *  @param    Integer The integer type (ring of integers) defaults to \c long
 *  @param    Id The multiple singleton identifier. The default value is \c 0.
 *  @note     This class uses the Multiple Myers Singleton paradigm (The Myers
 *            Singleton paradigm, but multiple non-compatible instances can be
 *            created via different Id parameters.)
 *  @author   Brian Tyler
 *  @version  2.1
 *  @date     2008-03-17
 */
template <class _Float = double, class _Integer = long, std::size_t _Id = 0>
    class iq_field
{
  public:
    //! Object type
    typedef iq_field<_Float,_Integer,_Id> self_type;
    //! The type of the field of fractions.
    typedef _Float float_type;
    //! The type of the ring of integers.
    typedef _Integer integer_type;
    //! The _Id of this singleton
    static const std::size_t identifier = _Id;
    
  private:
    //! The (integer) parameter type
    typedef typename boost::call_traits<integer_type>::param_type param_type;
    
    
  public:
    /**
     *  @brief  Provides access to the iq_field
     *  @return The iq_field instance.
     */
    static self_type& instance( void ) {
      static self_type iq;  /* Myers Singleton */
      return iq;
    }
    
    
    /**
     *  @brief  Sets the generator and initializes the numberfield for use.
     */
    void initialize( param_type __generator = -1 ){
      // Since the field is imaginary the generator is negative
      g_ = ( __generator < 0 ? static_cast<integer_type>(__generator)
                             : static_cast<integer_type>(-__generator)
           );
      
      /* The structure of the numberfield is different if the generator is
       * congruent to 1 mod 4, than if it is congruent to 2 or 3. Is is often
       * useful to know which situation we are in.
       */
      integer_type congruence( g_ % 4 );
      isCongruent_ = static_cast<bool>( congruence == 1 || congruence == -3 );
      
      // Square root of the generator
      gSqrt_ = utility::math::sqrt<float_type,integer_type>()( -g_ );
      
      /* Independent of the structure of the numberfield, it is generated as
       * a + b*w, where w = (D + sqrt(D))/2 and D is the discriminant.
       */
      D_ = ( isCongruent_ ? g_: 4 * g_ );
      
      // This constant is helps in computing ideals.
      mfactor_ = (D_ - D_*D_)/4;
      
      // Calculate the class number
      this->calculate_h();
      
      // Determine if the field is a ufd
      isUFD_ = ( h_ == 1 );
    }
    
    
    /**
     *  @brief  Returns the generator of the numberfield.
     */
    const integer_type& generator( void ) const {
      return g_;
    }
    
    
    /**
     *  @brief  Returns the square root ot the generator of the numberfield.
     */
    const float_type& sqrt_generator( void ) const {
      return gSqrt_;
    }
    
    
    /**
     *  @brief  Returns the class number of the numberfield.
     *  @note   Precomputed during initialization
     */
    const integer_type& class_number( void ) const {
      return h_;
    }
    
    
    /**
     *  @brief  Returns the discriminant of the numberfield.
     *  @note   Precomputed during initialization
     */
    const integer_type& discriminant( void ) const {
      return D_;
    }
    
    
    /**
     *  @brief  Returns the multiplication factor of the ring of integers
     *  @note   Precomputed during initialization \$ = \frac{D - D^2}{2} \$
     */
    const integer_type& mfactor( void ) const {
      return mfactor_;
    }
    
    
    /**
     *  @brief  Returns true if the generator is congruent to 1 mod 4.
     *  @note   Precomputed during initialization
     */
    bool is_congruent( void ) const {
      return isCongruent_;
    }
    
    
    /**
     *  @brief  Returns true if the ring of integers has unique factorisation.
     *  @note   Precomputed during initialization
     */
    const bool is_ufd( void ) const {
      return isUFD_;
    }
    
    
    /**
     *  @brief  Get a latex representation of the generator
     *  @return A Latex representation of the generator.
     */
    std::string tex_generator( void ) const {
      std::stringstream ss;
      if( g_ == -1 ) {
        ss << "\\imath{}";
      } 
      else if (is_congruent()) {
        ss << "\\sqrt{\\frac{1+" << g_ << "}{2}}";
      }
      else{
        ss << "\\sqrt{" << g_ << "}";
      }
      return ss.str();
    }
    
    /**
     *  @brief  Defines the output stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::ostream& operator<<( std::ostream& __os, const self_type& __t ){
      __os << '[' << __t.g_ << ',' << __t.D_ << ',' << __t.h_ << ']';
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
      
      string str;
      __is >> str;
      
      typename string_parser::string_vector_type vs = string_parser( )( str );
      
      // Don't require the discriminant or class number to initialize field.
      if( vs.size() == 0 ) throw( exceptions::bad_vector_input(1,str) );
      
      integer_type generator;
      stream_cast<string,integer_type>( )( vs.at(0), generator );
      __t.initialize( generator );
      
      return __is;
    }
    
    
  private:
    /**
     *  @brief  Computes the class number of the field.
     *  @note Method taken from Henri Cohen: A Course in Computational
     *  Number Theory; Algorithm 5.3.5 (h(D) counting reduced forms).
     * 
     *  Definition: A positive definite quadratic form (a, b, c) of
     *  discriminant D (\f$ = b^2 - 4ac \f$) is said to be reduced if
     *  \f$ |b| <= a <= c \f$ and if, in addition, when one of the two
     *  inequalities is an equality, then b >= 0.
     * 
     *  Proposition 5.3.3 In every class of positive definite quadratic forms
     *  of discriminant \f$ D < 0 \f$ there exists exactly one reduced form.
     *  In particular \f$ h(D) \f$ is equal to the number of primitive
     *  reduced forms of discriminant \f$ D \f$.
     */
    void calculate_h( void ) {
      //  The steps denote the steps in Cohen's algorithm.
      
      /* ### Step 1 ### */
      h_ = 1;

      integer_type b( ( D_ % 2 == 0 ) ? 0 : 1 ), a, aSq, q;
      float_type bound(utility::math::sqrt<float_type,integer_type>( )( -D_ ));
      bound /= 1.732050807568877; // ~= sqrt(3)
      
      while( b <= bound ) {
        /* ### Step 2 ### */
        q = ( b*b - D_ ) / 4;                                     /* = a*c */
        a = ( b == 0 ) ? 1 : b;            /* due to |b| <= a, then a >= 1 */
  
        /* ### Step 3 ### */
        aSq = a*a;
        while( aSq <= q ) {
          if( a != 1 ) {
            if( ( q % a ) == 0 ) {      /* c is integral; found a solution */
             /* Note that the condition comes from the fact the if any of
              * these is satisfied then b is of fixed sign and as such there
              * is one solution rather than two.
              */
              h_ += ( a == b || aSq == q || b == 0 ) ? 1 : 2;
            }
          }
    
          /* ### Step 4 ### */
          aSq += 2*a + 1; // == (a+1)^2
          ++a;
        }
  
        /* ### Step 5 ### */
        b += 2;
      }
    }
    
    
    /**
     *  @brief  Default constructor
     *  @note   Private enforces singleton paradigm
     */
    iq_field( )
  : g_(-1), gSqrt_(1), D_(-4), h_(1), isCongruent_(false), isUFD_(true) { }
    
    /**
     *  @brief  Copy constructor
     *  @note   Private enforces singleton paradigm
     */
    iq_field( const self_type& __that );
    
    /**
     *  @brief  Assignment
     *  @note   Private enforces singleton paradigm
     */
    self_type& operator=( const self_type& __that );
    
    /**
     *  @brief  Destructor
     *  @note   Private enforces singleton paradigm
     */
    ~iq_field( ){ }
    
    
    //! Generator of the field
    integer_type g_;
    
    //! Square root of the generator
    float_type gSqrt_;
    
    //! Discriminant of the field
    integer_type D_;
    
    //! Multiplication factor of the ring of integers
    integer_type mfactor_;
    
    //! Class number (h) of the field
    integer_type h_;
    
    //! True if \f$ g \equiv 1 \mod 4 \f$
    bool isCongruent_;
    
    //! True if \f$ h = 1 \f$
    bool isUFD_;
};

}// numerical
}// structure
}// sg
#endif
