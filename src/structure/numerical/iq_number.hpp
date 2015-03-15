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
 *  @file     iq_number.hpp
 *  @brief    An inline header file for the \c iq_number class.
 *  @note     A container for representing imaginary quadratic integers.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-29
 *  @note     In a previous incarnation of this software I developed a fairly
 *            full blown imaginary quadratic number class which interacted with
 *            an imaginary quadratic field class. In some sense this previous
 *            class was written in the "correct way," because it used an
 *            integral basis to represent it's numbers. However, I found that
 *            this actually created unecessary complications and overhead for
 *            my needs because it was much more natural to write
 *            \f$ x = a + b\sqrt{g} \f$ (with the standard qualifaction about
 *            factors of a half). As such this class is looks like it has been
 *            written in a clumsy way, but it hasn't (although it is imporant
 *            to remember that if \f$ D \eqiv 1 \mod 4 \f$ then
 *            <tt>iq_number(1) == 0</tt> and <tt>iq_number(2) == 1</tt>.
 */


#ifndef _SG_IQ_NUMBER_H
#define _SG_IQ_NUMBER_H 1

// Global includes
#include <cstddef>
#include <complex>

#include <boost/tuple/tuple_comparison.hpp>
#include <boost/operators.hpp>
#include <boost/call_traits.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/mpl/if.hpp>

// Local includes
#include "structure/numerical/iq_field.hpp"
#include "utility/functors/stream_cast_tuple.hpp"
#include "utility/math/square.hpp"
#include "utility/math/is_even.hpp"


namespace sg
{
namespace structure
{
namespace numerical
{
// Forward Declaeffectns
template <class _Float = double, class _Integer = long, std::size_t _Id = 0>
    class iq_number;
template <class _Float = double, class _Integer = long, std::size_t _Id = 0>
    class iq_canonical_form;

/**
 *  @class    iq_canonical_form
 *  @brief    Represents an imaginary quadratic integer in canonical form.
 *  @param    _Float The floating point type defaults to \c double
 *  @param    _Integer The integer type defaults to \c long
 *  @param    _Id The multiple singleton identifier of the iq_field
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-05-04
 *  @note     The canonical form is \f$ a + b( \frac{D + \sqrt{D}}{2} ) \f$
 *            where \f$ D \f$ is the discriminant of the field. This is the
 *            'right' way to represent algebraic integers computationally (ie
 *            with respect to some integral basis) however it is probably not
 *            the 'right' way to work with them in this case. Having said this
 *            the canonical form is necessary for the construction of ideals,
 *            and it is for this reason that this interface is provided.
 */
template <class _Float, class _Integer, std::size_t _Id>
    class iq_canonical_form
  : private boost::equality_comparable< iq_canonical_form<_Float,_Integer,_Id> >
{
  //! Object type
  typedef iq_canonical_form<_Float,_Integer,_Id> self_type;
  //! The type of the field of fractions.
  typedef _Float float_type;
  //! The type of the ring of integers.
  typedef _Integer integer_type;
  //! The iq_number type
  typedef iq_number<float_type,integer_type,_Id> normal_type;
  //! Field type
  typedef iq_field<float_type,integer_type,_Id> field_type;
  
  
  private:
    //! Type of the underlying tuple
    typedef boost::tuple<integer_type,integer_type> tuple_type;
    //! Integer param_type
    typedef typename boost::call_traits<integer_type>::param_type param_type;
  
  
  public:
    /**
     *  @brief  Default constructor
     */
    iq_canonical_form() : number_(0,0) { }
    
    /**
     *  @brief  Alpha / Beta constructor
     *  @param  alpha The alpha component
     *  @param  beta The beta component
     */
    iq_canonical_form( param_type __alpha, param_type __beta )
  : number_( __alpha, __beta ) { }
    
    /**
     *  @brief  Construct from normal form
     *  @param  normal The normal form to construct the canonical form from.
     */
    iq_canonical_form( normal_type __normal ) { from_normal_form( __normal ); }
    
    /**
     *  @brief  Assign from normal form
     *  @param  normal The normal form to assign the canonical form from.
     */
    iq_canonical_form operator=( normal_type __normal ) {
      from_normal_form( __normal );
      
      return *this;
    }
    
    /**
     *  @brief  Get a constant reference to the alpha part
     *  @return A constant reference to the real part of the number
     */
    const integer_type& alpha( void ) const {return boost::get<0>( number_ );}
    
    /**
     *  @brief  Get a reference to the alpha part
     *  @return A reference to the real part of the number
     */
    integer_type& alpha( void ) {
      return const_cast<integer_type&>(static_cast<const self_type&>(*this).alpha());
    }
    
    
    /**
     *  @brief  Get a constant reference to the beta part
     *  @return A constant reference to the beta part of the number
     */
    const integer_type& beta( void ) const {return boost::get<1>( number_ );}
    
    /**
     *  @brief  Get a reference to the beta part
     *  @return A reference to the beta part of the number
     */
    integer_type& beta( void ) {
      return const_cast<integer_type&>(static_cast<const self_type&>(*this).beta());
    }
    
    
    /**
     *  @brief  Get a representation of the number in normal form.
     *  @return The normal representation of the number
     */
    normal_type to_normal_form( ) const {
      return normal_type( *this );
    }
    
    
    /**
     *  @brief  Construct the canonical form from its normal form.
     *  @param  normal The normal form to construct the canonical form from
     */
    void from_normal_form( const normal_type& __normal ) {
      *this = __normal.to_canonical_form();
    }
    
    
    /**
     *  @brief  Compares two \c iq_canonical_form objects for equality.
     *  @param  rhs An \c iq_canonical_form
     */
    bool operator==( const self_type& __rhs ) const {
      return bool( number_ == __rhs.number_ );
    }
    
    /**
     *  @brief  Defines the output stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::ostream& operator<<( std::ostream& __os, const self_type& __t ) {
      using std::ostream;
      typedef utility::functors::stream_cast_tuple<tuple_type,ostream> _sc;
      
      _sc( )( __t.number_, __os );
      return __os;
    }
    
    /**
     *  @brief  Defines the input stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::istream& operator>>( std::istream& __is, self_type& __t ) {
      using std::istream;
      typedef utility::functors::stream_cast_tuple<istream,tuple_type> _sc;
      
      _sc( )( __is, __t.number_ );
      return __is;
    }
    
    
  private:
    tuple_type number_;
};


/**
 *  @class    iq_number
 *  @brief    Represents an imaginary quadratic integer.
 *  @param    _Float The floating point type defaults to \c double
 *  @param    _Integer The integer type defaults to \c long
 *  @param    _Id The multiple singleton identifier of the iq_field
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-29
 */
template <class _Float, class _Integer, std::size_t _Id>
    class iq_number
  : private boost::totally_ordered< iq_number<_Float,_Integer,_Id>,
            boost::ring_operators< iq_number<_Float,_Integer,_Id> > >
{
  public:
    //! Object type
    typedef iq_number<_Float,_Integer,_Id> self_type;
    //! The type of the field of fractions.
    typedef _Float float_type;
    //! The type of the ring of integers.
    typedef _Integer integer_type;
    //! The complex type (for calculations)
    typedef std::complex<float_type> complex_type;
    //! Field type
    typedef iq_field<float_type,integer_type,_Id> field_type;
    //! The type of the canonical form
    typedef iq_canonical_form<_Float,_Integer,_Id> canonical_type;
    //! The _Id of this field
    static const std::size_t identifier = _Id;
    
    
  private:
    //! The tuple type which holds the real and imaginary parts
    typedef boost::tuple<integer_type,integer_type> tuple_type;
    //! The (integer) parameter type
    typedef typename boost::call_traits<integer_type>::param_type param_type;
    
  public:
    /**
     *  @brief  Default constructor
     */
    iq_number( )
  : number_(0,0) {  }
    
    
    /**
     *  @brief  Real and imaginary constructor
     *  @param  real The real part
     *  @param  imag The imaginary part
     */
    iq_number( param_type __re, param_type __im )
  : number_( __re, __im ) { }
    
    
    /**
     *  @brief  Canonical form constructor
     *  @param  canonical The number in canonical form
     */
    iq_number( const canonical_type& __canonical )
  : number_() { from_canonical_form(__canonical); }
    
    
    /**
     *  @brief  Assignment from canonical form
     *  @param  canonical The number in canonical form
     */
    iq_number operator=( const canonical_type& __canonical ) {
      from_canonical_form( __canonical );
      
      return *this;
    }
    
    
    /**
     *  @brief  Get a constant reference to the real part
     *  @return A constant reference to the real part of the number
     */
    const integer_type& real( void ) const { return boost::get<0>( number_ ); }
    
    /**
     *  @brief  Get a reference to the real part
     *  @return A reference to the real part of the number
     */
    integer_type& real( void ) {
      return const_cast<integer_type&>(static_cast<const self_type&>(*this).real());
    }
    
    
    /**
     *  @brief  Get a constant reference to the imag part
     *  @return A constant reference to the imag part of the number
     */
    const integer_type& imag( void ) const { return boost::get<1>( number_ ); }
    
    /**
     *  @brief  Get a reference to the imag part
     *  @return A reference to the imag part of the number
     */
    integer_type& imag( void ) {
      return const_cast<integer_type&>(static_cast<const self_type&>(*this).imag());
    }
    
    
    /**
     *  @brief  Set the real and imaginary parts
     *  @param  real The real part
     *  @param  imag The imaginary part
     */
    void set_reim( param_type real, param_type imag ) {
      this->real() = real;
      this->imag() = imag;
    }
    
    
    /**
     *  @brief  Get the complex representation of this number
     *  @return A complex representation of the number.
     */
    complex_type to_complex( void ) const {
      using boost::mpl::if_;
      using boost::is_convertible;
      
      typedef typename if_< is_convertible<integer_type,float_type>,
                            prv_fast_to_cx, prv_safe_to_cx
                          >::type _op;
      
      complex_type output( _op()( real(), imag() ) );
      
      // Account for the congruent to 1 (4) case.
      if( field_type::instance().is_congruent( ) ) output /= 2.0;
      return output;
    }
    
  private:
    /**
     *  @brief  Helper struct: prefer native compile-time conversion.
     */
    struct prv_fast_to_cx {
      complex_type operator()( param_type __re, param_type __im ) {
        return complex_type( __re,
                             __im * field_type::instance().sqrt_generator()
                           );
      }
    };
    
    /**
     *  @brief  Helper struct: degrade to runtime conversion.
     *  @note   If \c stream_cast has been optimised for the integer and float
     *          types, this should still be fast.
     */
    struct prv_safe_to_cx {
      complex_type operator()( param_type __re, param_type __im ) {
        typedef utility::functors::stream_cast<integer_type,float_type> _sc;
      
        float_type re, im;
        _sc()( __re, re );
        _sc()( __im, im );
      
        return complex_type( re, im * field_type::instance().sqrt_generator() );
      }
    };
    
    
  public:
    /**
     *  @brief  Return a copy of the conjugate
     *  @return A copy of the complex conjugate of the integer.
     */
    self_type conj( void ) const { return self_type( real(), -imag() ); }
    
    
    /**
     *  @brief  Return the norm of the number
     *  @return The norm of the iq_number
     */
    integer_type norm( void ) const {
      typedef utility::math::square<integer_type> _sq;
      
      integer_type n( _sq()( real() ) - _sq()( imag() )
                      * field_type::instance().generator()
                    );
      
       // Account for the congruent to 1 (4) case.
      if( field_type::instance().is_congruent() ) n /= 4;
      
      return n;
    }
    
    
    /**
     *  @brief  Define conversion to the complex type
     *  @return A complex representation of the number.
     */
    operator complex_type( ) const {
      return complex_type( to_complex() );
    }
    
    
    /**
     *  @brief  Get a representation of the number in canonical form.
     *  @return The canonical form of the number
     *  @note   The canonical form is \f$ a + b( \frac{D + \sqrt{D}}{2} ) \f$
     *          where \f$ D \f$ is the discriminant of the field. This is needed
     *          to compute ideals.
     */
    canonical_type to_canonical_form( ) const {
      canonical_type can(0,imag());
      if( field_type::instance().is_congruent() ) {
        can.alpha() = (real() - imag() * field_type::instance().generator())/2;
      }
      else {
        can.alpha() = real() - 2*(imag() * field_type::instance().generator());
      }
      return can;
    }
    
    
    /**
     *  @brief  Construct the number from its canonical form.
     *  @param  canonical The canonical form to construct the number from
     */
    void from_canonical_form( const canonical_type& __canonical ) {
      imag() = __canonical.beta();
      if( field_type::instance().is_congruent() ) {
        real() =   2*__canonical.alpha()
                 + field_type::instance().generator()*__canonical.beta();
      } else {
        real() =   __canonical.alpha()
                 + 2*field_type::instance().generator()*__canonical.beta();
      }
    }
    
    
    /**
     *  @brief  Compares two \c iq_number objects for equality.
     *  @param  rhs An \c iq_number
     */
    bool operator==( const self_type& __rhs ) const {
      return bool( number_ == __rhs.number_ );
    }
    
    
    /**
     *  @brief  Less than operator supplied to use object in sets.
     *  @param  rhs An \c iq_number
     */
    bool operator<( const self_type& __rhs ) const
    {
      return bool(      real() < __rhs.real()
                   || ( real() == __rhs.real() && imag() < __rhs.imag() )
                 );
    }
    
    
    /**
     *  @brief  Adds two \c iq_number objects.
     *  @param  rhs An \c iq_number
     */
    self_type operator-( void ) const {
      return self_type( -real(), -imag() );
    }
    
    
    /**
     *  @brief  Adds two \c iq_number objects.
     *  @param  rhs An \c iq_number
     */
    self_type& operator+=( const self_type& __rhs ) {
      real() += __rhs.real();
      imag() += __rhs.imag();
      return *this;
    }
    
    
    /**
     *  @brief  Subtracts two \c iq_number objects.
     *  @param  rhs An \c iq_number
     */
    self_type& operator-=( const self_type& __rhs ) {
      real() -= __rhs.real();
      imag() -= __rhs.imag();
      return *this;
    }
    
    
    /**
     *  @brief  Multiplies two \c iq_number objects.
     *  @param  rhs An \c iq_number
     */
    self_type& operator*=( const self_type& __rhs ) {
      integer_type rl = real();
      
      real() *= __rhs.real();
      real() += imag() * __rhs.imag() * field_type::instance().generator();
      
      imag() *= __rhs.real();
      imag() += rl * __rhs.imag();
      
      if( field_type::instance().is_congruent() ) {
        real() /= 2;
        imag() /= 2;
      }
      return *this;
    }
    
    /**
     *  @brief  Get a latex representation of the number
     *  @return A Latex representation of the number.
     */
    std::string tex( void ) const {
      typedef utility::math::is_even<integer_type> _iev;
      std::stringstream ss;
      
      if( field_type::instance().is_congruent() ) {
        bool isEven = _iev()(real()) && _iev()(imag());
        integer_type re = isEven ? real() / 2 : real();
        integer_type im = isEven ? imag() / 2 : imag();
        
        if (!isEven ) { ss << "\\frac{"; }
        ss << sys_process_tex (re, im);
        if (!isEven ) { ss << "}{2}"; }
        
      }
      else {
        ss << sys_process_tex (real(), imag());
      }
      
      return ss.str();
    }
    
    
    std::string sys_process_tex( integer_type __re, integer_type __im ) const {
      std::stringstream ss;
      if( __re != 0 ) {
        ss << __re;
      }
      
      if( __im != 0 ) {
        if( __re != 0 && __im > 0 ) {
          ss << "+";
        }
        if (__im == 1) {}
        else if (__im == -1) { ss << "-"; }
        else { ss << __im; }
        ss << field_type::instance().tex_generator();
      }
       
      if( __re == 0 && __im == 0 ) {
        ss << "0";
      }
      
      return ss.str();
    }
    
    /**
     *  @brief  Defines the output stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::ostream& operator<<( std::ostream& __os, const self_type& __t ) {
      using std::ostream;
      typedef utility::functors::stream_cast_tuple<tuple_type,ostream> _sc;
      
      _sc( )( __t.number_, __os );
      return __os;
    }
    
    /**
     *  @brief  Defines the input stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::istream& operator>>( std::istream& __is, self_type& __t ) {
      using std::istream;
      typedef utility::functors::stream_cast_tuple<istream,tuple_type> _sc;
      
      _sc( )( __is, __t.number_ );
      return __is;
    }
    
    
  private:
    //! The storage tuple which holds the real and imaginary parts
    tuple_type number_;
};

}// numerical
}// structure
}// sg
#endif
