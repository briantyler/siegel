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
 *  @file     iq_ideal.hpp
 *  @brief    An inline header file for the \c iq_ideal class.
 *  @note     A class representing integral ideals in an imaginary quadratic
 *            field
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-05-02
 */


#ifndef _SG_IQ_IDEAL_H
#define _SG_IQ_IDEAL_H 1

// Global includes
#include <cassert>

#include <boost/call_traits.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/operators.hpp>
#include <boost/math/common_factor_rt.hpp>

// Local includes
#include "utility/math/gcd.hpp"
#include "utility/math/sgn.hpp"
#include "utility/functors/stream_cast_tuple.hpp"
#include "structure/numerical/iq_number.hpp"


namespace sg
{
namespace structure
{
namespace numerical
{
/**
 *  @struct   quadratic_form
 *  @brief    The reduced quadratic form representation of the ideal.
 *  @param    _Integer The integer type defaults to \c long
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-29
 *  @note     Equivalence of reduced quadratic forms indicates equivalence
 *            in the class group.
 */
template <class _Integer = long> struct quadratic_form
  : private boost::equality_comparable< quadratic_form<_Integer> >
{
  //! The object type
  typedef quadratic_form<_Integer> self_type;
  //! The integer type
  typedef _Integer integer_type;
  //! Integer param_type
  typedef typename boost::call_traits<integer_type>::param_type param_type;
  
  
  /**
   *  @brief  Default constructor
   */
  quadratic_form() : cXX(0), cXY(0), cYY(0) { }
  
  /**
   *  @brief  Component constructor
   *  @param  cXX The coefficient of \f$ X^2 \f$.
   *  @param  cXY The coefficient of \f$ XY \f$.
   *  @param  cYY The coefficient of \f$ Y^2 \f$.
   *  @note   This constructor is necessary because \c equality_comparable has
   *          a custom constructor, so a braces initializer is not possible.
   */
  quadratic_form( param_type __cXX, param_type __cXY, param_type __cYY )
  : cXX(__cXX), cXY(__cXY), cYY(__cYY) { }
  
  
  /**
   *  @brief  Defines the equality operator.
   *  @param  lhs A \c quadratic_form
   *  @param  rhs A \c quadratic_form
   *  @note   \c lhs and \c rhs are equal if they are equal componentwise.
   */
  bool operator==( const self_type& __rhs ) const {
    return bool( cXX == __rhs.cXX && cXY == __rhs.cXY && cYY == __rhs.cYY );
  }
  
  
  /**
   *  @brief  Defines the output stream operator
   *  @param  os An ostream.
   *  @param  t The object to stream out.
   */
  friend std::ostream& operator<< ( std::ostream& __os, const self_type& __t ) {
    __os << "[" << __t.cXX << "," << __t.cXY << "," << __t.cYY << "]";
    return __os;
  }
  
   //! The coefficient of \f$ X^2 \f$
  integer_type cXX;
  //! The coefficient of \f$ XY \f$
  integer_type cXY;
  //! The coefficient of \f$ Y^2 \f$
  integer_type cYY;
};


/**
 *  @class    iq_ideal
 *  @brief    Represents an integral ideal in an imaginary quadratic field.
 *  @param    _Float The floating point type defaults to \c double
 *  @param    _Integer The integer type defaults to \c long
 *  @param    _Id The multiple singleton identifier of the iq_field
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-29
 *  @note     See Cohen - A Course in Computational Algebraic Number Theory
 *            Section 5.2 for some limited information about the ideals of
 *            imaginary quadratic rings. Cohen does not seem to give an
 *            algorithm for adding ideals, or constructing ideals in HNF from a
 *            principal ideal generator, so I have had to write my own.
 */
template <class _Float = double, class _Integer = long, std::size_t _Id = 0>
    class iq_ideal
  : private boost::equality_comparable< iq_ideal<_Float,_Integer,_Id>,
            boost::addable< iq_ideal<_Float,_Integer,_Id>,
            boost::multipliable< iq_ideal<_Float,_Integer,_Id>
            > > >
{
  public:
    //! Object type
    typedef iq_ideal<_Float,_Integer,_Id> self_type;
    //! The type of the field of fractions.
    typedef _Float float_type;
    //! The type of the ring of integers.
    typedef _Integer integer_type;
    //! The imaginary quadratic number type
    typedef iq_number<_Float,_Integer,_Id> iq_number_type;
    //! The canonical iq number type
    typedef typename iq_number_type::canonical_type canonical_type;
    //! The imaginary quadratic field type
    typedef typename iq_number_type::field_type field_type;
    //! The quadratic form type
    typedef quadratic_form<integer_type> form_type;
    
  private:
    //! The tuple type which represents the ideal.
    typedef boost::tuple<integer_type,integer_type,integer_type> tuple_type;
    //! Sign functor
    typedef utility::math::sgn<integer_type> _sgn_fnc;
    //! Greatest common divisor functor
    typedef utility::math::gcd<integer_type> _gcd_fnc;
    //! GCD solution
    typedef typename _gcd_fnc::solution_type _gcdsol_typ;
    
  public:
    /**
     *  @brief  Default constructor
     *  @note   Constructs the zero ideal.
     */
    iq_ideal( )
  : ideal_(0,0,0), norm_(0), computeForm_(true) {  }
    
    
    /**
     *  @brief  Principal ideal constructor
     *  @param  iqNumber the number from which to construct the principal
     *          ideal.
     */
    iq_ideal( const iq_number_type& __iqNumber )
  : ideal_(0,0,0), norm_(0), computeForm_(true) { make_principal(__iqNumber); }
    
    
    /**
     *  @brief  Get the \c a factor of the ideal
     *  @return \c a
     *  @note   The ideal is internally represented as aZ + (b + cw)Z, where
     *          \$ w = \frac{D + \sqrt{D} }{2} \$ and \$ D \$ is the
     *          discriminant of the field.
     */
    const integer_type& a( void ) const { return boost::get<0>( ideal_ ); }
    
    /**
     *  @brief  Get the \c b factor of the ideal
     *  @return \c b
     */
    const integer_type& b( void ) const { return boost::get<1>( ideal_ ); }
    
    /**
     *  @brief  Get the \c c factor of the ideal
     *  @return \c c
     */
    const integer_type& c( void ) const { return boost::get<2>( ideal_ ); }
    
    
    /**
     *  @brief  Returns the canonical representation of the first generator of
     *          the ideal.
     *  @return The real generator of the ideal; which is the smallest, totally
     *          real, strictly positive element in the ideal.
     *  @note   This is primarily for testing purposes.
     */
    canonical_type first_generator_can( void ) const {
      return canonical_type( a(),0 );
    }
    
    /**
     *  @brief  Returns the iq_number representation of the first generator of
     *          the ideal.
     *  @note   This is primarily for testing purposes.
     */
    iq_number_type first_generator_num( void ) const {
      return iq_number_type( first_generator_can() );
    }
    
    
    /**
     *  @brief  Returns the canonical representation of the second generator of
     *          the ideal.
     *  @return The complex generator of the ideal; which is the element that
     *          has the smallest, strictly positive, beta part in the ideal
     *          such that the alpha part is non-negative and smaller than
     *          the first generator.
     *  @note   This is primarily for testing purposes.
     */
    canonical_type second_generator_can( void ) const {
      return canonical_type(b(),c());
    }
    
    /**
     *  @brief  Returns the iq_number representation of the second generator of
     *          the ideal.
     *  @note   This is primarily for testing purposes.
     */
    iq_number_type second_generator_num( void ) const {
      iq_number_type n( second_generator_can() );
      n.real() %= a();
      if( n.real() < 0 ) n.real() += a();
      
      return n;
    }
    
    
    /**
     *  @brief  Make the ideal into a principal ideal
     *  @param  iqNumber the number from which to make the principal ideal.
     */
    void make_principal( const iq_number_type& __iqNumber ) {
      if( __iqNumber.real() == 0 && __iqNumber.imag() == 0 ) {
        // This is the zero ideal
        make_zero_ideal();
        return;
      }
      // Norm of a principal ideal is equal to the norm of the element that
      // generates it.
      norm_ = __iqNumber.norm();
      
      // Let alpha + beta * omega be the canonical form of the generator.
      canonical_type can( __iqNumber );
      
      // a is the smallest positive integer in the ideal, so;
      // a = norm / gcd(alpha, beta)
      // and norm = a * c, thus
      // c = gcd(alpha, beta)
      _gcdsol_typ sol = _gcd_fnc()( can.alpha(), can.beta() );
      sys_c() = sol.gcd;
      sys_a() = norm_ / sys_c();
      
      
      // b + c*omega =   alpha*x + beta*y * (D - D^2)/4
      //               + ((alpha + beta*D)y + beta*x)omega
      // where x and y are integers,
      // since c | alpha, beta, then x and y can be found using the extended
      // Euclidean algorithm.
      sol = _gcd_fnc()( can.beta(),
                        can.alpha()
                        + field_type::instance().discriminant()*can.beta()
                      );
      
      // c divides both a and b
      assert( sys_c() == sol.gcd );
      sys_b() =   can.alpha()*sol.a
                + field_type::instance().mfactor()*can.beta()*sol.b;
      
      // 0 <= b < a, so reduce b modulo a and make the remainder positive
      sys_b() %= sys_a();
      if( sys_b() < 0 ) sys_b() += sys_a();
      
      computeForm_ = true;
    }
    
    
    /**
     *  @brief  Get the norm of the ideal
     *  @return The norm of the ideal
     */
    const integer_type& norm( void ) const { return norm_; }
    
    
    /**
     *  @brief  Get the reduced quadratic form representation of the ideal.
     *  @return The reduced quadratic form representation of the ideal.
     */
    const form_type& form( void ) const {
      // Do we need to recompute the form?
      sys_compute_form();
      
      return form_;
    }
    
    
    /**
     *  @brief  Makes the ideal into the maximal order (ring of integers)
     *  @note   It might be very slightly more efficient to check first to see
     *          if the ideal is a maximal order before calling this, but it is
     *          probably negligible.
     */
    void make_maximal_order( void ) {
      norm_ = sys_a() = sys_c() = 1;
      sys_b() = 0;
      computeForm_ = true;
    }
    
    
    /**
     *  @brief  Determines if the ideal is in fact the entire ring of integers
     *  @return True if the ideal generates the entire ring of integers.
     */
    bool is_maximal_order( void ) const {
      // If the norm is a unit then the ideal generates the ring
      return bool( norm_ == 1 );
    }
    
    
    /**
     *  @brief  Makes the ideal into the zero ideal
     */
    void make_zero_ideal( void ) {
      norm_ = sys_a() = sys_b() = sys_c() = 0;
      computeForm_ = true;
    }
    
    /**
     *  @brief  Determines if the ideal is the zero ideal
     *  @return True if the ideal is the zero ideal.
     */
    bool is_zero_ideal( void ) const {
      // Since 0 <= b < a this condition implies that b == 0.
      return bool( a() == 0 && c() == 0 );
    }
    
    
     /**
     *  @brief  Determines if the ideal is a principal ideal
     *  @return True if the ideal is a principal ideal.
      */
    bool is_principal_ideal( void ) const {
      // Note: use 2 instead of 1, since this deals with both congruence
      // classes.
      static integer_type generator = field_type::instance().generator();
      static self_type principal( iq_number_type( 2, 0 ) );
      
      // Account for the case where the generator changes
      if( generator != field_type::instance().generator() ) {
        generator = field_type::instance().generator();
        principal.make_maximal_order();
      }
      
      return bool( same_class( principal ) );
    }
    
    
    /**
     *  @brief  Compares the ideal classes of two ideals for equality.
     *  @param  rhs An \c iq_ideal
     *  @return True if *\c this and \c rhs are in the same ideal class
     */
    bool same_class( const self_type& __rhs ) const {
      return bool( form() == __rhs.form() );
    }
    
    
    /**
     *  @brief  Compares two \c iq_ideal objects for equality.
     *  @param  rhs An \c iq_ideal
     *  @return True if the ideals are equal
     */
    bool operator==( const self_type& __rhs ) const {
      return bool( ideal_ == __rhs.ideal_ );
    }
    
    
    /**
     *  @brief  Adds two \c iq_ideal objects.
     *  @param  rhs An \c iq_ideal
     *  @return The sum of \c this and \c rhs
     */
    self_type& operator+=( const self_type& __rhs ) {
      typedef boost::math::lcm_evaluator<integer_type> _lcm;
      
      // In the new ideal, the smallest integer will be the gcd of a(), rhs.a()
      // and any integer made by getting rid of the omega factors
      // (== (D + sqrt(D))/2 ) in the second generators.
      // ie: if this = a_1Z + (b_1 + c_1w)Z
      //     and rhs = a_2Z + (b_2 + c_2w)Z
      //     Then a_3 = gcd(a_1, a_2, b_1 * c_2 - b_2*c_1 )
      
      // Adding the zero ideal has no effect and adding the maximal order
      // gives the maximal order
      if( is_maximal_order() || __rhs.is_zero_ideal()  ) {
        // nothing to do
      }
      else if( is_zero_ideal() || __rhs.is_maximal_order() ) {
        *this = __rhs;
      }
      else{
        // Compute the initial a() value.
        integer_type lcm = _lcm()( c(), __rhs.c() );
        sys_a() = _gcd_fnc()( a(), __rhs.a() ).gcd;
        sys_a() = _gcd_fnc()( a(), b()*(lcm/c()) - __rhs.b()*(lcm/__rhs.c()) ).gcd;
        
        // The new c() value is given by gcd( c(), rhs.c() )
        _gcdsol_typ sol = _gcd_fnc()( c(), __rhs.c() );
        sys_c() = sol.gcd;
        
        // In order to construct c() it is necessary to employ Euclid's algorithm
        // so the value of b() is totally determined.
        sys_b() = sol.a * b() + sol.b * __rhs.b();
        
        // Computing the new value for a() means finding the smallest integer
        // in the new ideal.
        canonical_type can( b(), c() );
        
        // Refine a() by taking into accout the smallest integer in the ideal
        // generated by the new complex generator.
        sys_a() =
            _gcd_fnc()( a(),
                        (   iq_number_type( can ).norm()
                          / _gcd_fnc()( can.alpha(), can.beta() ).gcd
                        )
                      ).gcd;
        
        // Refine b();
        sys_b() %=  a();
        if( sys_b() < 0 ) sys_b() +=  a();
        
        // The norm of the ideal is just the determinant of the matrix
        //   [a b]
        //   [0 c]
        norm_ = a() * c();
        
        computeForm_ = true;
      }
      
      return *this;
    }
    
    
    /**
     *  @brief  Mutliplies two \c iq_ideal objects.
     *  @param  rhs An \c iq_ideal
     *  @return The product of \c *this and \c rhs
     *  @note   Even for ideals generated by small integers, the products get
     *          so large that even with 64-bit long integers the calculation
     *          overflows very easily; an example is with I generated by
     *          [34,-93] and J generated by [76,-82] in the field of
     *          discriminant -66*4, then the calculation of the b() coordinate
     *          of I*J overflows. As such it is necessary to use arbitrary
     *          precision integers (ie \c gmp++::mpz_class integers) for any
     *          kind of stability in this algorithm. Luckily I don't require
     *          ideal multiplication anywhere in the sg program so I don't need
     *          to slow things down with gmp integers.
     */
    self_type& operator*=( const self_type& __rhs ) {
      // This returns the product of *this and rhs.
      // This is the ideal generated by:
      //  #  a1.a2
      //  #  a1(b2 + c2w) == a1.b2 + a1.c2w
      //  #  a2(b1 + c1w) == a2.b1 + a2.c1w
      //  #  (b1 + c1w)(b2 + c2w) ==   (b1.b2 + c1.c2(D-D^2)/2)
      //                             + (c1.b2 + c2.b1 + D.c1.c2)w
      //  c3 is the smallest coefficient of w in the ideal so:
      //  c3 = gcd(a1.c2, a2.c1, c1.b2 + c2.b1 + D.c1.c2)
      //  
      //  The coefficients resulting from the extended eucliean algorithm
      //  now determine b3 up to reduction by a3.
      //  
      //  Since a3 * c3 = Norm(*this) * Norm(rhs) it follows that
      //  a3 = (Norm(*this) * Norm(rhs)) / c3
      //  
      //  And finally put b3 %= a3 (where b3 is the positive remainder)
      //  
      //  Note: I have kept the temporaries to make the code more readable as
      //  it was getting really hard to work out what was going on.
      
      // Multiplication by the zero ideal gives the zero ideal
      // Multiplication by the maximal order gives the original ideal
      if( is_zero_ideal() || __rhs.is_maximal_order() ) {
        // nothing to do
      }
      else if( is_maximal_order() || __rhs.is_zero_ideal() ) {
        *this = __rhs;
      }
      else {
        // Now it is necessary to compute the product since both ideals are
        // non-zero.
        integer_type gw1 = a()*__rhs.c();
        integer_type gw2 = __rhs.a()*c();
        integer_type gw3 = c()*__rhs.b() + __rhs.c()*b()
                           + field_type::instance().discriminant()*__rhs.c()*c();
        
        _gcdsol_typ sol = _gcd_fnc()( gw1, gw2 );
        integer_type v1 = sol.a;
        integer_type v2 = sol.b;
        
        sol = _gcd_fnc()( sol.gcd, gw3 );
        v1 *= sol.a;
        v2 *= sol.a;
        
        integer_type gu1 = a()*__rhs.b();
        integer_type gu2 = __rhs.a()*b();
        integer_type gu3 = b()*__rhs.b()
                           + field_type::instance().mfactor()*__rhs.c()*c();
        
        sys_c() = sol.gcd;
        sys_b() = v1*gu1 + v2*gu2 + sol.b*gu3;
        
        norm_ *= __rhs.norm_;
        sys_a() = norm_ / c();
        
        sys_b() %= a();
        if( b() < 0 ) sys_b() += a();
        
        computeForm_ = true;
      }
      
      return *this;
    }
    
    
    /**
     *  @brief  Defines the output stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::ostream& operator<<( std::ostream& __os, const self_type& __t ) {
      typedef utility::functors::stream_cast_tuple<tuple_type,std::ostream> sc;
      
      sc( )( __t.ideal_, __os );
      return __os;
    }
    
    /**
     *  @brief  Defines the input stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::istream& operator>>( std::istream& __is, self_type& __t ) {
      typedef utility::functors::stream_cast_tuple<std::istream,tuple_type> sc;
      
      sc( )( __is, __t.ideal_ );
      
      __t.norm_ = __t.a() * __t.c();
      __t.computeForm_ = true;
      
      return __is;
    }
    
    
  private:
    /**
     *  @brief  Computes the unique reduced quadratic form associated to the
     *          ideal. This quadratic form determines the ideal class.
     */
    void sys_compute_form( void ) const {
      if( computeForm_ ) {
        if( is_zero_ideal() ) {
          form_.cXX = form_.cXY = form_.cYY = 0;
        }
        else {
          // Construction
          sys_construct_form();
          
          // Reduction
          sys_reduce_form();
          
          // At debug time ensure that the quadratic form has the same
          // discriminant as the number field
          assert(    field_type::instance().discriminant()
                  == (form_.cXY * form_.cXY) - (4 * form_.cXX * form_.cYY)
                );
        }
        computeForm_ = false;
      }
    }
    
    /**
     *  @brief  Constructs the non-reduced quadratic form.
     */
    void sys_construct_form( void ) const {
      // Construction
      iq_number_type w1( canonical_type(a(), 0) ), w2( canonical_type(b(), c()) );
      
      assert( ( w1.norm() % norm_ == 0 ) && ( w2.norm() % norm_ == 0 ) );
      
      form_.cXX = w1.norm() / norm_;
      form_.cYY = w2.norm() / norm_;
      
      form_.cXY = - iq_number_type(w1 * w2.conj()).real();
      if( !field_type::instance().is_congruent() ) form_.cXY *= 2;
      
      // Ensure that the divisions are exact at debug time
      assert( form_.cXY % norm_ == 0 );
      
      form_.cXY /= norm_;
    }
    
    /**
     *  @brief  Reduces the quadratic form to its unique representative.
     */
    void sys_reduce_form( void ) const {
      
      if( -form_.cXX < form_.cXY && form_.cXY <= form_.cXX  ) {
        if( !sys_normalize_form() ) return;
      }
      
      // Euclidean step
      while( true ) {
        // Ensure division by zero never happens at debug time
        assert( ( form_.cXX != 0 ) );
        
        integer_type q = form_.cXY / (2*form_.cXX);
        integer_type r = form_.cXY % (2*form_.cXX);
        
        if( r <= -form_.cXX ) {
          r += 2*form_.cXX;
          --q;
        }
        else if( r > form_.cXX ) {
          r -= 2*form_.cXX;
          ++q;
        }
        
        form_.cYY -= ( (form_.cXY + r) * q )/2;
        form_.cXY = r;
        
        if( !sys_normalize_form() ) break;
      }
    }
    
    /**
     *  @brief  Puts the quadratic form in its normal form.
     */
    bool sys_normalize_form( void ) const {
      bool output = true;
      if( form_.cXX > form_.cYY ) {
        
        form_.cXY = -form_.cXY;
          // Swap XX and YY
        integer_type tmp(form_.cXX);
        form_.cXX = form_.cYY;
        form_.cYY = tmp;
      }
      else if (form_.cXX == form_.cYY) {
        if( form_.cXY < 0 ) form_.cXY = -form_.cXY;
        output = false;
      }
      else {
        output = false;
      }
      
      return output;
    }
    
    /**
     *  @brief  Get a writeable reference to a
     */
    integer_type& sys_a() {
      return const_cast<integer_type&>(static_cast<const self_type&>(*this).a());
    }
    
    /**
     *  @brief  Get a writeable reference to b
     */
    integer_type& sys_b() {
      return const_cast<integer_type&>(static_cast<const self_type&>(*this).b());
    }
    
    /**
     *  @brief  Get a writeable reference to c
     */
    integer_type& sys_c() {
      return const_cast<integer_type&>(static_cast<const self_type&>(*this).c());
    }
    
    //! The storage tuple which holds the ideal
    tuple_type ideal_;
    //! The norm of the ideal
    integer_type norm_;
    // Since the quadratic form is only needed to determine the class of the
    // ideal there is no need to compute it after every time the ideal changes,
    // so it makes sense to make it mutable.
    //! The associated quadratic form
    mutable form_type form_;
    //! Should we compute the form?
    mutable bool computeForm_;
};

}// numerical
}// structure
}// sg

#endif
