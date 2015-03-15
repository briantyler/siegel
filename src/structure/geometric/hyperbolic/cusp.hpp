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
 *  @file     cusp.hpp
 *  @brief    An inline header file for \c cusp.
 *  @note     Include this file to create cusps.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-01
 */

#ifndef _SG_CUSP_H
#define _SG_CUSP_H 1

// Global includes
#include <cassert>
#include <algorithm>
#include <iterator>

#include <boost/call_traits.hpp>
#include <boost/operators.hpp>

// Local includes
#include "structure/numerical/iq_ideal.hpp"
#include "structure/numerical/accessors/iq_number_accessor.hpp"
#include "structure/geometric/euclidean/accessors/real_accessor.hpp"
#include "structure/geometric/detail/heisenberg_structure.hpp"
#include "structure/geometric/hyperbolic/hyperbolic_point.hpp"
#include "utility/math/sqrt.hpp"
#include "utility/math/is_equal_cx.hpp"
#include "utility/functors/stream_cast.hpp"
#include "utility/io/string_parser.hpp"


namespace sg
{
namespace structure
{
namespace geometric
{
namespace hyperbolic
{
/**
 *  @class    cusp
 *  @brief    A class representing a cusp on the boundary of hyperbolic space
 *  @author   Brian Tyler
 *  @version  1.2
 *  @date     2008-04-01
 *  @param    N The dimension of the point
 *  @param    _Float A floating point type; defaults to \c double
 *  @param    _Integer An integral type; defaults to \c long
 *  @param    _Id The identifier of the field
 */
  template < std::size_t N, class _Float = double, class _Integer = long,
             std::size_t _Id = 0>
    class cusp
  : public detail::heisenberg_structure<
                 N,
                 numerical::accessor::iq_number_accessor<_Float,_Integer,_Id>,
                 euclidean::accessor::real_accessor<_Integer>
           >
   , private boost::equality_comparable< cusp<N,_Float,_Integer,_Id>,
               boost::equality_comparable< cusp<N,_Float,_Integer,_Id>,
                                           hyperbolic_point<N,_Float>,
                 boost::less_than_comparable< cusp<N,_Float,_Integer,_Id>
             > > >
{
  public:
    //! Object type
    typedef cusp<N,_Float,_Integer,_Id> self_type;
    //! Float type
    typedef _Float float_type;
    //! Integer type
    typedef _Integer integer_type;
    //! The IQ Ideal type
    typedef numerical::iq_ideal<_Float,_Integer,_Id> iq_ideal_type;
    //! IQ Number type
    typedef typename iq_ideal_type::iq_number_type iq_number_type;
    //! IQ Field type
    typedef typename iq_number_type::field_type field_type;
    //! Complex type
    typedef typename iq_number_type::complex_type complex_type;
    //! Point type
    typedef hyperbolic_point<N,float_type> point_type;
    
    //! The _Id of this field
    static const std::size_t identifier = _Id;
    
    
  private:
    //! Integer parameter type
    typedef typename boost::call_traits<integer_type>::param_type param_type;
    
    
  public:
    /**
     *  @brief  Default constructor
     */
    cusp()
  : dilation_(1),
    rotation_( ( field_type::instance().is_congruent() ? 2 : 1 ), 0 ),
    point_()
    { }
    
    
    /**
     *  @brief  Get a constant reference to the dilation factor
     *  @return A constant reference to the dilation factor
     */
    const integer_type& dilation( void ) const { return dilation_; }
    
    /**
     *  @brief  Get a reference to the dilation factor
     *  @return A reference to the dilation factor
     */
    integer_type& dilation( void ) { return dilation_; }
    
    
    /**
     *  @brief  Get a constant reference to the rotation factor
     *  @return A constant reference to the rotation factor
     */
    const iq_number_type& rotation( void ) const { return rotation_; }
    
    /**
     *  @brief  Get a reference to the rotation factor
     *  @return A reference to the rotation factor
     */
    iq_number_type& rotation( void ) { return rotation_; }
    
    
    /**
     *  @brief  Get a constant reference to the cusp threshold
     *  @return A constant reference to the threshold
     */
    const float_type& threshold( void ) const { return threshold_; }
    
    /**
     *  @brief  Get a reference to the cusp threshold
     *  @return A reference to the threshold
     */
    float_type& threshold( void ) { return threshold_; }
    
    
    /**
     *  @brief  Get a constant reference to the point representation
     *  @return A constant reference to the point representation
     *  @note   Only a constant representation is supplied as the point is
     *          constructed from the other factors of the cusp.
     */
    const point_type& point( void ) const { return point_; }
    
    /**
     *  @brief  Get a constant reference to the ideal generated by the cusp
     *  @return A constant reference to the ideal generated by the integral
     *          representation of the cusp
     */
    const iq_ideal_type& ideal( void ) const { return ideal_; }
    
    
    /**
     *  @brief Get the first coordinate of the projective point.
     */
    const iq_number_type& first_coordinate( void ) const { return rotation_; }
    
    
    /**
     *  @brief Compute the final coordinate of the projective point.
     */
    const iq_number_type& final_coordinate( void ) const { return final_; }
    
    
    /**
     *  @brief  Compute the inner quadratic form of the zeta component
     */
    integer_type compute_inner_qf( void ) const {
      integer_type qf = 0;
      
      typename self_type::const_zeta_ref_iterator it = self_type::zeta_ref_begin();
      for( ; it != self_type::zeta_ref_end(); ++it ) { qf += it->norm(); }
      
      return qf;
    }
    
    
    /**
     *  @brief Initialize the cusp.
     *  @note  Computes the point representation
     */
    void initialize( void ) {
      typedef utility::functors::stream_cast<integer_type,float_type> _sc;
      typedef utility::math::is_equal_cx<float_type> _iecx;
      
      // Compute the final integral coordinate of the cusp
      final_.set_reim( -compute_inner_qf(), self_type::r_ref() );
      if( !field_type::instance().is_congruent() ) final_.real() /= 2;
      final_ *= rotation_;
      final_.real() /= dilation_;
      final_.imag() /= dilation_;
      
      // This assert ensures that the cusp is really a cusp (ie is on the
      // boundary)
      assert(      ( field_type::instance().is_congruent() ? 1 : 2 )
                 * ((first_coordinate().conj()*final_coordinate()).real())
              == -compute_inner_qf()
            );
      
      // Setting the hight as zero probably isn't necessary, but it doesn't hurt
      point_.height_ref() = 0.0;
      
      // To projectivise divide the zeta component by the first coordinate
      complex_type proj(   complex_type(1.0,0.0)
                         / first_coordinate().to_complex()
                       );
      std::transform( self_type::zeta_ref_begin(), self_type::zeta_ref_end(),
                      point_.zeta_ref_begin(),
                      std::bind2nd( std::multiplies<complex_type>(), proj )
                    );
      
      // Divide out r by a factor of dilation_ (and account for the generator)
      float_type dilationFlt;
      _sc()( dilation_, dilationFlt );
      threshold_ = 1.0/dilationFlt;
      
      _sc()( self_type::r_ref(), point_.r_ref() );
      point_.r_ref() *= field_type::instance().sqrt_generator() / dilationFlt;
      
      // Account for the factor of two in the 1 mod 4 case.
      if( field_type::instance().is_congruent() ) point_.r_ref() *= 0.5;
      
      // Initialize the point
      point_.initialize();
      
      
      assert( _iecx()( point_.dependent(),
                         final_coordinate().to_complex()
                       / first_coordinate().to_complex()
                     )
            );
      
      // Build the ideal - this will become useful when someone works out how
      // to deal with the case when the the class number of the field is
      // non-trivial, currently it only offers a rather elaborate way to remove
      // multiple instances of the same cusp.
      ideal_.make_principal( rotation_ );
      
      typename self_type::const_zeta_ref_iterator it = self_type::zeta_ref_begin();
      for( ; it != this->zeta_ref_end(); ++it ) {
        ideal_ += iq_ideal_type( *it );
      }
      
      ideal_ += iq_ideal_type( final_coordinate() );
    }
    
    
    /**
     *  @brief   Print the cusp to a stream in a human friendly format.
     *  @param   os The out stream to print to, defaults to \c std::cout.
     *  @note    Why didn't I just output a string?
     */
    void pretty_print( std::ostream& __os = std::cout ) const {
      using utility::functors::stream_cast;
      
      typedef stream_cast<iq_number_type,std::ostream> _scq;
      typedef stream_cast<integer_type,float_type> _sci;
      typedef utility::math::sqrt<float_type,integer_type> _sqrt;
      typedef typename self_type::const_zeta_ref_iterator _citer;
      
      __os << "data:\n"
           << "dilation factor internal: " << dilation_
           << "\ndilation factor: " << _sqrt()(dilation_)
           << "\nrotation factor internal: " << rotation_
           << "\nrotation factor: "
           << (first_coordinate().to_complex()/_sqrt()(dilation_))
           << "\nzeta component internal:\n";
      std::copy( self_type::zeta_ref_begin(), self_type::zeta_ref_end(),
                 std::ostream_iterator<iq_number_type>(__os, "\n")
               );
      __os << "zeta component:\n";
      
      for( _citer it = self_type::zeta_ref_begin();
           it != self_type::zeta_ref_end(); ++it )
      { __os << it->to_complex() << "\n"; }
      
      float_type rFloat;
      _sci()( self_type::r_ref(), rFloat );
      rFloat *= field_type::instance().sqrt_generator();
      if( field_type::instance().is_congruent() ) rFloat *= 0.5;
      __os << "r component internal: " << this->r_ref()
           << "\nr component: " << rFloat
           << "\nis primitive: "
           << (   ideal().is_maximal_order()
                ? "yes"
                : (   field_type::instance().is_ufd()
                    ? "no"
                    : "unsure"
                  )
              );
      
      __os << "\n--\nintegral cusp:\n";
      _scq()( first_coordinate(), __os );
      __os << '\n';
      std::copy( self_type::zeta_ref_begin(), self_type::zeta_ref_end(),
                 std::ostream_iterator<iq_number_type>(__os, "\n")
               );
      _scq()( final_coordinate(), __os );
      __os << "\n--\nprojective point:\n";
      point_.pretty_print( __os );
    }
    
    
    /**
     *  @brief  Get a latex representation of the cusp
     *  @return A Latex representation of the cusp.
     */
    std::string tex( void ) const {
      std::stringstream ss;
      ss << "\\left[\\begin{smallmatrix}";
      ss << first_coordinate().tex() << " & ";
      
      typename self_type::const_zeta_ref_iterator it = self_type::zeta_ref_begin();
      for( ; it != this->zeta_ref_end(); ++it ) {
        ss << it->tex() << " & ";
      }
      
      ss << final_coordinate().tex();
      ss << "\\end{smallmatrix}\\right]^{\\mathrm{t}}";
      
      return ss.str();
    }
    
    /**
     *  @brief   Compare a point to this cusp for equality
     *  @param   rhs A hyperbolic point
     *  @return  True if the point of the cusp and the point are equal.
     */
    bool operator==( const point_type& __rhs ) const {
      return bool( point_ == __rhs );
    }
    
    
    /**
     *  @brief   Compare another cusp with this one for equality
     *  @param   rhs A cusp
     *  @return  True if the points of both cusps are equal
     */
    bool operator==( const self_type& __rhs ) const {
      return bool( *this == __rhs.point_ );
    }
    
    
    /**
     *  @brief   Weak ordering for \c cusps
     *  @param   rhs A cusp
     *  @return  True if this point is less than that on the rhs.
     *  @note    This is a totally arbitrary ordering, and has no actual
     *           meaning. However, it is a requirement of the \c std::set
     *           container that the object be ordered, it seems natural to
     *           put cusps into sets as we require them to be unique, hence
     *           the need for an ordering.
     */
    bool operator<( const self_type& __rhs ) const {
      return bool( point_ < __rhs.point_ );
    }
    
    
    /**
     *  @brief  Defines the output stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::ostream& operator<<( std::ostream& __os, const self_type& __t ){
      using std::ostream;
      using utility::functors::stream_cast;
      
      typedef stream_cast<typename self_type::heisenberg_structure_type,ostream> _scb;
      typedef stream_cast<integer_type,ostream> _sci;
      typedef stream_cast<float_type,ostream> _scf;
      typedef stream_cast<iq_number_type,ostream> _scq;
      
      __os << '[';
      _sci( )( __t.dilation_, __os );
      __os << ',';
      _scq( )( __t.rotation_, __os );
      __os << ',';
      _scb( )( __t, __os );
      __os << ',';
      _scf( )( __t.threshold_, __os );
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
      
      typedef stream_cast<string,typename self_type::heisenberg_structure_type> _scb;
      typedef stream_cast<string,integer_type> _sci;
      typedef stream_cast<string,float_type> _scf;
      typedef stream_cast<string,iq_number_type> _scq;
      
      string s;
      __is >> s;
      
      string_parser::string_vector_type vs = string_parser( )(s);
      
      _sci( )( vs.at(0), __t.dilation_ );
      _scq( )( vs.at(1), __t.rotation_ );
      _scb( )( vs.at(2), __t );
      _scf( )( vs.at(3), __t.threshold_ );
      
      __t.initialize();
      
      return __is;
    }
    
    
  private:
    //! The dilation factor of the cusp (delta)
    integer_type dilation_;
    //! The rotation factor type of the cusp (beta)
    iq_number_type rotation_;
    //! The final coordinate of the integral representation of the cusp
    iq_number_type final_;
    //! The threshold for effectiveness of the effect function
    float_type threshold_;
    //! The hyperbolic point representation
    point_type point_;
    //! The ideal generated by the integral coordinates of the cusp
    iq_ideal_type ideal_;
};

}// hyperbolic
}// geometric
}// structure
}// sg
#endif
