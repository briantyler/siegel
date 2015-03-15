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
 *  @file     region_lattice.hpp
 *  @brief    An inline header file for the \c region_lattice class.
 *  @note     Increments over imaginary quadratic lattice points in a region.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-05
 */


#ifndef _SG_REGION_LATTICE_H
#define _SG_REGION_LATTICE_H 1

// Global includes
#include <iostream>

#include <boost/operators.hpp>
#include <boost/call_traits.hpp>

// Local includes
#include "utility/math/is_even.hpp"
#include "utility/math/is_less.hpp"
#include "utility/math/comparison.hpp"
#include "utility/functors/location_builder.hpp"
#include "utility/functors/stream_cast.hpp"
#include "utility/io/string_parser.hpp"

#include "structure/numerical/iq_field.hpp"
#include "structure/numerical/iq_number.hpp"

#include "structure/geometric/euclidean/complex_region.hpp"
#include "structure/geometric/euclidean/rectangle.hpp"
#include "structure/geometric/lattice/interval_lattice.hpp"
#include "structure/geometric/lattice/iterator/region_lattice_iterator.hpp"


namespace sg
{
namespace structure
{
namespace geometric
{
namespace lattice
{
/**
 *  @class    region_lattice
 *  @brief    A class for iterating ofer the points in an imaginary quadratic
 *            lattice.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-05
 *  @param    _Float A floating point type; defaults to \c double
 *  @param    _Integer The integer type which iterates over the interval
 *  @param    _Id is the imaginary quadratic field the lattice exists in
 *  @note     Because the region gets bound to the lattice in an unsafe way,
 *            it is the responsibility of the programmer to manage the lifetime
 *            of the region in such a way that it exists at least as long as
 *            the lattice.
 */
template <class _Float = double, class _Integer = long, std::size_t _Id = 0>
    class region_lattice
  : private boost::equality_comparable< region_lattice<_Float,_Integer,_Id> >
{
  public:
    //! Object type
    typedef region_lattice<_Float,_Integer,_Id> self_type;
     //! Float type
    typedef _Float float_type;
    //! Float type
    typedef _Integer integer_type;
    //! The complex region type
    typedef euclidean::complex_region<float_type> region_type;
    //! The complex type
    typedef typename region_type::complex_type complex_type;
    //! Interval lattice type
    typedef interval_lattice<float_type, integer_type> interval_lattice_type;
    //! The interval type
    typedef typename interval_lattice_type::interval_type interval_type;
    //! Field type
    typedef numerical::iq_field<_Float,_Integer,_Id> field_type;
    //! Iq Number type
    typedef numerical::iq_number<_Float,_Integer,_Id> iq_number_type;
    
    // Add friendship to the iterator to allow access to the private enums
    friend class iterator::region_lattice_iterator<self_type>;
    
    
    /**
     *  @brief The iterator which is used for iterating over the lattice
     *  @note  A reverse iterator is not provided as there I have no need for
     *        one, but it could be done easily with a reverse iterator adaptor
     *        as the iterator is bi-directional. I have not made it random
     *        access as this seems like a load of un-necessary hassle for a
     *        virtual object that does not need random access.
     */
    typedef iterator::region_lattice_iterator<self_type> const_iterator;
    
    
  private:
    //! Integer parameter type
    typedef typename boost::call_traits<integer_type>::param_type
        integer_param_type;
    //! Float parameter type
    typedef typename boost::call_traits<float_type>::param_type
        float_param_type;
    //! The rectangle type which performs the transformation
    typedef euclidean::rectangle<float_type> rectangle_type;
    //! Even functor
    typedef utility::math::is_even<integer_type> even_func;
    //! Even functor
    typedef utility::math::is_odd<integer_type> odd_func;
    
    
  public:
    /**
     *  @brief  Default constructor
     */
    region_lattice( )
  : real_( tfRegion_.real(), sys_real_stride() ),
    imag_( tfRegion_.imag(), sys_imag_stride() ),
    bound_(0.0), transform_(1.0,0.0), invTransform_(1.0)
    { }
    
    
    /**
     *  @brief  Get a constant reference to the bound
     *  @return A constant reference to the bound
     */
    const float_type& bound( void ) const { return bound_; }
    
    /**
     *  @brief  Get a reference to the bound
     *  @return A reference to the bound
     *  @note   After assignement, initialize should be called.
     */
    float_type& bound( void ) { return bound_; }
    
    
    /**
     *  @brief  Get a constant reference to the transformation factor
     *  @return A constant reference to the transformation factor
     */
    const complex_type& transform( void ) const { return transform_; }
    
    /**
     *  @brief  Get a reference to the transformation factor
     *  @return A reference to the transformation factor
     *  @note   After assignement, initialize should be called.
     */
    complex_type& transform( void ) { return transform_; }
    
    
    /**
     *  @brief  Get a constant reference to the original region
     *  @return A constant reference to the original region
     */
    const region_type& original_region( void ) const { return orRegion_; }
    
    /**
     *  @brief  Get a reference to the original region
     *  @return A reference to the original region
     *  @note   After assignement, initialize should be called.
     */
    region_type& original_region( void ) { return orRegion_; }
    
    
    /**
     *  @brief  Get a constant reference to the transformed region
     *  @return A constant reference to the transformed region
     */
    const region_type& transformed_region( void ) const { return tfRegion_; }
    
    /**
     *  @brief  Get a constant reference to the start
     *  @return A constant reference to the start
     */
    const iq_number_type& start( void ) const { return start_; }
    
    /**
     *  @brief  Get a constant reference to the start
     *  @return A constant reference to the start
     */
    const iq_number_type& stop( void ) const { return stop_; }
    
    /**
     *  @brief  Get a constant reference to the size
     *  @return A constant reference to the size
     */
    const integer_type& size( void ) const { return size_; }
    
    
    /**
     *  @brief  Sets the stride back to the stride required by the field.
     *  @note   Mainly for testing, should be called after changing the
     *          field generator.
     */
    void reset_stride( void ) {
      real().stride() = sys_real_stride();
      imag().stride() = sys_imag_stride();
    }
    
    
    /**
     *  @brief  Get a constant reference to the real interval lattice
     */
    const interval_lattice_type& real( void ) const { return real_; }
    
    /**
     *  @brief  Get a reference to the real interval lattice
     */
    interval_lattice_type& real( void ) {
      return const_cast<interval_lattice_type&>
          ( static_cast<const self_type&>( *this ).real() );
    }
    
    
    /**
     *  @brief  Get a constant reference to the imaginary interval lattice
     */
    const interval_lattice_type& imag( void ) const { return imag_; }
    
    /**
     *  @brief  Get a reference to the imaginary interval lattice
     */
    interval_lattice_type& imag( void ) {
      return const_cast<interval_lattice_type&>
          ( static_cast<const self_type&>( *this ).imag() );
    }
    
    
    /**
     *  @brief Determine if the lattice is empty
     */
    bool empty ( void ) const { return bool( size_ == 0 ); }
    
    
    /**
     *  @brief  Access the iq_number at a given location
     *  @param  loc The location
     *  @return The iq_number at \c loc in the lattice
     */
    iq_number_type operator( )( integer_param_type __loc ) const {
      iq_number_type iqNum;
      
      switch ( latticeType_ ) {
        case prv_lattice_enum::not_congruent: {
           // In this case there is nothing odd going on, so we can just use
           // the __location_builder to get the incremental indices
          utility::functors::location_builder<integer_type> _lb( __loc );
          
          iqNum.real() = real()( _lb( real().size() ) );
          iqNum.imag() = imag()( _lb( imag().size() ) );
          
          break;
        }
        default:
          
          integer_type imagLoc( (2 * __loc) / real().size() );
          integer_type realLoc( (2 * __loc) % real().size() );
          
          switch( sizeType_ ) {
            case prv_size_enum::even_even:
            case prv_size_enum::even_odd:
              // The lattice has same the number of points in even / odd rows
              switch( startType_ ) {
                case prv_start_enum::even_even:
                case prv_start_enum::odd_odd:
                  if( odd_func()(imagLoc) ) ++realLoc;
                  break;
                default:
                  if( even_func()(imagLoc) ) ++realLoc;
              }
              break;
            default:
            // The lattice has a different number of points in even / odd rows
              switch( startType_ ) {
                case prv_start_enum::even_even:
                case prv_start_enum::odd_odd:
                // In   this case the lattice has one more point in the even rows
                // counting the first row as the zeroth.
                  break;
                default:
                  ++realLoc;
                  if( realLoc == real().size() ) {
                    realLoc = 0;
                    ++imagLoc;
                  }
                  break;
              }
          }
          iqNum.real() = real()( realLoc );
          iqNum.imag() = imag()( imagLoc );
      }

      return iqNum;
    }
    
    /**
     *  @brief  Provide an iterator to move through the lattice
     *  @param  begin The iq_number which plays the part of the lattice
     *          point during the iteeffectn.
     *  @note   Because the lattice is virtual, we need to supply a point
     *          to act as the iteratand.
     */
    const_iterator begin( iq_number_type& __begin ) const {
      __begin = start_;
      return const_iterator( *this, __begin );
    }
    
    
    /**
     *  @brief  Sentinal iterator marking the end of the lattice
     *  @param  end The iq_number which plays the part of the one-past the
     *          end point in the lattice
     *  @note   Because the lattice is virtual, we need to supply a point
     *          to act as the end point.
     *  @note   Due to there being some overhead in the construction of this
     *          sentinal iterator, it is probably best to construct it once
     *          outside of the iteeffectn loop, then test against this. There
     *          are no real problems with this, because if the lattice changes
     *          mid-iteeffectn it is invalidated anyway.
     */
    const_iterator end( iq_number_type& __end ) const {
      __end = stop_;
      return const_iterator( *this, __end );
    }
    
    
    /**
     *  @brief  Intitialize the interval latice
     */
    void initialize( void ) {
      sys_compute_transformed();
      
       // If either lattice is empty then there is nothing to do.
      if( real().empty() || imag().empty() ) size_ = 0;
      
      // Although this might look silly, it is in fact much easier to compute
      // the congruence classes of the field, the sizes and the starts and then
      // compute the starts and the sizes, rather than try to do both at once
      // with a load of if's and switches.
      latticeType_ =   field_type::instance().is_congruent()
          ? prv_lattice_enum::congruent
        : prv_lattice_enum::not_congruent;
      
      if( even_func()(real().size()) ) {
        if( even_func()(imag().size()) ) sizeType_ = prv_size_enum::even_even;
        else sizeType_ = prv_size_enum::even_odd;
      }
      else {
        if( even_func()(imag().size()) ) sizeType_ = prv_size_enum::odd_even;
        else sizeType_ = prv_size_enum::odd_odd;
      }
      
      if( even_func()(real().start()) ) {
        if( even_func()(imag().start()) ) startType_ = prv_start_enum::even_even;
        else startType_ = prv_start_enum::even_odd;
      }
      else {
        if( even_func()(imag().start()) ) startType_ = prv_start_enum::odd_even;
        else startType_ = prv_start_enum::odd_odd;
      }
      
      sys_compute_size();
      sys_compute_start();
      sys_compute_stop();
    }
    
    
    /**
     *  @brief  Validate a point to see if, under transformation it is close
     *          enough to the original region
     *  @param  value The \c iq_number to check
     *  @return True if the value is closer than the bound to the original
     *          region under the inverse transformation.
     */
    bool validate( const iq_number_type& value ) {
      typedef utility::math::is_less<float_type> _il;
      
      return bool( _il()( orRegion_.distance(   value.to_complex()
                                              * invTransform_
                                            ), bound_
                        )
                 );
    }
    
    
  private:
    /**
     *  @brief  Compute the transformed region
     */
    void sys_compute_transformed( void ) {
      // Set the inverse transform
      invTransform_ = complex_type(1.0,0.0) / transform_;
      
      // Initialize the transformation region
      tfRegion_ = orRegion_;
      tfRegion_.extend( bound_ );
      
      // create the transforming rectangle
      rectangle_type r( tfRegion_.bl(), tfRegion_.tr() );
      // Transform the rectangle
      r.transform_contain( transform_ );
      // Update the transformed region
      tfRegion_.from_rectangle( r );
      
      // Update the interval lattices
      real().initialize();
      imag().initialize();
    }
    
    
    /**
     *  @brief  Compute the number of points in the lattice
     */
    void sys_compute_size( void ) {
      size_ = ( real().size() * imag().size() );
      if( latticeType_ == prv_lattice_enum::not_congruent ) return;
      
      size_ /= 2;
      if( sizeType_ != prv_size_enum::odd_odd ) return;
      
      if(    startType_ == prv_start_enum::even_even
          || startType_ == prv_start_enum::odd_odd ) ++size_;
    }
    
    
    /**
     *  @brief  Compute the first point in the lattice
     */
    void sys_compute_start( void ) {
      if(   latticeType_ == prv_lattice_enum::not_congruent
          || startType_ == prv_start_enum::even_even
          || startType_ == prv_start_enum::odd_odd ) {
        start_.set_reim( real().start(), imag().start() );
        
      }
      else {
        if( real().size() == 1 ) {
          start_.set_reim( real().start(), imag().start() + 1 );
        }
        else {
          start_.set_reim( real().start() + 1, imag().start() );
        }
      }
    }
    
    
    /**
     *  @brief  Compute the first point after the end of the lattice
     *  @note   Must be called after computing the start value.
     */
    void sys_compute_stop( void ) {
      if( size() == 0 ) {
        stop_ = start_;
        return;
      }
      
      stop_.set_reim( real().start(), imag()( imag().size() - 1 ) + 1 );
      
      if( latticeType_ == prv_lattice_enum::congruent ) {
        if( even_func()(stop_.real()) == odd_func()(stop_.imag()) ) {
          if( real().size() == 1 ) {
            ++(stop_.imag());
          }
          else {
            ++(stop_.real());
          }
        }
      }
    }
    
    
    /**
     *  @brief  Get the standard real stride for the current field.
     */
    static float_type sys_real_stride() {
      return float_type( field_type::instance().is_congruent() ? 0.5 : 1.0 );
    }
    
    /**
     *  @brief  Get the standard real stride for the current field.
     */
    static float_type sys_imag_stride() {
      return float_type(   (field_type::instance().is_congruent() ? 0.5 : 1.0)
                         * field_type::instance().sqrt_generator()
                       );
    }
    
  public:
    /**
     *  @brief  Compares two region_lattice objects for equality.
     *  @param  lhs A \c region_lattice
     *  @param  rhs A \c region_lattice
     *  @note   Two lattices are equal if they have the same orginal region,
     *          same bound and transformation factors and the same interval
     *          strides.
     */
    bool operator==( const self_type& __rhs ) const {
      typedef utility::math::is_equal<float_type> _ie;
      typedef utility::math::is_equal_cx<float_type> _iecx;
      
      return bool(    orRegion_ == __rhs.orRegion_
                   && _ie()(bound_, __rhs.bound_)
                   && _iecx()(transform_, __rhs.transform_)
                   && _ie()(real().stride(), __rhs.real().stride())
                   && _ie()(imag().stride(), __rhs.imag().stride())
                 );
    }
    
    
    /**
     *  @brief  Defines the output stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     *  @note   The output stream takes the form:
     *          [original region, bound, transform, real stride, imag stride]
     */
    friend std::ostream& operator<<( std::ostream& __os, const self_type& __t ) {
      using utility::functors::stream_cast;
      
      typedef stream_cast<region_type,std::ostream> _scr;
      typedef stream_cast<float_type,std::ostream> _scf;
      typedef stream_cast<complex_type,std::ostream> _scc;
      
      __os << '[';
      _scr()( __t.orRegion_, __os );
      __os << ',';
      _scf()( __t.bound_, __os );
      __os << ',';
      _scc()( __t.transform_, __os );
      __os << ',';
      _scf( )( __t.real().stride(), __os );
      __os << ',';
      _scf( )( __t.imag().stride(), __os );
      __os << ']';
      
      return __os;
    }
    
    /**
     *  @brief  Defines the input stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream in.
     */
    friend std::istream& operator>>( std::istream& __is, self_type& __t ) {
      using utility::functors::stream_cast;
      
      typedef stream_cast<std::string,region_type> _scr;
      typedef stream_cast<std::string,float_type> _scf;
      typedef stream_cast<std::string,complex_type> _scc;
      typedef utility::io::string_parser _sp;
      
      std::string s;
      __is >> s;
      _sp::string_vector_type vs = _sp()( s );
      
      _scr( )( vs.at(0), __t.orRegion_ );
      _scf( )( vs.at(1), __t.bound_ );
      _scc( )( vs.at(2), __t.transform_ );
      _scf( )( vs.at(3), __t.real().stride() );
      _scf( )( vs.at(4), __t.imag().stride() );
      __t.initialize();
      
      return __is;
    }
    
    
  private:
    // Encode information about the structure of the lattice
    struct prv_lattice_enum {
      enum lattice {
        not_congruent,
        congruent
      };
    };
    
    struct prv_start_enum {
      enum start {
        even_even,
        even_odd,
        odd_even,
        odd_odd
      };
    };
    
    struct prv_size_enum {
      enum size {
        even_even,
        even_odd,
        odd_even,
        odd_odd
      };
    };
    
    //! The original region
    region_type orRegion_;
    //! The transformed region
    region_type tfRegion_;
    //! The real interval lattice
    interval_lattice_type real_;
    //! The imaginary interval lattice
    interval_lattice_type imag_;
    
    //! The first point in the lattice
    iq_number_type start_;
    //! The first point after the end of the lattice_
    iq_number_type stop_;
    //! The number of lattice points in the region
    integer_type size_;
    
    //! A bound on the inner region
    float_type bound_;
    //! The transformation factor original region -> transformed region
    complex_type transform_;
    //! The inverse transformation factor transformed region -> original region
    complex_type invTransform_;
    
    //! A classification of the lattice type
    typename prv_lattice_enum::lattice latticeType_;
    //! A classification of the starting coordinate
    typename prv_start_enum::start startType_;
    //! A classification of the  stopping coordinate
    typename prv_size_enum::size sizeType_;
};

}// lattice
}// geometric
}// structure
}// sg
#endif
