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
 *  @file     interval_lattice.hpp
 *  @brief    An inline header file for the \c interval_lattice class.
 *  @note     Increments over lattice points in an interval.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-05
 */


#ifndef _SG_INTERVAL_LATTICE_H
#define _SG_INTERVAL_LATTICE_H 1

// Global includes
#include <cassert>

#include <boost/operators.hpp>
#include <boost/call_traits.hpp>

// Local includes
#include "structure/geometric/euclidean/real_interval.hpp"
#include "structure/geometric/lattice/detail/basic_range_fix.hpp"
#include "utility/math/is_equal.hpp"
#include "utility/math/rounding.hpp"
#include "utility/functors/stream_cast.hpp"
#include "utility/io/string_parser.hpp"


namespace sg
{
namespace structure
{
namespace geometric
{
namespace lattice
{
/**
 *  @class    interval_lattice
 *  @brief    A class for slicing an interval.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-05
 *  @param    _Float A floating point type; defaults to \c double
 *  @param    _Integer The integer type which iterates over the interval
 *  @note     Because the interval gets bound to the lattice in an unsafe way,
 *            it is the responsibility of the programmer to manage the lifetime
 *            of the interval in such a way that it exists at least as long as
 *            the lattice.
 *  @note     Due to a bug in the boost \c basic_range class preventing reverse
 *            iteeffectn through negative ranges, I have made a local copy and
 *            corrected the bug in that.
 */
template <class _Float = double, class _Integer = long>
    class interval_lattice
  : private boost::equality_comparable< interval_lattice<_Float,_Integer> >
{
  public:
    //! Object type
    typedef interval_lattice<_Float,_Integer> self_type;
     //! Float type
    typedef _Float float_type;
    //! Float type
    typedef _Integer integer_type;
    //! Interval type
    typedef euclidean::real_interval<float_type> interval_type;
    
    
  private:
    //! Integer parameter type
    typedef typename boost::call_traits<integer_type>::param_type
        integer_param_type;
    //! Float parameter type
    typedef typename boost::call_traits<float_type>::param_type
        float_param_type;
    //! Range type
    typedef boost::numeric::ublas::basic_range_fix<integer_type,integer_type>
        range_type;
    
  public:
    //! Forward iterator
    typedef typename range_type::const_iterator const_iterator;
    //! Reverse iterator
    typedef typename range_type::const_reverse_iterator const_reverse_iterator;
    
    
    /**
     *  @brief  Default constructor
     */
    interval_lattice( ) : interval_(0), stride_(1.0) { }
    
    
    /**
     *  @brief  Interval constructor
     */
    interval_lattice
        ( const interval_type& __interval, float_param_type __stride )
  : interval_( &__interval ), stride_( __stride ) { initialize(); }
    
    
    /**
     *  @brief  Bind an interval to the lattice
     *  @param  interval The interval to bind.
     */
    void bind( const interval_type& __interval ) {
      interval_ = &__interval;
    }
    
    /**
     *  @brief  Release the interval from the lattice
     */
    void release( void ) { interval_ = 0; }
    
    
    /**
     *  @brief  Get the stride of the interval.
     *  @return The stride of the interval
     */
    float_type stride( void ) const { return stride_; }
    
    /**
     *  @brief  Get a reference to the stride of the interval.
     *  @return A reference to the stride of the interval.
     */
    float_type& stride( void ) { return stride_; }
    
    /**
     *  @brief  Get the start of the interval.
     *  @return The start of the interval.
     */
    integer_type start( void ) const { return range_.start(); }
    
    /**
     *  @brief  Get the value one after the end of the interval.
     *  @return The value one after the end of the interval
     */
    integer_type stop( void ) const { return stop_; }
    
    /**
     *  @brief  Get the number of elements in the interval.
     *  @return The number of elements in the interval
     */
    integer_type size( void ) const { return range_.size(); }
    
    
    /**
     *  @brief Determine if the lattice is empty
     */
    bool empty ( void ) const { return bool( range_.empty() ); }
    
    
    /**
     *  @brief  Access the value at a given location
     *  @param  loc The location
     *  @return start + loc
     */
    integer_type operator( )( integer_param_type __loc ) const {
      return integer_type( range_(__loc) );
    }
    
    
    /**
     *  @brief Get the iterator pointing at start
     */
    const_iterator begin() const { return range_.begin(); }
    
    /**
     *  @brief Get the iterator pointing one past the end
     */
    const_iterator end() const { return range_.end(); }
    
    /**
     *  @brief Get the iterator pointing at the reverse start
     */
    const_reverse_iterator rbegin() const { return range_.rbegin(); }
    
    /**
     *  @brief Get the iterator pointing one before the reverse end
     */
    const_reverse_iterator rend() const { return range_.rend(); }
    
    
    /**
     *  @brief  Intitialize the interval latice
     *  @note   Computes the start and size of the lattice
     */
    void initialize( void ) {
      typedef utility::math::ceil<integer_type,float_type> _ceil;
      typedef utility::math::floor<integer_type,float_type> _floor;
      // At debug time check that the lattice is bound to an interval and that
      // the stride is strictly positive.
      assert( interval_ != 0 && stride_ >= 0 );
      
      interval_type tmp( *interval_ / stride_ );
      integer_type start = _ceil()( tmp.lower() );
      stop_ = _floor()( tmp.upper() ) + 1;
      
      range_ = range_type(start,stop_);
    }
    
    
    /**
     *  @brief  Compares two interval_lattice objects for equality.
     *  @param  lhs A \c interval_lattice
     *  @param  rhs A \c interval_lattice
     *  @note   Two ranges are equal if they both point at the same interval
     *          and have the same stride
     */
    bool operator==( const self_type& __rhs ) const {
      typedef utility::math::is_equal<float_type> _ie;
      
      return bool(    interval_ == __rhs.interval_
                   && _ie()( stride_, __rhs.stride_ )
                 );
    }
    
    /**
     *  @brief  Defines the output stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::ostream& operator<<( std::ostream& __os, const self_type& __t ) {
      using utility::functors::stream_cast;
      typedef stream_cast<float_type,std::ostream> _scf;
      typedef stream_cast<interval_type,std::ostream> _sci;
      
      __os << '[';
      _scf( )( __t.stride_, __os );
      __os << ',';
      _sci( )( *(__t.interval_), __os );
      __os << ']';
      
      return __os;
    }
    
    /**
     *  @brief  Defines the input stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream in.
     */
    friend std::istream& operator>>( std::istream& __is, self_type& __t ) {
      typedef utility::functors::stream_cast<std::string,float_type> _sc;
      typedef utility::io::string_parser _sp;
      
      std::string s;
      __is >> s;
      _sp::string_vector_type vs = _sp()( s );
      
      _sc( )( vs.at(0), __t.stride_ );
      __t.initialize();
      
      return __is;
    }
    
    
  private:
    //! The underlying interval
    const interval_type* interval_;
    //! The stride
    float_type stride_;
    //! The range which handles most of the work
    range_type range_;
    //! One past the final element
    integer_type stop_;
};

}// lattice
}// geometric
}// structure
}// sg
#endif
