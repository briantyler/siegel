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
 *  @file     heisenberg_slice.hpp
 *  @brief    An inline header file for the \c heisenberg_slice class.
 *  @note     Heisenberg Slice specifies a slice of indices in a horosphere of
 *            \c hyperbolic_space
 *  @author   Brian Tyler
 *  @version  1.1
 *  @date     2008-03-18
 */

#ifndef _SG_HEISENBERG_SLICE_H
#define _SG_HEISENBERG_SLICE_H 1

// Global includes
#include <algorithm>
#include <functional>

#include <boost/call_traits.hpp>

// Local includes
#include "structure/geometric/euclidean/accessors/interval_data_accessor.hpp"
#include "structure/geometric/euclidean/accessors/region_data_accessor.hpp"
#include "structure/geometric/detail/heisenberg_structure.hpp"
#include "structure/geometric/hyperbolic/hyperbolic_space.hpp"
#include "structure/geometric/hyperbolic/hyperbolic_point.hpp"
#include "structure/geometric/hyperbolic/heisenberg_cube.hpp"
#include "structure/geometric/hyperbolic/incrementor/heisenberg_incrementor.hpp"
#include "structure/geometric/hyperbolic/iterators/cube_iterator.hpp"
#include "structure/geometric/hyperbolic/iterators/slice_iterator.hpp"
#include "utility/math/pow.hpp"
#include "utility/math/is_even.hpp"
#include "utility/math/floor.hpp"
#include "utility/algorithms/mutate.hpp"
#include "utility/functors/location_builder.hpp"
#include "geometry/measure/measure.hpp"
#include "geometry/measure/interval_measure.hpp"


namespace sg
{
namespace structure
{
namespace geometric
{
namespace hyperbolic
{
/**
 *  @class    heisenberg_slice
 *  @brief    A class for slicing a compact region of hyperbolic space into
 *            iterable subdivisions.
 *  @author   Brian Tyler
 *  @version  1.1
 *  @date     2008-03-20
 *  @param    N The dimension of the space to slice.
 *  @param    _Float The \c float_type of the \c hyperbolic_space class to
 *            slice. Defaults to \c double.
 *  @param    _Integer The integer type used for indexing the slice. Defaults
 *            to \c unsigned_long
 */
template < std::size_t N, class _Float = double, class _Integer = std::size_t >
    class heisenberg_slice
  : public geometric::detail::heisenberg_structure<
             N,
             euclidean::accessor::region_data_accessor<_Float,_Integer>,
             euclidean::accessor::interval_data_accessor<_Float,_Integer>
           >
{
  public:
    //! Object type
    typedef heisenberg_slice<N,_Float,_Integer> self_type;
    //! The size type
    typedef std::size_t size_type;
    //! Floating point type
    typedef _Float float_type;
    //! Integer type
    typedef _Integer integer_type;
    //! Space type
    typedef hyperbolic_space<N,_Float> space_type;
    //! Point type
    typedef hyperbolic_point<N,_Float> point_type;
    //! Cube type
    typedef heisenberg_cube<N,_Float> cube_type;
    //! Interval type
    typedef typename space_type::real_type interval_type;
    //! Interval data type
    typedef typename self_type::real_type interval_data_type;
    //! Heisenberg incrementor type
    typedef incrementor::heisenberg_incrementor<point_type> incrementor_type;
    //! Hypercube iterator type
    typedef iterators::cube_iterator<self_type> cube_iterator;
    //! Slice iterator type
    typedef iterators::slice_iterator<self_type> slice_iterator;
    
    
  private:
    //! Parameter type of the integer type
    typedef typename boost::call_traits<integer_type>::param_type int_param_type;
    //! Parameter type of the float type
    typedef typename boost::call_traits<float_type>::param_type flt_param_type;
    
  public:
    /**
     *  @brief  Default constructor
     */
    heisenberg_slice( )
    : space_(), resolution_(1) { }
    
    
    /**
     *  @brief  Space constructor
     *  @param  space The \c hyperbolic_space to slice
     *  @note   The slice must have the resolution set using set_resolution
     *          before it is usable.
     */
    heisenberg_slice( const space_type& __space )
    : space_( __space ), resolution_(1) { }
    
    
    /**
     *  @brief  Space and resolution constructor
     *  @param  space The \c hyperbolic_space to slice
     *  @param  resolution The resolution at which to slice.
     */
    heisenberg_slice( const space_type& __space, int_param_type __resolution )
    : space_( __space ), resolution_(1)
    { this->sys_set_resolution( __resolution ); }
    
    
    /**
     *  @brief  Initialize the slice after changing the space or resolution.
     */
    void initialize( void ) {
      integer_type resolution = resolution_;
      sys_set_resolution( resolution );
    }
    
    
    /**
     *  @brief  Get a constant reference to the resolution of the slice
     *  @return A constant reference to the resolution of the space
     */
    const integer_type& resolution( void ) const { return resolution_; }
    
    /**
     *  @brief  Get a refernce to the resolution of the slice
     *  @return A reference to the resolution of the space
     *  @note   After changing the resolution initialize must be called before
     *          attempting to use the slice.
     */
    integer_type& resolution( void ) { return resolution_; }
    
    /**
     *  @brief  Get constant reference to the \c hyperbolic_space
     *  @return A constant reference to the \c hyperbolic_space
     */
    const space_type& space( void ) const { return space_; }
    
    /**
     *  @brief  Get reference to the \c hyperbolic_space
     *  @return A reference to the \c hyperbolic_space
     *  @note   After changing the space initialize must be called before
     *          attempting to use the slice.
     */
    space_type& space( void ) { return space_; }
    
    
    const incrementor_type& forward( void ) const { return forward_; }
    const incrementor_type& reverse( void ) const { return reverse_; }
    const incrementor_type& reset( void ) const { return reset_; }
    const incrementor_type& rreset( void ) const { return rreset_; }
    
    
    /**
     *  @brief  Return the space at a given index
     *  @param  loc The lindex ocation of the subspace to return.
     *  @return A copy of the subspace at \c loc
     *  @note   Due to the virtual nature of the slice, it is legal to
     *          construct subspaces which lie outside of the slice. There
     *          seems to be no good reason to dissallow this possibly useful
     *          behaviour.
     */
    space_type subspace_at( int_param_type __loc ) const {
      using std::transform;
      
      // Create the subspace
      space_type subspace;
      
      // Set the height to be the same as the referenced space
      subspace.height() = space_.height();
      
      // Construct the subspace at loc
      transform( self_type::begin(), self_type::heisenberg_end(),
                 space_.begin(), subspace.begin(),
                 prv_subspace_from_data( __loc )
               );
      
      return subspace;
    }
    
    
    /**
     *  @brief  Return the point at a given index
     *  @param  loc The index location of the point to return.
     *  @return A copy of the point at \c loc
     *  @note   The point returned is equal to the lower corner of the
     *          subspace at \c loc
     *  @note   Due to the virtual nature of the slice, it is legal to
     *          construct points which lie outside of the slice.
     */
    point_type point_at( int_param_type __loc ) const {
      using std::transform;
      
      // Create the point
      point_type point;
      
      // Set the height to be the same as the referenced space
      point.height() = space_.height_ref().lower();
      
      // Construct the point at loc
      transform( self_type::begin(), self_type::heisenberg_end(),
                 space_.begin(), point.begin(),
                 prv_point_from_data( __loc )
               );
      
      // Deal with points outside the space
      if( __loc >= resolution_ ) {
        typedef utility::functors::stream_cast<integer_type, float_type> _sc;
        float_type factor;
        _sc()( (__loc / resolution_), factor );
        point.r_ref() += factor * space_.r_ref().length();
      }
      
      // Initialize (this means compute the dependent variable)
      point.initialize();
      
      return point;
    }
    
    
    /**
     *  @brief  Return the hypercube at a given index
     *  @param  loc The index location of the hypercube to return.
     *  @return A copy of the hypercube at \c loc
     *  @note   The first vertex of the returned hypercube is equal to the
     *          point returned by \c this->point_at(loc) and is also equal
     *          to the subspace of minimal height returned by
     *          \c this->subspace_at(loc)
     *  @note   Due to the virtual nature of the slice, it is legal to
     *          construct hypercubes which lie outside of the slice.
     */
    cube_type cube_at( int_param_type __loc ) const {
      typedef utility::math::is_odd<integer_type> _od;
      typedef typename cube_type::iterator _iter;
      
      cube_type cube;
      
      // Take the base point of the cube to be the point at loc.
      std::fill( cube.begin(), cube.end(), point_at(__loc) );
      
      // Iterate through all points in the cube.
      integer_type index = 0;
      for( _iter it = cube.begin(); it != cube.end(); ++it, ++index ) {
        // Record the current index of the point in the hypercube
        integer_type indexTmp = index;
        // Record the current index of the coordinate in the point
        size_type i = 0;
        
        /* The idea is that the incrementor at "i" is used if and only if the
         * i'th binary digit of index is non-zero.
         * This is an easy way to ensure all points propogate correctly from
         * the base point, since it will give all possible combinations of
         * incrementors; it works because the cardinatlity of the cube is
         * 2^(heisenberg_size).
         */
        while( indexTmp != 0 ) {
          // Is the i'th digit of index 1?
          if( _od()(indexTmp) ) {
            // YES, apply the i'th incrementor to the ith coordinate.
            forward_.common_at(i)( *it, it->common_at(i) );
          }
          // Strip off the last digit and increment i.
          indexTmp /= 2;
          ++i;
        }
      }
      
      return cube;
    }
    
    
    /**
     *  @brief  Return the index of a hypercube where \c point is the leading
     *          vertex (the back() vertex).
     *  @param  point The point which determines the leading vertex
     *  @return The index of a hypercube which is guaranteed have leading
     *          vertex which is less than \c point in each dimension if point
     *          is inside the original space.
     */
    integer_type location_at( const point_type& __point ) const {
      typedef utility::math::floor<integer_type,float_type> _floor;
      typedef typename point_type::const_iterator _piter;
      typedef typename space_type::const_iterator _spiter;
      typedef typename self_type::const_iterator _sliter;
      
      _piter pit = __point.begin();
      _spiter spit = space_.begin();
      _sliter slit = self_type::begin();
      
      integer_type loc = 0;
      integer_type res = 1;
      for( ; slit != self_type::heisenberg_end(); ++slit, ++spit, ++pit ) {
        // If the point is outside of the space then make the index zero
        if( *pit < spit->lower() + slit->stride() || *pit > spit->upper() ) {
          loc = 0;
          break;
        }
        
        // Compute the index of the point below the current point.
        float_type x = ( *pit - spit->lower() ) / slit->stride() - 1.0;
        if( x < 0.0 ) x = 0.0;
        
        // Compute the new location and resolution.
        loc += _floor()(x) * res;
        res *= slit->resolution();
      }
      
      return loc;
    }
    
    
  private:
    //! Pointer to a heisenberg cube
    typedef typename cube_iterator::cube_ptr_type cube_ptr_type;
    
    
  public:
    /**
     *  @brief  Return a slice iterator pointing to the first point in the slice
     *  @return A starting \c slice_iterator iterator.
     */
    slice_iterator slice_begin( void ) const {
      return slice_iterator( this, 0 );
    }
    
    
    /**
     *  @brief  Return a cube iterator pointing to the first index
     *  @return A starting \c heisenberg_cube iterator.
     */
    cube_iterator cube_begin( void ) const {
      cube_ptr_type ptr( new cube_type( cube_at(0) ) );
      return cube_iterator( *this, 0, ptr );
    }
    
    
    /**
     *  @brief  Return a cube iterator pointing one past the last index
     *  @return A one past the end \c cube_iterator.
     */
    cube_iterator cube_end( void ) const {
      cube_ptr_type ptr( new cube_type( cube_at(resolution_) ) );
      return cube_iterator( *this, resolution_, ptr );
    }
    
    
    /**
     *  @brief  Return a sentinal end cube iterator.
     *  @return A sentinal ending \c heisenberg_cube iterator.
     *  @note   This iterator is not valid for iteeffectn, but is much more
     *          efficient than \c cube_end for the purposes of an iteeffectn
     *          sentinal since no \c heisenberg_cube is actually constructed.
     */
    cube_iterator cube_shallow_end( void ) const {
      return cube_iterator( *this, resolution_ );
    }
    
    
  private:
    /**
     *  @brief  Set the resolution of the slice
     *  @return A constant reference to the new resolution of the space
     *  @note   The input resolution is very unlikely to be the same as the
     *          output resolution because the final resolution is chosen to
     *          exactly and evenly divide the space, whilst being as close to
     *          the original input as possible. For example; in a 2D-space
     *          there is no way of evenly dividing the 1-D subspaces to give
     *          a resolution of 97, because \f$ \sqrt{97} \not \in \mathbf{Z}
     *          \f$.
     */
    const integer_type& sys_set_resolution( int_param_type __resolution ) {
      using geometry::measure::measure;
      using utility::algorithms::mutate;
      
      typedef interval_data_type _idt;
      typedef typename _idt::template compute_data<interval_type> _ci;
      typedef geometry::measure::interval_measure<interval_type> _im;
      typedef utility::math::powf<float_type> _pow;
      
      /*  len ~ length of the side of a hypercube with the same measure
       *        as the space.
       *  res ~ the resolution of each interval if the space was a hypercube.
       *  exp ~ the root exponent.
       */
      typedef utility::functors::stream_cast<integer_type, float_type> _sc;
      float_type resolution;
      _sc()( __resolution, resolution );
      
      float_type exp( 1.0/static_cast<float_type>(self_type::heisenberg_size) );
      float_type res( _pow()( resolution, exp ) );
      float_type len( _pow()( measure( space_.begin(),
                                       space_.heisenberg_end(), _im( )
                                     ), exp
                            )
                    );
      
      resolution_ = mutate( space_.begin(), space_.heisenberg_end(),
                            self_type::begin(), _ci( res, len )
                          );
      
      sys_initialize_incrementors();
      
      return resolution_;
    }
    
    
    //! The location builder type
    typedef utility::functors::location_builder<integer_type> lb_func;
    
    
    /**
     *  @struct   prv_subspace_from_data
     *  @brief    A functor constructing a subspace from interval data, an
     *            interval and a location
     *  @author   Brian Tyler
     *  @version  1.0
     */
    struct prv_subspace_from_data
    : public std::binary_function< const interval_type&,
                                   const interval_data_type&,
                                   interval_type>
     {
       /**
        *  @brief  Location constructor.
        *  @param  loc the index of the subspace to construct
        */
       prv_subspace_from_data( int_param_type __loc )
         : lb_(__loc) { }
       
       /**
        *  @brief  Get the required subinterval.
        *  @param  data The interval data describing how this interval has
        *          been sliced
        *  @param  interval The interval to take the sub interval from.
        *  @return The subinterval for this component of the Heisenberg
        *          structure.
        */
       interval_type operator( )
           ( const interval_data_type& data, const interval_type& interval )
       {
         /* Each time it is called lb_ returns the index of the subspace in
          * this interval.
          */
         return interval_type( data.subinterval_at( interval,
                                                    lb_( data.resolution() )
                                                  )
                             );
       }
       
       private:
         //! Location builder
         lb_func lb_;
     };
     
     
     /**
      *  @struct   prv_point_from_data
      *  @brief    A functor constructing a point from interval data, an
      *            interval and a location
      *  @author   Brian Tyler
      *  @version  1.0
      */
    struct prv_point_from_data
    : public std::binary_function< const interval_type&,
                                   const interval_data_type&,
                                   float_type
                                 >
     {
       /**
        *  @brief  Location constructor.
        *  @param  loc the index of the subspace to construct
        */
       prv_point_from_data( int_param_type __loc )
         : lb_(__loc) { }
       
       /**
        *  @brief  Get the required point.
        *  @param  data The interval data describing how this interval has
        *          been sliced
        *  @param  interval The interval to take the sub interval from.
        *  @return The value for this coordinate of the Heisenberg
        *          structure.
        */
       float_type operator( )
           ( const interval_data_type& data, const interval_type& interval )
       {
         typedef utility::functors::stream_cast<integer_type, float_type> _sc;
         float_type loc;
         _sc()( lb_( data.resolution() ), loc );
         
         return float_type( interval.lower() + data.stride()*loc );
       }
       
       private:
         //! Location builder
         lb_func lb_;
     };
     
     
     /**
      *  @brief  Initialize the incrementors
      */
     void sys_initialize_incrementors( void ) {
       /* because the r_incrementor and the zeta_incrementors have different
        * types it is not possible to increment through the
        * heisenberg_incrementor to set the incrementor strides. Therefore
        * the iteeffectn must be done "by hand."
        */
       for( size_type i = 0; i != self_type::zeta_size; ++i ) {
         // Set the stride of the real component
         sys_set_incrementors
             ( forward_.zeta_at(i).real_ref(),
               reverse_.zeta_at(i).real_ref(),
               reset_.zeta_at(i).real_ref(),
               rreset_.zeta_at(i).real_ref(),
               self_type::zeta_at(i).real_ref().stride(),
               self_type::zeta_at(i).real_ref().resolution()
             );
  
         // Set the stride of the imag component
         sys_set_incrementors
             ( forward_.zeta_at(i).imag_ref(),
               reverse_.zeta_at(i).imag_ref(),
               reset_.zeta_at(i).imag_ref(),
               rreset_.zeta_at(i).imag_ref(),
               self_type::zeta_at(i).imag_ref().stride(),
               self_type::zeta_at(i).imag_ref().resolution()
             );
       }
       // Zeta initialized!
       
       
       // Set the stride of the r component
       sys_set_incrementors
           ( forward_.r_ref(),reverse_.r_ref(),
             reset_.r_ref(), rreset_.r_ref(),
             self_type::r_ref().stride(), self_type::r_ref().resolution()
           );
       // r initialized!
     }
     
     
     /**
      *  @brief  Set the forward, reverse and reset iterators
      *  @param  forward The forward incrementor
      *  @param  reverse The reverse incrementor
      *  @param  reset The forward reset incrementor
      *  @param  rreset The reverse reset incrementor
      *  @param  stride The stride of the interval
      *  @param  resolution The resolution of the interval
      */
     template <class _Incrementor>
         void sys_set_incrementors
         ( _Incrementor& __forward, _Incrementor& __reverse,
           _Incrementor& __reset, _Incrementor& __rreset,
           flt_param_type __stride, int_param_type __resolution )
     {
       __forward.set_stride( __stride );
       __reverse.set_stride( -__stride );
       
       typedef utility::functors::stream_cast<integer_type, float_type> _sc;
       float_type tmp;
       _sc()(  __resolution - 1, tmp );
       
       tmp *= __stride;
       __reset.set_stride( -tmp );
       __rreset.set_stride( tmp );
     }
    
    
    //! The \c hyperbolic_space to slice
    space_type space_;
    //! The resolution of the slice.
    integer_type resolution_;
    
    //! The forward incrementor
    incrementor_type forward_;
    //! The reverse incrementor
    incrementor_type reverse_;
    //! The reset incrementor
    incrementor_type reset_;
    //! The reverse reset incrementor
    incrementor_type rreset_;
};

}// hyperbolic
}// geometric
}// structure
}// sg
#endif
