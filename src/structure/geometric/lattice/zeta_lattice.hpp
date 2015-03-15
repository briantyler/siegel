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
 *  @file     zeta_lattice.hpp
 *  @brief    An inline header file for zeta_lattice.
 *  @note     Include this file to create zeta lattices
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-08
 */

#ifndef _SG_ZETA_LATTICE_H
#define _SG_ZETA_LATTICE_H 1

// Global includes
#include <cassert>
#include <algorithm>
#include <functional>

#include <boost/call_traits.hpp>

// Local includes
#include "structure/geometric/detail/zeta_array.hpp"
#include "structure/geometric/lattice/accessors/region_lattice_accessor.hpp"
#include "structure/geometric/lattice/iterator/zeta_lattice_iterator.hpp"
#include "structure/geometric/lattice/detail/lattice_bound.hpp"
#include "structure/geometric/lattice/zeta_lattice_point.hpp"
#include "structure/geometric/hyperbolic/hyperbolic_space.hpp"

#include "utility/functors/location_builder.hpp"
#include "utility/algorithms/copy_member.hpp"
#include "utility/algorithms/fill_member.hpp"
#include "utility/math/is_greater.hpp"
#include "utility/math/is_zero.hpp"
#include "utility/math/square.hpp"


namespace sg
{
namespace structure
{
namespace geometric
{
namespace lattice
{
/**
 *  @class    zeta_lattice
 *  @brief    A class representing a lattice in complex space
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-08
 *  @param    N The hyperbolic dimension of the lattice
 *  @param    _Float A floating point type; defaults to \b double
 *  @param    _Integer The integer type which iterates over the interval
 *  @param    _Id is the id of the imaginary quadratic field the lattice exists
 *            in.
 *  @note     This object is specific to the SU(n,1) problem.
 */
  template <std::size_t N, class _Float = double, class _Integer = long,
            std::size_t _Id = 0>
  class zeta_lattice
  : public geometric::detail::zeta_array< N,
               accessor::region_lattice_accessor<_Float,_Integer,_Id>
           >
{
  public:
    //! Object type
    typedef zeta_lattice<N,_Float,_Integer,_Id> self_type;
    //! Float type
    typedef _Float float_type;
    //! The integer type
    typedef _Integer integer_type;
    //! The Region Lattice type;
    typedef typename self_type::zeta_type region_lattice_type;
    //! The Interval Lattice type
    typedef typename region_lattice_type::interval_lattice_type
        interval_lattice_type;
    //! The region type
    typedef typename region_lattice_type::region_type region_type;
    //! The imaginary quadratic number type
    typedef typename region_lattice_type::iq_number_type iq_number_type;
    //! The imaginary quadratic field type
    typedef typename iq_number_type::field_type field_type;
    //! The complex type
    typedef typename region_type::complex_type complex_type;
    //! The interval type
    typedef typename interval_lattice_type::interval_type interval_type;
    //! The space type to bind to.
    typedef hyperbolic::hyperbolic_space<N,float_type> space_type;
    //! The type of points in the lattice
    typedef zeta_lattice_point<N,float_type,integer_type,_Id>
        lattice_point_type;
    //! The type of the zeta lattice underlying the point
    typedef typename lattice_point_type::zeta_array_type point_array_type;
    //! The lattice iterator
    typedef iterators::zeta_lattice_iterator<self_type> lattice_iterator;
    //! The object representing a bound on the lattice
    typedef detail::lattice_bound<self_type::zeta_size,float_type,integer_type,_Id>
        bound_type;
    
    // Pull up the zeta typedef
    typedef typename self_type::zeta_type zeta_type;
    
  private:
    //! Integer parameter type;
    typedef typename boost::call_traits<integer_type>::param_type int_param_type;
    
  public:
    /**
     *  @brief  Default constructor
     */
    zeta_lattice( ): untransform_(1.0,0.0) { }
    
    
    /**
     *  @brief  Space constructor
     *  @param  space The space to set the lattice to
     */
    zeta_lattice( const space_type& __space ) : untransform_(1.0,0.0)
    { set_space(__space); }
    
    
    /**
     *  @brief  Set the lattice space
     *  @param  space The space to set
     */
    void set_space( const space_type& __space ) {
      // Copy the space's regions into the region lattices
      utility::algorithms::copy_member<region_type&, zeta_type>
          ( __space.zeta_ref_begin(), __space.zeta_ref_end(),
            self_type::zeta_ref_begin(), &zeta_type::original_region
          );
    }
    
    
    /**
     *  @brief  Set the lattice transform
     *  @param  transform The transform value to set
     */
    void set_transform( const complex_type& __transform ) {
      typedef utility::math::is_zero<float_type> _iz;
      assert( !(_iz()(__transform.real()) && _iz()(__transform.imag() )) );
      
      // Copy the transform into the region lattices
      untransform_ = complex_type(1.0,0.0) / __transform;
      utility::algorithms::fill_member<complex_type&, zeta_type>
          ( self_type::zeta_ref_begin(), self_type::zeta_ref_end(),
            &zeta_type::transform, __transform );
    }
    
    
    /**
     *  @brief  Intitialize the interval latice
     *  @note   Call this after changing the space, bound or transform to
     *          revalidate the lattice.
     */
    void initialize( void ) {
      // If the lattice is empty don't do anything
      if ( self_type::zeta_size == 0 ) return;
      
      lbound_.initialize();
      
      // Copy the bound into the region lattices
      utility::algorithms::fill_member<float_type&, zeta_type>
          ( self_type::zeta_ref_begin(), self_type::zeta_ref_end(), &zeta_type::bound,
            lbound_( self_type::zeta_size - 1 )
          );
      
      // Begin by initializing the region lattices
      typedef region_lattice_type _rlt;
      std::for_each( self_type::zeta_ref_begin(), self_type::zeta_ref_end(),
                     std::mem_fun_ref<void, _rlt>( &_rlt::initialize )
                   );
      
      // Now that the region lattices are ready the size, start and end points
      // can be calculated.
      typedef typename self_type::const_zeta_ref_iterator _iter;
      typedef typename lattice_point_type::zeta_ref_iterator _lpiter;
      
      // Set all the iterators to the beginning of their respective containers
      _iter it = self_type::zeta_ref_begin();
      _lpiter startIt = start_.zeta_ref_begin();
      _lpiter stopIt = stop_.zeta_ref_begin();
      
      // start_ Is pointwise a copy of all the starting values
      // stop_  Is the first point outside the lattice, so the first n-1
      //        coordinates are start values, and the nth is the stop value
      for( ; it != self_type::zeta_ref_end() - 1; ++it, ++startIt, ++stopIt ) {
        *stopIt = *startIt = it->start();
      }
      
      *startIt = it->start();
      *stopIt = it->stop();
    }
    
    
    /**
     *  @brief  Get a constant reference to the lattice bound.
     *  @return A constant reference to the lattice bound.
     */
    const bound_type& bound( void ) const { return lbound_; }
    
    /**
     *  @brief  Get a constant reference to the lattice bound.
     *  @return A constant reference to the lattice bound.
     */
    bound_type& bound( void ) { return lbound_; }
    
    
    /**
     *  @brief  Get a constant reference to the untransform factor
     *  @return A constant reference to the untransform factor
     */
    const complex_type& untransform( void ) const { return untransform_; }
    
    
    /**
     *  @brief  Compute the size of the lattice
     *  @return The size of the lattice
     *  @note   This actually computes the size of the lattice when it is
     *          called. This is necessary as the size can change during
     *          iteeffectn, but makes calling this function slightly more
     *          expensive than one might like.
     */
    integer_type size( void ) const {
      typedef typename self_type::const_zeta_ref_iterator _iter;
      
      integer_type s(1);
      
      if ( self_type::zeta_size == 0 ) {
        s = 0;
      }
      else {
        for( _iter it = self_type::zeta_ref_begin(); it != self_type::zeta_ref_end();
             ++it)
        { s *= it->size(); }
      }
      return s;
    }
    
    /**
     *  @brief  Get a constant reference to the start point
     *  @return A constant reference to the start point
     */
    const lattice_point_type& start( void ) const { return start_; }
    
    /**
     *  @brief  Get a constant reference to the stop point
     *  @return A constant reference to the stop point
     */
    const lattice_point_type& stop( void ) const { return stop_; }
    
    
    /**
     *  @brief Determine if the lattice is empty
     */
    bool empty ( void ) const { return bool( size() == 0 ); }
    
    
    /**
     *  @brief   Get an iterator pointing to the begining of the lattice
     *  @param   point The point which will act as the virtual lattice point
     *  @note    To iterate through the virtual space in the most efficient
     *           way it is neccessary to modify the bounds of the underlying
     *           region lattices during the iteeffectn process. This means that
     *           at the end of the iteeffectn the stop point is no longer valid.
     *  @warning Only one point can be iterated at a time because to
     *           efficiently iterate, the lattice itself must be modified (in a
     *           non-destructive way). Concurrently iterating more than one
     *           point leads to undefined behaviour.
     */
    lattice_iterator lattice_begin( point_array_type& __point ) {
      // Since the space is virtual we must assign a starting value
      __point = start_;
      lattice_iterator it( *this, __point );
      // The iterator needs to do some computation before it can be iterated
      // this means binding the coordinates of the point to the sub-iterators.
      it.initialize();
      return it;
    }
    
    
    /**
     *  @brief  Determine if the point which is currently being iterated is
     *          contained in the original space + bound;
     */
    bool validate( void ) const { return lbound_.validate(); }
    
    
  private:
    //! The transformation factor
    complex_type untransform_;
    //! The first point in the lattice
    lattice_point_type start_;
    //! One after the last point in the lattice
    lattice_point_type stop_;
    //! The bounds on the reference regions
    bound_type lbound_;
};

}// lattice
}// geometric
}// structure
}// sg
#endif
