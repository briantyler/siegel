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
 *  @file     make_siegel.hpp
 *  @brief    This is a header implementation file for make_siegel.
 *  @note     Include this file to make a hyperbolic space a Siegel container.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-13
 */


#ifndef _SG_MAKE_SIEGEL_H
#define _SG_MAKE_SIEGEL_H 1

// GLOBAL INCLUDES
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <functional>

#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_same.hpp>


namespace sg
{
namespace geometry
{
namespace algorithms
{
namespace make_siegel_helpers
{
/**
 *  @struct  siegel_1_init
 *  @brief   Initializes the zeta subspace when \f$ d = -1 \f$
 */
template <class _Region>
    struct siegel_1_init
  : public std::unary_function< _Region&, void >
{
  void operator( ) ( _Region& __region ) const {
    __region.real().lower() = __region.imag().lower() = 0.0;
    __region.real().upper() = __region.imag().upper() = 0.5;
  }
};


/**
 *  @struct  siegel_3_init
 *  @brief   Initializes the zeta subspace when \f$ d = -3 \f$
 */
template <class _Region, class _IQField>
    struct siegel_3_init
  : public std::unary_function< _Region&, void >
{
  void operator( ) ( _Region& __region ) const {
    const _IQField& field = _IQField::instance();
    
    __region.real().lower() = __region.imag().lower() = 0.0;
    __region.real().upper() = 0.5;
    __region.imag().upper() = field.sqrt_generator() * 0.25;
  }
};


/**
 *  @struct  siegel_congruent_init
 *  @brief   Initializes the zeta subspace when \f$ d \equiv 1 \mod 4 \f$
 */
template <class _Region, class _IQField>
    struct siegel_congruent_init
  : public std::unary_function< _Region&, void >
{
  void operator( ) ( _Region& __region ) const {
    const _IQField& field = _IQField::instance();
    
    __region.real().lower()= -( __region.real().upper() = 0.5 );
    __region.imag().lower() = 0.0;
    __region.imag().upper() = field.sqrt_generator() * 0.25;
  }
};


/**
 *  @struct  siegel_non_congruent_init
 *  @brief   Initializes the zeta subspace when \f$ d \not \equiv 1
 *           \mod 4 \f$
 */
template <class _Region, class _IQField>
    struct siegel_non_congruent_init
  : public std::unary_function< _Region&, void >
{
  void operator( ) ( _Region& __region ) const {
    const _IQField& field = _IQField::instance();
    
    __region.real().lower()= -( __region.real().upper() = 0.5 );
    __region.imag().lower() = 0.0;
    __region.imag().upper() = field.sqrt_generator() * 0.5;
  }
};

}// make_siegel_helpers


/**
 *  @brief   Initializes the hyperbolic space to the canonical Siegel container
 *  @param   _IQField is an imaginary quardatic field type. The real type
 *           should be the same as that used to generate the \c _HyperbolicSpace.
 *           Note that this <em>is</em> the right way to do this because there
 *           is no way of knowing the integer type that is in use.
 *  @param   _HyperbolicSpace The hyperbolic space type, this will be deduced
 *           by the compiler.
 *  @param   space the Hyperbolic space.
 */
template <class _IQField, class _HyperbolicSpace>
    void make_siegel( _HyperbolicSpace& __space ) {
  
  using std::for_each;
  using namespace make_siegel_helpers;
  typedef typename _HyperbolicSpace::zeta_type region_type;
  typedef _IQField field_type;
  
  // Make sure that the field and the __space have compatible types
  BOOST_MPL_ASSERT( (
                        boost::is_same<
                          typename _HyperbolicSpace::float_type
                        ,
                          typename field_type::float_type
                        >
                      )
                    );
  
  const field_type& field = field_type::instance();
  
  // Set the vertical sub-space
  __space.r().ref().lower() = -(   __space.r().ref().upper()
                                 = field.sqrt_generator() * 0.5
                               );
  
  // Set the height sub-space.
  __space.height().ref().lower() = 0.0;
  __space.height().ref().upper() = 2.0;
  
  // If the dimension is one then there is nothing more to do.
  if( _HyperbolicSpace::dimension_size == 1 ) return;
  
  if( field.generator( ) == -1 ) { // Special Case
    __space.real_ref_at(0).lower() = -( __space.real_ref_at(0).upper() = 1.0 );
    __space.imag_ref_at(0).lower() = 0.0;
    __space.imag_ref_at(0).upper() = 0.5;
    
    if( _HyperbolicSpace::dimension_size == 3 ) {
      __space.real_ref_at(1).lower() = -0.5;
      __space.real_ref_at(1).upper() = __space.imag_ref_at(1).upper() = 0.5;
      __space.imag_ref_at(1).lower() = 0.0;
    }
    else if( _HyperbolicSpace::dimension_size > 3 ) {
      __space.real_ref_at(1).lower() = __space.imag_ref_at(1).lower() = -0.5;
      __space.real_ref_at(1).upper() = __space.imag_ref_at(1).upper() = 0.5;
      
      for_each( __space.zeta_ref_begin( ) + 2, __space.zeta_ref_end( ),
                siegel_1_init<region_type>( )
              );
    }
  }
  else if ( field.is_congruent( ) ) {
    if( _HyperbolicSpace::dimension_size == 2 ) {
      __space.real_ref_at(0).lower() = -0.5;
      __space.real_ref_at(0).upper() = 0.5;
      __space.imag_ref_at(0).lower() = 0.0;
      __space.imag_ref_at(0).upper() = field.sqrt_generator() * 0.25;
    }
    else {
      if( field.generator( ) == -3 ) { // Special Case
        __space.real_ref_at(0).lower() = -0.5;
        __space.real_ref_at(0).upper() = 0.5;
        __space.imag_ref_at(0).lower() = -field.sqrt_generator() * 0.25;
        __space.imag_ref_at(0).upper() = field.sqrt_generator() * 0.25;
        
        for_each( __space.zeta_ref_begin( ) + 1, __space.zeta_ref_end( ),
                  siegel_3_init<region_type,field_type>( )
                );
      }
      else {
        for_each( __space.zeta_ref_begin( ),
                  __space.zeta_ref_end( ),
                  siegel_congruent_init<region_type,field_type>( )
                );
      }
    }
  }
  else { // Not congruent, not == -1
    __space.real_ref_at(0).lower() = -( __space.real_ref_at(0).upper() = 1.0 );
    
    if( _HyperbolicSpace::dimension_size == 3 ) { // Special Case
      __space.real_ref_at(1).lower() = -(   __space.real_ref_at(1).upper()
                                          = 0.5
                                        );
      __space.imag_ref_at(1).lower() = -(   __space.imag_ref_at(0).upper()
                                          = __space.imag_ref_at(1).upper()
                                          = field.sqrt_generator() * 0.5
                                        );
      __space.imag_ref_at(0).lower() = 0.0;
    }
    else {
      __space.imag_ref_at(0).lower() = 0.0;
      __space.imag_ref_at(0).upper() = field.sqrt_generator() * 0.5;
      
      for_each( __space.zeta_ref_begin( ) + 1, __space.zeta_ref_end( ),
                siegel_non_congruent_init<region_type,field_type>( )
              );
    }
  }
  
  __space.initialize();
}

}// algorithms
}// geometry
}// sg
#endif
