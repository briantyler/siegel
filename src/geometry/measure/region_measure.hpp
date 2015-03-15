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
 *  @file     region_measure.hpp
 *  @brief    An inline header file for the \c region_measure class.
 *  @note     Include this file to compute the Euclidean measure of a
 *            complex region, or a product of regions.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-12
 */


#ifndef _SG_REGION_MEASURE_H
#define _SG_REGION_MEASURE_H 1

// Global includes
#include <functional>

// Local includes
#include "geometry/measure/interval_measure.hpp"


namespace sg
{
namespace geometry
{
namespace measure
{
/**
 *  @class    region_measure
 *  @brief    A function object that computes the euclidean measure of a
 *            region in complex space.
 *  @see      euclidean_measure
 *  @author   Brian Tyler
 *  @version  1.1
 *  @date     2008-03-20
 *  @param    Region The type of the region to measure.
 */
template <class _Region> struct region_measure
  : public std::unary_function<const _Region&,typename _Region::float_type>
{
  //! Object type
  typedef region_measure<_Region> self_type;
  //! Region type
  typedef _Region region_type;
  //! Floating point type
  typedef typename _Region::float_type float_type;
  //! Interval type
  typedef typename _Region::interval_type interval_type;
  
  
  /**
   *  @brief  Computes the area of a complex region
   *  @param  region A \c complex_region.
   *  @return The area defined by the \c complex_region
   */
  float_type operator( )( const region_type& __region ) {
    measure_( __region.real() );
    return float_type( measure_( __region.imag() ) );
  }
  
  
  /**
   *  @brief  Gets the cumulative value of measurements.
   *  @return Cumulative value of measurements since last reset.
   */
  operator float_type ( ) const { return static_cast<float_type>(measure_); }
  
  
  /**
   *  @brief  Resets the cumulative value of measurements.
   */
  void reset( void ) { measure_.reset(); }
  
  
  private:
    //! The cumulative measure of the intervals parsed to the function
    interval_measure<interval_type> measure_;
};

}// measure
}// geometry
}// sg
#endif
