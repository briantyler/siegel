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
 *  @file     interval_measure.hpp
 *  @brief    An inline header file for the \c interval_measure class.
 *  @note     Include this file to compute the Euclidean measure of a
 *            \c real_interval or a product of intervals
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-20
 */


#ifndef _SG_INTERVAL_MEASURE_H
#define _SG_INTERVAL_MEASURE_H 1

// Global includes
#include <functional>

// Local includes
#include "geometry/measure/euclidean_measure.hpp"


namespace sg
{
namespace geometry
{
namespace measure
{
/**
 *  @class    interval_measure
 *  @brief    A function object that computes the euclidean measure of an
 *            interval on the real line.
 *  @see      euclidean_measure
 *  @author   Brian Tyler
 *  @version  1.1
 *  @param    Interval The type of the interval to measure.
 *  @date     2008-03-20
 */
template <class _Interval> struct interval_measure
  :public std::unary_function<const _Interval&,typename _Interval::float_type>
{
  //! Object type
  typedef interval_measure<_Interval> self_type;
  //! Interval type
  typedef _Interval interval_type;
  //! The floating point type
  typedef typename _Interval::float_type float_type;
  
  
  /**
   *  @brief  Computes the length of an interval
   *  @param  interval An interval on the real line
   *  @return The measure of the interval.
   */
  float_type operator( )( const interval_type& __interval ) {
    return float_type( measure_( __interval.lower( ), __interval.upper( ) ) );
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
    euclidean_measure<float_type> measure_;
};

}// measure
}// geometry
}// sg
#endif
