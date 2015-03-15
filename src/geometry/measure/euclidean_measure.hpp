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
 *  @file     euclidean_measure.hpp
 *  @brief    An inline header file for the \c euclidean_measure class.
 *  @note     Include this file to compute the Euclidean measure of a product
 *            of intervals on the real line.
 *  @author   Brian Tyler
 *  @version  1.1
 *  @date     2008-02-14
 */


#ifndef _SG_EUCLIDEAN_MEASURE_H
#define _SG_EUCLIDEAN_MEASURE_H 1

// Global includes
#include <functional>
#include <boost/call_traits.hpp>


namespace sg
{
namespace geometry
{
namespace measure
{
/**
 *  @class    euclidean_measure
 *  @brief    A function object that computes the euclidean measure of an
 *            interval on the real line.
 *  @note     The cumulative result of multiple applications is recorded
 *            internally so by applying the functor repeatedly to a number of
 *            intervals one computes the measure of the prouct of these
 *            intervals.
 *  @author   Brian Tyler
 *  @version  1.1
 *  @date     2008-02-14
 *  @param    Float The floating point type of the interval to measure.
 *  @note     v1.1 Split out the region and interval operators to make a true
 *            binary functor
 */
template <class _Float = double> struct euclidean_measure
  : public std::binary_function< typename boost::call_traits<_Float>::param_type,
                                 typename boost::call_traits<_Float>::param_type,
                                 _Float
                               >
{
  //! Object type
  typedef euclidean_measure<_Float> self_type;
  //! Real type
  typedef _Float float_type;
  //! The parameter type
  typedef typename boost::call_traits<_Float>::param_type param_type;
  
  /**
   *  @brief  Default constructor
   */
  euclidean_measure( ) : measure_(1) { };
  
  
  /**
   *  @brief  Computes the distance between two real numbers
   *  @param  lhs The lower bound of an interval.
   *  @param  rhs The upper bound of an interval.
   *  @return The distance between \c lhs and \c rhs.
   */
  float_type operator( )( param_type __lhs, param_type __rhs ) {
    float_type output( __lhs > __rhs ? __lhs - __rhs : __rhs - __lhs );
    measure_ *= output ;
    return output;
  }
  
  /**
   *  @brief  Gets the cumulative value of measurements.
   *  @return Cumulative value of measurements since last reset.
   */
  operator float_type ( ) const { return measure_; }
  
  
  /**
   *  @brief  Resets the cumulative value of measurements.
   */
  void reset( void ) { measure_ = 1; }
  
  
  private:
    //! The cumulative measure of the intervals parsed to the function
    float_type measure_;
};

}// measure
}// geometry
}// sg
#endif
