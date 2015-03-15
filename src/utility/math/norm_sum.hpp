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
 *  @file     norm_sum.hpp
 *  @brief    This is a header implementation file for \c norm_sum.
 *  @note     Include this file to sum the norms of values in a sequence
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-11
 */

#ifndef _SG_NORM_SUM_H
#define _SG_NORM_SUM_H 1

// Global includes
#include<functional>


namespace sg
{
namespace utility
{
namespace math
{
/**
 *  @struct   norm_sum
 *  @brief    A unary function object that sums norms of values in a sequence
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-11
 *  @param    _Cx A complex type.
 *
 *  Suppose \f$ z, y \in \mathbb{C} \f$. then this function object returns
 *  \f$ z\overline{y} \f$.
 */
template <class _Cx>
    struct norm_sum
  : public std::unary_function<const _Cx&, void>
{
  //! Object type
  typedef norm_sum<_Cx> self_type;
  //! Complex Type
  typedef _Cx complex_type;
  //! Value type
  typedef typename _Cx::value_type float_type;
  
  
  /**
   *  @brief  Default constructor
   *  @note   Initialises the total sum to zero
   */
  norm_sum( ): sum_(0.0) { }
  
  
  /**
   *  @brief  Increments the current sum by the norm of \c cx
   *  @param  cx A complex number
   */
  void operator( ) ( const complex_type& __cx ) {
    sum_ += norm( __cx );
  }
  
  
  /**
   *  @brief  Get the norm sum of the sequence so far
   *  @return A copy of the total sum.
   */
  operator float_type ( ) { return sum_; }
  
  
  /**
   *  @brief  Resets the cumulative value of measurements.
   */
  void reset( void ) { sum_ = 0.0;}
  
  
  private:
    //! The norm sum of the sequence
    float_type sum_;
};

}// math
}// utility
}// sg
#endif
