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
 *  @file     measure.hpp
 *  @brief    This is a header implementation file for measure.
 *  @note     Include this file to measure a container.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-13
 *  
 *  This is basically \c std::for_each with a wrapper around it
 */


#ifndef _SG_MEASURE_H
#define _SG_MEASURE_H 1

// Global includes
#include <algorithm>


namespace sg
{
namespace geometry
{
namespace measure
{
/**
 *  @brief  Measures a container using a measuring functor
 *  @param  first1     An input iterator.
 *  @param  last1      An input iterator.
 *  @param  measurer   A unary operator which measures its inputs.
 *
 */
template<class _InputIterator, class _UnaryOperation> inline
    typename _UnaryOperation::result_type measure
    (_InputIterator __first, _InputIterator __last, _UnaryOperation __measurer)
{
  typedef typename _UnaryOperation::result_type result_type;
  return result_type( std::for_each( __first, __last, __measurer ) );
}

}// measure
}// utility
}// sg
#endif
