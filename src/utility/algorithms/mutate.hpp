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
 *  @file     mutate.hpp
 *  @brief    This is a header implementation file for mutate.
 *  @note     Include this file to mutate a container.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-05
 *  
 *  \c std::for_each can only mutate one container and \c std::transform
 *  works by assignment, so is not particularly efficient for arrays.
 *  \c mutate can mutate up to two input ranges simultaneously.
 */


#ifndef _SG_MUTATE_H
#define _SG_MUTATE_H 1


namespace sg
{
namespace utility
{
namespace algorithms
{
/**
 *  @brief  Perform an operation on corresponding elements of two sequences.
 *  @param  first1     An input iterator.
 *  @param  last1      An input iterator.
 *  @param  first2     An input iterator.
 *  @param  mutator    A binary operator which mutates its inputs.
 *  @return The return value of the binary operation
 *
 *  Applies the operator to the corresponding elements in the two
 *  input ranges.
 *  Evaluates @p binary_op(*(first1+N),*(first2+N)) for each
 *  @c N in the range @p [0,last1-first1).
 *
 */
template<typename _InputIterator1, typename _InputIterator2,
         typename _BinaryOperation> inline
  _BinaryOperation
  mutate( _InputIterator1 __first1, _InputIterator1 __last1,
          _InputIterator2 __first2, _BinaryOperation __mutator)
{
  for ( ; __first1 != __last1; ++__first1, ++__first2 ) {
    __mutator(*__first1, *__first2);
  }
  
  return __mutator;
}

}// algorithms
}// utility
}// sg
#endif
