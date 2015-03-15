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
 *  @file     copy_member.hpp
 *  @brief    An inline header file for the \c copy_member algorithm.
 *  @note     Allows a container to be copied to a data member (accessed
 *            via a function call) of another container's coordinates.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-25
 */


#ifndef _SG_COPY_MEMBER_H
#define _SG_COPY_MEMBER_H 1

namespace sg
{
namespace utility
{
namespace algorithms
{
/**
 *  @brief  Copies the range [first,last) to result.*(pf)()
 *  @param  _Ret Return type of a member function.
 *  @param  _Tp The object type that the member function is called on.
 *  @param  _Ret(_Tp::*__pf)() A function pointer to a member function
 *  @param  first  A forward iterator.
 *  @param  last   A forward iterator.
 *  @param  result An output iterator.
 *  @return The result output iterator.
 */
template< class _Ret, class _Tp, typename _InputIterator,
          typename _OutputIterator> inline
  _OutputIterator
  copy_member( _InputIterator __first, _InputIterator __last,
               _OutputIterator __result, _Ret (_Tp::*__pf)() )
{
  for ( ; __first != __last; ++__first, ++__result ) {
    ( (*__result).*__pf )() = *__first;
  }
  
  return __result;
}

}// algorithms
}// utility
}// sg
#endif
