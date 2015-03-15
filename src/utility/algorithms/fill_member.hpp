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
 *  @file     fill_member.hpp
 *  @brief    An inline header file for the \c fill_member algorithm.
 *  @note     Allows a member function to be filled with a value
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-26
 */


#ifndef _SG_FILL_MEMBER_H
#define _SG_FILL_MEMBER_H 1

// Global includes
#include <boost/call_traits.hpp>
#include <boost/type_traits/remove_reference.hpp>

namespace sg
{
namespace utility
{
namespace algorithms
{
/**
 *  @brief Fills the range [first.(*_pf),last(*_pf)) with copies of value.
 *  @param  _Ret Return type of a member function.
 *  @param  _Tp The object type that the member function is called on.
 *  @param  first  A forward iterator.
 *  @param  last   A forward iterator.
 *  @param  _Ret(_Tp::*__pf)() A function pointer to a member function
 *  @param  value  A reference-to-const of arbitrary type.
 *  @return   Nothing.
 *
 *  This function fills a range with copies of the same value after making
 *  a call to a non-constant member function returning a non-constant
 *  reference. This allows the setting of data members.
 */
template<class _Ret, class _Tp,typename _InputIterator> inline
  void
  fill_member( _InputIterator __first, _InputIterator __last,
               _Ret (_Tp::*__pf)(),
               typename boost::call_traits<
                   typename boost::remove_reference<_Ret>::type
                   >::param_type __value )
{
  for ( ; __first != __last; ++__first ) { ( (*__first).*__pf )() = __value; }
}

}// algorithms
}// utility
}// sg
#endif
