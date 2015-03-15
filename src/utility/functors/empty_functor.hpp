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
 *  @file     empty_functor.hpp
 *  @brief    An inline header file for the \c empty_functor class.
 *  @note     A dummy functor which does nothing that can be used to terminate
 *            recursive meta-algorithms.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-05-21
 */


#ifndef _SG_EMPTY_FUNCTOR_H
#define _SG_EMPTY_FUNCTOR_H 1


namespace sg
{
namespace utility
{
namespace functors
{
/**
 *  @struct   empty_functor
 *  @brief    A functor which does nothing
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-05-21
 *  @param  ... Matches any number of POD types.
 *  @note   This functor can be used to terminate a recursive meta-algorithm.
 *          The calling functor must take only POD types; this means that
 *          arguments which would normally be passed by refernce must be passed
 *          as pointers, except for references to POD types.
 */
struct empty_functor { void operator() ( ... ) const { } };

}// functors
}// utility
}// sg
#endif
