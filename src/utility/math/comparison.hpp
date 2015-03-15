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
 *  @file     comparison.hpp
 *  @brief    This is a header implementation file for is_less
 *  @note     Safe(er) floating point comparison
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-17
 *  @note     Comparing floating point numbers is a dodgy business because small
 *            errors in computation lead to things not being equal that should
 *            be (for instance 1 * sqrt(2) / sqrt(2) != 1, which is a problem).
 *            These classes make these sorts of comparison much safer.
 *            For conversion to integers, use the rounding classes.
 */


#ifndef _SG_COMPARISON_H
#define _SG_COMPARISON_H 1

#include "utility/math/is_zero.hpp"
#include "utility/math/is_equal.hpp"
#include "utility/math/is_equal_cx.hpp"
#include "utility/math/is_less.hpp"
#include "utility/math/is_less_equal.hpp"
#include "utility/math/is_greater.hpp"
#include "utility/math/is_greater_equal.hpp"

#endif
