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
 *  @file     heisenberg_incrementor.hpp
 *  @brief    An inline header file for the \c heisenberg_incrementor class.
 *  @note     Heisenberg incrementor efficiently increments the Heisenberg
 *            component of a \c hyperbolic_point
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-24
 */

#ifndef _SG_HEISENBERG_INCREMENTOR_H
#define _SG_HEISENBERG_INCREMENTOR_H 1

// Local includes
#include "structure/geometric/hyperbolic/incrementor/accessors/r_incrementor_accessor.hpp"
#include "structure/geometric/hyperbolic/incrementor/accessors/zeta_incrementor_accessor.hpp"
#include "structure/geometric/detail/heisenberg_structure.hpp"


namespace sg
{
namespace structure
{
namespace geometric
{
namespace hyperbolic
{
namespace incrementor
{
/**
 *  @class    heisenberg_incrementor
 *  @brief    A class for efficiently incrementing the heisenberg component of
 *            a hyperbolic point.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-24
 *  @param    _Point The \c hyperbolic_point class to increment
 */
template <class _Point> class heisenberg_incrementor
  : public geometric::detail::heisenberg_structure< _Point::dimension_size,
                               accessor::zeta_incrementor_accessor<_Point>,
                               accessor::r_incrementor_accessor<_Point>
           >
{
  public:
    //! Object type
    typedef heisenberg_incrementor<_Point> self_type;
    //! Point type
    typedef _Point point_type;
};

}// incrementor
}// hyperbolic
}// geometric
}// structure
}// sg
#endif
