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
 *  @file     zeta_lattice_point.hpp
 *  @brief    An inline header file for zeta_lattice_point.
 *  @note     Include this file to create points in a zeta lattice
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-09
 */

#ifndef _SG_ZETA_LATTICE_POINT_H
#define _SG_ZETA_LATTICE_POINT_H 1

// Local includes
#include "structure/geometric/detail/zeta_array.hpp"
#include "structure/numerical/accessors/iq_number_accessor.hpp"


namespace sg
{
namespace structure
{
namespace geometric
{
namespace lattice
{
/**
 *  @class    zeta_lattice_point
 *  @brief    A class representing a point in an integral lattice
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-09
 *  @param    N The hyperbolic dimension of the lattice
 *  @param    _Float A floating point type; defaults to \b double
 *  @param    _Integer The integer type which iterates over the interval
 *  @param    _Id is the imaginary quadratic field the lattice exists in
 */
template <std::size_t N, class _Float = double, class _Integer = long,
          std::size_t _Id = 0>
  class zeta_lattice_point
  : public geometric::detail::zeta_array< N,
                numerical::accessor::iq_number_accessor< _Float,_Integer,_Id >
           >
{
  public:
    //! Object type
    typedef zeta_lattice_point<N,_Float,_Integer,_Id> self_type;
    //! Float type
    typedef _Float float_type;
    //! The integer type
    typedef _Integer integer_type;
};

}// lattice
}// geometric
}// structure
}// sg
#endif
