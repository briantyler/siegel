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
 *  @file     hyperbolic_base.hpp
 *  @brief    An inline header file for hyperbolic_base
 *  @note     Include this file and derive from \c hyperbolic_base create
 *            hyperbolic objects.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-19
 */


#ifndef _SG_HYPERBOLIC_BASE_H
#define _SG_HYPERBOLIC_BASE_H 1

// Global includes
#include <cstddef>

#include <boost/mpl/assert.hpp>
#include <boost/mpl/greater.hpp>
#include <boost/mpl/long.hpp>

// Local includes
#include "utility/math/pow.hpp"


namespace sg
{
namespace structure
{
namespace geometric
{
namespace detail
{
/**
 *  @struct   hyperbolic_base
 *  @brief    A struct containing dimension values for hyperbolic objects
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-19
 *  @param    N The dimension of the hyperbolic object
 *  @note     Hyperbolic objects of dimension < 1 cannot be constructed,
 *            attempting to do so will lead to a compile-time error.
 */
template <std::size_t N> struct hyperbolic_base
{
  // Ensure that the dimension is greater than zero.
  BOOST_MPL_ASSERT( ( boost::mpl::greater< boost::mpl::long_<N>,
                                           boost::mpl::long_<0>
                                          >
                    )
                  );
  
  //! The type of the object
  typedef hyperbolic_base<N> self_type;
  //! The size type for the hyperbolic object
  typedef std::size_t size_type;
  
  //! The value of the input parameter
  static const size_type value = N;
  
  
  /**
   *  @enum sizes
   *  @note This \b enum contains important constants for hyperbolic objects
   *        <ul>
   *          <li>\c dimension_size = the complex hyperbolic dimension of
   *              the space.</li>
   *          <li>\c zeta_size = the dimension of the free complex zeta
   *               subspace.</li>
   *          <li>\c heisenberg_size = the real dimension of the Heisenberg
   *              space at each height</li>
   *        </ul>
   */
  enum sizes {
    dimension_size = N,
    zeta_size = N - 1,
    zeta_real_size = 2*(N - 1),
    heisenberg_size = 2*N - 1,
    hyperbolic_size = 2*N,
    hypercube_size = utility::math::powi_ct<2,heisenberg_size>::value
  };

};

}// detail
}// geometric
}// structure
}// sg
#endif
