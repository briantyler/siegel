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
 *  @file     hyperbolic_space.hpp
 *  @brief    An inline header file for \c hyperbolic_space.
 *  @note     Include this file to create (hyper)cubic compact regions of
 *            hyperbolic space.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-12
 */


#ifndef _SG_HYPERBOLIC_SPACE_H
#define _SG_HYPERBOLIC_SPACE_H 1

// Global includes
#include <algorithm>
#include <functional>

// Local includes
#include "structure/geometric/detail/hyperbolic_structure.hpp"
#include "structure/geometric/euclidean/accessors/interval_accessor.hpp"
#include "structure/geometric/euclidean/accessors/region_accessor.hpp"
#include "structure/geometric/hyperbolic/hyperbolic_point.hpp"


namespace sg
{
namespace structure
{
namespace geometric
{
namespace hyperbolic
{
/**
 *  @class    hyperbolic_space
 *  @brief    A class representing a compact region of hyperbolic space.
 *  @author   Brian Tyler
 *  @version  1.2
 *  @date     2008-03-12
 *  @param    N The dimension of the space
 *  @param    _Float A floating point type; defaults to \c double
 *  @note     v1.2 Moved over to general hyperbolic structure
 */
template <std::size_t N, class _Float = double>
    class hyperbolic_space
  :public geometric::detail::hyperbolic_structure<
                  N,
                  euclidean::accessor::region_accessor<_Float>,
                  euclidean::accessor::interval_accessor<_Float>,
                  euclidean::accessor::interval_accessor<_Float>
          >
{
  public:
    //! Object type
    typedef hyperbolic_space<N,_Float> self_type;
    //! Float type
    typedef _Float float_type;
    //! The type of points in the space
    typedef hyperbolic_point<N,float_type> point_type;
    
    
    /**
     *  @brief  Initialize the space so that all regions and intervals are valid
     */
    void initialize( void ) {
      typename self_type::iterator it = self_type::begin();
      for( ; it != self_type::hyperbolic_end(); ++it ) { it->initialize(); }
    }
    
    
    /**
     *  @brief  Determines if the space contains a point
     *  @param  value The point in hyperbolic space
     *  @return True if \c value is contained in the space
     *          ( i.e. each coordinate is in its respective interval )
     */
    bool contains( const point_type& __point ) const {
      typedef typename self_type::const_iterator _spc_iter;
      typedef typename point_type::const_iterator _pnt_iter;
      
      // RVO
      bool b = false;
      
      _spc_iter spcIt = self_type::begin();
      _pnt_iter pntIt = __point.begin();
      
      // Check all the intervals until one is found that does not contain
      // the point's coordinate. If this happens then the point is not in the
      // space.
      for( ; spcIt != self_type::hyperbolic_end(); ++spcIt, ++pntIt ) {
        if( !spcIt->contains( *pntIt ) ) { return b; }
      }
      
      b = true;
      return b;
    }
    
    
    /**
     *  @brief  Tex string representation of the complex hyperbolic space that
     *          this space lives in
     */
    std::string tex_complex_hyperbolic( void ) const {
      std::stringstream ss;
      ss << "\\mathbb{H}_\\mathbb{C}^" << self_type::dimension_size;
      return ss.str();
    }
};

}// hyperbolic
}// geometric
}// structure
}// sg
#endif
