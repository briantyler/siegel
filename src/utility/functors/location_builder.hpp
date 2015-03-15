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
 *  @file     location_builder.hpp
 *  @brief    This is a header implementation file for location_builder.
 *  @note     Include this file to sequentially construct indices from a
 *            location and a sequence of resolutions
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-25
 */


#ifndef _SG_LOCATION_BUILDER_H
#define _SG_LOCATION_BUILDER_H 1

// Global includes
#include <functional>

#include <boost/call_traits.hpp>


namespace sg
{
namespace utility
{
namespace functors
{
/**
 *  @class    location_builder
 *  @brief    A function object to sequentially construct indices from a
 *            location and a sequence of resolutions
 *  @param    _Integer an integer type
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-25
 */
template <class _Integer> class location_builder
  : public std::unary_function<
                typename boost::call_traits<_Integer>::param_type,
                _Integer
           >
{
  public:
    //! Object type
    typedef location_builder<_Integer> self_type;
    //! Value type
    typedef _Integer integer_type;
    
    
  private:
    //! Parameter type
    typedef typename boost::call_traits<integer_type>::param_type param_type;
    
    
  public:
    /**
     *  @brief  Location constructor
     *  @param  loc The location of the object we want to construct
     */
    location_builder( param_type __loc ) : loc_(__loc) { }
    
    
    /**
     *  @brief  Returns the sequential indices required to construct the object.
     *  @param  resolution The resolution of the current dimension
     *  @return The corresponding index.
     */
    integer_type operator( )( param_type __resolution ) {
      integer_type loci( loc_ % __resolution );
      loc_ -= loci;
      loc_ /= __resolution;
      
      return loci;
    }
    
  private:
    //! The current location
    integer_type loc_;
};

}// functors
}// utility
}// sg
#endif
