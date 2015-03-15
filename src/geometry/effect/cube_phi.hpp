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
 *  @file     cube_phi.hpp
 *  @brief    An inline header file for the \c cube_phi class.
 *  @note     Computes the maximum height attainable by the action of a cusp on
 *            a hypercube.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-06-03
 */


#ifndef _SG_CUBE_PHI_H
#define _SG_CUBE_PHI_H 1

// Global includes
#include <functional>
#include <cmath>

#include <boost/call_traits.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_same.hpp>

// Local includes
#include "utility/math/hermitian_inner_product.hpp"
#include "utility/math/is_zero.hpp"


namespace sg
{
namespace geometry
{
namespace effect
{
/**
 *  @class    cube_phi
 *  @brief    Computes the maximum height attainable by the action of a cusp on
 *            a hypercube.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-20
 *  @param    _Cube The hypercube type
 *  @param    _Cusp The cusp type
 */
template <class _Cube, class _Cusp> class cube_phi
: public std::binary_function< const _Cube&,
                               const _Cusp&,
                               typename _Cusp::float_type>
{
  // Ensure the floating point types are the same at compile time
  BOOST_MPL_ASSERT(( boost::is_same< typename _Cube::float_type,
                                     typename _Cusp::float_type >
                   ));
  
  // Ensure the hyperbolic dimensions are equal at compile time
  BOOST_MPL_ASSERT_RELATION(  _Cube::dimension_size,
                              ==,
                              _Cusp::dimension_size
                            );
  
  public:
    //! The object type
    typedef cube_phi<_Cube,_Cusp> self_type;
    //! The Hypercube type
    typedef _Cube cube_type;
    //! The Cusp type
    typedef _Cusp cusp_type;
    //! The floating point type
    typedef typename cusp_type::float_type float_type;
    //! The complex type
    typedef typename cusp_type::complex_type complex_type;
    //! The Hypercube vertex type
    typedef typename cube_type::vertex_type vertex_type;
    
  private:
    //! Floating point param type
    typedef typename boost::call_traits<float_type>::param_type param_type;
    
  public:
    /**
     *  @brief  Compute the maximum height of a hypercube which can be raised
     *          by a cusp.
     *  @param  cube The hypercube to attempt to raise.
     *  @param  cusp The cusp to raise the hypercube with.
     *  @return The height at which \c cusp no longer has an upward effect on
     *          \c cube. If no upward action is possible at any height then
     *          \c -1.0 is returned.
     *  @note   The hypercube should be in a heisenberg slice at height zero.
     *          It would be easy to make a small modification to the algorithm
     *          so the height of the hypercube is unimportant, but there seems
     *          to be no benefit in introducing this inefficiency.
     *  @note   The algorithm fails early, in the sense that at soon as the
     *          cusp fails to raise one vertex, the algorithm returns the fail
     *          value.
     */
    float_type operator()
        ( const cube_type& __cube, const cusp_type& __cusp ) const
    {
      typename cube_type::const_iterator it = __cube.begin();
      
      // The heisenberg cube has at least one vertex (even in dimension zero)
      // so this is always valid.
      float_type output = sys_cusp_phi( *it, __cusp.point(), __cusp.threshold() );
      
      //std::cout << "cusp: " << output << "\n";
      // Negative value means fail
      if( output < 0.0 ) return output;
      
      // Iterate through all of the vertices to find the minimum
      ++it;
      for( ; it != __cube.end(); ++it ) {
        float_type tmp = sys_cusp_phi( *it, __cusp.point(), __cusp.threshold() );
        
        //std::cout << "cusp: " << tmp << "\n";
        if( tmp < output ) output = tmp;
        if( output < 0.0 ) break;
      }
      
      // At this stage the cusp has been shown to have an upward effect on the
      // entire cube.
      return output;
    }
    
    
    /**
     *  @brief  Bound the maximum height attainable by any possible cusp for a
     *          given hypercube.
     */
    float_type operator()
        ( const cube_type& __cube, param_type __threshold ) const
    {
      vertex_type midpoint( __cube.midpoint() );
      
      typename cube_type::const_iterator it = __cube.begin();
      float_type output = sys_cusp_phi( *it, midpoint, __threshold );
      if( output < 0.0 ) return output;
      
      ++it;
      for( ; it != __cube.end(); ++it ) {
        float_type tmp = sys_cusp_phi( *it, midpoint, __threshold );
        
        if( tmp < output ) output = tmp;
        if( output < 0.0 ) break;
      }
      
      // At this stage the cusp has been shown to have an upward effect on the
      // entire cube.
      return output;
    }
    
    
  private:
    /**
     *  @brief  Compute the maximum height of a vertex such that it can be
     *          raised by a cusp.
     *  @param  vertex The vertex to try and raise.
     *  @param  cusp The cusp to raise the vertex with.
     *  @return The height at which \c cusp no longer has an upward effect on
     *          \c vertex. If no upward action is possible at any height then
     *          \c -1.0 is returned.
     */
    const float_type sys_cusp_phi
        ( const vertex_type& __vertex, const vertex_type& __point,
          const param_type __threshold ) const
    {
      using utility::math::hermitian_inner_product;
      typedef utility::math::is_zero<float_type> _iz;
      
      static const float_type fail = -1.0;
      
      // Compute the standard hermitian inner product of the zeta components
      // of the vertex and cusp point.
      complex_type hip
          = hermitian_inner_product( __vertex.zeta_ref_begin(),
                                     __vertex.zeta_ref_end(),
                                     __point.zeta_ref_begin()
                                   );
      
      // This is the imaginary part of the effect function
      float_type output = __vertex.r_common() - __point.r_common()
                          + hip.imag();
      output *= output;
      
      
      // Since we are taking square roots then __threshold - output must be
      // positive, if not then the cusp cannot increase the height of the cube
      // so fail.
      if( output > __threshold ) {
        output = fail;
        return output;
      }
      
      // Now add on the real contribution of the effect function.
      // Note that the real dependent coordinates are half the negative
      // quadratic form of the zeta component.
      output = 2.0 * (   ::sqrt( __threshold - output ) + hip.real()
                       + __vertex.dependent().real_common()
                       + __vertex.height_ref() * 0.5
                       + __point.dependent().real_common()
                     );
      
      // A negative value means failure
      // We also take a value of very close to zero as a failure as we have no
      // confidence in it really being positive.
      if( output < 0.0 || _iz()( output ) ) {
        output = fail;
      }
      
      return output;
    }
};

}// effect
}// geometry
}// sg
#endif
