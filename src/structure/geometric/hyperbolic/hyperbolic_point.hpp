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
 *  @file     hyperbolic_point.hpp
 *  @brief    An inline header file for hyperbolic_point.
 *  @note     Include this file to create hyperbolic points
 *  @author   Brian Tyler
 *  @version  1.2
 *  @note     In version 1.2 moved over to \c hyperbolic_structure base
 *  @date     2008-01-31
 */

#ifndef _SG_HYPERBOLIC_POINT_H
#define _SG_HYPERBOLIC_POINT_H 1

// Global includes
#include <algorithm>
#include <functional>

#include <boost/operators.hpp>
#include <boost/call_traits.hpp>
#include <boost/iterator/reverse_iterator.hpp>

// Local includes
#include "structure/geometric/detail/hyperbolic_structure.hpp"
#include "structure/geometric/euclidean/accessors/complex_accessor.hpp"
#include "structure/geometric/euclidean/accessors/real_accessor.hpp"
#include "structure/geometric/hyperbolic/iterators/hyperbolic_point_iterator.hpp"
#include "structure/geometric/hyperbolic/incrementor/r_incrementor.hpp"
#include "structure/geometric/hyperbolic/incrementor/zeta_incrementor.hpp"
#include "utility/math/hermitian_quadratic_product.hpp"
#include "utility/math/is_equal.hpp"


namespace sg
{
namespace structure
{
namespace geometric
{
namespace hyperbolic
{
/**
 *  @class    hyperbolic_point
 *  @brief    A class representing a point in hyperbolic space
 *  @author   Brian Tyler
 *  @version  1.2
 *  @date     2008-01-11
 *  @param    N The dimension of the point
 *  @param    _Float A floating point type; defaults to \b double
 *  @note     v1.1: modified the type of zeta: \c vector -> \c boost::array
 *            This should be more efficient as the extra overhead
 *            required for dynamically sizing is removed. It also creates
 *            greater compile time safety since the dimension of the point is
 *            known. The slight disadvantage is that the application needs to
 *            be compiled for each dimension.
 *  @note     v1.2 Is a complete rewrite using the generalised hyperbolic
 *            structure as a base to build on. This dramatically simplifies the
 *            code needed for the point itself.
 */
template <std::size_t N, class _Float = double>
  class hyperbolic_point
  :public detail::hyperbolic_structure<
                  N,
                  euclidean::accessor::complex_accessor<_Float>,
                  euclidean::accessor::real_accessor<_Float>,
                  euclidean::accessor::real_accessor<_Float>
          >
  ,private boost::less_than_comparable< hyperbolic_point<N,_Float>,
             boost::additive< hyperbolic_point<N,_Float>,
               boost::multiplicative< hyperbolic_point<N,_Float>, _Float
           > > >
{
  public:
    //! Object type
    typedef hyperbolic_point<N,_Float> self_type;
    //! Float type
    typedef _Float float_type;
    
    // Incrementor friendship; gives write access to dependent variable
    friend class incrementor::r_incrementor<self_type>;
    friend class incrementor::zeta_incrementor<self_type>;
    
    // Pull up the base types which are needed
    typedef typename self_type::zeta_type zeta_type;
    typedef typename self_type::zeta_coordinate_type zeta_coordinate_type;
    typedef typename self_type::iterator iterator;
    typedef typename self_type::const_iterator const_iterator;
    
    //! A constant point iterator
    typedef iterators::hyperbolic_point_iterator<self_type> const_point_iterator;
    //! A reverse constant point iterator
    typedef boost::reverse_iterator<const_point_iterator>
        const_reverse_point_iterator;
    
  private:
    //! Floating point parameter type
    typedef typename boost::call_traits<float_type>::param_type param_type;
    
    
  public:
    /**
     *  @brief  Initializes the point by computing the dependent variable.
     */
    void initialize( void ) {
      dependent_.real_ref( ) = sys_compute_qf( ) - (0.5 * self_type::height_ref( ));
      dependent_.imag_ref( ) = self_type::r_ref( );
    }
    
    
    /**
     *  @brief  Sets the height of the point in an efficient way
     *  @param  height The new height of the point
     *  @note   There is no need to initialize the point after calling this
     *          function.
     */
    void set_height( param_type __height ) {
      dependent_.real_ref( ) += 0.5*( self_type::height_ref( ) - __height );
      self_type::height( ) = __height;
    }
    
    
    /**
     *  @brief  Sets the height of the point in an efficient way
     *  @param  height The new height of the point
     *  @note   There is no need to initialize the point after calling this
     *          function.
     */
    void set_r( param_type __r ) {
      dependent_.imag_ref( ) = __r;
      self_type::r( ) = __r;
    }
    
    
    /**
     *  @brief  Get the first coordinate under the projective representation.
     *  @return The first coordinate under the projective representation.
     *  @note   This is provided so that the point iterator has something
     *          concrete to return for the first coordinate.
     */
    static const zeta_type& first( void ) {
      static const zeta_type z(1.0,0.0);
      return z;
    }
    
    
    /**
     *  @brief  Constant access to the \c dependent component.
     *  @return Const reference to \c dependent component.
     *  @note   Full access is not given to the dependent variable as it
     *          should be modified via the \c zeta, \c r and \c height
     *          components.
     */
    const zeta_coordinate_type& dependent( void ) const { return dependent_; }
    
    
    /**
     *  @brief  Negates the heisenberg component
     *  @return A copy of this after incrementing
     */
    self_type operator-( void ) {
      self_type point;
      
      const_iterator meIt = self_type::begin();
      iterator pointIt = point.begin();
      
      for( ; meIt != self_type::heisenberg_end(); ++meIt, ++pointIt ) {
        *pointIt = -(*meIt);
      }
      // This isn't particularly clever, but seems necessary. If this function
      // is called a lot it would be a good idea to implement a "lazy" negation
      // function which didn't initialize, so that computing the dependent
      // variable could be held off until necessary.
      point.initialize();
      return point;
    }
    
    /**
     *  @brief  Performs addition in the heisenberg component
     *  @param  rhs The right hand point
     *  @return A copy of this after incrementing
     */
    self_type& operator+=( const self_type& __rhs ) {
      iterator meIt = self_type::begin();
      const_iterator rhsIt = __rhs.begin();
      
      for( ; meIt != self_type::heisenberg_end(); ++meIt, ++rhsIt ) {
        *meIt += *rhsIt;
      }
      initialize();
      
      return *this;
    }
    
    /**
     *  @brief  Performs subtraction in the heisenberg component
     *  @param  rhs The right hand point
     *  @return A copy of this after decrementing
     */
    self_type& operator-=( const self_type& __rhs ) {
      iterator meIt = self_type::begin();
      const_iterator rhsIt = __rhs.begin();
      
      for( ; meIt != self_type::heisenberg_end(); ++meIt, ++rhsIt ) {
        *meIt -= *rhsIt;
      }
      initialize();
      
      return *this;
    }
    
    
    /**
     *  @brief  Performs scalar multiplication in the heisenberg component
     *  @param  rhs A scalar
     *  @return A copy of this after multiplying
     */
    self_type& operator*=( param_type __rhs ) {
      iterator meIt = self_type::begin();
      for( ; meIt != self_type::heisenberg_end(); ++meIt ) {
        *meIt *= __rhs;
      }
      
      // The dependent variable can be updated quite efficiently
      float_type rhssq(__rhs*__rhs);
      dependent_.real_ref() *= rhssq;
      dependent_.real_ref() += 0.5 * ( (rhssq - 1.0) * self_type::height_ref() );
      dependent_.imag_ref() = self_type::r_ref();
      
      return *this;
    }
    
    
    /**
     *  @brief  Performs scalar division in the heisenberg component
     *  @param  rhs A scalar
     *  @return A copy of this after dividing
     */
    self_type& operator/=( param_type __rhs ) {
      float_type rhsinv(1.0 / __rhs);
      *this *= rhsinv;
      return *this;
    }
    
    
    /**
     *  @brief   A stric weak ordering on \c hyperbolic_point's 
     *  @param   rhs A Hyperbolic point
     *  @return  A weak ordering on the \c *this point and \c rhs
     */
    bool operator<( const self_type& __rhs ) const {
      typedef typename self_type::const_iterator _citer;
      
      //RVO
      bool b = false;
      
      _citer meIt = self_type::begin();
      _citer rhsIt = __rhs.begin();
      
      for( ; meIt != self_type::hyperbolic_end() ; ++meIt, ++rhsIt ) {
        
        typedef utility::math::is_equal<float_type> _ie;
        
        if( _ie( )( *meIt, *rhsIt ) ) continue; // indeterminate
        else if( *meIt < *rhsIt ) b = true;
        break;
      }
      return b;
    }
    
    
    /**
     * @brief  Forward starting point iterator.
     * @return A constant iterator pointing to the first point coordinate.
     */
    const_point_iterator point_begin() const {
      return const_point_iterator( *this, 0 );
    }
    
    /**
     * @brief  Forward end point iterator.
     * @return A constant iterator pointing one past the final point coordinate.
     * @note   The final point coordinate is the dependent coordinate.
     */
    const_point_iterator point_end() const {
      return const_point_iterator( *this, self_type::dimension_size + 1 );
    }
    
    
    /**
     * @brief  Reverse starting point iterator.
     * @return A reverse constant iterator pointing to the reverse first point
     *         coordinate, this being the dependent coordinate.
     */
    const_reverse_point_iterator point_rbegin() const {
      return const_reverse_point_iterator( point_end() );
    }
    
    
    /**
     * @brief  Reverse end point iterator.
     * @return A reverse constant iterator pointing one past the reverse final
     *         point coordinate; ie the normal first coordinate.
     */
    const_reverse_point_iterator point_rend() const {
      return const_reverse_point_iterator( point_begin() );
    }

    /**
     *  @brief   Print the point to a stream in a human friendly format.
     *  @param   os The out stream to print to, defaults to \c std::cout.
     */
    void pretty_print( std::ostream& __os = std::cout ) const {
      std::copy( point_begin(), point_end(),
                 std::ostream_iterator<zeta_type>(__os, "\n")
               );
    }
    
    
    /**
     *  @brief  Defines the input stream operator
     *  @param  __is An istream.
     *  @param  __t The object to stream out.
     */
    friend std::istream& operator>>( std::istream& __is, self_type& __t ) {
      __is >> static_cast<typename self_type::hyperbolic_structure_type&>(__t);
      __t.initialize();
      return __is;
    }
    
  private:
    /**
     *  @brief  Computes the quadratic form of \b zeta:
     *          \f$ \mathrm{Q}(\zeta) = \Sigma |\zeta_i|^2 \f$
     */
    float_type sys_compute_qf( void ) const {
      using utility::math::hermitian_quadratic_product;
      return float_type( hermitian_quadratic_product( self_type::zeta_ref_rbegin(),
                                                      self_type::zeta_ref_rend()
                                                    ) * (-0.5)
                       );
    }

    //! The dependent coordinate
    zeta_coordinate_type dependent_;
};

}// hyperbolic
}// geometric
}// structure
}// sg
#endif
