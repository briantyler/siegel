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
 *  @file     heisenberg_cube.hpp
 *  @brief    An inline header file for the \c heisenberg_cube class.
 *  @note     A hypercube in Heisenberg space made from hyperbolic points.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-26-03
 */


#ifndef _SG_HEISENBERG_CUBE_H
#define _SG_HEISENBERG_CUBE_H 1

// Global includes
#include <boost/array.hpp>
#include <boost/operators.hpp>

// Local includes
#include "structure/geometric/detail/hyperbolic_base.hpp"
#include "structure/geometric/hyperbolic/hyperbolic_point.hpp"
#include "utility/io/container_to_string.hpp"
#include "utility/io/string_to_array.hpp"


namespace sg
{
namespace structure
{
namespace geometric
{
namespace hyperbolic
{
  /**
   *  @class    heisenberg_cube
   *  @brief    A hypercube in Heisenberg space made from hyperbolic points.
   *  @author   Brian Tyler
   *  @version  1.1
   *  @date     2008-03-20
   *  @param    N The hyperbolic dimension of the space the cube is in
   *  @param    _Float The floating point type
   */
template <std::size_t N, class _Float = double> class heisenberg_cube
  : public geometric::detail::hyperbolic_base<N>,
    private boost::equality_comparable< heisenberg_cube<N, _Float> >
{
  public:
    //! Object type
    typedef heisenberg_cube<N,_Float> self_type;
    //! Real type
    typedef std::size_t size_type;
    //! Real type
    typedef _Float float_type;
    //! Vertex type: hyperbolic points
    typedef hyperbolic_point<N,float_type> vertex_type;
    
  private:
    //! Array type
    typedef boost::array< vertex_type, self_type::hypercube_size > array_type;
    
  public:
    //! The vertex iterator type
    typedef typename array_type::iterator iterator;
    //! The constant vertex iterator type
    typedef typename array_type::const_iterator const_iterator;
    //! The reverse vertex iterator type
    typedef typename array_type::reverse_iterator reverse_iterator;
    //! The constant reverse vertex iterator type
    typedef typename array_type::const_reverse_iterator const_reverse_iterator;
    
    //! The static size of the cube
    enum { static_size = array_type::static_size };
    
    
    /**
     *  @brief   Get a constant reference to the vertex at \c loc
     *  @param   loc The index of the vertex to return
     *  @return  A constant reference to the vertex at \c loc
     */
    const vertex_type& vertex_at( size_type __loc ) const {
      return vertices_.at(__loc);
    }
    
    /**
     *  @brief   Get a reference to the vertex at \c loc
     *  @param   loc The index of the vertex to return
     *  @return  A reference to the vertex at \c loc
     */
    vertex_type& vertex_at( size_type __loc ) {
      return vertices_.at(__loc);
    }
    
    
    /**
     *  @brief   Get a constant reference to the first vertex.
     *  @return  A constant reference to the first vertex.
     */
    const vertex_type& front( void ) const { return vertices_.front();}
    
    /**
     *  @brief   Get a reference to the first vertex.
     *  @return  A reference to the first vertex.
     */
    vertex_type& front( void ) { return vertices_.front(); }
    
    
    /**
     *  @brief   Get a constant reference to the last vertex.
     *  @return  A constant reference to the last vertex.
     */
    const vertex_type& back( void ) const { return vertices_.back(); }
    
    /**
     *  @brief   Get a reference to the last vertex.
     *  @return  A reference to the last vertex.
     */
    vertex_type& back( void ) { return vertices_.back(); }
    
    
    /**
     *  @brief   Get the Euclidean midpoint of the cube
     *  @return  The euclidean midpoint of the cube.
     */
    vertex_type midpoint( void ) const {
      typedef typename vertex_type::iterator _iter;
      typedef typename vertex_type::const_iterator _citer;
      
      // Iterate through the first and last points of the hypercube and take
      // the mid-point of each.
      vertex_type v;
      _iter  vit = v.begin();
      _citer fit = front().begin();
      _citer bit = back().begin();
      
      for( ; vit != v.heisenberg_end(); ++vit, ++fit, ++bit ) {
        *vit = (*fit + *bit) * 0.5;
      }
      v.height() = front().height();
      v.initialize();
      
      return v;
    }
    
    
    /**
     *  @def     __SG_ITERATOR_DEF__(ITERATOR,DIRECTION,CONST)
     *  @brief   Defines the iterator access functions
     *  @note    Defines:
     *           <ul>
     *              <li>\c begin
     *              <li>\c end
     *              <li>\c rbegin
     *              <li>\c rend
     *           </ul>
     *           Including \c const versions.
     */
#define __SG_ITERATOR_DEF__(ITERATOR,DIRECTION,CONST)                         \
        ITERATOR DIRECTION( void ) CONST {                                    \
          return vertices_.DIRECTION( );                                      \
        }
    __SG_ITERATOR_DEF__(iterator,begin,);
    __SG_ITERATOR_DEF__(iterator,end,);
    __SG_ITERATOR_DEF__(const_iterator,begin,const);
    __SG_ITERATOR_DEF__(const_iterator,end,const);
    __SG_ITERATOR_DEF__(reverse_iterator,rbegin,);
    __SG_ITERATOR_DEF__(reverse_iterator,rend,);
    __SG_ITERATOR_DEF__(const_reverse_iterator,rbegin,const);
    __SG_ITERATOR_DEF__(const_reverse_iterator,rend,const);
#undef __SG_ITERATOR_DEF__
    
    
    /**
     *  @brief  Compares \c heisenberg_cube objects for equality.
     *  @param  lhs A \c heisenberg_cube
     *  @return \c true if all the vertices are equal.
     *  @note   Only the front and back vertices are checked for equality since
     *          these completely determine all the others due to the geometry
     *          of the object.
     */
    bool operator==( const self_type& __rhs ) const {
      return bool( front() == __rhs.front() && back() == __rhs.back() );
    }
    
    
    /**
     *  @brief  Defines the output stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::ostream& operator<<( std::ostream& __os, const self_type& __t ){
      using utility::io::container_to_string;
      typedef container_to_string<const_iterator> _cts;
      
      __os << _cts( )( __t.begin(), __t.end() );
      return __os;
    }
    
    /**
     *  @brief  Defines the input stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::istream& operator>>( std::istream& __is, self_type& __t ){
      using utility::io::string_to_array;
      typedef string_to_array<iterator> _sta;
      
      std::string s;
      __is >> s;
      _sta( )( s, __t.begin(), __t.end() );
      return __is;
    }
    
  private:
    //! Array of vertices which describe the hypercube.
    array_type vertices_;
};

}// hyperbolic
}// geometric
}// structure
}// sg
#endif
