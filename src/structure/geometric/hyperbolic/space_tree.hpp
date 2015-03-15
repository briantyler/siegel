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
 *  @file     space_tree.hpp
 *  @brief    An inline header file for the \c space_tree class.
 *  @note     A static binary tree of cascading spaces.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-05-20
 */


#ifndef _SG_SPACE_TREE_H
#define _SG_SPACE_TREE_H 1

// Global includes
#include <algorithm>

// Local includes
#include "structure/geometric/hyperbolic/detail/static_tree.hpp"
#include "structure/geometric/hyperbolic/hyperbolic_space.hpp"
#include "utility/functors/empty_functor.hpp"


namespace sg
{
namespace structure
{
namespace geometric
{
namespace hyperbolic
{
/**
 *  @class   space_tree
 *  @brief   A binary tree of fixed sized arrays of spaces
 *  @param   N The dimension of the space
 *  @param   Depth The depth of the tree.
 *  @param   Float The floating point type
 *  @author  Brian Tyler
 *  @version 1.0
 *  @date    2008-05-20
 *  @note    If Depth is too large (depending on the dimension), then the
 *           object might segment fault, this will happen immediately on
 *           propagation. If this happens reduce the depth.
 */
template <std::size_t _Depth, std::size_t N, class _Float = double>
    class space_tree
  : public detail::static_tree< hyperbolic_space<N,_Float>, _Depth >
{
  public:
    //! The object type
    typedef space_tree<_Depth, N,_Float> self_type;
    //! The floating point type
    typedef _Float float_type;
    //! The size type
    typedef std::size_t size_type;
    //! The hyperbolic space type
    typedef hyperbolic_space<N,float_type> space_type;
    //! The hyperbolic point type
    typedef typename space_type::point_type point_type;
    //! The dimension of the space
    static const size_type dimension = N;
    
    
  public:
    /**
     *  @brief  Default constructor.
     */
    space_tree( ) { }
    
    
    /**
     *  @brief  Root constructor.
     *  @param  space The space to set as the root.
     */
    space_tree( const space_type& __space ) { propagate(__space); }
    
    
    /**
     *  @brief  Propagate the root space.
     *  @note   This function subdivides the space down the tree, dividing the
     *          first interval by two, then the next, until it loops and starts
     *          again or reaches the final level.
     */
    void propagate( void ) {
      using namespace boost::mpl;
      
      // If the starting depth is zero, then there is nothing to do
      typedef typename if_< equal_to< int_<0>, int_<_Depth> >,
                            utility::functors::empty_functor,
                            prv_propagate<0>
                          >::type _fn;
      
      // Propagate
      _fn()( this );
    }
    
    
    /**
     *  @brief  Set and propagate the root space.
     *  @param  space  The space to propagate.
     */
    void propagate( const space_type& __space ) {
      // Set the new root space
      this->root() = __space;
      
      propagate();
    }
    
    
    /**
     *  @brief  Compute the index of the top level space containing the point.
     *  @param  point The point to to compute the index of.
     *  @return The index of the subspace that contains \c point
     *  @note   This assumes that the point is in the root space, if this is
     *          not the case then the index returned is that of the closest
     *          subspace to the point.
     */
    size_type get_index( const point_type& __point ) const {
      using namespace boost::mpl;
      
      // Since the structure of the tree does not result in the standard index
      // structure, the easiest way to compute the index is to filter the point
      // through the same process used to construct the tree.
      size_type index = 0;
      
      // If the starting depth is zero, then there is nothing to do
      typedef typename if_< equal_to< int_<0>, int_<_Depth> >,
                            utility::functors::empty_functor,
                            prv_get_index<0>
                          >::type _fn;
      
      _fn()( &__point, this, &index );
      
      return index;
    }
    
    
  private:
    /**
     *  @struct prv_get_index
     *  @brief  A helper function object that computes the index of a point
     *  @param  Level The level to compute the index at.
     */
    template <size_type _Level> struct prv_get_index {
      /**
       *  @brief  The functor which carries out the recursion.
       *  @param  p The point to compute the index of
       *  @param  t The space tree to propagate.
       *  @param  loc The index computed so far.
       */
      void operator()
          ( const point_type* __p, const self_type* __t, size_type* __loc ) const
      {
        using boost::mpl::if_;
	using boost::mpl::int_;
	using boost::mpl::equal_to;
        
        // index is the index of the space interval to examine on this round
        static const size_type index = _Level % space_type::heisenberg_size;
        
        *__loc *= 2;
        // We are interested in filtering the point to the closest index rather
        // than determining if the point is or is not contained in a given
        // subspace, so it is best to filter on the basis of the upper bound.
        if( !( __t->template at<_Level+1>(*__loc)
                  .common_at(index).upper() > __p->common_at(index) )
          )
        {
          // In the second sub interval
          ++(*__loc);
        }
        
        // Define mpl integer types from the template arguments
        typedef int_<_Level + 1> next_level_int;
        typedef int_<_Depth> depth_int;
        
        // Repeat until the maximum depth is reached
        typedef if_< equal_to<next_level_int,depth_int>,
                     utility::functors::empty_functor,
                     prv_get_index<next_level_int::value>
                   > if_type;
        
        typename if_type::type()( __p, __t, __loc );
      }
    };
    
    
    /**
     *  @struct prv_propagate
     *  @brief  A helper function object which recursively propagates the
     *          space threough each level
     *  @param  Level The level to propagate
     */
    template <size_type _Level> struct prv_propagate {
      /**
       *  @brief  The functor which carries out the recursive propagation.
       *  @param  t The space tree to propagate.
       */
      void operator() ( self_type* __t ) const {
        using boost::mpl::if_;
	using boost::mpl::int_;
	using boost::mpl::equal_to;
        
        // index is the index of the space interval that is going to be halfed.
        // on this round
        static const size_type index = _Level % space_type::heisenberg_size;
        
        // All intervals have the same length so only compute this once
        float_type len
            = __t->template at<_Level>(0).common_at(index).length() * 0.5;
        
        for( size_type i = 0; i != self_type::template size<_Level>(); ++i ) {
          // Get references to the spaces to save writing the ugly syntax
          // more than once
          space_type& s0 = __t->template at<_Level + 1>(2*i);
          space_type& s1 = __t->template at<_Level + 1>(2*i + 1);
          
          // Set the two lower spaces equal to the one above them
          s0 = s1 = __t->template at<_Level>(i);
          
          // Now we need to modify the intervals for this loop:
          
          // Half both intervals
          s0.common_at(index).upper() -= len;
          s1.common_at(index).lower() += len;
        }
        // Now decide what happens on the next loop:
        
        // Define mpl integer types from the template arguments
        typedef int_<_Level + 1> next_level_int;
        typedef int_<_Depth> depth_int;
        
        // If on the next loop the level and depth will be equal then the space
        // has been propagated through all levels, so terminate: call
        // empty_functor(), otherwise we need to propagate to the next level.
        typedef if_< equal_to<next_level_int,depth_int>,
                     utility::functors::empty_functor,
                     prv_propagate<next_level_int::value>
                   > if_type;
        
        // Both functors need the same signature, hence the dots "..." in the
        // empty_functor. Note that "..." matches any number of POD types
        typename if_type::type()( __t );
      }
    };
};

}// hyperbolic
}// geometric
}// structure
}// sg
#endif
