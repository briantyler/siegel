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
 *  @file     cusp_sieve.hpp
 *  @brief    An inline header file for the \c cusp_sieve class.
 *  @note     Sieves a cusp to find which subspaces of a compact hyperbolic
 *            space are affected by a cusp.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-05-21
 */


#ifndef _SG_CUSP_SIEVE_H
#define _SG_CUSP_SIEVE_H 1

// Global includes
#include <cassert>
#include <cstddef>
#include <utility>
#include <set>

#include <boost/array.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/tuple/tuple.hpp>

// Local includes
#include "structure/geometric/hyperbolic/space_tree.hpp"
#include "engine/cusp_validator.hpp"
#include "utility/functors/empty_functor.hpp"


namespace sg
{
namespace engine
{
/**
 *  @class    cusp_sieve
 *  @brief    An engine for sieving cusps.
 *  @param    N The dimension of the space
 *  @param    _Depth The depth of the sieve
 *  @param    _Float A floating point type; defaults to \c double
 *  @param    _Integer An integral type; defaults to \c long
 *  @param    _Id The identifier of the field
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-05-25
 */
template < std::size_t _Depth, std::size_t N, class _Float = double,
           class _Integer = long, std::size_t _Id = 0>
    class cusp_sieve
{
  public:
    //! The object type
    typedef cusp_sieve<_Depth,N,_Float,_Integer,_Id> self_type;
    //! The float type
    typedef _Float float_type;
    //! The integral type
    typedef _Integer integer_type;
    //! The validation engine type
    typedef cusp_validator<N,float_type,integer_type,_Id> validator_type;
    //! The cusp type
    typedef typename validator_type::cusp_type cusp_type;
    //! The point type
    typedef typename validator_type::point_type point_type;
    //! The field type
    typedef typename validator_type::field_type field_type;
    //! The space tree type
    typedef structure::geometric::hyperbolic::space_tree<_Depth,N,float_type>
        tree_type;
    //! The shared cusp pointer type
    typedef boost::shared_ptr<cusp_type> cusp_sptr_type;
    //! Weak pointer to a shared cusp pointer type
    typedef boost::weak_ptr<cusp_type> cusp_wptr_type;
    //! The data type
    typedef boost::tuple<float_type, cusp_wptr_type> data_type;
    //! The space type
    typedef typename tree_type::space_type space_type;
    
  private:
    //! The data comparison functor
    struct data_compare {
      bool operator()( const data_type& lhs, const data_type& rhs ) {
        return bool( lhs.template get<0>() > rhs.template get<0>() );
      }
    };
    
  public:
    //! The data container type
    // Note: when iterating through the multiset we want to access the larger
    // elements first, as these are the most effective. We could just use a
    // reverse iterator to do this, but it is more natural to order the data
    // sets in descending order. A set type is not appropriate, since two cusps
    // could have an identical effect on a region in space.
    typedef std::multiset<data_type, data_compare> data_container_type;
    
    //! The type of the array of data containers
    typedef boost::array<data_container_type,tree_type::static_size>
        data_array_type;
    
    /**
     *  @brief  Default constructor
     *  @note   Before sieving a cusp it is necessary to propagate the engine
     *          from a space.
     */
    cusp_sieve( ) : tree_(), validator_(), data_(0) { }
    
    
    /**
     *  @brief  Bind the engine to a data store
     *  @param  data  The data store to bind to
     */
    void bind_data( data_array_type& __data ) { data_ = &__data; }
    
    
    /**
     *  @brief  Set and propagate the root space.
     *  @param  space The space to propagate.
     */
    void propagate( const space_type& __space ) { tree_.propagate( __space ); }
    
    
    /**
     *  @brief  Compute the index of the top level space containing the point.
     *  @param  point The point to compute the index of.
     *  @return The index of the subspace that contains \c point
     */
    std::size_t get_index( const point_type& __point ) const {
      return std::size_t( tree_.get_index( __point ) );
    }
    
    
    /**
     *  @brief  Add a cusp to the bound data array of sieved cusps
     *  @param  cuspptr A shared pointer to the cusp which is to be sieved.
     *  @return \c true if the input cusp is effective for some point in the
     *          root space.
     */
    bool operator() ( const cusp_sptr_type& __cuspptr )
    {
      using namespace boost::mpl;
      
      // RVO
      bool b = false;
      
      // Ensure that there is somewhere to put the data
      assert( data_ != 0 );
      
      // Make a weak copy of the shared pointer. The reason for using a weak
      // pointer (which acts like an observor) is so that at some future point,
      // when extending this algorithm to deal with the full problem of
      // finding the fundamental domain, one can remove non-contributing cusps
      // by deleting the original smart pointer and then get rid of all the
      // weak pointers which point to null pointers.
      cuspptr_ = __cuspptr;
      
      // Bind the cusp to the validator
      validator_.bind_cusp( *__cuspptr );
      
      
       // Define mpl integer types from the template arguments
      typedef int_<0> zero_int;
      typedef int_<_Depth> depth_int;
      
      // Determine what to do on the next loop
      typedef typename if_< equal_to<zero_int,depth_int>,
                            sys_sieve_final,
                            sys_sieve_recursive<0>
                          >::type _fn;
      
      validator_.bind_space( tree_.root() );
      if( validator_() ) {
        // Start sieving
        std::size_t loc = 0;
        _fn()( this, loc );
        b = true;
      }
      
      // Otherwise the cusp has no effect on the space and we ignore it.
      return b;
    }
    
    
    /**
     *  @brief  Access to the tree for testing purposes.
     */
    const tree_type& tree( void ) const { return tree_; }
    
    
  private:
    /**
     *  @struct sys_sieve_final
     *  @brief  A helper function object that terminates a cusp sieve
     */
    struct sys_sieve_final {
      /**
       *  @brief  The functor which terminates the recursion and records data.
       *  @param  t The cusp sieve engine.
       *  @param  loc The current index
       */
      void operator() ( self_type* __t, std::size_t __loc ) const
      {
        // It is generally more efficient to insert a chunk of data into a
        // vector and then sort it, rather than insert objects into sets,
        // which automatically sorts them. However, here the insertion time is
        // negligible compared to the time it takes to sieve a cusp, so the
        // inefficiency introduced is negigible, but the benefit is that the
        // sieve is always correctly ordered.
        __t->data_->at(__loc).insert(
                                data_type( __t->validator_.effect(),
                                           __t->cuspptr_
                                         )
                              );
        assert( !( __t->cuspptr_.expired() ) );
      }
    };
    
    /**
     *  @struct sys_sieve_recursive
     *  @brief  A helper function object that sieves a cusp
     *  @param  Level The level to compute the index at.
     */
    template <std::size_t _Level> struct sys_sieve_recursive {
      /**
       *  @brief  The functor which carries out the recursion.
       *  @param  t The cusp sieve engine.
       *  @param  loc The current index
       */
      void operator() ( self_type* __t, std::size_t __loc ) const
      {
        using boost::mpl::if_;
	using boost::mpl::int_;
	using boost::mpl::equal_to;
        
        // Define mpl integer types from the template arguments
        typedef int_<_Level + 1> next_level_int;
        typedef int_<_Depth> depth_int;
        
        // Repeat until the maximum depth is reached
        typedef typename if_< equal_to<next_level_int,depth_int>,
                              sys_sieve_final,
                              sys_sieve_recursive<next_level_int::value>
                            >::type _fn;
        
        __loc *= 2;
        __t->validator_.bind_space( __t->tree_.template at<_Level+1>(__loc) );
        if( __t->validator_() ) _fn()( __t, __loc );
        
        ++(__loc);
        __t->validator_.bind_space( __t->tree_.template at<_Level+1>(__loc) );
        if( __t->validator_() ) _fn()( __t, __loc );
      }
    };
    
    
    //! The tree of spaces
    tree_type tree_;
    //! The validation engine
    validator_type validator_;
    //! Pointer to the data array that holds the sieving data
    data_array_type* data_;
    //! The cusp pointer
    cusp_wptr_type cuspptr_;
};

}// engine
}// sg
#endif
