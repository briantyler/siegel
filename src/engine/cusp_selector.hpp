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
 *  @file     cusp_selector.hpp
 *  @brief    An inline header file for the \c cusp_selector class.
 *  @note     Selects cusp candiates that might be close to a given point
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-05-25
 */


#ifndef _SG_CUSP_SELECTOR_H
#define _SG_CUSP_SELECTOR_H 1

// Global includes
#include <cassert>
#include <cstddef>
#include <set>

// Local includes
#include "engine/cusp_sieve.hpp"
#include "engine/iterators/cusp_selector_iterator.hpp"


namespace sg
{
namespace engine
{
/**
 *  @class    cusp_selector
 *  @brief    An engine for selecting cusps that might be close to a point.
 *  @param    N The dimension of the space
 *  @param    _Depth The depth of the sieve to use for selection
 *  @param    _Float A floating point type; defaults to \c double
 *  @param    _Integer An integral type; defaults to \c long
 *  @param    _Id The identifier of the field
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-05-25
 */
template < std::size_t _Depth, std::size_t N, class _Float = double,
           class _Integer = long, std::size_t _Id = 0>
    class cusp_selector
{
  public:
    //! The object type
    typedef cusp_selector<_Depth,N,_Float,_Integer,_Id> self_type;
    //! The float type
    typedef _Float float_type;
    //! The integral type
    typedef _Integer integer_type;
    //! The cusp sieve type
    typedef cusp_sieve<_Depth,N,float_type,integer_type,_Id> sieve_type;
    //! The cusp type
    typedef typename sieve_type::cusp_type cusp_type;
    //! The point type
    typedef typename sieve_type::point_type point_type;
    //! The field type
    typedef typename sieve_type::field_type field_type;
    //! The space type
    typedef typename sieve_type::space_type space_type;
    //! The iterator that iterates through all the cusp candidates
    typedef iterators::cusp_selector_iterator<self_type> const_iterator;
    //! The shared cusp pointer type
    typedef typename sieve_type::cusp_sptr_type cusp_sptr_type;
    //! Weak pointer to a shared cusp pointer type
    typedef typename sieve_type::cusp_wptr_type cusp_wptr_type;
    
  private:
    //! Defines comparison operator for the set of cusp shared pointers.
    struct cusp_compare {
      bool operator( )
          ( const cusp_sptr_type& __lhs, const cusp_sptr_type& __rhs ) const
      { return bool( *__lhs < *__rhs ); }
    };
    
  public:
    //! The type of the array of data containers
    typedef typename sieve_type::data_array_type data_array_type;
    //! The shared cusp container type
    typedef std::set<cusp_sptr_type,cusp_compare> cusp_container_type;
    
    
    /**
     *  @brief  Default constructor
     */
    cusp_selector( ) : sieve_(), dilation_(0) { sieve_.bind_data( data_ ); }
    
    
    /**
     *  @brief  Get read access to the cusps container.
     *  @return A const reference to the cusps container.
     *  @note   The container is a unique sorted associative container.
     */
    const cusp_container_type& cusps( void ) const { return cusps_; }
    
    
    /**
     *  @brief  Get access to the array of sieving data.
     *  @note   This is supplied primarily for testing purposes.
     */
    data_array_type& data( void ) { return data_; }
    
    
    /**
     *  @brief  Get the maximum dilation factor over all cusps in the engine.
     *  @return The maximum dilation factor over all cusps in the engine
     */
    const integer_type& dilation( void ) const { return dilation_; }
    
    
    /**
     *  @brief  Propagate the root space..
     *  @param  space The space which candiate cusps must affect.
     *  @note   It is safe to re-propagate a space if only the height bound
     *          changes as the tree only propagates through the heisenberg part
     *          of the space.
     */
    void propagate( const space_type& __space ){ sieve_.propagate( __space ); }
    
    
    /**
     *  @brief  Add a cusp to the list of candidate cusps
     *  @param  cusp A cusp to add to the list of candiates.
     *  @param  check If true check to see if the cusp is already contained in
     *          the data set before sieving. When the field is a UFD this is
     *          unnecessary as we can filter out duplicate cusps based on the
     *          ideal generated by the cusp. Defaults to \c false
     *  @return \c true if the cusp was added to the selector. The only reasons
     *          for not adding a cusp are that: 1) it is a duplicate cusp of
     *          one already added, 2) it has no effect on the root space.
     *  @note   This is a time consuming operation, since to ensure rapid
     *          extraction of good candidates the cusps must be ordered in a
     *          specific way.
     *  @todo   Prior to sieving it may be a good idea to try pruning out cusps
     *          that clearly have no effect. However, this is probably
     *          expensive so should be something to look at after the
     *          \c siegel_engine works.
     */
    bool add_cusp( const cusp_type& __cusp, bool check = false ) {
      bool added = true;
      // If some kind of non-effective cusp pruning needs to be done then this
      // is the place to do it.
      
      // Create a smart pointer to a copy of the cusp.
      cusp_sptr_type ptr( new cusp_type( __cusp ) );
      
      if( check ) {
        // In this case we first ask if the cusp already exists in the set ...
        if( cusps_.count( ptr ) ){
          // ... it doesn't, so sieve the cusp for rapid extraction
          if( added = sieve_( ptr ) ) {
            // The cusp affects the root space, so by implication some subspace
            // at the maximum depth, so we can insert the cusp into our set of
            // valid cusps.
            cusps_.insert( ptr );
            if( dilation_ < __cusp.dilation() ) dilation_ = __cusp.dilation();
          }
        }
        else {
          // ... it already exists, so there is nothing to do.
          added = false;
        }
      }
      else {
        // We assume (well, it shouldn't be an assumption, we should know!)
        // that the cusp is not a duplicate of one that already exists, so
        // sieve the cusp for rapid extraction.
        assert( cusps_.count( ptr ) == 0 );
        
        if( added = sieve_( ptr ) ) {
          // Cusp affects the root space.
          cusps_.insert( ptr );
          if( dilation_ < __cusp.dilation() ) dilation_ = __cusp.dilation();
        }
      }
      
      return added;
    }
    
    
    /**
     *  @brief  Set the point to find close candidate cusps for.
     *  @param  point The point in hyperbolic space that we are seeking to
     *          raise by a cusp action.
     */
    const_iterator begin( const point_type& __point ) const {
      std::size_t loc = sieve_.get_index( __point );
      
      return const_iterator( data_.at(loc).begin(), loc );
    }
    
    
    /**
     *  @brief  Get an end sentinal iterator for the input iterator
     *  @param  iterator The iterator to get the sentinal for.
     *  @note   Since there are many possible end iterators it is necessary to
     *          know what is being iterated over to get a sentinal, hence the
     *          need for an input iterator.
     */
    const_iterator end( const const_iterator& __iterator ) const {
      return const_iterator( data_.at( __iterator.loc() ).end(),
                             __iterator.loc()
                           );
    }
    
    
    /**
     *  @brief  Given an iterator, asks if the container it points to is empty
     *  @param  iterator The iterator to ask if it container contains anything
     *  @return True if the underlying container is empty
     */
    bool empty( const const_iterator& __iterator ) const {
      return bool( data_.at( __iterator.loc() ).empty() );
    }
    
    
    /**
     *  @brief Clear the cusps and sieving data.
     */
    void clear( void ) {
      typedef typename data_array_type::iterator _iter;
      
      // Clear all of the cusps
      cusps_.clear();
      
      // Remove all of the sieving data
      for( _iter it = data_.begin(); it != data_.end(); ++it ) { it->clear(); }
      
      // Set the dilation factor back to 0
      dilation_ = 0;
    }
    
    
    /**
     *  @brief  Asks if there are any slots in the sieve which are empty
     *  @return True if all of the slots contain at least one cusp.
     */
    bool complete( void ) const {
      typedef typename data_array_type::const_iterator _iter;
      // RVO
      bool b = true;
      for( _iter it = data_.begin(); it != data_.end(); ++it ) {
        if( it->empty() ) {
          b = false;
          break;
        }
      }
      return b;
    }
    
    
  private:
    //! The set of all cusps currently under consideeffectn
    cusp_container_type cusps_;
    //! The sieved data to aid rapid extraction of candidates
    data_array_type data_;
    //! The sieve which processes the cusp information
    sieve_type sieve_;
    //! The current maximum dilation factor
    integer_type dilation_;
};

}// engine
}// sg
#endif
