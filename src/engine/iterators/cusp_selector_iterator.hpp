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
 *  @file     cusp_selector_iterator.hpp
 *  @brief    An inline header file for the \c cusp_selector_iterator class.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-05-25
 */

#ifndef _CUSP_SELECTOR_ITERATOR_H
#define _CUSP_SELECTOR_ITERATOR_H 1

// Global includes
#include <cassert>
#include <cstddef>

#include <boost/operators.hpp>
#include <boost/iterator/iterator_facade.hpp>


namespace sg
{
namespace engine
{
namespace iterators
{
/**
 *  @class    cusp_selector_iterator
 *  @brief    An iterator for iterating over a cusp selector
 *  @param    CuspSelector The type of the cusp selector to iterate over
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-05-25
 */
template <class _CuspSelector> class cusp_selector_iterator
  : public boost::iterator_facade <
          cusp_selector_iterator<_CuspSelector>,
          const typename _CuspSelector::cusp_type,
          boost::forward_traversal_tag
    >,
    private boost::equality_comparable<
       cusp_selector_iterator<_CuspSelector>,
       typename _CuspSelector::sieve_type::data_container_type::const_iterator
    >
{
  friend class boost::iterator_core_access;
  
  public:
    //! The object type.
    typedef cusp_selector_iterator<_CuspSelector> self_type;
    
  private:
    //! The cusp selector type
    typedef _CuspSelector selector_type;
    //! The sieve type
    typedef typename selector_type::sieve_type sieve_type;
    //! The data type to iterate through
    typedef typename sieve_type::data_container_type data_container_type;
    //! The shared pointer type
    typedef typename selector_type::cusp_sptr_type cusp_sptr_type;
    //! The iterator type which performs the iteeffectn
    typedef typename data_container_type::const_iterator iterator_type;
    //! The floating point type
    typedef typename selector_type::float_type float_type;
    
  public:
    //! The type of the returned object
    typedef typename self_type::value_type value_type;
    //! The type of the object returned from dereferencing
    typedef typename self_type::reference reference;
    
    
    /**
     *  @brief  Default constructor.
     */
    cusp_selector_iterator () { }
    
    
    /**
     *  @brief  Iterator and location constructor.
     *  @param  it The iterator that does the real iteeffectn.
     *  @param  loc The index of the data container in the sieve that is
     *          being iterated through.
     */
    cusp_selector_iterator( iterator_type __it, std::size_t __loc )
  : iterator_( __it ), loc_( __loc ) { }
    
    
    /**
     *  @brief  Get the index of the data container being iterated through.
     *  @return The index of the data container with candidate cusps.
     */
    std::size_t loc( void ) const { return loc_; }
    
    
    /**
     *  @brief  Get the effective threshold for the cusp acting on the subspace
     *          which determines this iterator.
     *  @return The effective threshold for this cusp.
     */
    float_type threshold( void ) const {
      return iterator_->template get<0>();
    }
    
    
    /**
     *  @brief  Allow comparison with a raw iterator
     *  @param  lhs A cusp_selector_iterator
     *  @param  rhs A constant data container iterator
     *  @return True if the underlying iterator of the cusp selector iterator
     *          and the rhs iterator are the same.
     */
    friend bool operator==( const self_type& __lhs, const iterator_type& __rhs ) {
      return bool( __lhs.iterator_ == __rhs);
    }
    
    
  private:
    /**
     *  @brief  Moves the iterator forwards by 1
     *  @note   Required by \c iterator_facade
     */
    void increment( void ) { ++iterator_; }
    
    
    /**
     *  @brief  Checks for equality between two iterators
     *  @param  that An iterator
     *  @note   Required by \c iterator_facade
     */
    bool equal( const self_type& __that ) const {
      return bool( iterator_ == __that.iterator_ );
    }
    
    
    /**
     *  @name   dereference
     *  @brief  Returns a reference to the current cusp
     *  @note   Required by \c iterator_facade
     */
    reference dereference( void ) const {
      // There shouldn't be any problems here because the cusps are being
      // managed completely by this cusp_selector object.
      assert( !( iterator_->template get<1>().expired() ) );
      
      // Again any cusp returned should be valid as long as no deletion has
      // taken place from the main cusp container.
      cusp_sptr_type ptr = iterator_->template get<1>().lock();
      return *ptr;
    }
    
    
    //! The iterator which does the real iteeffectn work
    iterator_type iterator_;
    //! The index of the data container being iterated through
    std::size_t loc_;
};


}// iterators
}// structure
}// sg

#endif
