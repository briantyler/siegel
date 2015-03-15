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
 *  @file     common_iterator.hpp
 *  @brief    An inline header file for the common type iterator.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-26
 */

#ifndef _SG_COMMON_ITERATOR_H
#define _SG_COMMON_ITERATOR_H 1

// Global includes
#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/iterator_traits.hpp>
#include <boost/type_traits.hpp>
#include <boost/mpl/if.hpp>
#include <boost/utility/enable_if.hpp>

namespace sg
{
namespace structure
{
namespace geometric
{
namespace detail
{
namespace iterators
{
/**
 *  @struct   iterator_fu
 *  @brief    A helper struct for getting the correct value_type from an input
 *            iterator.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-07-29
 *  @param    Iterator The pointer iterator to get the value_type from.
 */
template <class _Iterator> struct iterator_fu {
  typedef typename boost::remove_pointer<_Iterator>::type pointer_type;
  typedef typename boost::remove_pointer<pointer_type>::type type;
};


/**
 *  @class    common_iterator
 *  @brief    An iterator for iterating over the real and imaginary parts
 *            of the zeta component of a heisenberg structure
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-26
 *  @param    Iterator The base iterator to use
 */
template <class _Iterator> class common_iterator
  : public boost::iterator_facade < common_iterator<_Iterator>,
                                    typename iterator_fu<_Iterator>::type,
                                    boost::random_access_traversal_tag
                                  >
{
  friend class boost::iterator_core_access;
  // Self friendship makes const - non-const compatible
  template <class> friend class common_iterator;
  
  /**
   *  @struct enabler
   *  @brief  Needed for boost kung fu.
   */
  struct enabler {};
  
  
  public:
    //! The object type
    typedef common_iterator<_Iterator> self_type;
    
  private:
    //! Iterator type
    typedef _Iterator iterator_type;
    
  public:
    //! The difference type
    typedef typename self_type::difference_type difference_type;
    //! The type of the returned object
    typedef typename self_type::value_type value_type;
    //! the type of the object returned from dereferencing
    typedef typename self_type::reference reference;
    
    
    /**
     *  @brief Default constructor
     */
    common_iterator() : iterator_() { }
    
    /**
     *  @brief  Generic copy constructor
     *  @param  that Another \c common_iterator
     */
    template <class _OtherIterator>
    common_iterator
        ( const common_iterator<_OtherIterator>& __that,
          typename boost::enable_if<
                           boost::is_convertible<_OtherIterator,_Iterator>
                         , enabler
          >::type = enabler()
        )
    : iterator_(__that.iterator_) { }
    
    
    /**
     *  @brief  Iterator constructor
     *  @param  iterator an iterator pointing to a zeta array
     */
    common_iterator( const iterator_type __iterator )
    : iterator_(__iterator) { }
    
    
  private:
    /**
     *  @brief  Moves the iterator index forward by 1
     *  @note   Required by @c iterator_facade
     */
    void increment( void ) { ++iterator_; }
    
    
    /**
     *  @brief  Moves the iterator index backward by 1
     *  @note   Required by @c iterator_facade
     */
    void decrement( void ) { --iterator_; }
    
    
    /**
     *  @brief  Checks for equality between two iterators
     *  @param  that An iterator
     *  @note   Required by @c iterator_facade
     */
    void advance( difference_type __n ) { iterator_ += __n; }
    
    
    /**
     *  @brief  Checks for distance between two iterators
     *  @param  that An iterator
     *  @note   Required by @c iterator_facade
     */
    difference_type distance_to( const self_type& __that ) {
      return difference_type( __that.iterator_ - iterator_ );
    }
    
    
    /**
     *  @brief  Checks for equality between two iterators
     *  @param  that An iterator
     *  @note   Required by @c iterator_facade
     */
    template <class _OtherIterator> bool equal
        ( const common_iterator<_OtherIterator>& __that ) const
    {
      return bool( iterator_ == __that.iterator_ );
    }
    
    
    /**
     *  @name   dereference
     *  @brief  Returns a reference to the current bound
     *  @note   Required by @c iterator_facade
     */
    reference dereference( void ) const { return **iterator_; }
    
    
    //! Base iterator object
    iterator_type iterator_;
};

}// iterators
}// detail
}// geometric
}// structure
}// sg
#endif
