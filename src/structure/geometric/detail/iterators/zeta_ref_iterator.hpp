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
 *  @file     zeta_ref_iterator.hpp
 *  @brief    An inline header file for the zeta reference iterator.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-25
 */

#ifndef _SG_ZETA_REF_ITERATOR_H
#define _SG_ZETA_REF_ITERATOR_H 1

// Global includes
#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/iterator_traits.hpp>
#include <boost/type_traits.hpp>
#include <boost/mpl/if.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/assert.hpp>

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
 *  @struct   ref_kung_fu
 *  @brief    A helper struct for getting the correct value_type from an input
 *            iterator.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-07-29
 *  @param    Iterator The iterator to get the value_type from. This can be a
 *            pointer iterator, or an stl compatible iterator.
 */
template <class _Iterator> struct ref_kung_fu {
  typedef typename boost::iterator_value<_Iterator>::type zeta_coordinate_type;
  typedef typename zeta_coordinate_type::zeta_type zeta_type;
  typedef typename boost::remove_reference<
                        typename boost::iterator_reference<_Iterator>::type
                   >::type const_type;
  
  typedef boost::mpl::if_< boost::is_const<const_type>,
                           typename boost::add_const<zeta_type>::type,
                           zeta_type
                         > if_type;
  typedef typename if_type::type type;
};


/**
 *  @class    zeta_ref_iterator
 *  @brief    An iterator for iterating over the zeta component of a
 *            heisenberg structure
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-25
 *  @param    _Iterator The base iterator to use
 *  @param    _Value The return value type
 *  @note     \c iterator_facade is a base class template that implements
 *            the interface of standard iterators in terms of a few
 *            core functions and associated types, to be supplied by a
 *            derived iterator class. The template takes the following
 *            parameters:
 *            \verbatim
 *            template <
 *                class Derived
 *              , class Value
 *              , class CategoryOrTraversal
 *              , class Reference  = _Value&
 *              , class Difference = ptrdiff_t
 *            >
 *            \endverbatim
 *            <ul>
 *              <li>\c Derived - The type deriving from \c iterator_facade
 *              <li>\c Value - The return type.
 *              <li>\c CategoryOrTraversal - The iterator type see:
 *                  \link http://tinyurl.com/2swl53
 *              <li>\c Reference - The type of the object returned when the
 *                  iterator is dereferenced.
 *              <li>\c Difference - determines how the distance between two
 *                  iterators will be measured.
 *            </ul>
 *  @link     http://www.boost.org/libs/iterator/doc/index.html
 */
template <class _Iterator> class zeta_ref_iterator
  : public boost::iterator_facade< zeta_ref_iterator<_Iterator>,
                                   typename ref_kung_fu<_Iterator>::type,
                                   boost::random_access_traversal_tag
                                 >
{
  friend class boost::iterator_core_access;
  // Self friendship makes const - non-const compatible
  template <class> friend class zeta_ref_iterator;
  
  public:
    //! The object type.
    typedef zeta_ref_iterator<_Iterator> self_type;
    
  private:
    //! Iterator type
    typedef _Iterator iterator_type;
    
    /**
     *  @struct enabler
     *  @brief  Needed for more boost kung fu.
     */
    struct enabler {};
        
  public:
    //! The difference type
    typedef typename self_type::difference_type difference_type;
    //! The type of the returned object
    typedef typename self_type::value_type value_type;
    //! the type of the object returned from dereferencing
    typedef typename self_type::reference reference;
    
    
  public:
    /**
     *  @brief Default constructor
     */
    zeta_ref_iterator() : iterator_() { }
    
    
    /**
     *  @brief  Generic copy constructor
     *  @param  that Another \c zeta_ref_iterator
     */
    template <class _OtherIterator> zeta_ref_iterator
        ( const zeta_ref_iterator<_OtherIterator>& __that,
                     /* What follows is some boost kung fu which stops the code
                      * compiling if the iterator types are not convertible in
                      * a given direction:
                      * Specifically const -convert-> non-const, can't happen
                      */
         typename boost::enable_if<
                        boost::is_convertible<_OtherIterator,_Iterator>,
                        enabler
                  >::type = enabler()
        )
    : iterator_(__that.iterator_) { }
    
    
    /**
     *  @brief  _Iterator constructor
     *  @param  iterator An iterator pointing to a zeta array
     */
    zeta_ref_iterator( const _Iterator& __iterator )
    : iterator_(__iterator) { }
    
    
  private:
    /**
     *  @brief  Moves the iterator index forward by 1
     *  @note   Required by @c iterator_facade
     */
    void increment( void ) {
      ++iterator_;
    }
    
    
    /**
     *  @brief  Moves the iterator index backward by 1
     *  @note   Required by @c iterator_facade
     */
    void decrement( void ) {
      --iterator_;
    }
    
    
    /**
     *  @brief  Advance an iterator by a given amount
     *  @param  n the amount to advance by
     *  @note   Required by @c iterator_facade
     */
    void advance( difference_type __n ) {
      iterator_ += __n;
    }
    
    
    /**
     *  @brief  Compute the distance between two iterators
     *  @param  iterator An iterator
     *  @return The distance between the iterators (compatible with advance)
     *  @note   Required by @c iterator_facade
     */
    difference_type distance_to( const self_type& __iterator ) const {
      return difference_type( __iterator.iterator_ - iterator_ );
    }
    
    
    /**
     *  @brief  Checks for equality between two iterators
     *  @param  that An iterator
     *  @note   Required by @c iterator_facade
     */
    template <class _OtherIterator>
    bool equal( const zeta_ref_iterator<_OtherIterator>& __that ) const {
      return bool( iterator_ == __that.iterator_ );
    }
    
    
    /**
     *  @name   dereference
     *  @brief  Returns a reference to the current zeta coordinate
     *  @note   Required by @c iterator_facade
     */
    reference dereference( void ) const {
      return (*iterator_)();
    }
    
    
    //! Base iterator object
    iterator_type iterator_;
};

}// iterators
}// detail
}// geometric
}// structure
}// sg
#endif
