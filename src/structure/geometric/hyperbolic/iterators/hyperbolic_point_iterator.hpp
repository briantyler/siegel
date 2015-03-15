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
 *  @file     hyperbolic_point_iterator.hpp
 *  @brief    An inline header file for \c hyperbolic_point_iterator.
 *  @note     Include this file to iterate a through a hyperbolic point.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-17
 */

#ifndef _SG_HYPERBOLIC_POINT_ITERATOR_H
#define _SG_HYPERBOLIC_POINT_ITERATOR_H 1

// Global includes
#include <boost/iterator/iterator_facade.hpp>


namespace sg
{
namespace structure
{
namespace geometric
{
namespace hyperbolic
{
namespace iterators
{
/**
 *  @class    hyperbolic_point_iterator
 *  @brief    A class for iterating through a hyperbolic_point.
 *  @author   Brian Tyler
 *  @version  1.1
 *  @date     2008-04-17
 *  @param    _Point The structure type to iterate over.
 */

template <class _Point> class hyperbolic_point_iterator
: public boost::iterator_facade< hyperbolic_point_iterator<_Point>,
                                 const typename _Point::zeta_type,
                                 boost::random_access_traversal_tag
                               >
{
  friend class boost::iterator_core_access;
  
  public:
    //! Object type
    typedef hyperbolic_point_iterator<_Point> self_type;
    //! The difference type
    typedef typename self_type::difference_type difference_type;
    //! The type of the returned object
    typedef typename self_type::value_type value_type;
    //! the type of the object returned from dereferencing
    typedef typename self_type::reference reference;
    
  private:
    //! The point type
    typedef _Point point_type;
    //! The size type
    typedef typename point_type::size_type size_type;
    //! The underlying iterator type
    typedef typename point_type::const_zeta_ref_iterator iterator_type;
    
    
  public:
    /**
     *  @brief Default constructor
     */
    hyperbolic_point_iterator() : index_(0), point_(0), iterator_() {}
    
    /**
     *  @brief  Index constructor
     *  @param  point The point to iterate through.
     *  @param  index The index to construct at.
     */
    hyperbolic_point_iterator( const point_type& __point, size_type __index )
  : index_( __index ), point_( &__point ),
    iterator_( __point.zeta_ref_begin() + (__index - 1)) { }
    
    
  private:
    /**
     *  @brief  Moves the iterator index forward by 1
     *  @note   Required by @c iterator_facade
     */
    void increment( void ) {
      ++iterator_;
      ++index_;
    }
    
    
    /**
     *  @brief  Moves the iterator index backward by 1
     *  @note   Required by @c iterator_facade
     */
    void decrement( void ) {
      --iterator_;
      --index_;
    }
    
    
    /**
     *  @brief  Checks for equality between two iterators
     *  @param  that An iterator
     *  @note   Required by @c iterator_facade, although for our purposes we
     *          check for equality against an integer rather than an iterator
     */
    bool equal( const self_type& __that ) const {
      return bool( point_ == __that.point_ && index_ == __that.index_ );
    }
    
    
    /**
     *  @brief  Advance the iterator by \c n places
     *  @param  n The amount to advance by
     *  @note   Required by @c iterator_facade
     */
    void advance( difference_type __n ) {
      index_ += __n;
      iterator_ += __n;
    }
    
    
    /**
     *  @brief  Checks for distance between two iterators
     *  @param  that An iterator
     *  @note   Required by @c iterator_facade
     */
    difference_type distance_to( const self_type& __that ) const {
      return difference_type( __that.index_ - index_ );
    }
    
    
    /**
     *  @brief  Returns a reference to the current object
     *  @note   Required by @c iterator_facade
     */
    reference dereference( void ) const {
      switch( index_ ) {
        case 0:
          return point_type::first();
        case point_type::dimension_size:
          return point_->dependent()();
        default:
          return *iterator_;
      }
    }
    
    
  private:
    //! The index of the iterator
    size_type index_;
    //! Pointer to the Hyperbolic point
    const point_type* point_;
    //! The iterator which iterates through the real and imaginary parts
    iterator_type iterator_;
};

}// iterators
}// hyperbolic
}// geometric
}// structure
}// sg
#endif
