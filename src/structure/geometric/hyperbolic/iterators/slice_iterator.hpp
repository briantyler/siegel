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
 *  @file     slice_iterator.hpp
 *  @brief    An inline header file for \c slice_iterator.
 *  @note     Include this file to iterate a hyperbolic point through a
 *            heisenberg slice
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-25
 */

#ifndef _SG_SLICE_ITERATOR_H
#define _SG_SLICE_ITERATOR_H 1

// Global includes
#include <algorithm>

#include <boost/array.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/call_traits.hpp>

// Local includes
#include "utility/functors/location_builder.hpp"
#include "utility/functors/mem_fun_adaptor.hpp"


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
 *  @class    slice_iterator
 *  @brief    A class for iterating a point through a heisenberg space.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-09-07
 *  @param    _Slice The heisenberg slice type to iterate over.
 */
template <class _Slice> class slice_iterator
  : public boost::iterator_facade< slice_iterator<_Slice>,
                                   const typename _Slice::point_type,
                                   boost::bidirectional_traversal_tag
                                 >
{
  friend class boost::iterator_core_access;
  
  public:
    //! Object type
    typedef slice_iterator<_Slice> self_type;
    //! Slice type
    typedef _Slice slice_type;
    //! Point type
    typedef typename slice_type::point_type point_type;
     //! Integer type
    typedef typename slice_type::integer_type integer_type;
    //! Float type
    typedef typename slice_type::float_type float_type;
    //! Size type
    typedef typename slice_type::size_type size_type;
    
    
  private:
    //! Type of the index array
    typedef boost::array<integer_type, slice_type::heisenberg_size>
        index_array_type;
    //! Incrementor type
    typedef typename slice_type::incrementor_type incrementor_type;
    //! Parameter type of the integer type
    typedef typename boost::call_traits<integer_type>::param_type int_param_type;
    
    
  public:
    //! The difference type
    typedef typename self_type::difference_type difference_type;
    //! The type of the returned object
    typedef typename self_type::value_type value_type;
    //! the type of the object returned from dereferencing
    typedef typename self_type::reference reference;
    
    
    /**
     *  @brief  Default constructor
     */
    slice_iterator( )
  : slice_(0), index_(0), point_(), indices_(), count_(0), notFinished_(true)
    { }
    
    
    /**
     *  @brief  Index constructor
     *  @param  slice The slice to iterate through
     *  @param  index The index to construct at
     */
    slice_iterator( const slice_type* __slice, int_param_type __index )
  : slice_(__slice), index_(__index), point_( slice_->point_at(__index) ),
    indices_( ), count_(0), notFinished_(true)
    { sys_build_index_array(); }
    
    
    /**
     *  @brief  Get the current index of the point in the slice.
     *  @return The current index of the point in the slice.
     */
    const integer_type& index( void ) const { return index_; }
    
    
    /**
     *  @brief  Conversion to a bool indicates the end of the iteeffectn
     *  @return True if the slice has been iterated through.
     */
    operator bool ( ) const { return notFinished_; }
    
    
  private:
    /**
     *  @brief  Build the array of indices
     */
    void sys_build_index_array( void ) {
      using namespace utility::functors;
      typedef typename slice_type::return_type _rt;
      typedef location_builder<integer_type> _lb;
      typedef const_mem_fun_ref_adaptor< const integer_type&, _rt, _lb,
                                         typename _lb::result_type
                                       > _lb_ad;

      std::transform( slice_->begin(), slice_->heisenberg_end(),
                      indices_.begin(),
                      _lb_ad( &_rt::resolution, _lb(index_) )
                    );
    }
    
    
    /**
     *  @brief  Moves the iterator forward by 1
     *  @note   Required by @c iterator_facade
     */
    void increment( void ) {
      // Iterate over each interval
      size_type i = 0;
      ++count_;
      
      for( ; i != index_array_type::static_size - 1; ++i ){
        // Check to see if we are at the end of the interval
        if( ++( indices_[i] ) == slice_->common_at(i).resolution() ) {
          // We have reached the end of this interval, so reset the point
          movement_helper( i, slice_->reset(), false );
          // Return the index to zero
          indices_[i] = 0;
        }
        else {
          // Increase the global index
          ++index_;
          
          // We are not at the end of the interval so increment the cube here
          movement_helper( i, slice_->forward() );
          return;
        }
      }
      
      /* If we get here then we are at the final coordinate, we want to
       * continually increase this.
       */
      ++index_;
      if( ++( indices_[i] ) == slice_->common_at(i).resolution() ) {
        notFinished_ = false;
      }
      movement_helper( i, slice_->forward() );
    }
    
    
    /**
     *  @brief  Moves the iterator backward by 1
     *  @note   Required by @c iterator_facade
     */
    void decrement( void ) {
      // Iterate over each interval
      size_type i = 0;
      ++count_;
      
      for( ; i != index_array_type::static_size - 1; ++i ){
        // Check to see if we are at the start of the interval
        if( --( indices_[i] ) == static_cast<integer_type>(-1) ) {
          // We have gone one past the start of this interval, so reset the point
          movement_helper( i, slice_->rreset(), false );
          // Return the index to one less than the resolution
          indices_[i] = slice_->common_at(i).resolution() - 1;;
        }
        else {
          // Increase the global index
          --index_;
          
          // We are not at the end of the interval so increment the cube here
          movement_helper( i, slice_->reverse() );
          return;
        }
      }
      
      --index_;
      movement_helper( i, slice_->reverse() );
    }
    
    
    /**
     *  @brief  Increment the point
     *  @param  loc The index to increment at.
     *  @param  incrementor The incrementor which does the incrementing
     *  @param  rebuild Should the point be rebuilt
     */
    void movement_helper
        ( size_type __loc, const incrementor_type& __incrementor,
          bool __rebuild = true
        )
    {
      // Every so many iteeffectns the point is refreshed.
      if( count_ == refresh_ && __rebuild ) {
        count_ = 0;
        point_ = slice_->point_at(index_);
      }
      else {
        __incrementor.common_at(__loc)( point_,
                                        point_.common_at(__loc)
                                      );
      }
    }
    
    
    /**
     *  @brief  Checks for equality between two iterators
     *  @param  that An iterator
     *  @note   Required by @c iterator_facade, although there is no mechanism
     *          to check against an end iterator.
     */
    bool equal( const self_type& __that ) const {
      return bool( index_ == __that.index_ );
    }
    
    
    /**
     *  @brief  Returns a reference to the current bound
     *  @note   Required by @c iterator_facade
     */
    reference dereference( void ) const { return point_; }
    
    
  public:
    /**
     *  @brief  Prepare the iterator to be updated.
     *  @note   Call before changing the resolution of the slice.
     *  @see    update()
     */
    void prepare_update( void ) {
      point_ = slice_->cube_at( index_ ).back();
    }
    
    
    /**
     *  @brief  Update the iterator after changing the resolution of the space.
     *  @note   After updating, the new iterator is either logically less than
     *          the original iterator in all dimensions, or is the starting
     *          iterator.
     *  @note   HOW TO UPDATE:
     *          1) Decrement the iterator (and any others connected to the slice)
     *             REASON: This takes the point to the last validated point.
     *          2) Call pre_update()
     *          3) Modify the resolution of the slice.
     *          4) Call update()
     */
    void update( void ) {
      // Compute the new index
      index_ = slice_->location_at( point_ );
      
      // Move the point to the new index
      point_ = slice_->point_at( index_ );
      
      // Build the array of indices at this index.
      sys_build_index_array();
      
      // This shouldn't be necessary, but you never know...
      if( index_ < slice_->resolution() ) {
        notFinished_ = true;
      }
      else {
        notFinished_ = false;
      }
      
      count_ = 0;
    }
    
    
  private:
    //! Constant reference to the slice
    const slice_type* slice_;
    //! Current index
    integer_type index_;
    //! The point to iterate
    point_type point_;
    //! Pointer to the array of indices
    index_array_type indices_;
    //! The refresh interval
    static const std::size_t refresh_ = 65536UL;
    //! The refresh counter
    std::size_t count_;
    //! Fast end checking
    bool notFinished_;
};

}// iterators
}// hyperbolic
}// geometric
}// structure
}// sg
#endif
