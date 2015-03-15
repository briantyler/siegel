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
 *  @file     cube_iterator.hpp
 *  @brief    An inline header file for \c cube_iterator.
 *  @note     Include this file to iterate a hyperbolic cube through a
 *            heisenberg slice
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-25
 */

#ifndef _SG_CUBE_ITERATOR_H
#define _SG_CUBE_ITERATOR_H 1

// Global includes
#include <memory>
#include <algorithm>

#include <boost/scoped_ptr.hpp>
#include <boost/array.hpp>
#include <boost/operators.hpp>
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
 *  @class    cube_iterator
 *  @brief    A class for iterating through a heisenberg space.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-18
 *  @param    _Slice The heisenberg slice type to iterate over.
 *  @note     When writing this class I had the idea of partitioning the
 *            vertices of the heisenberg cube in such a way that only half the
 *            vertices needed to be incremented on each iteeffectn. The
 *            partitioning can be done very efficiently using bitwise operations
 *            (although it is a bit fiddly); given this I expected a decent
 *            reduction in in the overall iteeffectn time. However the exact
 *            opposite happens, iteeffectn time seems to increase by about 30% in
 *            -O3 optimised code. I suppose that this is testament to the
 *            efficiency of the incrementor objects; it is better to call an
 *            incrementor twice than absorb the overhead required to implement
 *            the logic needed for partitioning.
 *            
 *            Either way iterating is surprisingly rapid and is unlikely to be
 *            a bottleneck in the final algorithm. A speed increase is likely
 *            to be seen most dramatically by re-thinking how the cube_phi
 *            function is evaluated on a cube; instinct tells me that the
 *            thing to do is to create a virtual cube, based on just one point
 *            and then calculate the cube_phi function using that point and an
 *            incrementor style of thinking. However, since iteeffectn speed is
 *            likely to become a problem when dim > 5, and since there are
 *            significantly greater problems to be dealt with in dimensions
 *            of this size, I have not deviated from what I consider to be an
 *            object theoretically understandable way of doing things.
 *            (After the next can of Special Brew I may think differently)
 */
template <class _Slice> class cube_iterator
  : public boost::iterator_facade< cube_iterator<_Slice>,
                                   const typename _Slice::cube_type,
                                   boost::bidirectional_traversal_tag
                                 >
  , private boost::equality_comparable< cube_iterator<_Slice>,
                                        typename _Slice::integer_type
                                      >
{
  friend class boost::iterator_core_access;
  
  public:
    //! Object type
    typedef cube_iterator<_Slice> self_type;
    //! Slice type
    typedef _Slice slice_type;
    //! Cube type
    typedef typename slice_type::cube_type cube_type;
     //! Integer type
    typedef typename slice_type::integer_type integer_type;
    //! Float type
    typedef typename slice_type::float_type float_type;
    //! Size type
    typedef typename slice_type::size_type size_type;
    //! Pointer to the cube to iterate
    typedef std::auto_ptr<cube_type> cube_ptr_type;
    
  private:
    //! Type of the index array
    typedef boost::array<integer_type, slice_type::heisenberg_size>
        index_array_type;
    //! Pointer to the index array type
    typedef boost::scoped_ptr<index_array_type> index_array_ptr_type;
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
    cube_iterator( )
  : slice_(0), index_(0), cubePtr_(0), indexArrPtr_(0)
    { }
    
    
    /* This probably requires a little bit of explanation: The idea is that a
     * pointer to the hypercube is supplied, and this pointer is automatically
     * deleted when the iterator goes out of scope. The reason for doing this
     * is that a cube_type is an expensive object to construct, so when just
     * testing for the end of an iteeffectn sequence it is preferable not to
     * construct this object, in this case a null pointer is passed in. The
     * iterator is not actually functional, but this shouldn't be a problem.
     */
    /**
     *  @brief  Index constructor
     *  @param  slice The slice to iterate through
     *  @param  index The index to construct at
     *  @param  cubePtr an auto_ptr to the cube at index in slice.
     */
    cube_iterator( const slice_type& __slice, int_param_type __index,
                   cube_ptr_type __cubePtr = cube_ptr_type(0) )
        
  : slice_(&__slice), index_(__index), cubePtr_(__cubePtr),
    indexArrPtr_( cubePtr_.get() ? new index_array_type() : 0 )
    {
      /* Here is an important point to note about the constructor list above;
       * at the point we set "cubePtr_( cubePtr )" ownership of the pointer is
       * transferred from cubePtr to cubePtr_, as such
       * "indexArrPtr_( cubePtr.get() .. "
       * is wrong, as it always evaluates to 0.
       */
      
      build_index_array();
    }
    
    
    /**
     *  @brief  Copy constructor
     *  @param  that the cube_iterator to copy
     *  @note   Deep copies the pointers if they are non-null
     */
    cube_iterator( const self_type& __that )
  : slice_( __that.slice_ ), index_( __that.index_ ),
    cubePtr_(   __that.cubePtr_.get()
              ? new cube_type( *(__that.cubePtr_) )
              : 0
            ),
    indexArrPtr_(   __that.indexArrPtr_.get()
                  ? new index_array_type( *(__that.indexArrPtr_) )
                  : 0
                )
    /* It makes sense to deep copy the pointers since otherwise the original
     * object would become invalidated; this may not cause problems most of the
     * time, but is probably wrong and liable to lead to subtle and hard to
     * find bugs.
     */
    { }
    
    
    /**
     *  @brief  Assignment operator
     *  @param  that the cube_iterator to copy
     *  @note   Deep copies the pointers if they are non-null.
     */
    cube_iterator operator=( const self_type& __that ) {
      slice_ = __that.slice_;
      index_ = __that.index_;
      
      if( __that.cubePtr_.get() ) {
        if( cubePtr_.get() ) {
          *cubePtr_ = *(__that.cubePtr_);
        }
        else {
          cube_ptr_type ptr( new cube_type( *(__that.cubePtr_) ) );
          cubePtr_ = ptr;
        }
      }
      else {
        cubePtr_.reset();
      }
      
      if( __that.indexArrPtr_.get() ) {
        if( indexArrPtr_.get() ) {
          *indexArrPtr_ = *(__that.indexArrPtr_);
        }
        else {
          index_array_ptr_type
              ptr( new index_array_type( *( __that.indexArrPtr_ ) ) );
          indexArrPtr_.swap( ptr );
        }
      }
      else {
        indexArrPtr_.reset();
      }
      
      return *this;
    }
    
    
    /**
     *  @brief  Get the current index of the cube in the slice.
     *  @return The current index of the cube in the slice.
     */
    const integer_type& index( void ) const { return index_; }
    
    
  private:
    /**
     *  @brief  Build the array of indices
     */
    void build_index_array( void ) {
      using namespace utility::functors;
      typedef typename slice_type::return_type _rt;
      typedef location_builder<integer_type> _lb;
      typedef const_mem_fun_ref_adaptor< const integer_type&, _rt, _lb,
                                         typename _lb::result_type
                                       > _lb_ad;
      
      /* Translation: If the index array is not null, then build the index
       * array. This is done via a location builder object which, given a
       * global index and sequence of local resolutions, returns a sequence
       * of local indices, the "const_mem_fun_adaptor" just puts a wrapper
       * around this functor so that we can pass the resolution component of
       * each interval, rather than the interval itself.
       */
      if( indexArrPtr_.get() ) {
        std::transform( slice_->begin(), slice_->heisenberg_end(),
                        indexArrPtr_->begin(),
                        _lb_ad( &_rt::resolution, _lb(index_) )
                      );
      }
    }
    
    
    /**
     *  @brief  Moves the iterator index forward by 1
     *  @note   Required by @c iterator_facade
     */
    void increment( void ) {
      // Iterate over each interval
      size_type i = 0;
      
      // This is efficient because we are just creating a pointer
      typename index_array_type::iterator indexIt = indexArrPtr_->begin();
      
      
      /* It is wrong to create a iterator to iterate through the
       * slice as there is some overhead which is negated by one call to
       * common_at
       */
      for( ; i != index_array_type::static_size - 1; ++i, ++indexIt ){
        // Check to see if we are at the end of the interval
        if( ++( *indexIt ) == slice_->common_at(i).resolution() ) {
          // We have reached the end of this interval, so reset the cube
          movement_helper( i, slice_->reset() );
          // Return the index to zero
          *indexIt = 0;
        }
        else {
          // Increase the global index
          ++index_;
          
          // We are not at the end of the interval so increment the cube here
          movement_helper( i, slice_->forward(), true );
          return;
        }
      }
      /* If we get here then we are at the final coordinate. we want to
       * continually increase this.
       */
      ++index_;
      ++( *indexIt );
      movement_helper( i, slice_->forward(), true );
    }
    
    
    /**
     *  @brief  Moves the iterator index backward by 1
     *  @note   Required by @c iterator_facade
     */
    void decrement( void ) {
      size_type i = 0;
      typename index_array_type::iterator indexIt = indexArrPtr_->begin();
      
      for(; i != index_array_type::static_size - 1; ++i, ++indexIt ){
        if( --( *indexIt ) == static_cast<integer_type>(-1) ) {
          movement_helper( i, slice_->rreset() );
          *indexIt = slice_->common_at(i).resolution() - 1;
        }
        else {
          --index_;
          movement_helper( i, slice_->reverse(), true );
          return;
        }
      }
      --index_;
      --( *indexIt );
      movement_helper( i, slice_->reverse(), true );
    }
    
    
    /**
     *  @brief  Increment the cube
     *  @param  loc The index to increment at.
     *  @param  incrementor The incrementor which does the incrementing
     *  @param  rebuild Should the entire cube be rebuilt
     */
    void movement_helper
        ( size_type __loc, const incrementor_type& __incrementor,
          const bool __rebuild = false
        )
    {
      /* Every so many iteeffectns (65536 = 2^16 - which is chosen
       * so the compiler can optimise the modulus operation) the cube is
       * refreshed. This stops small errors from propagating. Adds about
       * 2-5% to the iteeffectn time, but this is a small price to pay for
       * accuracy.
       */
      if( index_ % refresh_ == 0 || !__rebuild ) {
        for( typename cube_type::iterator it = cubePtr_->begin();
             it != cubePtr_->end(); ++it )
        {
          __incrementor.common_at(__loc)( *it, it->common_at(__loc) );
        }
      }
      else {
        *cubePtr_ = slice_->cube_at(index_);
      }
    }
    
    
    /**
     *  @brief  Checks for equality between two iterators
     *  @param  that An iterator
     *  @note   Required by @c iterator_facade, although for our purposes we
     *          check for equality against an integer rather than an iterator
     */
    bool equal( const self_type& __that ) const {
      return bool( index_ == __that.index_ );
    }
    
    
    /**
     *  @brief  Returns a reference to the current bound
     *  @note   Required by @c iterator_facade
     */
    reference dereference( void ) const { return *cubePtr_; }
    
    
  public:
    /**
     *  @brief  Update the iterator after changing the resolution of the space.
     *  @note   After updating, the new iterator is either logically less than
     *          the original iterator in all dimensions, or is the starting
     *          iterator.
     *  @note   HOW TO UPDATE:
     *          1) Decrement the iterator (and any others connected to the slice)
     *             REASON: This takes the cube to the last validated hypercube.
     *          2) Modify the resolution of the slice.
     *          3) Call the update method on each iterator connected to the
     *             slice.
     */
    void update( void ) {
      // Compute the new index
      index_ = slice_->location_at( cubePtr_->back() );
      // Move the cube to the new index
      *cubePtr_ = slice_->cube_at(index_);
      // Build the array of indices at this index.
      build_index_array();
    }
    
    
    /**
     *  @brief  Compares a \c cube_iterator with an integer equality.
     *  @param  lhs A \c cube_iterator
     *  @param  rhs A positive integer
     *  @note   This is supplied for speed; rather than constructing a new
     *          iterator on each loop to check for the end of the sequence,
     *          the index can be checked directly.
     */
    friend bool operator==( const self_type& __lhs, int_param_type __rhs ) {
      return bool(  __lhs.index_ == __rhs );
    }
    
    
  private:
    //! Constant reference to the slice
    const slice_type* slice_;
    //! Current index
    integer_type index_;
    //! Pointer to the hypercube
    cube_ptr_type cubePtr_;
    //! Pointer to the array of indices
    index_array_ptr_type indexArrPtr_;
    //! The refresh interval
    static const std::size_t refresh_ = 65536UL;
};

}// iterators
}// hyperbolic
}// geometric
}// structure
}// sg
#endif
