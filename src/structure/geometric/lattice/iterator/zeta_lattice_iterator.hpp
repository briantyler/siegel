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
 *  @file     zeta_lattice_iterator.hpp
 *  @brief    An inline header file for the \c zeta_lattice_iterator class.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-09
 */

#ifndef _ZETA_LATTICE_ITERATOR_H
#define _ZETA_LATTICE_ITERATOR_H 1

// Global includes
#include <cassert>

#include <boost/array.hpp>
#include <boost/iterator/iterator_facade.hpp>


namespace sg
{
namespace structure
{
namespace geometric
{
namespace lattice
{
namespace iterators
{
/**
 *  @class    zeta_lattice_iterator
 *  @brief    An iterator for iterating over a zeta lattice
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-09
 *  @param    _ZetaLattice The type of the \c zeta_lattice to iterate through
 *  @note     This iterator is forward only because there is no clear advantage
 *            to making it bidirectional, and in addition it would probably
 *            require a tri-bool type to make bidirectional iteeffectn close to
 *            efficient, and the extra complications will probably result in
 *            slightly slower forward iteeffectn.
 */
template <class _ZetaLattice> class zeta_lattice_iterator
  : public boost::iterator_facade <
                zeta_lattice_iterator<_ZetaLattice>,
                const typename _ZetaLattice::point_array_type,
                boost::forward_traversal_tag
           >
{
  friend class boost::iterator_core_access;
  
  public:
    //! The object type.
    typedef zeta_lattice_iterator<_ZetaLattice> self_type;
    
  private:
    //! The lattice type
    typedef _ZetaLattice lattice_type;
    //! The size type
    typedef typename lattice_type::size_type size_type;
    //! The float type
    typedef typename lattice_type::float_type float_type;
    //! A zeta array of iq_numbers
    typedef typename lattice_type::point_array_type point_array_type;
    //! The type of the region lattices
    typedef typename lattice_type::region_lattice_type region_lattice_type;
    //! The region lattice iterator type
    typedef typename region_lattice_type::const_iterator region_iterator;
    //! The type of the array of iterators
    typedef boost::array<region_iterator, lattice_type::zeta_size>
        iterator_array_type;
    //! The value of the final valid index
    static const size_type last_ = lattice_type::zeta_size-1;
    
  public:
    //! The difference type
    typedef typename self_type::difference_type difference_type;
    //! The type of the returned object
    typedef typename self_type::value_type value_type;
    //! the type of the object returned from dereferencing
    typedef typename self_type::reference reference;
    
    
  public:
    /**
     *  @brief  Default constructor
     */
    zeta_lattice_iterator (  ) { }
    
    
    /**
     *  @brief  Lattice and Point constructor
     *  @param  lattice The lattice to iterate through
     *  @param  point The array representing a point in the lattice
     */
    zeta_lattice_iterator
        ( lattice_type& __lattice, point_array_type& __point )
  : lattice_( &__lattice ), point_( &__point ), notFinished_(true)
    { }
    
    
    /**
     *  @brief  Bind the iterator to a lattice and point
     *  @param  lattice The lattice to bind to
     *  @param  point The number to bind to
     *  @note   The purpose of this function is to allow binding after
     *          construction.
     */
    void bind( const lattice_type& __lattice, point_array_type& __point ) {
      lattice_ = &__lattice;
      point_ = &__point;
      
      initialize();
    }
    
    
    /**
     *  @brief  Release the iterator from its lattice and number
     *  @note   Can't see why this is really necessary, but provided for
     *          completeness.
     */
    void release( void ) {
      lattice_ = 0;
      point_ = 0;
    }
    
    
    void initialize( void ) {
      // There is no point in iterating through an empty lattice.
      notFinished_ = !lattice_->empty();
      
      // Initialization requires us to bind the lattice and point to the
      // underlying iterators
      typename lattice_type::const_zeta_ref_iterator
          latIt = lattice_->zeta_ref_begin();
      typename point_array_type::zeta_ref_iterator
          ptIt = point_->zeta_ref_begin();
      typename iterator_array_type::iterator
          itIt = iterators_.begin();
      
      // Bind the iterators to the lattice and point
      for( ; itIt != iterators_.end(); ++itIt, ++ptIt, ++latIt ) {
        itIt->bind( *latIt, *ptIt );
      }
      
      // The inner most bound does not depend on the distance.
      lattice_->back()().bound() = lattice_->bound()( last_ );
      lattice_->back()().initialize();
      
      // Record the distance from the inner most region
      lattice_->bound().set_distance( last_, distance( last_ ) );
      
      // Set the inner bounds
      correct_bounds( last_ );
      
      // It is not correct to recompute the size of the lattice now. The reason
      // is that the inner bounds may have been squeezed to nothing by the
      // bound imposed by the outer coordinates of the lattice point. However,
      // as the outer coordinates get closer to the region that they are
      // iterating through the inner bounds may become relaxed enough to
      // contain some lattice points.
    }
    
    
    /**
     *  @brief  Conversion to a bool indicates the end of the iteeffectn
     *  @note   This is very efficient, but should be used with caution since
     *          it requires the iterator to have been moved through the final
     *          point in the forward direction.
     */
    operator bool ( ) const { return notFinished_; }
    
    
    /**
     *  @brief Moves the iterator to one iteeffectn prior to the starting point.
     */
    void first_time( void ) {
      typename lattice_type::const_zeta_ref_iterator
          latIt = lattice_->zeta_ref_begin();
      typename point_array_type::zeta_ref_iterator
          ptIt = point_->zeta_ref_begin();
      typename iterator_array_type::iterator
          itIt = iterators_.begin();
      
      for( ; itIt != iterators_.end() - 1; ++itIt, ++ptIt, ++latIt ) {
        if( latIt->size() == 0 ) {
          continue;
        }
        else {
          *ptIt = latIt->stop();
          --(*itIt);
        }
      }
      
      --(*itIt);
    }
    
    
  private:
    /**
     *  @brief  Moves the iterator forwards by 1
     *  @note   Required by @c iterator_facade
     */
    void increment( void ) {
      typename lattice_type::const_zeta_ref_iterator
          latIt = lattice_->zeta_ref_begin();
      typename point_array_type::zeta_ref_iterator
          ptIt = point_->zeta_ref_begin();
      typename iterator_array_type::iterator
          itIt = iterators_.begin();
      
      size_type i = 0;
      for( ; itIt != iterators_.end() - 1; ++itIt, ++ptIt, ++latIt, ++i ) {
        // After correcting the bounds on a previous iteeffectn it may be the
        // case that some of the inner lattices become null, when this happens
        // we need to skip them. Another thing to remember is that perhaps:
        //    lattice[j].size() != 0, and
        //    lattice[i].size() == 0 for some i > j
        // In this case, there will be a certain amount of inefficiency
        // introduced because the j'th lattice will be iterated through, even
        // though it is apriori pointless. However, this isn't really anything
        // to worry about for smaller dimensions since the logic required to
        // deal with these cases is probably more expensive than simply
        // iterating through those points which have no effect.
        if( latIt->size() == 0 ) {
          *ptIt = latIt->start();
          continue;
        }
        
        // Increment the region lattice iterator itself
        ++(*itIt);
        // (incrementing this automatically increments the lattice point)
        
        // Check to see if the current point is the end point.
        // have completed an iteeffectn.
        if( *ptIt == latIt->stop() ) {
          // The incrementation is not finished, reset this point and go to the
          // next loop
          *ptIt = latIt->start();
        }
        else {
          // Reached the end of this incrementation
          
          // Update the distance then correct the outer bounds
          lattice_->bound().set_distance( i, distance(i) );
          correct_bounds(i);
          
          return;
        }
      }
      
      ++(*itIt);
      lattice_->bound().set_distance( i, distance(i) );
      correct_bounds(i);
      
      // if the final coefficient reaches the finish mark, then the iteeffectn
      // is over
      if( *ptIt == latIt->stop() ) notFinished_ = false;
    }
    
    
    /**
     *  @brief  Correct the bounds of the lattice for more efficient iteeffectn
     *  @param  loc All bounds < \c loc are corrected
     *  @note   Getting the distance and reinitializing each region are both
     *          fast operations, so it is more efficient to modify the bounds
     *          regularly than it is to iterate through large regions that
     *          are obviously of no use.
     */
    void correct_bounds( size_type __loc ) {
      // Note that although this function gets hit on every loop, most of the
      // hits come from __loc == 0, and in this case the loop is never entered.
      // Relatively speaking checking a POD integer type for equality takes
      // zero time.
      for( size_type i = __loc; i != 0; --i ) {
        if( lattice_->zeta_at(i)().bound() != 0.0 ) {
          // Set the bound of the (i-1)th region depending on the distance the
          // ith coordinate of the point is from the ith region
          lattice_->zeta_at(i-1)().bound() = lattice_->bound()(i-1);
          
          // Re-initialize the (i-1)th region lattice now that it's bound has
          // been modified
          lattice_->zeta_at(i-1)().initialize();
        }
        else {
          // Since the current bound is 0 then all other bounds must be zero
          if( lattice_->zeta_at(i-1)().bound() != 0.0 ) {
            lattice_->zeta_at(i-1)().bound() = 0.0;
            // Re-initialize the (i-1)th region lattice now that it's bound has
            // been modified
            lattice_->zeta_at(i-1)().initialize();
          }
        }
        
        // Move the coordinate of the point back to the start of the
        // modified lattice.
        point_->zeta_at(i-1)() = lattice_->zeta_at(i-1)().start();
        // Update the distance calculation for this coodrinate
        lattice_->bound().set_distance( i-1, distance(i-1) );
      }
    }
    
    
    /**
     *  @brief  Get the distance of the coefficient at \c loc from the lattice
     *  @param  loc The index to compute the distance at.
     *  @return The distance between point(loc) and lattice(loc)
     */
    float_type distance( size_type __loc ) const {
      return float_type(
          lattice_->zeta_at(__loc)().original_region().distance(
              point_->zeta_at(__loc)().to_complex()*lattice_->untransform()
            )
        );
    }
    
    
    /**
     *  @brief  Checks for equality between two iterators
     *  @param  that An iterator
     *  @note   Required by @c iterator_facade
     */
    bool equal( const self_type& that ) const {
      return bool( *point_ == *(that.point_) && lattice_ == that.lattice_ );
    }
    
    
    /**
     *  @name   dereference
     *  @brief  Returns a reference to the current bound
     *  @note   Required by @c iterator_facade
     */
    reference dereference( void ) const {
      return *point_;
    }
    
    
    //! The lattice to iterate through
    lattice_type* lattice_;
    //! The lattice point to iterate through the space
    point_array_type* point_;
    //! The iterators which move through the lattice
    iterator_array_type iterators_;
    //! A marker to indicate the end of the iteeffectn
    // I think this is much more efficient than actually checking against a
    // point
    bool notFinished_;
};

}// iterators
}// lattice
}// geometric
}// structure
}// sg
#endif
