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
 *  @file     region_lattice_iterator.hpp
 *  @brief    An inline header file for the \c region_lattice_iterator class.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-08
 */

#ifndef _SG_REGION_LATTICE_ITERATOR_H
#define _SG_REGION_LATTICE_ITERATOR_H 1

// Global includes
#include <cassert>

#include <boost/iterator/iterator_facade.hpp>

// Local includes
#include "utility/math/is_even.hpp"


namespace sg
{
namespace structure
{
namespace geometric
{
namespace lattice
{
namespace iterator
{
/**
 *  @class    region_lattice_iterator
 *  @brief    An iterator for iterating over a region lattice
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-08
 *  @param    _RegionLattice The type of the region_lattice to iterate through
 *  @link     http://www.boost.org/libs/iterator/doc/index.html
 */
template <class _RegionLattice> class region_lattice_iterator
  : public boost::iterator_facade <
          region_lattice_iterator<_RegionLattice>,
          const typename _RegionLattice::iq_number_type,
          boost::bidirectional_traversal_tag
    >
{
  friend class boost::iterator_core_access;
  
  public:
    //! The object type.
    typedef region_lattice_iterator<_RegionLattice> self_type;
    //! The lattice type
    typedef _RegionLattice lattice_type;
    //! The iq_number type
    typedef typename lattice_type::iq_number_type iq_number_type;
    //! The integer type
    typedef typename lattice_type::integer_type integer_type;
    //! The float type
    typedef typename lattice_type::float_type float_type;
    
  private:
    //! Even functor
    typedef utility::math::is_even<integer_type> even_func;
    //! Even functor
    typedef utility::math::is_odd<integer_type> odd_func;
        
  public:
    //! The difference type
    typedef typename self_type::difference_type difference_type;
    //! The type of the returned object
    typedef typename self_type::value_type value_type;
    //! The type of the object returned from dereferencing
    typedef typename self_type::reference reference;
    
    
  public:
    /**
     *  @brief Default constructor
     */
    region_lattice_iterator( ) : lattice_(0), iqNumber_(0) { }
    
    /**
     *  @brief  Lattice and Point constructor
     *  @param  lattice The lattice to iterate through
     *  @param  iqNumber The number that gets iterated
     */
    region_lattice_iterator
        ( const lattice_type& __lattice, iq_number_type& __iqNumber )
    : lattice_( &__lattice ), iqNumber_( &__iqNumber )
    {
      // During testing it is a good idea to ensure that the iterator is
      // properly initialised.
      assert(      lattice_->latticeType_
                == lattice_type::prv_lattice_enum::congruent
              ?
                   even_func()(iqNumber_->real())
                ==
                   even_func()(iqNumber_->imag()) : true
              );
    }
    
    
    /**
     *  @brief  Bind the iterator to a lattice and number
     *  @param  lattice The lattice to bind to
     *  @param  iqNumber The number to bind to
     *  @note   The purpose of this function is to allow binding after
     *          construction which is necessary to build the zeta lattice
     *          iterator.
     */
    void bind( const lattice_type& __lattice, iq_number_type& __iqNumber ) {
      lattice_ = &__lattice;
      iqNumber_ = &__iqNumber;
    }
    
    
    /**
     *  @brief  Release the iterator from its lattice and number
     *  @note   Can't see why this is really necessary, but provided for
     *          completeness.
     */
    void release( void ) {
      lattice_ = 0;
      iqNumber_ = 0;
    }
    
    
  private:
    /**
     *  @brief  Moves the iterator index forwards by 1
     *  @note   Required by @c iterator_facade
     */
    void increment( void ) {
      if(    lattice_->latticeType_
          == lattice_type::prv_lattice_enum::congruent
        )
      {
        if( (iqNumber_->real() += 2) >= lattice_->real().stop() ) {
          // This might look a bit lazy, but in all likelihood, the overhead
          // incurred from trying to be clever would probably all but negate
          // any speed increase.
          iqNumber_->real() = lattice_->real().start();
          ++iqNumber_->imag();
          if(    even_func()(iqNumber_->real())
              == odd_func()(iqNumber_->imag()) )
          {
            // The only strange case is when the real interval is only wide
            // enough to allow one integral value (this never happens when
            // searching on a full Siegel container, but may do if it is
            // subdivided quite finely). In this case only the imaginary
            // part gets iterated (it's slightly inefficient, but nothing
            // close to a bottleneck, and does not need optimising).
            if( lattice_->real().size() == 1 ) {
              ++iqNumber_->imag();
            }
            else {
              ++iqNumber_->real();
            }
          }
        }
      }
      else {
        if( ++(iqNumber_->real()) == lattice_->real().stop() ) {
          iqNumber_->real() = lattice_->real().start();
          ++(iqNumber_->imag());
        }
      }
    }
    
    
    /**
     *  @brief  Moves the iterator index backwards by 1
     *  @note   Required by @c iterator_facade
     */
    void decrement( void ) {
      if(    lattice_->latticeType_
          == lattice_type::prv_lattice_enum::congruent
        )
      {
        if( (iqNumber_->real() -= 2) < lattice_->real().start() ) {
          iqNumber_->real() = lattice_->real().stop() - 1;
          --(iqNumber_->imag());
          if(    even_func()(iqNumber_->real())
              == odd_func()(iqNumber_->imag()) )
          {
            if( lattice_->real().size() == 1 ) {
              --iqNumber_->imag();
            }
            else {
              --iqNumber_->real();
            }
          }
        }
      }
      else {
        if( --(iqNumber_->real()) < lattice_->real().start() ) {
          iqNumber_->real() = lattice_->real().stop() - 1;
          --(iqNumber_->imag());
        }
      }
    }
    
    
    /**
     *  @brief  Checks for equality between two iterators
     *  @param  that An iterator
     *  @note   Required by @c iterator_facade
     */
    bool equal( const self_type& that ) const {
      // The lattice points and the lattices are the same
      return bool(    ( *iqNumber_ == *(that.iqNumber_) )
                   && ( lattice_ == that.lattice_ )
                 );
    }
    
    
    /**
     *  @name   dereference
     *  @brief  Returns a reference to the current bound
     *  @note   Required by @c iterator_facade
     */
    reference dereference( void ) const {
     
      // Make sure at debug time that the lattice points are valid
      assert(      lattice_->latticeType_
                == lattice_type::prv_lattice_enum::congruent
                    ?    even_func()(iqNumber_->real())
                      ==
                         even_func()(iqNumber_->imag()) : true
            );
             
      return *iqNumber_;
    }
    
    
    //! Pointer to the lattice to iterate through
    const lattice_type* lattice_;
    //! Pointer to the iq number that iterates through the lattice
    iq_number_type* iqNumber_;
    
};

}// iterator
}// lattice
}// geometric
}// structure
}// sg
#endif
