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
 *  @file     cusp_engine.hpp
 *  @brief    An inline header file for the \c cusp_engine class.
 *  @note     The engine for building cusps
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-02
 */


#ifndef _SG_CUSP_ENGINE_H
#define _SG_CUSP_ENGINE_H 1

// Global includes
#include <cassert>
#include <cmath>
#include <climits>
#include <algorithm>
#include <functional>
#include <set>

#include <boost/array.hpp>
#include <boost/numeric/ublas/storage.hpp>

// Local includes
#include "structure/geometric/hyperbolic/cusp.hpp"
#include "structure/geometric/hyperbolic/hyperbolic_space.hpp"
#include "structure/geometric/lattice/zeta_lattice.hpp"

#include "utility/math/is_square.hpp"
#include "utility/math/is_even.hpp"
#include "utility/math/congruence_system.hpp"
#include "utility/math/rounding.hpp"
#include "utility/math/square.hpp"
#include "utility/math/abs.hpp"
#include "utility/functors/stream_cast.hpp"
#include "utility/math/is_zero.hpp"



namespace sg
{
namespace engine
{
/**
 *  @class    cusp_engine
 *  @brief    An engine for generating cusps.
 *  @param    N The dimension of the space
 *  @param    _Float A floating point type; defaults to \c double
 *  @param    _Integer An integral type; defaults to \c long
 *  @param    _Id The identifier of the field
 *  @author   Brian Tyler
 *  @version  1.1
 *  @date     2008-04-02
 *  @note     \c v1.1 \c 2008-06-05 Modified the engine so that the threshold
 *            and bound does not depend on the norm even when the class number
 *            is non-trivial.
 *  @note     \c 2009-04-30 I discovered an issue with the code; if the lower
 *            bound of the height component of the space is non-zero, then when
 *            the dilation factor of the cusp gets bigger than (2/h)^2, the
 *            bounds can no longer be satisfied and the cusp generator sits in
 *            an infinite loop. As I am about two days away from submitting my
 *            phd I kind of don't want to get involved with this as it will
 *            mean thinking about how the siegel_engine copes with this; the
 *            thing to do is probably to set a flag which says that no more
 *            cusps can be created.
 */
template < std::size_t N, class _Float = double, class _Integer = long,
           std::size_t _Id = 0> class cusp_engine
{
  public:
    //! The object type
    typedef cusp_engine<N,_Float,_Integer,_Id> self_type;
    //! The float type
    typedef _Float float_type;
    //! The integral type
    typedef _Integer integer_type;
    //! The cusp type
    typedef structure::geometric::hyperbolic::cusp<N,float_type,integer_type,_Id>
        cusp_type;
    //! The complex
    typedef typename cusp_type::complex_type complex_type;
    //! The space type
    typedef structure::geometric::hyperbolic::hyperbolic_space<N,float_type>
        space_type;
    //! The field type
    typedef typename cusp_type::field_type field_type;
    
    // Keep public for testing
    //! The zeta lattice type
    typedef structure::geometric::lattice::zeta_lattice< N,float_type,
                                                         integer_type,_Id >
        lattice_type;
    //! The iterator which moves through the lattice
    typedef typename lattice_type::lattice_iterator lattice_iterator_type;
    //! Slice type for generating the r component
    typedef boost::numeric::ublas::basic_slice<integer_type,integer_type>
        slice_type;
    typedef typename slice_type::const_iterator slice_iterator_type;
    //! The imaginary quadratic number type
    typedef typename cusp_type::zeta_type iq_number_type;
    
  private:
    //! The type of a set of imaginary quadratic integers
    typedef std::set<iq_number_type> iqn_set_type;
    
    
  public:
    /**
     *  @brief  Default constructor
     */
    cusp_engine( void )
  : space_(0), cusp_(), lattice_(), latticeIterator_()
    { cusp_.dilation() = 1; }
    
    
    /**
     *  @brief  Bind the a space to the engine
     *  @param  space The space to bind
     */
    void bind_space( const space_type& __space ) { space_ = &__space; }
    
    
    /**
     *  @brief  Is the engine finished building cusps?
     */
    bool finished() const { return finished_; }
    
    
    /**
     *  @brief  Initialize the space to begin generating cusps
     *  @param  dilation The starting dilation factor.
     */
    void initialize( integer_type __dilation = integer_type(1) ) {
      typedef utility::math::square<float_type> _sq;
      typedef utility::math::is_zero<float_type> _isz;
      typedef utility::math::floor<integer_type,float_type> _flr;
      
      // Ensure that the starting value is non-negative
      assert( __dilation >= 1 );
      
      finished_ = false;
      cusp_.dilation() = __dilation - 1;
      maxDilation_ =   _isz()(space_->height_ref().lower())
                     ? LONG_MAX
                     : _flr()(_sq()(2.0 / space_->height_ref().lower()));

      // Copy the space into the lattice
      lattice_.set_space( *space_ );
      
      // Move to the next valid dilation factor
      if( next_dilation() ) {
        // Move the r iterator back by one so that next time operator() is called
        // the engine moves to the first valid cusp.
        --rIterator_;
      }
    }
    
    
    /**
     *  @brief  Return the next logical cusp
     *  @return A reference to the cusp which has been generated
     */
    const cusp_type& operator( ) ( void ) {
      do {
        // Try to increment the r value, then zeta, then rotation, then dilation.
        // Note that since an "if statement" is guaranteed to be evaluated from
        // left to right, terminating as soon as it is satisifed, then this
        // performs as required (cheeky eh?)
        if( next_r() || next_zeta() || next_rotation() || next_dilation() ){ }
        else {
          // Oh dear! there are no more cusps that can be created from this
          // space
          break;
        }
        
        // When the class number is trivial, the ideal generated by the cusp
        // can be assumed to be the maximal order, so if this is not the case,
        // then reject the cusp.
        // In fact furthermore, if the ideal is principal then it can be assume
        // to be a maximal order, irrespective of the class number of the field.
        if(    (    field_type::instance().is_ufd()
                 || cusp_.ideal().is_principal_ideal()
               )
            && cusp_.ideal().is_maximal_order()
          )
        { break; }
        
      } while( true );
      
      return cusp_;
    }
    
    
    /**
     *  @brief  Return a reference to the most recently found cusp.
     *  @return The most recently found cusp.
     */
    const cusp_type& cusp ( void ) const { return cusp_; }
    
    
    /**
     *  @brief  Tex string representation of the group that the cusps are
     *          generated to proxy.
     */
    std::string tex_group( void ) const {
      std::stringstream ss;
      ss << "\\mathrm{SU}\\left(" << space_type::dimension_size
         << ",1;\\mathbb{Z}\\left[" << field_type::instance().tex_generator()
         << "\\right]\\right)";
      return ss.str();
    }
    
  private:
    /**
     *  @brief  Advance the dilation factor to the next logical value.
     *  @return Always evaluates to true.
     */
    bool next_dilation( void ) {
      // Increment the dilation factor
      do {
        next_logical_dilation();
        if( cusp_.dilation() > maxDilation_ ) {
          finished_ = true;
          return false;
        }
        post_dilation();
        
      // If there is no valid rotation for this dilation factor then
      // increment the dilation factor.
      } while( !next_rotation() );
      
      return true;
    }
    
    
    /**
     *  @brief  Advance the rotation factor to the next logical value.
     *  @return Returns true if the increment was successful and false if the
     *          final rotation factor in this loop has been reached.
     */
    bool next_rotation( void ) {
      // RVO
      bool b = true;
      
      do {
        if( !next_logical_rotation() ) {
          b = false;
          return b;
        }
        post_rotation();
      
        // If there is no valid zeta for this rotation factor then
        // increment the rotation factor.
      } while( !next_zeta() );
      
      return b;
    }
    
    
    /**
     *  @brief  Advance the zeta factor to the next logical value.
     *  @return Returns true if the increment was successful and false if the
     *          final zeta factor in this loop has been reached.
     */
    bool next_zeta( void ) {
      // RVO
      bool b = true;
      
      do {
        if( !next_logical_zeta() ) {
          b = false;
          return b;
        }
        post_zeta();
        
        // If there is no valid r for this zeta factor then
        // recursively increment the zeta factor.
      } while( !next_r() );
      
      return b;
    }
    
    
    /**
     *  @brief  Advance the r factor to the next logical value.
     *  @return Returns true if the increment was successful and false if the
     *          final r factor in this loop has been reached.
     */
    bool next_r( void ) {
      // RVO
      bool b = true;
      
      if( !next_logical_r() ) {
        b = false;
        return b;
      }
      post_r();
      
      return b;
    }
    
    
  public:
    //! Get a reference to the lattice for testing
    lattice_type& lattice( void ) { return lattice_; }
    
    
    // Leave these public for testing
    /**
     *  @brief  Reset the dilation.
     *  @note   For testing.
     */
    void reset_dilation( void ) { cusp_.dilation() = 1; }
    
    /**
     *  @brief  Advance the dilation factor to the next logical value.
     */
    void next_logical_dilation( void ) { ++( cusp_.dilation() ); }
    
    
    /**
     *  @brief Perform post dilation calculations
     */
    void post_dilation( void ) {
      typedef utility::math::sqrt<float_type,integer_type> _sqrt;
      typedef utility::math::square<integer_type> _sq;
      typedef utility::functors::stream_cast<integer_type,float_type> _sc;
      typedef utility::math::is_square<integer_type> __issq;
      
      // Set the rotation variables
      rotation_.square
          =   field_type::instance().is_congruent()
            ? 4 * cusp_.dilation()
            : cusp_.dilation();
      
      _sc()( cusp_.dilation(), rotation_.fltDilation );
      
      // Set the zeta lattice bound variables for this loop
      lattice_.bound().dilation() = cusp_.dilation();
      lattice_.bound().height() = space_->height_common().lower();
      
      
      // Compute all rotations possible for this dilation factor
      rotation_.rotations_.clear();
      rotation_.deduplicated_.clear();
      iq_number_type rotationTmp(0,0);
      while( rotation_.square >= 0 ) {
        if( __issq( )( rotation_.square, rotationTmp.real() ) ) {
          rotation_.rotations_.insert( rotationTmp );
          rotation_.rotations_.insert( rotationTmp.conj() );
        }
        rotation_.square +=   field_type::instance().generator()
                            * ( 2 * ( rotationTmp.imag() ) + 1 );
        ++( rotationTmp.imag() );
      }
      // Now we have a complete set of rotations
      
      static const iq_number_type unitsnc[] = {
        iq_number_type(-1,0)
      };
      
      static const iq_number_type unitsc[] = {
        iq_number_type(-2,0)
      };
      
      static const iq_number_type units1[] = {
        iq_number_type(-1,0), iq_number_type(0,1),iq_number_type(0,-1)
      };
      
      static const iq_number_type units3[] = {
        iq_number_type(-2,0), iq_number_type(1,1), iq_number_type(-1,1),
        iq_number_type(1,-1), iq_number_type(-1,-1)
      };
      
      static const field_type& field = field_type::instance();
      // Get the right set of units for the ring under investigation
      const iq_number_type* units
          = ( field.generator() == -1 ? units1
              : ( field.generator() == -3 ? units3
                  : ( field.is_congruent() ? unitsc : unitsnc )
                )
            );
      
      const int numUnits
          = field.generator() == -1 ? 3 : ( field.generator() == -3 ? 5 : 1 );
      
      // Go through all of the rotations and store a representative for each
      // class (after quotienting out by the group of units)
      typename iqn_set_type::const_iterator rIt = rotation_.rotations_.begin();
      for( ; rIt != rotation_.rotations_.end(); ++rIt ) {
        bool insert = true;
        for( int i = 0; i != numUnits; ++i ) {
          if( rotation_.deduplicated_.count( units[i] * (*rIt) ) ) {
            insert = false;
            break;
          }
        }
        if( insert ) {
          rotation_.deduplicated_.insert( *rIt );
        }
      }
      
      rotation_.iterator_ = rotation_.deduplicated_.begin();
    }
    
    
    /**
     *  @brief  Advance the rotation factor to the next logical value.
     *  @return Returns true if the increment was successful and false if the
     *          final rotation in this loop has already been reached.
     */
    bool next_logical_rotation( void ) {
      // RVO
      bool b = true;
      
      if( rotation_.iterator_ == rotation_.deduplicated_.end() ) {
        b = false;
      }
      else {
        cusp_.rotation() = *(rotation_.iterator_);
        ++( rotation_.iterator_ );
      }
      
      return b;
    }
    
    
    /**
     *  @brief Perform post rotation calculations
     */
    void post_rotation( void ) {
      r_.c10 = field_type::instance().generator() * cusp_.rotation().imag();
      r_.c11 = cusp_.rotation().real();
      r_.zetaMod = std::conj( cusp_.rotation().to_complex() );
      
      // The space and bound are already set, so the only thing left to do is
      // set the transform
      lattice_.set_transform( cusp_.rotation().to_complex() );
      
      // Initialize the lattice
      lattice_.initialize();
      
      // Put the iterator back to the begining of the lattice
      latticeIterator_ = lattice_.lattice_begin( cusp_ );
      latticeIterator_.first_time();
    }
    
    
    /**
     *  @brief  Advance the zeta factor to the next logical value.
     *  @return Returns true if the increment was successful and false if the
     *          final zeta in this loop has been reached.
     */
    bool next_logical_zeta( void ) {
      // Increment the iterator
      ++latticeIterator_;
      
      // Move to the next valid zeta value
      // Keep incrementing the iterator until either the end of the iteeffectn
      // has been reached, or until the point is valid
      while( latticeIterator_ && !lattice_.validate() ) { ++latticeIterator_; }
      
      return bool( latticeIterator_ );
    }
    
    
    /**
     *  @brief  Perform post zeta calculations
     *  @note   This prepares the engine to solve the congruence system for r.
     */
      void post_zeta( void ) {
      typedef utility::math::congruence_system_solver<integer_type> _solver;
      typedef utility::math::ceil<integer_type,float_type> _ceil;
      typedef utility::math::floor<integer_type,float_type> _floor;
      typedef utility::math::is_even<integer_type> _even;
      typedef utility::math::abs<float_type> _abs;
      typedef utility::functors::stream_cast<integer_type,float_type> _sc;
      
      typename _solver::system_type system;
      
      integer_type qf = cusp_.compute_inner_qf();
      
      // The quadratic from must be even in the non-congruent case as it is
      // divisible by 2.
      if( !field_type::instance().is_congruent() && !_even()(qf) ) {
        r_.solution.xN = 0;
        return;
      }
      
      if( !field_type::instance().is_congruent() ) qf /= 2;
      
      integer_type c00 = cusp_.rotation().real() * qf;
      integer_type c01 = cusp_.rotation().imag() * qf;
      
      if( field_type::instance().is_congruent() ) {
        system.add_equation( r_.c10, c00, 2 * cusp_.dilation() );
        system.add_equation( r_.c11, c01, 2 * cusp_.dilation() );
        system.add_equation( r_.c10 + r_.c11, c00 + c01, 4 * cusp_.dilation() );
      }
      else {
        system.add_equation( r_.c10, c00, cusp_.dilation() );
        system.add_equation( r_.c11, c01, cusp_.dilation() );
      }
      
      r_.solution = _solver()( system );
      
      // There is no solution so fail
      if ( r_.solution.xN == 0 ) return;
      
      // There is a solution
      compute_r_bound();
      
      float_type fltxN, fltx0;
      _sc()( r_.solution.xN, fltxN );
      _sc()( r_.solution.x0, fltx0 );
      
      // Account for a bit of error at each end.
      integer_type rMin = _floor()( ( r_.boundMin - _abs()(fltx0) ) / fltxN );
      integer_type rMax = _ceil()( ( r_.boundMax + _abs()(fltx0) ) / fltxN );
      integer_type size = rMax - rMin + 1;
      
      r_.slice_ = slice_type( rMin * r_.solution.xN + r_.solution.x0,
                              r_.solution.xN, size );
      
      rIterator_ = r_.slice_.begin();
    }
    
    
    /**
     *  @brief  Compute the bounds on the r component
     */
    void compute_r_bound( void ) {
      typedef utility::math::square<float_type> _sq;
      typedef utility::math::abs<float_type> _abs;
      typedef typename cusp_type::const_zeta_ref_iterator _cusp_it ;
      typedef typename space_type::const_zeta_ref_iterator _space_it ;
      typedef typename space_type::zeta_type::const_iterator _corner_it;
      
      r_.boundMin = r_.boundMax = 0.0;
      
      // Compute the amount of bound remaining after accounting for the real
      // component of the height effect function
      r_.boundMin = -( r_.boundMax = lattice_.bound().r_bound() );
      
      // Account for the bounds on r itself
      r_.boundMax += rotation_.fltDilation * space_->r_common().upper();
      r_.boundMin += rotation_.fltDilation * space_->r_common().lower();
      
      
      // Now account for the bound imposed by the imaginary part of height
      // effect function.
      float_type zetaMin(0.0), zetaMax(0.0);
      
      // Iterate through each zeta component and complex region simultaneously
      _cusp_it cuspIt = cusp_.zeta_ref_begin();
      _space_it spaceIt = space_->zeta_ref_begin();
      for( ; cuspIt != cusp_.zeta_ref_end(); ++cuspIt, ++spaceIt ) {
        // For each region we need to find the maximum and minimum over all corners
        complex_type zeta( cuspIt->to_complex() * r_.zetaMod );
        bool firstTime = true;
        
        for( _corner_it cornerIt = spaceIt->begin();
             cornerIt != spaceIt->end(); ++cornerIt )
        {
          // = Im <z, conj( cusp::rotation ) * cusp::zeta )>
          float_type value( zeta.real() * cornerIt->imag()
                            - zeta.imag() * cornerIt->real() );
          if( firstTime ) {
            // Rather than messing around precomputing one variable, it is much
            // easier to have a "first time" variable. Checking a bool is
            // probably the fastest thing that can be done, so this is not a
            // bottleneck since it is insignificant in comparison to the other
            // operations performed.
            zetaMin = zetaMax = value;
            firstTime = false;
          }
          else if( value < zetaMin ) {
            zetaMin = value;
          }
          else if ( zetaMax < value ) {
            zetaMax = value;
          }
        }
        
        r_.boundMax += zetaMax;
        r_.boundMin += zetaMin;
      }
      
      r_.boundMin /= field_type::instance().sqrt_generator();
      r_.boundMax /= field_type::instance().sqrt_generator();
      
      if( field_type::instance().is_congruent() ) {
        r_.boundMin *= 2.0;
        r_.boundMax *= 2.0;
      }
    }
    
    
    /**
     *  @brief  Advance the r factor to the next logical value.
     *  @return Returns true if the increment was successful and false if the
     *          final r in this loop has been reached.
     */
    bool next_logical_r( void ) {
      // RVO
      bool b = true;
      
      if ( r_.solution.xN == 0 || rIterator_ == r_.slice_.end() ) {
        b = false;
      }
      else {
        cusp_.r_ref() = *rIterator_;
        ++rIterator_;
      }
      
      return b;
    }
    
    
    /**
     *  @brief  Perform post r calculations
     *  @note   This finalises the cusp. After calling this function the cusp
     *          is valid.
     */
    void post_r( void ) { cusp_.initialize(); }
    
    
  private:
    //! variables required to caclulate rotation
    struct prv_rotation {
      integer_type square;
      float_type fltDilation;
      iqn_set_type rotations_;
      iqn_set_type deduplicated_;
      typename iqn_set_type::const_iterator iterator_;
    };
    
    
    //! variables required to calculate r
    struct prv_r {
      integer_type c10;
      integer_type c11;
      float_type boundMin;
      float_type boundMax;
      complex_type zetaMod;
      utility::math::congruence_solution<integer_type> solution;
      slice_type slice_;
    };
    
    
    //! The space which the cusp is required to affect
    const space_type* space_;
    //! The current generated cusp
    cusp_type cusp_;
    //! The maximum dilation factor which is attainable
    integer_type maxDilation_;
    //! Have all possible cusps that could affect the space have been created?
    bool finished_;
    //! The lattice which builds the zeta components
    lattice_type lattice_;
    //! The iterator which traverses the lattice
    lattice_iterator_type latticeIterator_;
    //! The slice iterator for generating the r values
    slice_iterator_type rIterator_;
    //! The variables required for calculating the rotation component
    prv_rotation rotation_;
    //! The variables required for calculation the r component
    prv_r r_;
};

}// engines
}// sg
#endif
