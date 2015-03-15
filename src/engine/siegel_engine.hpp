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
 *  @file     siegel_engine.hpp
 *  @brief    An inline header file for the \c siegel_engine class.
 *  @note     The engine which constructs siegel sets.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-05-26
 */


#ifndef _SG_SIEGEL_ENGINE_H
#define _SG_SIEGEL_ENGINE_H 1

// Global includes
#include <exception>
#include <algorithm>
#include <iostream>

#include <boost/progress.hpp>
#include <boost/call_traits.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>

// Local includes
#include "engine/cusp_engine.hpp"
#include "engine/cusp_selector.hpp"
#include "structure/geometric/hyperbolic/heisenberg_slice.hpp"
#include "utility/math/round.hpp"
#include "geometry/effect/point_phi.hpp"



namespace sg
{
namespace engine
{
/**
 *  @class    siegel_engine
 *  @note     The engine which constructs siegel sets.
 *  @param    N The dimension of the space
 *  @param    Depth The depth of the sieve to use for cusp selection
 *  @param    Float A floating point type; defaults to \c double
 *  @param    Integer An integral type; defaults to \c long
 *  @param    Id The identifier of the field
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-05-26
 */
template < std::size_t _Depth, std::size_t N, class _Float = double,
           class _Integer = long, std::size_t _Id = 0 >
    class siegel_engine
{
  public:
    //! The object type
    typedef siegel_engine<_Depth,N,_Float,_Integer,_Id> self_type;
    //! The float type
    typedef _Float float_type;
    //! The integral type
    typedef _Integer integer_type;
    //! The size type of the slice
    typedef std::size_t size_type;
    //! The type of the cusp engine
    typedef cusp_engine<N,float_type,integer_type,_Id> cusp_engine_type;
    //! The cusp selector engine type
    typedef cusp_selector<_Depth,N,float_type,integer_type,_Id> selector_type;
    //! The heisenberg slice type
    typedef structure::geometric::hyperbolic::heisenberg_slice
            <N,float_type,integer_type> slice_type;
    //! The cusp type
    typedef typename cusp_engine_type::cusp_type cusp_type;
    //! The space type
    typedef typename cusp_engine_type::space_type space_type;
    //! The point type
    typedef typename space_type::point_type point_type;
    //! The field type
    typedef typename cusp_engine_type::field_type field_type;
    
    //! Failure class
    class engine_fail: public std::exception
    {
      virtual const char* what() const throw()
      {
        return "The resolution needed to satisfy the maxmal height condition "
               "is too great to computationally find a valid height bound for "
               "this space";
      }
    };
    
    //! Template constants
    enum constants {
      depth = _Depth,
      dimension = N
    };
    
    
  private:
    //! The point_phi functor finds the maximum height of a hypercubecube
    typedef geometry::effect::point_phi<point_type,slice_type,cusp_type> phi_func;
    
    //! Type of the mersenne twister
    typedef boost::mt19937 base_generator_type;
    //! Type of the probability distribution (uniform)
    typedef boost::uniform_int<std::size_t> distribution_type;
    //! Type of the random number generator
    typedef boost::variate_generator<base_generator_type&, distribution_type>
        generator_type;
    
    //! Type of the cusp selector iterator
    typedef typename selector_type::const_iterator selector_iterator;
    //! Type of the heisenberg slice iterator
    typedef typename slice_type::slice_iterator slice_iterator;
    
    //! Floating point param type
    typedef typename boost::call_traits<float_type>::param_type flt_param_type;
    //! Integer param type
    typedef typename boost::call_traits<integer_type>::param_type int_param_type;
    
    //! The progress bar type
    typedef boost::progress_display progress_type;
    //! The timer type
    typedef boost::timer timer_type;
    
    
  public:
    /**
     *  @brief  Default constructor
     */
    siegel_engine( )
  : space_(), engine_(), selector_(), slice_(), phi_( slice_ ), hthresh_(2.0),
    hcur_(2.0), match_(1.0), balance_(1.0), sample_( 1000000ul ),
    generator_( sys_base_generator(), distribution_type( 0ul, 1ul ) ), output_(0)
    { engine_.bind_space(space_); }
    
    
    /**
     *  @brief  Get a constant reference to the underlying space.
     *  @return A constant reference to the underlying space.
     */
    const space_type& space( void ) const { return space_; }
    
    /**
     *  @brief  Get a reference to the underlying space.
     *  @return A reference to the underlying space.
     */
    space_type& space( void ) { return space_; }
    
    
    /**
     *  @brief  Set the imaginary quadratic field generator.
     */
    void set_field_generator( int_param_type __generator ) {
      field_type::instance().initialize( __generator );
    }
    
    
    /**
     *  @brief  Set the output stream.
     */
    void set_output( std::ostream& __output ) { output_ = &__output; }
    
    
    /**
     *  @brief  Get a reference to the sample variable.
     *  @return A reference to sample.
     *  @note   \c sample determines how many random points to sample when
     *          estimating a minimum and maximum height, clearly the larger
     *          sample is, the more accurate the result, but the slower the
     *          algorithm during the sampling phase. Something like 10^6 gives
     *          a good balance between accuracy and speed.
     */
    size_type& sample( void ) { return sample_; }
    
    /**
     *  @brief  Get a copy of the sample variable.
     *  @return The sample.
     */
    size_type sample( void ) const { return sample_; }
    
    
    /**
     *  @brief  Get a reference to the match variable.
     *  @return A reference to match.
     *  @note   \c match should be in the set \f$ \[1,0) \f$. It describes how
     *          accurately the minimum height is matched to the maximum
     *          theoretical height for the current set of cusps. If it is too
     *          large the algorithm could be VERY slow and will likely fail
     *          due to running out of memory, convversely, if it is too small,
     *          then the result returned will probably be far from optimal.
     *          The best thing to do is to do a run with a small value to get
     *          some kind of minimum and then to increase this in subsequent
     *          runs to improve the approximation.
     */
    float_type& match( void ) { return match_; }
    
    /**
     *  @brief  Get a copy of the match variable.
     *  @return The match.
     */
    float_type match( void ) const { return match_; }
    
    
    /**
     *  @brief  Get a reference to the balance variable.
     *  @return A reference to balance.
     *  @note   \c balance is a bit like \c match except that it is only used
     *          during initialization. During this phase of the algorithm, an
     *          initial list of cusps is generated based on an estimated
     *          minimum height \c hmin and a maximum theoretical height \c hmax
     *          The \c balance determines the effect of these two numbers and
     *          should be in the set  \f$ \[1,0) \f$, the closer to 1 it is the
     *          more likely it is that the initial cusp set will generate a
     *          Siegel container (note that a effect of exactly 1 is not
     *          guaranteed to give a complete set. However, the closer it is to
     *          1 the slower the algorithm is going to be. It is not necessary
     *          to generate a complete list of cusps during initialization as
     *          more cusps can be added during the search phase. During cusp
     *          geneeffectn, \c hmin will stabilise at a certain value, \c hmax
     *          will continually decrease as \f$ 1 / \sqrt{\delta} \f$ where
     *          \f$ \delta \f$ is the dilation factor of the smallest cusp (i.e.
     *          it is the square root of the integer value in the cusp object).
     *          As such at some stage \f$ hmax < balance * hmin \f$ for any
     *          value of \c balance. However, if balance is close to 1 this
     *          may take a painfully long time, and most of the cusps that are
     *          generated will be totally unnecessary. A good plan seems to be
     *          to set balance ~ 0.1 - 0.2 and then let the search phase take
     *          over cusp geneeffectn. For higher dimensions and dilation
     *          factors one may want to settle on a lower balance than that.
     */
    float_type& balance( void ) { return balance_; }
    
    /**
     *  @brief  Get a copy of the balance variable.
     *  @return The balance.
     */
    float_type balance( void ) const { return balance_; }
    
    
    /**
     *  @brief  Get a copy minimum height variable.
     *  @return The minimum height
     */
    float_type min_height( void ) const { return hcur_; }
    
    
    /**
     *  @brief  Initialize the engine ready for use
     */
    void initialize( ) {
      using std::cout;
      using std::endl;
      using std::flush;
      
      // Start the clock for timing initialization
      timer_.restart();
      
      // Set the precision of the output stream
      output_->precision( utility::precision::stream() );
      
      cout << "Initialization" << endl;
      cout << "^^^^^^^^^^^^^^" << endl;
      cout << flush;
      
      *output_ << "Initialization" << endl;
      *output_ << "^^^^^^^^^^^^^^" << endl;
      *output_ << flush;
      
      
      // Initialize the space in case it has not already been done.
      space_.initialize();
      
      cout << "Space initialized!" << endl;
      *output_ << "Space initialized!" << endl;
      
      // Initialize the engine and move to the first valid cusp.
      engine_.initialize(1);
      engine_();
      
      cout << "Engine initialized!" << endl;
      *output_ << "Engine initialized!" << endl;
      
      // Clear all cusp data
      selector_.clear();
      
      // Set the cusp selector's space
      selector_.propagate( space_ );
      
      cout << "Selector initialized!\n" << endl;
      *output_ << "Selector initialized!\n" << endl;
      
      // Build the first lot of cusps, and keep building until at least all of
      // the slots in the selector are filled. Note that this does not mean
      // that we are guaranteed to have enough cusps to get a non-zero height
      // bound, but it does mean that we will always have at least one cusp
      // candidate for every hypercube in the slice and due to this there will
      // not be any iterator dereferencing errors.
      cout << "Building initial cusp set..." << endl;
      *output_ << "Building initial cusp set..." << endl;
      
      while( !selector_.complete() ) {
        sys_build_next_dilation();
      }
      
      // Set the slice space to the generator's space
      slice_.space() = space_;
      
      // Initialize the heisenberg slice with a very large resolution
      // Having a large resolution mitigates the error induced by discretising
      // the space.
      // 2^40 ~ 10^12 is about the computational limit for the resolution
      slice_.resolution() = static_cast<integer_type>(2ul << 40);
      slice_.initialize();
      
      
      // Put an estimated lower bound on the height.
      hcur_ = sys_estimate_height();
      // Put an upper bound on the height
      hthresh_ = sys_thresh_height( threshold_ );
      
      // Display the estimates
      output_estimates();
      
      if( hthresh_ <= 0.0 ) throw engine_fail();
      
      // Now try to find the right dilation factor to start validating the
      // height with. By balancing the maximum attainable height using the
      // current set of cusps with an actual estimate taken from a random
      // sampling, we hope to find roughly the right dilation factor and
      // roughly the right minimum height.
      while( hcur_ < balance_ * hthresh_ ) {
        sys_build_next_dilation();
        
        hcur_ = sys_estimate_height();
        hthresh_ = sys_thresh_height( threshold_ );
        
        // Display the estimates
        output_estimates();
        
        // Display the time taken so far
        sys_short_precision();
        cout << "Time elapsed: " << timer_.elapsed() << "s" << endl;
        *output_ << "Time elapsed: " << timer_.elapsed() << "s" << endl;
        sys_long_precision();
      }
      
      // Now that an estimate for the dilation factor and height has been
      // reached we put the resolution of the slice back to something more
      // managable.
      slice_.resolution() = static_cast<integer_type>(2ul << 24);
      slice_.initialize();
      
      cout << "Resolution set to: " << slice_.resolution() << "\n" << endl;
      *output_ << "Resolution set to: " << slice_.resolution() << "\n" << endl;
      
      
      // Move the silce iterator to its starting position
      sliceIterator_ = slice_.slice_begin();
      
      // Match the resolution to the minimum cusp threshold
      sys_match_resolution();
      hthresh_ = sys_thresh_height( threshold_ );
      
      // Reset the current height to the maximum value, it will then be
      // set on the first loop, or it can be set manually to a minimum value
      // before searching.
      hcur_ = 2.0;
      
      // Get the first cusp candidate for the starting cube
      cuspIterator_ = selector_.begin( *sliceIterator_ );
      
      sys_short_precision();
      cout << "Initialization complete.\n"
           << "Time elapsed: " << timer_.elapsed() << "s" << endl
           << "--\n\n" << endl;
      
      *output_ << "Initialization complete.\n"
               << "Time elapsed: " << timer_.elapsed() << "s" << endl
               << "--\n\n" << endl;
      sys_long_precision();
      
      init_ = timer_.elapsed();
    }
    
    
    /**
     *  @brief  Start the engine to find a Siegel Set.
     *  @return A lower bound on the height of a Siegel Set.
     */
    float_type operator()( void ) {
      using std::cout;
      using std::endl;
      
      sys_loop();
      
      std::cout << "Minimum height: " << hcur_ << endl;
      *output_ << "Minimum height: " << hcur_ << endl;
      
      // Display the time taken so far
      sys_short_precision();
      cout << "Total time elapsed: " << (init_ + timer_.elapsed())
           << "s" << endl;
      *output_ << "Total time elapsed: " << (init_ + timer_.elapsed())
               << "s" << endl;
      sys_long_precision();
      
      return hcur_;
    }
    
    
  private:
    /**
     *  @brief  Display and log the current height estimates
     */
    void output_estimates( void ) {
      using std::cout;
      using std::endl;
      
      cout << "Current height estimate: " << hcur_ << endl;
      cout << "Threshold height estimate: " << hthresh_ << endl;
        
      *output_ << "Current height estimate: " << hcur_ << endl;
      *output_ << "Threshold height estimate: " << hthresh_ << endl;
    }
    
    
    /**
     *  @brief  Begin the search loop.
     */
    void sys_loop( void ) {
      typedef utility::functors::stream_cast<integer_type,std::size_t> _sc;
      
      using std::cout;
      using std::endl;
      
      // Restart the timer
      timer_.restart();
      
      cout << "Search" << endl;
      cout << "^^^^^^" << endl;
      
      *output_ << "Search" << endl;
      *output_ << "^^^^^^" << endl;
      
      std::size_t tmp;
      _sc()( slice_.resolution(), tmp );
      progress_type progress( tmp );
      
      for( ;sliceIterator_ ;++progress )
      {
        
        float_type height = sys_compute_height( *sliceIterator_, cuspIterator_,
                                                hcur_
                                              );
        
        // The current cube can be raised above our current minimum so we can
        // move on to the next point.
        if( height >= hcur_ ) {
          ++sliceIterator_;
          continue;
        }
        
        // Now we have a problem: the maximum height attainable by this cube is
        // less than our current minimum height. There are two options:
        // 
        // 1) Accept the new height as the minimum and continue
        // 2) Build more cusps in the expectation that a new cusp of smaller
        //    dilation factor will be successful in raising the cube further.
        // 
        // So the question is how do we decide on which option to choose?
        
        // 1) If the height is greater than the maximum height attainable at
        // the largest dilation factor currently generated then we can reduce
        // the minimum height bound and continue (this maximum height is not an
        // absolute maximum, but it can be thought of as a good approximation
        // to it).
        
        
        // This is now something of a balancing act: Computing the cube_phi at
        // the current cube and threshold tells us its theoretical maximum
        // height of effect. If we have enough cusps built, then the cube_phi
        // and the height we have just computed should be fairly close. If they
        // are not close it suggests that we need to increase the dilation
        // factor of the smallest cusp.
        
        float_type phi = phi_.threshold_height( *sliceIterator_, threshold_ );
        
        // 2) Assume that there is a cusp of smaller dilation factor that can
        //    increase the height of the cusp. Calling sys_build_next_dilation:
        //      # Decreases threshold_
        //      # Decreases hthresh_
        //    Thus on the next loop there are lower bounds to satisfy. Calling
        //    sys_increase_resolution matches the resolution of the space to
        //    the new minimum cusp threshold, this will increase hthresh_
        //    slightly.
        if( height < match_ * phi ) {
          cout << "\n\nFail at this dilation factor, increasing...\n"
               << "Fail point:\n" << *sliceIterator_ << endl;
          
          *output_ << "\n\nFail at this dilation factor, increasing...\n"
                   << "Fail point:\n" << *sliceIterator_ << endl;
          
          sys_build_next_dilation();
          sys_match_resolution();
          
          // Display the estimates
          output_estimates();
          
          // Display the time taken so far
          sys_short_precision();
          cout << "Time elapsed: " << (init_ + timer_.elapsed())
               << " s" << endl;
          *output_ << "Time elapsed: " << (init_ + timer_.elapsed())
                   << " s" << endl;
          
          double percentage =   100.0 * static_cast<double>( progress.count() )
                              / static_cast<double>( progress.expected_count() );
          
          cout << "Percentage complete: " << percentage << "%\n" << endl;
          
          _sc()( slice_.resolution(), tmp );
          progress.restart( tmp );
          
          _sc()( sliceIterator_.index(), tmp );
          progress += tmp;
          
          percentage =   100.0 * static_cast<double>( progress.count() )
                       / static_cast<double>( progress.expected_count() );
          
          *output_ << "Percentage complete: " << percentage << "%\n" << endl;
          sys_long_precision();
        }
        // 1) Although the height is not as large as we would like, we know it
        //    is not being influenced too much by the error induced by
        //    discretising the space and it is close to the maximum we could
        //    expect from our smallest possible cusp, so we reduce our height
        //    expectations to this new value.
        else {
          hcur_ = height;
          ++sliceIterator_;
          continue;
        }
      }
      
      cout << "Search complete." << endl;
      *output_ << "Search complete." << endl;
    }
    
    
    /**
     *  @brief  Build the set of cusps at the next logical dilation factor.
     */
    void sys_build_next_dilation( void ) {
      using std::cout;
      using std::endl;
      
      // Record the dilation and the threshold
      integer_type dilation = engine_.cusp().dilation();
      threshold_ = engine_.cusp().threshold();
      
      do {
        // Note the cusp is valid even before calling engine_() so we just
        // leave the cusp with new dilation factor until the next time this
        // function is called.
        selector_.add_cusp( engine_.cusp() );
        
        // Move to the next cusp
        engine_();
        
        // Keep looping until a new dilation factor is reached
      } while( dilation == engine_.cusp().dilation() );
      
      cout << "Dilation factor increased to: " << selector_.dilation() << endl;
      cout << "Total cusps: " << selector_.cusps().size() << endl;
      
      *output_ << "Dilation factor increased to: " << selector_.dilation()
               << endl;
      *output_ << "Total cusps: " << selector_.cusps().size() << endl;
    }
    
    
    /**
     *  @brief  Match the resolution of the slice with the smallest threshold.
     */
    void sys_match_resolution( void ) {
      using std::cout;
      using std::endl;
      
      typedef utility::math::round<integer_type,float_type> _round;
      typedef utility::functors::stream_cast<integer_type,float_type> _sc;
      
      // The factor by which to increase the resolution of the space
      static const float_type factor = ::sqrt(10.0);
      static const float_type tfactor = 0.90;
      
      const float_type tollerance = tfactor * threshold_;
      
      hthresh_ = sys_thresh_height( threshold_ );
      
      // The resolution of the space is adequate for the current set of cusps.
      if( hthresh_ >= tollerance ) return;
      
      // Changing the resolution is not so straight forward as upping the
      // dilation factor because the cube iterator actually gets moved
      // backwards to a place in the slice where we have satisfied all previous
      // cubes. Increasing the resolution means that an individual cusp needs
      // to do less work, but it means that there is more work to do to get
      // through the space. As such it is better not to increase the resolution
      // if we can help it.
      
      --sliceIterator_;
      sliceIterator_.prepare_update();
      
      
      while( hthresh_ < tollerance ) {
        float_type tmpResolution;
        _sc()( slice_.resolution(), tmpResolution );
        
        integer_type r = _round()( tmpResolution * factor );
        
        slice_.resolution() = r;
        slice_.initialize();
        sliceIterator_.update();
        
        hthresh_ = sys_thresh_height( threshold_ );
      }
      
      cout << "Resolution increased to: " << slice_.resolution() << endl;
      *output_ << "Resolution increased to: " << slice_.resolution() << endl;
    }
    
    
    /**
     *  @brief  Estimates the maximum theoretical height attainable at the
     *          current resolution given a bound on the dilation factor.
     *  @param  threshold The threshold of the minimum dilation factor assumed
     *          to give the maximum height.
     *  @return The estimate to the maximum height.
     *  @note   The function samples a fairly large number of uniformly
     *          distributed points from the slice and takes the minimum height
     *          over all of them.
     */
    float_type sys_thresh_height( flt_param_type __threshold ) {
      point_type point = sys_random_point();
      float_type height = phi_.threshold_height( point, __threshold );
      float_type tmpheight = 0.0;
      for( size_type i = 0; i != sample_; ++i ) {
        point = sys_random_point();
        tmpheight = phi_.threshold_height( point, __threshold );
        
        if( tmpheight < height ) height = tmpheight;
        if( height < 0.0 ) return height;
      }
      
      return height;
    }
    
    
    /**
     *  @brief  Estimates the minimum height attainable with the
     *          current set of cusps at the current resolution.
     *  @return The estimate to the minimum height.
     *  @note   This function acts as a starting value for a height bound
     *          estimate, it does not give any kind of mathematically valid
     *          lower bound and the lower bound is very likely to be lower than
     *          the one returned.
     */
    float_type sys_estimate_height( void ) {
      static const float_type bound = 2.0;
      
      // Find a starting value for the minimum
      point_type point( sys_random_point() );
      
      selector_iterator iterator = selector_.begin( point );
      float_type height = sys_compute_height( point, iterator, bound );
      
      // Find a minimum over the sample number of cubes
      for( size_type i = 0; i != sample_; ++i ) {
        point = sys_random_point();
        iterator = selector_.begin( point );
        float_type tmp = sys_compute_height( point, iterator, bound );
        
        if( tmp < height ) height = tmp;
        // Fail: there is no need to calculate any further, the resolution must
        // be decreased, or the number of cusps increased in order to find a
        // valid height.
        if( height < 0.0 ) break;
      }
      
      return height;
    }
    
    
    /**
     *  @brief  Compute the maximum height attainable for a hypercube with the
     *          current set of cusps.
     *  @param  point The front vertex of the hypercube.
     *  @param  iterator The iterator pointing to candidate cusps
     *  @param  bound The bound at which to accept the height attained
     *  @return The maximum height of the cube such that there is a cusp that
     *          can move it upwards.
     */
    float_type sys_compute_height
        ( const point_type& __point, selector_iterator& __iterator,
          flt_param_type __bound )
    {
      // Our best guess for the cusp that is most likely to do the job is the
      // one that worked for the last cube.
      float_type height = phi_( __point, *__iterator );
      
      // We are looking to satisfy a bound on the height.
      if( height >= __bound  ) {
        // This cusp raises the cube sufficiently here, so we don't need to
        // look any further.
        return height;
      }
      
      selector_iterator tmpIterator = __iterator
                                    = selector_.begin( __point );
      selector_iterator endIterator = selector_.end( __iterator );
      
      // The search now begins for a cusp that can raise the cube beyond the
      // current minimum bound.
      for( ; __iterator != endIterator; ++__iterator )
      {
        float_type tmp = phi_( __point, *__iterator );
        if( tmp > height ) {
          // Record the cusp which gave the maximum height and record that
          // height.
          tmpIterator = __iterator;
          height = tmp;
        }
        
        // NOTE This is really important: __iterator.threshold() <-- note the
        // dot, this gives the effectiveness of the cusp in terms of it's
        // ability to raise the height, dereferenceing the iterator gives the
        // threshold with respect to the effect function and should not be used.
        if( height >= __iterator.threshold() || height >= __bound ) {
          // The cube has been raised to the maximum attainable height with
          // these cusps, or the height it has been raised to is good enough.
          break;
        }
      }
      
      // The end of the iteeffectn has been reached and there are no cusps
      // capable of raising the cube above either the bound, or the individual
      // thresholds.
      
      // Set the iterator to the position which gave the maximum height
      __iterator = tmpIterator;
      
      return height;
    }
    
    
    /**
     *  @brief  Get a random hypercube from the slice.
     *  @return A random hypercube from the heisenberg slice.
     */
    point_type sys_random_point( void ) const {
      // If the resolution changes then it is necessary to alter the bounds on
      // the random number generator.
      if(    static_cast<integer_type>(generator_.max())
          != slice_.resolution() - 1 )
      {
        typedef utility::functors::stream_cast<integer_type,std::size_t> _sc;
        std::size_t max;
        _sc()( slice_.resolution()-1, max );
        generator_.distribution() = distribution_type( 0, max );
      }
      
      return point_type( slice_.point_at( generator_() ) );
    }
    
    
    /**
     *  @brief  Get a constant reference to a static random number generator
     *  @return A constant reference to a static random number generator
     */
    base_generator_type& sys_base_generator( void ) const {
      static base_generator_type generator;
      return generator;
    }
    
    
    void sys_short_precision( void ) {
      std::cout.precision( 2 );
      output_->precision( 2 );
    }
    
    void sys_long_precision( void ) {
      std::cout.precision( utility::precision::stream() );
      output_->precision( utility::precision::stream() );
    }
    
    
    
    //! The space that is to become the Siegel set
    space_type space_;
    //! The engine which builds the cusps
    cusp_engine_type engine_;
    //! The current cusp threshold
    float_type threshold_;
    
    //! The cusp selector which stores all candiate cusps
    selector_type selector_;
    //! An iterator pointing to the current cusp candidate
    selector_iterator cuspIterator_;
    
    //! The heisenberg slice which performs the iteeffectn
    slice_type slice_;
    //! The slice iterator
    slice_iterator sliceIterator_;
    //! The point phi functor
    phi_func phi_;
    
    //! The height bound imposed by the threshold in the worst case scenario
    float_type hthresh_;
    //! The current lower bound on the height
    float_type hcur_;
    //! The factor with which to match hypercubes and cusps
    float_type match_;
    //! The initial factor used to the minimum and maximum heights
    float_type balance_;
    //! The amount of points to sample when computing a maximum/minimum height
    size_type sample_;
    
    //! Random number generator
    mutable generator_type generator_;
    //! The output stream to log information to
    std::ostream* output_;
    //! The timer for timing operations
    boost::timer timer_;
    //! Initialisation time
    double init_;
    
};

}// engine
}// sg
#endif
