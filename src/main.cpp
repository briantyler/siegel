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

#ifndef _SG_MAIN_H
#define _SG_MAIN_H 1

#define NDEBUG

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

// Global includes
#include <iostream>
#include <fstream>
#include <sstream>

#include <boost/algorithm/string.hpp>

#ifdef HAVE_GMP
#include <gmpxx.h>
#endif

// Local includes
#include "utility/functors/stream_cast.hpp"
#include "initialization.hpp"

namespace sg {
  template <class _Float, class _Integer> class settings {
    public:
      typedef settings<_Float,_Integer> self_type;
      typedef _Float float_type;
      typedef _Integer integer_type;
      typedef std::size_t size_type;
      
      /**
       * @brief Default Constructor
       * @param argc The size of \c argv
       * @param argv An array of inputs which set up the application
       */
      settings( int argc, char* argv[] )
      : dimension_(2), generator_(-1), balance_(0.9), match_(0.9),
        resolution_(1ul), first_(0ul), count_(0ul), sample_(100000),
        sieve_('m'), filename_("/dev/null"), help_(false),
        getResolution_(false)
        #ifdef __GMP_PLUSPLUS__
        , useGmp_(false)
        #endif
      {
        using boost::find_first;
        // Process the inputs
        for(int i = 1; i != argc; i++) {
          if(  parse_input( argv[i], "-g=", generator_ )
            || parse_input( argv[i], "-d=", dimension_ )
            || parse_input( argv[i], "-b=", balance_ )
            || parse_input( argv[i], "-m=", match_ )
            || parse_input( argv[i], "-r=", resolution_ )
            || parse_input( argv[i], "-f=", first_ )
            || parse_input( argv[i], "-c=", count_ )
            || parse_input( argv[i], "-s=", sample_ )
            || parse_input( argv[i], "-o=", filename_ )
            || parse_input( argv[i], "--sieve=", sieve_ )
          )
          { continue; }
          else if( find_first( argv[i], "--get-resolution" ) ) {
            getResolution_ = true;
            continue;
          }
          else if( find_first( argv[i], "--help" ) ) {
            help_ = true;
            continue;
          }
          #ifdef __GMP_PLUSPLUS__
          else if( find_first( argv[i], "--use-gmp" ) ) {
            useGmp_ = true;
            continue;
          }
          #endif
        }
      }
      
      integer_type generator( void ) const { return generator_; }
      size_type dimension( void ) const { return dimension_; }
      float_type balance( void ) const { return balance_; }
      float_type match( void ) const { return match_; }
      size_type resolution( void ) const { return resolution_; }
      size_type first( void ) const { return first_; }
      size_type count( void ) const { return count_; }
      size_type sample( void ) const { return sample_; }
      char sieve( void ) const { return sieve_; }
      std::string filename( void ) const { return filename_; }
      bool help( void ) const { return help_; }
      bool get_resolution( void ) const { return getResolution_; }
      bool use_gmp( void ) const { return useGmp_; }
      
    private:
      /**
       *  @brief  Parses the input variables.
       *  @param  arg The input string
       *  @param  param The input identifier
       *  @param  value A reference to the input variable.
       *  @return True if \c param was found in \c arg.
       */
      template <class _Type>
      bool parse_input ( const char* __arg, const char* __param, _Type& __value )
      {
        using namespace std;
        using namespace boost;
        typedef utility::functors::stream_cast<string,_Type> _sc;
        
        bool b = false;
        
        string arg( __arg );
        if( find_first( arg, __param ) ) {
          string value = erase_all_copy( arg, __param ) ;
          _sc()( value, __value );
          b = true;
        }
        
        return b;
      }
      
      size_type dimension_;
      integer_type generator_;
      float_type balance_;
      float_type match_;
      size_type resolution_;
      size_type first_;
      size_type count_;
      size_type sample_;
      char sieve_;
      std::string filename_;
      bool help_;
      bool getResolution_;
      #ifdef __GMP_PLUSPLUS__
      bool useGmp_;
      #else
      static const bool useGmp_ = false;
      #endif
  };
  
  /**
   *  @brief  Initiates the search.
   *  @param  Dimension the complex hyperbolic dimension to search in.
   */
  template <std::size_t _Dimension, std::size_t _Depth, class _Settings>
      void go( _Settings __settings )
  {
    using namespace std;
    using namespace boost;
    
    static const std::size_t dimension = _Dimension;
    static const std::size_t depth = _Depth;
    
    // Display the computed resolution
    if( __settings.get_resolution() ) {
      sg::get_resolution<dimension>
        ( __settings.generator(), __settings.resolution() );
    }
    // Compute the Siegel Set
    else {
      ofstream output;
      try {
        // Try to open the file.
        output.open( __settings.filename().c_str(), ios::out | ios::trunc );
        if( !output.is_open() ) {
          // The output file either has not been supplied, or is invalid,
          // so use the default value which is "/dev/null"
          output.open( "/dev/null", ios::out | ios::trunc );
        }
    
        // If we are unable to open an output file then abort.
        if( !output.is_open() ) {
          cout << "Unable to open output file, aborting..." << endl;
        }
        
        // Compute the Siegel set.
        sg::compute_siegel<depth,dimension>
        ( output, __settings.generator(), __settings.balance(), __settings.match(),
          __settings.resolution(), __settings.first(), __settings.count(),
          __settings.sample()
        );
        
        // Close the output file.
        output.close();
      }
      catch( std::exception& e ) {
        cout << "There was an error initializing the algorithm.\nDetails: "
             << e.what() << endl;
        // Close the output file.
        output.close();
      }
    }
  }
  
  
  /**
  *  @brief  Initiates the search.
  *  @param  Dimension the complex hyperbolic dimension to search in.
  */
  template <std::size_t _Depth, class _Settings>
  void go_sieve( _Settings __settings )
  {
    using namespace std;
    static const std::size_t depth = _Depth;
    
    switch( __settings.dimension() ) {
      case 2:
        go<2,depth,_Settings>(__settings);
        break;
      case 3:
        go<3,depth,_Settings>(__settings);
        break;
      case 4:
        go<4,depth,_Settings>(__settings);
        break;
      case 5:
        go<5,depth,_Settings>(__settings);
        break;
      case 6:
        go<6,depth,_Settings>(__settings);
        break;
      case 7:
        go<7,depth,_Settings>(__settings);
        break;
      case 8:
        go<8,depth,_Settings>(__settings);
        break;
      default:
        cout << "The dimension must be between 2 and 8 inclusive" << endl;
    }
  }
  
}

int main( int argc, char* argv[] ) {
  using namespace std;
  using sg::go_sieve;
  
  cout << "siegel  Copyright (C) 2009  Brian Tyler\n"
          "This program comes with ABSOLUTELY NO WARRANTY.\n"
          "This is free software, and you are welcome to redistribute it\n"
          "under certain conditions.\n"
          "For more information see <http://www.gnu.org/licenses/>"
       << endl;
  
  typedef sg::settings<double,long> native_settings;
  native_settings nativeSettings( argc, argv );
  
  #ifdef __GMP_PLUSPLUS__
  typedef sg::settings<mpf_class,mpz_class> gmp_settings;
  gmp_settings gmpSettings( argc, argv );
  #endif
  
  // Display help information
  if( nativeSettings.help() ) {
    sg::help_information();
    return 0;
  }
  
  if( !nativeSettings.use_gmp() ) {
    switch( nativeSettings.sieve() ) {
      case 'n':
        go_sieve<1,native_settings>(nativeSettings);
        break;
      case 't':
        go_sieve<2,native_settings>(nativeSettings);
        break;
      case 's':
        go_sieve<4,native_settings>(nativeSettings);
        break;
      case 'm':
        go_sieve<8,native_settings>(nativeSettings);
        break;
      case 'l':
        go_sieve<12,native_settings>(nativeSettings);
        break;
      default:
        cout << "The sieve must be one of; n, t, s, m, l" << endl;
    }
  }
#ifdef __GMP_PLUSPLUS__
  // It probably isn't a good idea to have many options for GMP types as they
  // are so slow.
  else {
    switch( gmpSettings.sieve() ) {
      case 'n':
        go_sieve<1,gmp_settings>(gmpSettings);
        break;
      case 't':
        go_sieve<2,gmp_settings>(gmpSettings);
        break;
      case 's':
        go_sieve<4,gmp_settings>(gmpSettings);
        break;
      case 'm':
        go_sieve<8,gmp_settings>(gmpSettings);
        break;
      default:
        cout << "The sieve must be either n, t, s, m" << endl;
    }
  }
#endif
  
  return 0;
}

#endif
