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
#include <set>

#include <boost/algorithm/string.hpp>


// Local includes
#include "engine/cusp_candidates_engine.hpp"
#include "geometry/algorithms/make_siegel.hpp"


namespace sg {
  class settings {
    public:
      typedef settings self_type;
      typedef double float_type;
      typedef long integer_type;
      typedef std::size_t size_type;
      
      /**
       * @brief Default Constructor
       * @param argc The size of \c argv
       * @param argv An array of inputs which set up the application
       */
      settings( int argc, char* argv[] )
      : dimension_(2), generator_(-1), height_(1.0), tex_(false), store_(false),
        filename_("/dev/null"), help_(false)
      {
        using boost::find_first;
        // Process the inputs
        for(int i = 1; i != argc; i++) {
          if(  parse_input( argv[i], "-g=", generator_ )
            || parse_input( argv[i], "-d=", dimension_ )
            || parse_input( argv[i], "-o=", filename_ )
            || parse_input( argv[i], "-h=", height_ )
          )
          { continue; }
          else if( find_first( argv[i], "--tex" ) ) {
            tex_ = true;
            continue;
          }
          else if( find_first( argv[i], "--store" ) ) {
            store_ = true;
            continue;
          }
          else if( find_first( argv[i], "--help" ) ) {
            help_ = true;
            continue;
          }
        }
      }
      
      integer_type generator( void ) const { return generator_; }
      size_type dimension( void ) const { return dimension_; }
      float_type height( void ) const { return height_; }
      std::string filename( void ) const { return filename_; }
      bool tex( void ) const { return tex_; }
      bool store( void ) const { return store_; }
      bool help( void ) const { return help_; }
      
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
      float_type height_;
      bool tex_;
      bool store_;
      std::string filename_;
      bool help_;
  };
  
  
 /**
  *  @brief  Print help information to the screen.
  */
void help_information( void ) {
  using std::cout;
  using std::endl;
  
  cout << "Usage: siegelcl [options]\n"
          "Options:\n"
          "  -d=<arg>\tThe complex hyperbolic dimension of the space.\n"
          "\t\tDefault is 2.\n"
          "  -g=<arg>\tGenerator of the number field. Default is -1.\n"
          "  -h=<arg>\tThe lower height bound of the Siegel set, Default is 1.\n"
          "  -o=<arg>\tName of the output file. Default is /dev/null\n"
          "  --tex\t\tGenerate LaTeX representation of the data.\n"
          "  --store\tStore the cusps rather than just count them. This\n"
          "\t\toption is only useful when used in conjunction with the --tex\n"
          "\t\tflag, as this will output all the generated cusps in a format\n"
          "\t\tthat can be pasted straight into LaTeX. This option consumes\n"
          "\t\tmore memory, so don't use it for complicated Siegel sets.\n"
          "  --help\tDisplays this information.\n"
        << endl;
}

/**
 *  @struct bad_sg_input
 *  @brief  A lazy exception object
 */
struct bad_sg_input: public std::exception
{
  std::string parameters;
  
  virtual const char* what() const throw()
  {
    return parameters.c_str();
  }
  
  virtual ~bad_sg_input() throw() { }
};


/**
 *  @brief  Determine if the input generator is valid.
 *  @param  generator The generator of the number field.
 *  @note   Throws if \c generator is not a Heegner number.
 */
void validate_generator( long __generator ) {
  std::set<long> generators;
  generators.insert( -1 );
  generators.insert( -2 );
  generators.insert( -3 );
  generators.insert( -7 );
  generators.insert( -11 );
  generators.insert( -19 );
  generators.insert( -43 );
  generators.insert( -67 );
  generators.insert( -163 );
  
  if( generators.count( __generator ) == 0 ) {
    std::stringstream ss;
    ss << "Bad generator input, generator given: " << __generator
       << ". Generator must be a Heegner number (-1, -2, -3, -7, -11, -19,"
        " -43, -67, -163).";
    
    bad_sg_input ex;
    ex.parameters = ss.str();
    throw ex;
  }
}


void validate_dimension( std::size_t __dimension ) {
  std::set<std::size_t> dimensions;
  dimensions.insert( 2 );
  dimensions.insert( 3 );
  dimensions.insert( 4 );
  dimensions.insert( 5 );
  dimensions.insert( 6 );
  dimensions.insert( 7 );
  dimensions.insert( 8 );
  
  if( dimensions.count( __dimension ) == 0 ) {
    std::stringstream ss;
    ss << "Bad dimension input, dimension given: " << __dimension
        << ". Dimension must be 2,3,4,5,6,7 or 8.";
    
    bad_sg_input ex;
    ex.parameters = ss.str();
    throw ex;
  }
}

void validate_height( double __height ) {
  if ( __height <= 0.0 || 2.0 <= __height ) {
    std::stringstream ss;
    ss << "Bad height input, height given: " << __height
        << ". Height must be between 0 and 2.0.";
    
    bad_sg_input ex;
    ex.parameters = ss.str();
    throw ex;
  }
}


template< std::size_t Dimension> void generate( settings& __settings ) {
  using namespace std;
  using sg::geometry::algorithms::make_siegel;
  
  typedef engine::cusp_candidates_engine<Dimension,double,long,0>
      candidates_engine_type;
  typedef typename candidates_engine_type::space_type space_type;
  typedef typename candidates_engine_type::field_type field_type;
  
  ofstream output;
  // Try to open the file.
  output.open( __settings.filename().c_str(), ios::out | ios::trunc );
  if( !output.is_open() ) {
        // The output file either has not been supplied, or is invalid,
        // so use the default value which is "/dev/null"
    output.open( "/dev/null", ios::out | ios::trunc );
    cout << "\n\n### Outputting to /dev/null ###\n\n" << endl;
  }
  // If we are unable to open an output file then abort.
  if( !output.is_open() ) {
    cout << "Unable to open output file, aborting..." << endl;
    return;
  }
  
  output << "Siegel Set Cusp List Generator v1.0 (written by Brian Tyler)" << endl;
  output << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n" << endl;
  
  // Output the setup information
  cout << "Dimension: " << __settings.dimension() << endl;
  cout << "Generator: " << __settings.generator() << endl;
  cout << "Height: " << __settings.height() << endl;
  
  // Log the setup information
  output << "Dimension: " << __settings.dimension() << endl;
  output << "Generator: " << __settings.generator() << endl;
  output << "Height: " << __settings.height() << endl;
  
  candidates_engine_type engine;
  engine.set_output( output );
  engine.set_field_generator( __settings.generator() );
  make_siegel<field_type>( engine.space() );
  engine.space().height_ref().lower() = __settings.height();
  
  cout << "Space: " << engine.space() << endl << endl;
  output << "Space: " << engine.space() << endl << endl;
  
  engine( __settings.store() && __settings.tex() );
  
  if( __settings.tex() ) {
    engine.tex();
  }
  
  output.close();
}


} // namespace sg


int main( int argc, char* argv[] ) {
  using std::cout;
  cout << "siegelcl  Copyright (C) 2009  Brian Tyler\n"
          "This program comes with ABSOLUTELY NO WARRANTY.\n"
          "This is free software, and you are welcome to redistribute it\n"
          "under certain conditions.\n"
          "For more information see <http://www.gnu.org/licenses/>"
       << std::endl;
  
  sg::settings settings(argc, argv);
  
  // Display help information
  if( settings.help() ) {
    sg::help_information();
    return 0;
  }
  
  try {
    sg::validate_generator( settings. generator() );
    sg::validate_dimension( settings.dimension() );
    sg::validate_height( settings.height() );
  }
  catch (std::exception& e){
    std::cout << e.what() << std::endl;
    return 0;
  }
  
  
#define _SGSWITCH(N) case N: sg::generate<N>( settings ); break;
  switch( settings.dimension() ) {
    _SGSWITCH(2)
    _SGSWITCH(3)
    _SGSWITCH(4)
    _SGSWITCH(5)
    _SGSWITCH(6)
    default:
      break;
  }
#undef _SGSWITCH
  
  return 0;
}

#endif
