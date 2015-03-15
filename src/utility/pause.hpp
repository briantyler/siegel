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
 *  @file     pause.hpp
 *  @brief    An inline header file for the \c pause class.
 *  @note     Pauses the program and intitiates a countdown timer to restart.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-01-01
 */


#ifndef _SG_PAUSE_H
#define _SG_PAUSE_H 1

// Global includes
#include <iostream>
#include <boost/timer.hpp>


namespace sg
{
namespace utility
{
/**
 *  @struct  pause
 *  @brief   Pauses the program and counts down to restart
 *  @author  Brian Tyler
 *  @version 1.0
 *  @date    2008-06-20
 */
struct pause {
  /**
   *  @brief  Pauses the program and counts down to restart.
   *  @param  countdown  The number of seconds to count down for.
   */
  void operator( ) ( unsigned int __countdown ) const {
    using std::cout;
    using std::flush;
    using std::endl;
    
    double total = __countdown, update = 0.0;
    boost::timer t;
    
    cout << endl << "Starting in:";
    // Loop until the timer is beyond the maximum time
    while( t.elapsed() < total ) {
      // After each second increase the update and print the next countdown
      // number to the screen.
      if( t.elapsed() > update ) {
        cout << " " << __countdown << flush ;
        update += 1.0;
        --__countdown;
      }
    }
    cout << endl;
  }
};

}// utility
}// sg
#endif
