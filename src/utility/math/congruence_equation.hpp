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
 *  @file     congruence_equation.hpp
 *  @brief    An inline header file for solving linear congruence equations
 *  @note     Include to solve linear congruences
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-14
 */


#ifndef _SG_LINEAR_CONGRUENCE_H
#define _SG_LINEAR_CONGRUENCE_H 1

// Global includes
#include <boost/array.hpp>
#include <boost/call_traits.hpp>

// Local includes
#include "utility/math/gcd.hpp"


namespace sg
{
namespace utility
{
namespace math
{
/**
 *  @class    congruence_equation
 *  @brief    Class representing a linear congruence equation
 *            \f$ c_1 X \equiv c_0 \mod m \f$
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-14
 *  @param    _Integer An integral type; defaults to \c long
 */
template <class _Integer = long> struct congruence_equation {
  //! Object type
  typedef congruence_equation<_Integer> self_type;
  //! Integer type
  typedef _Integer integer_type;
  
  private:
    //! Parameter type
    typedef typename boost::call_traits<integer_type>::param_type param_type;
  
  public:
    /**
     *  @brief  Default constructor
     */
    congruence_equation( )
  : c1(1), c0(1), m(1) { }
    
    /**
     *  @brief  Equation constructor
     *  @param  c1 The coefficient of X
     *  @param  c0 The constant coefficient
     *  @param  m The modulus
     */
    congruence_equation( param_type __c1, param_type __c0, param_type __m )
  : c1(__c1), c0(__c0), m(__m) { }
    
    /**
     *  @brief  General copy constructor
     *  @param  that A congruence equation.
     *  @note   This is mainly for testing.
     */
    template<class _OtherInt>
        congruence_equation ( const congruence_equation<_OtherInt>& __that )
  : c1(__that.c1), c0(__that.c0), m(__that.m) { }
    
    
    //! The coefficient of X
    integer_type c1;
    //! The constant coefficient
    integer_type c0;
    //! The modulus
    integer_type m;
    
    /**
     *  @brief  Defines the output stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::ostream& operator<<( std::ostream& __os, const self_type& __t ){
      __os << '[' << __t.c0 << ',' << __t.c1 << ',' << __t.m << ']';
      return __os;
    }
};


/**
 *  @class    congruence_solution
 *  @brief    Class representing the solution to a linear congruence equation
 *            \f$ X = n*x_N + x_0 \f$
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-14
 *  @param    _Integer An integral type; defaults to \c long
 */
template <class _Integer = long> struct congruence_solution {
  //! Object type
  typedef congruence_solution<_Integer> self_type;
  //! Integer type
  typedef _Integer integer_type;
  
  //! A base solution
  integer_type x0;
  //! The skip between solutions
  integer_type xN;
  
  /**
   *  @brief  Defines the output stream operator
   *  @param  __os An ostream.
   *  @param  __t The object to stream out.
   */
  friend std::ostream& operator<<( std::ostream& __os, const self_type& __t ){
    __os << '[' << __t.x0 << ',' << __t.xN << ']';
    return __os;
  }
};

/**
 *  @class    congruence_solver
 *  @brief    Functor for solving linear congruences
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-14
 *  @param    _Integer An integral type; defaults to \c long
 */
template <class _Integer = long> class congruence_solver
{
  public:
    //! The object type
    typedef congruence_solver<_Integer> self_type;
    //! The integer type
    typedef _Integer integer_type;
    //! The solution type
    typedef congruence_equation<integer_type> equation_type;
    //! The solution type
    typedef congruence_solution<integer_type> solution_type;
    //! The result type
    typedef solution_type result_type;
    
  private:
    //! Parameter type
    typedef typename boost::call_traits<integer_type>::param_type param_type;
    //! The GCD functor;
    typedef gcd<integer_type> gcd_func;
    
  public:
    /**
     *  @brief  Solve the linear congruence \f$ c_1 X \equiv c_0 \mod m \f$
     *  @param  c1 The coefficient of X
     *  @param  c0 The constant coefficient
     *  @param  m The modulus
     *  @return A congruence solution \f$ \{x_0,x_N\}\f$ such that
     *          \f$ X =  nx_1 + x_0 \f$ gives all possible solutions.
                If no solution exists then \f$ \{x,0\}\f$ is returned where
                \f$ x \f$ is arbitrary.
     *  @note   If a solution exists, then the base solution returned is always
     *          the smallest possible positive solution.
     */
    const result_type operator( )
        ( param_type __c1, param_type __c0, param_type __m ) const
    {
      // Avoids compilation warning and keeps congruence_solution as pod
      result_type result = {0,0};
      typename gcd_func::result_type gcd( gcd_func()( __c1, __m ) );
      
      if( __c0 % gcd.gcd != 0 ) {
        // No solution exists
        result.xN = 0;
        return result;
      }
      
      // A solution exists
      result.xN = __m / gcd.gcd;
      result.x0 = (gcd.a * __c0 / gcd.gcd) % result.xN;
      
      // Not necessary, but easier if the base solution is always positive
      if( result.x0 < 0 ) result.x0 += result.xN;
      
      return result;
    }
    
    
    /**
     *  @brief  Solve the linear congruence \f$ c_1 X \equiv c_0 \mod m \f$
     *  @param  equation The congruence equation to solve.
     *  @return A congruence solution \f$ \{x_0,x_N\}\f$ such that
     *          \f$ X =  nx_1 + x_0 \f$ gives all possible solutions.
     *  @note   Because of its custom constructor \c equation_type is not POD.
     */
    result_type operator( )( const equation_type& __equation ) {
      return result_type( operator()( __equation.c1, __equation.c0,
                                      __equation.m )
                        );
    }
};

}// math
}// utility
}// sg
#endif
