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
 *  @file     congruence_system.hpp
 *  @brief    An inline header file for the \c congruence_system class.
 *  @note     Include to solve systems linear congruences:
 *            \f$ a_i X = b_i (m_i) \f$
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-14
 */


#ifndef _SG_CONGRUENCE_SYSTEM_H
#define _SG_CONGRUENCE_SYSTEM_H 1

// Global includes
#include <cstddef>
#include <vector>

#include <boost/call_traits.hpp>

// Local includes
#include "utility/math/congruence_equation.hpp"


namespace sg
{
namespace utility
{
namespace math
{
/**
 *  @class    congruence_system
 *  @brief    Class representing a system of linear congruences
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-15
 *  @param    _Integer An integral type; defaults to \c long
 */
template <class _Integer = long> class congruence_system
{
  public:
    //! The object type
    typedef congruence_system<_Integer> self_type;
    //! The integer type
    typedef _Integer integer_type;
    //! The equation type
    typedef congruence_equation<integer_type> equation_type;
    
  private:
    //! The type of the structure which describes the system
    typedef std::vector<equation_type> system_type;
    //! Parameter type
    typedef typename boost::call_traits<integer_type>::param_type param_type;
    
  public:
    //! Iterator
    typedef typename system_type::iterator iterator;
    //! Constant iterator
    typedef typename system_type::const_iterator const_iterator;
    //! Reverse iterator
    typedef typename system_type::reverse_iterator reverse_iterator;
    //! Constant reverse iterator
    typedef typename system_type::const_reverse_iterator const_reverse_iterator;
    
    /**
     *  @brief  Add a new equation to the system
     *  @param  equation The equation to add to the system
     */
    void add_equation( equation_type __equation ) {
      system_.push_back( __equation );
    }
    
    /**
     *  @brief  Add a new equation to the system
     *  @param  c1 The coefficient of X
     *  @param  c0 The constant coefficient
     *  @param  m The modulus
     */
    void add_equation( param_type __c1, param_type __c0, param_type __m ) {
      system_.push_back( equation_type( __c1, __c0, __m ) );
    }
    
    
    /**
     *  @brief  Remove the last equation that was added to the system.
     */
    void remove_equation( void ) { system_.pop_back(); }
    
    
    /**
     *  @brief  Clear the system and remove all equations
     */
    void clear( void ) { system_.clear(); }
    
    /**
     *  @brief  Determine if the system is empty
     *  @return true - if the system contains no equations.
     */
    bool empty( void ) const { return bool( system_.empty() ); }
    
    // Iterator access
    //! Iterator pointing to the beginning of the system.
    iterator begin( void ) { return system_.begin(); }
    //! Iterator pointing one past the end of the system.
    iterator end( void ) { return system_.end(); }
    //! Const iterator pointing to the beginning of the system.
    const_iterator begin( void ) const { return system_.begin(); }
    //! Const iterator pointing one past the end of the system.
    const_iterator end( void ) const { return system_.end(); }
    //! Reverse iterator pointing to the r-beginning of the system.
    reverse_iterator rbegin( void ) { return system_.rbegin(); }
    //! Reverse iterator pointing one past the r-end of the system.
    reverse_iterator rend( void ) { return system_.rend(); }
    //! Const reverse iterator pointing to the r-beginning of the system.
    const_reverse_iterator rbegin( void ) const { return system_.rbegin(); }
    //! Const reverse iterator pointing one past the r-end of the system.
    const_reverse_iterator rend( void ) const { return system_.rend(); }
    
    /**
     *  @brief  Defines the output stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::ostream& operator<<( std::ostream& __os, const self_type& __t ){
      if( __t.empty() ) {
        __os << "[]";
      }
      else {
        __os << '[';
        const_iterator it = __t.begin();
        for( ; it != __t.end() - 1; ++it ) { __os << *it << ","; }
        __os << *it;
        __os << ']';
      }
      return __os;
    }
    
  private:
    //! The system data
    system_type system_;
};


/**
 *  @class    congruence_system_solver
 *  @brief    Functor for solving a system of linear congruences
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-15
 *  @param    _Integer An integral type; defaults to \c long
 */
template <class _Integer = long> class congruence_system_solver
{
  public:
    //! The object type
    typedef congruence_system_solver<_Integer> self_type;
    //! The integer type
    typedef _Integer integer_type;
    //! The type of the structure which describes the system
    typedef congruence_system<integer_type> system_type;
    //! The type of equations in the system
    typedef typename system_type::equation_type equation_type;
    
  private:
    //! Congruence solver type
    typedef congruence_solver<integer_type> solver_type;
    
  public:
    //! The solution type
    typedef typename solver_type::solution_type solution_type;
    //! The result type
    typedef solution_type result_type;
    
    /**
     *  @brief  Solve the system
     *  @param  system The system of equations to solve
     *  @return A \c congruence_solution {x0,xN} where X = x0 + n*xN for any
     *          integer n.
     *  @note   If \f$  xN = 0 \f$ then there is no solution.
     */
    result_type operator( ) ( const system_type& __system ) const {
      typedef typename system_type::const_iterator _iter;
      
      result_type result;
      if( __system.empty() ){ result.xN = 0; return result; }
      
      _iter it = __system.begin();
      result = solver_type()( *it );
      ++it;
      for(; it != __system.end(); ++it ) {
        
        // Fail there is no solution, so no point in continuing
        if( result.xN == 0 ) return result;
        
        // A solution exists
        
        // Note that if one is willing to invalidate the system in the
        // process of solving it then it can be solved slightly more
        // efficiently by writing:
        // it->c0 -= (it->c1) * result.x0;
        // it->c1 *= result.xN;
        
        // Substitution step:
        integer_type c0(it->c0 - (it->c1) * result.x0);
        integer_type c1(it->c1 * result.xN);
        
        // Solve and update the system
        result_type current = solver_type()( c1, c0, it->m );
        result.x0 += result.xN * current.x0;
        result.xN *= current.xN;
      }
      return result;
    }
    
    
    /**
     *  @brief  Validate a solution to a system
     *  @param  system The system of equations to solve
     *  @param  solution The solution to validate
     *  @return true if the solution satisfies the system
     *  @note   This is primarily for testing.
     */
    bool validate_solution
        ( const system_type& __system, const solution_type& __solution ) const
    {
      bool b = true;
      
      if( __solution.xN == 0 ) {
        return b;
      }
      
      typename system_type::const_iterator it = __system.begin();
      for( ; it != __system.end(); ++it ) {
        b &= (    (it->c1*__solution.x0 - it->c0) % it->m == 0
               &&
                  (it->c1*( __solution.x0 + __solution.xN ) - it->c0) % it->m
                  == 0
             );
      }
      
      return b;
    }
};

}// math
}// utility
}// sg
#endif
