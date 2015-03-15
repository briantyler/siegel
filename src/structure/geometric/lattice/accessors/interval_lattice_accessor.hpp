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
 *  @file     interval_lattice_accessor.hpp
 *  @brief    An inline header file for the \c interval_lattice_accessor class.
 *  @note     The \c interval_lattice_accessor class creates a common interface
 *            to interval lattices for hyperbolic structures.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-08
 */


#ifndef _SG_INTERVAL_LATTICE_ACCESSOR_H
#define _SG_INTERVAL_LATTICE_ACCESSOR_H 1

// Global includes
#include <boost/operators.hpp>

// Local includes
#include "structure/geometric/lattice/interval_lattice.hpp"


namespace sg
{
namespace structure
{
namespace geometric
{
namespace lattice
{
namespace accessor
{
/**
 *  @class    interval_lattice_accessor
 *  @brief    An accessor class for real numbers
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-04-08
 *  @param    _Float A floating point type; defaults to \c double
 *  @param    _Integer The integer type which iterates over the interval
 */
template <class _Float = double, class _Integer = long>
    class interval_lattice_accessor
  : private boost::equality_comparable<
                            interval_lattice_accessor<_Float,_Integer>
              , boost::equality_comparable<
                                interval_lattice_accessor<_Float,_Integer>,
                                interval_lattice<_Float,_Integer>
                >
            >
{
  public:
    //! Object's type
    typedef interval_lattice_accessor<_Float,_Integer> self_type;
    //! The floating point type
    typedef _Float float_type;
    //! The integer type
    typedef _Integer integer_type;
    //! The real type
    typedef interval_lattice<float_type,integer_type> value_type;
    //! The common return type of with zeta
    typedef value_type return_type;
    
    
  public:
    /**
     *  @brief  Empty constructor
     */
    interval_lattice_accessor( ) : value_() { }
    
    
    /**
     *  @brief  Value constructor
     *  @param  value The lattice to construct from
     */
    interval_lattice_accessor( const value_type& __value ) : value_(__value) { }
    
    
    /**
     *  @brief   Value assignment
     *  @param   value The interval lattice to assign from
     *  @return  A reference to this.
     */
    self_type& operator=( const value_type& __value ) {
      value_ = __value;
      return *this;
    }
    
    
    /**
     *  @brief   A common interface between r and zeta
     *  @return  An object which is constructed from r and shares a common
     *           type with the zeta return type.
     */
    const return_type& common( void ) const { return value_; }
    
    
    /**
     *  @brief  Compares two interval_lattice_accessor objects for equality.
     *  @param  rhs A \c interval_lattice_accessor
     */
    bool operator==( const self_type& __rhs ) const {
      return bool( *this == __rhs.value_ );
    }
    
    
    /**
     *  @brief  Compares a interval_lattice_accessor to an interval_lattice for
     *          equality.
     *  @param  rhs A \c value_type
     */
    bool operator==( const value_type& __rhs ) const {
      return bool( value_ == __rhs );
    }
    
    
    /**
     *  @brief  Defines the output stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::ostream& operator<<( std::ostream& __os, const self_type& __t ) {
      __os << __t.value_;
      return __os;
    }
    
    /**
     *  @brief  Defines the input stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::istream& operator>>( std::istream& __is, self_type& __t ) {
      __is >> __t.value_;
      return __is;
    }
    
    
  protected:
    //! An interval_lattice
    value_type value_;
};

}// accessor
}// lattice
}// geometric
}// structure
}// sg
#endif
