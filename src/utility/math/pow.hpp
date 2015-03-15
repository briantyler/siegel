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
 *  @file     pow.hpp
 *  @brief    This is a header implementation file for pow
 *  @note     Include this file to compute (floating point) powers.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-25
 */


#ifndef _SG_POW_H
#define _SG_POW_H 1

// Global includes
#include <cstddef>
#include <cmath>
#include <functional>

#include <boost/call_traits.hpp>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_same.hpp>

// Local includes
#include "utility/functors/stream_cast.hpp"


namespace sg
{
namespace utility
{
namespace math
{
/**
 *  @struct   powf
 *  @brief    Floating point powering functor.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-25
 *  @param    Type The type of the object to power
 */
template <class _Type> class powf
  :public std::binary_function< typename boost::call_traits<_Type>::param_type,
                                typename boost::call_traits<_Type>::param_type,
                                _Type>
{
  public:
    //! The object type
    typedef powf<_Type> self_type;
    //! The value type
    typedef _Type value_type;
    //! The parameter type
    typedef typename boost::call_traits<_Type>::param_type param_type;
    
    
    /**
     *  @brief Functor operator.
     *  @param   value The value to power
     *  @param   exp The exponent
     *  @return  \f$ value^{exp} \f$
     */
    value_type operator()( param_type __value, param_type __exp ) const {
      using boost::is_same;
      using boost::mpl::if_;
      typedef typename if_< is_same<value_type,double>, d_op, x_op >::type op;
      
      return value_type( op()( __value, __exp ) );
    }
    
    
  private:
    /**
     *  @struct  d_op
     *  @brief   Helper functor for double specialisation.
     */
    struct d_op
    {
      double operator()(const double __value, const double __exp) const {
        return double( ::pow( __value, __exp ) );
      }
    };
    
    /**
     *  @struct  x_op
     *  @brief   Helper functor for non double types.
     */
    struct x_op
    {
      value_type operator()(param_type __value, param_type __exp) const {
        typedef functors::stream_cast<value_type,double> _scvd;
        typedef functors::stream_cast<double,value_type> _scdv;
        
        double dblvalue;
        _scvd( )( __value, dblvalue );
  
        double dblexp;
        _scvd( )( __exp, dblexp );
  
        value_type output;
        _scdv( )(::pow( dblvalue, dblexp ),output);
  
        return output;
      }
    };
};


/**
 *  @struct   powi
 *  @brief    Integer powering functor.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-25
 *  @param    _Type The type of the object to power
 */
template <class _Type> struct powi
  :public std::binary_function<_Type, long, _Type>
{
  public:
    //! The object type
    typedef powi<_Type> self_type;
    //! The value type
    typedef _Type value_type;
    
    /**
     *  @brief Functor operator.
     *  @param   value The value to power
     *  @param   exp The exponent
     *  @return  \f$ value^{exp} \f$
     */
    value_type operator()( value_type __value, long __exp ) const {
      // Explicitly compute the first few powers, this is never going to be as
      // efficient as simply putting x*x (for example) in the code, but it is
      // better than needing to go through initialization and loop constructs.
      //
  
      _Type z( 1 );
  
      if( __exp < 0 ) {
        z /= operator()( __value, -__exp );
        return z;
      }
  
      switch( __exp ){
        case 7:
          z *= __value;
        case 6:
          z *= __value;
        case 5:
          z *= __value;
        case 4:
          z *= __value;
        case 3:
          z *= __value;
        case 2:
          z *= __value;
        case 1:
          z *= __value;
        case 0:
          return z;
      }
  
      do {
        if( __exp % 2 != 0 )
        {
          z *= __value;
        }
        // Repeatedly square on each bit.
        __value *= __value;
        // Iterate through all the bits.
        __exp /= 2;
      
      }while( __exp != 0 );     // If exp == 0 then we are done
    
      return z;
    }
};


/**
 *  @class    powi_ct
 *  @brief    Compile time unsigned integer powering metafunction
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-18
 *  @param    Base The base to power
 *  @param    Exponent The exponent to raise the base to.
 *  @return   \c value
 */
template <std::size_t _Base, std::size_t _Exponent> struct powi_ct
{
  //! The object type
  typedef powi_ct<_Base,_Exponent> self_type;
  //! The input type
  typedef std::size_t value_type;
  
  //! The base to power
  static const value_type base = _Base;
  //! The exponent to raise the base to
  static const value_type exponent = _Exponent;
  //! The value of the \c Base to the power \c Exponent
  static const value_type value = _Base * powi_ct<_Base,_Exponent-1>::value;
};


template <std::size_t _Base> struct powi_ct<_Base,0>
{
  //! The object type
  typedef powi_ct<_Base,0> self_type;
  //! The input type
  typedef std::size_t value_type;
  
  //! The base to power
  static const value_type base = _Base;
  //! The exponent to raise the base to
  static const value_type exponent = 0;
  //! The base case
  static const value_type value = 1;
};

}// math
}// utility
}// sg
#endif
