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
 *  @file     stream_cast_tuple.hpp
 *  @brief    This is a header implementation file for stream_cast_tuple.
 *  @note     Include this file to cast between tuple types using streams.
 *  @author   Brian Tyler
 *  @version  1.1
 *  @note     Made \c stream_cast_tuple into a function object
 *  @date     2008-03-12
 */


#ifndef _SG_STREAM_CAST_TUPLE_H
#define _SG_STREAM_CAST_TUPLE_H 1

// Global includes
#include <sstream>
#include <iostream>
#include <functional>

#include <boost/mpl/assert.hpp>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/tuple/tuple_io.hpp>

// Local includes
#include "utility/precision.hpp"
#include "utility/functors/stream_cast.hpp"


namespace sg
{
namespace utility
{
namespace functors
{
/**
 *  @struct   x_type_op
 *  @brief    A helper function object for casting between tuple types
 *  @param    _From The tuple type to convert from.
 *  @param    _To The tuple type to convert to.
 *  @note     Uses the stream cast operator for fast conversion between types
 *  @author   Brian Tyler
 *  @version  1.1
 *  @date     2008-03-12
 */
template <class _From, class _To, int _N> struct x_type_op {
  void operator() ( const _From& __from, _To& __to ) const {
    using boost::tuples::element;
    typedef typename element<_N,_From>::type from_type;
    typedef typename element<_N,_To>::type to_type;
    typedef stream_cast<from_type,to_type> _sc;
    
    _sc()( __from.template get<_N>(), __to.template get<_N>() );
    x_type_op<_From,_To,_N-1>()( __from, __to );
  }
};
  
  
template <class _From, class _To> struct x_type_op<_From,_To,0> {
  void operator() ( const _From& __from, _To& __to ) const {
    using boost::tuples::element;
    typedef typename element<0,_From>::type from_type;
    typedef typename element<0,_To>::type to_type;
    typedef stream_cast<from_type,to_type> _sc;
    
    _sc()( __from.template get<0>(), __to.template get<0>() );
  }
};
  
  
/**
 *  @class    stream_cast_tuple
 *  @brief    A function object for casting between tuple types
 *  @param    _From The tuple type to convert from.
 *  @param    _To The tuple type to convert to.
 *  @note     Attempts to use native casting, if no native cast exists, a
 *            stream cast is used.
 *  @author   Brian Tyler
 *  @version  1.2
 *  @date     2008-03-12
 */
template <class _From, class _To> class stream_cast_tuple
  :public std::binary_function< const _From&, _To&, void >
{
    struct s_type{
      void operator() ( const _From& __from, _To& __to ) const {__to = __from;}
    };
    
  public:
    //! Object type
    typedef stream_cast_tuple<_From,_To> self_type;
    //! The type to convert from
    typedef _From from_type;
    //! The type to convert to
    typedef _To to_type;
    
    /**
     * @brief Generic casting operator.
     * @param from The variable to convert from
     * @param to The variable to convert to
     * @note  Conversion takes place throught a stringstream, so conversion
     *        will happen if both types have compatible string representations.
     *        This function is inefficient, so if conversion between types is
     *        speed critical then a native function should be prefered.
     */
    void operator( )( const from_type& __from, to_type& __to ){
      using boost::is_same;
      using boost::mpl::if_;
      using boost::tuples::length;
      
      
      //BOOST_MPL_ASSERT(( length<_From>::value == length<_To>::value ));
      
      /* The point of the stream cast functor is that it should be able to
       * cast between types. However, if the types are the same it is a really
       * bad idea to push the conversion through a stream as that is incredibly
       * inefficient. This bit of template magic ensures that stream conversion
       * is only carried out when the types don't match, and this should be
       * optimised at compile time.
       */
      typedef x_type_op<_From,_To, length<_From>::value - 1> x_op;
      typedef typename if_< is_same<_From,_To>, s_type, x_op >::type op;
      op()( __from, __to );
    }
};


/**
 *  @brief Specialisation for ostreams
 */
template <class _From>
    class stream_cast_tuple<_From,std::ostream>
  :public std::binary_function< _From&, std::ostream&, void >
{
  public:
    //! Object type
    typedef stream_cast_tuple<_From,std::ostream> self_type;
    //! The type to convert from
    typedef _From from_type;
    //! The type to convert to
    typedef std::ostream to_type;
    
    void operator( )( const from_type& __from, to_type& __to ) const {
      __to.precision( precision::stream() );
      __to.flags( std::ios::fixed );
      
      __to << boost::tuples::set_open('[')
           << boost::tuples::set_close(']')
           << boost::tuples::set_delimiter(',')
           << __from;
    }
};

/** @brief Specialisation for istreams
 *  @bug   The istream operator (>>) for the \c gmp::mpq_class is not
 *         compatible with the ostream operator (<<). The upshot of this is
 *         that if a mpq_class object is streamed out and then in it is a
 *         completely different value. This is a GMP bug, but something to be
 *         aware of here. This could be fixed by rewriting the istream functor
 *         in a similar way as the general functor.
 */
template <class _To>
    class stream_cast_tuple<std::istream, _To>
  :public std::binary_function< std::istream&, _To&, void >
{
  public:
    //! Object type
    typedef stream_cast_tuple<std::istream,_To> self_type;
    //! The type to convert from
    typedef std::istream from_type;
    //! The type to convert to
    typedef _To to_type;
    
    void operator( )( from_type& __from, to_type& __to ) const {
      __from.precision( precision::stream() );
      __from.flags( std::ios::fixed );
      
      __from >> boost::tuples::set_open('[')
             >> boost::tuples::set_close(']')
             >> boost::tuples::set_delimiter(',')
             >> __to;
    }
};


}// functors
}// utility
}// sg

#endif
