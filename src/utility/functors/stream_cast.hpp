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
 *  @file     stream_cast.hpp
 *  @brief    This is a header implementation file for stream_cast.
 *  @note     Include this file to cast between types.
 *  @author   Brian Tyler
 *  @version  1.1
 *  @note     Made \c stream_cast into a function object
 *  @date     2008-03-03
 */


#ifndef _SG_STREAM_CAST_H
#define _SG_STREAM_CAST_H 1

// Local includes
#include <sstream>
#include <iostream>
#include <functional>

#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/call_traits.hpp>

// Global includes
#include "utility/precision.hpp"


namespace sg
{
namespace utility
{
namespace functors
{
/**
 *  @class    stream_cast
 *  @brief    A function object for casting between types
 *  @param    _From The type to convert from
 *  @param    _To The type to convert to
 *  @note     Attempts to use native casting, if no native cast exists, a
 *            stream cast is used.
 *  @author   Brian Tyler
 *  @version  1.1
 *  @date     2008-03-03
 */
template <class _From, class _To> struct stream_cast
  :public std::binary_function< typename boost::call_traits<_From>::param_type,
                                _To&, void >
{
  public:
    //! Object type
    typedef stream_cast<_From,_To> self_type;
    //! The type to convert from
    typedef _From from_type;
    //! The type to convert to
    typedef _To to_type;
    
    
  private:
    typedef typename boost::call_traits<_From>::param_type param_type;
    
    struct s_type{
      void operator() ( param_type __from, _To& __to ) const { __to = __from; }
    };
    struct d_type {
      void operator() ( param_type __from, _To& __to ) const {
        std::stringstream ss;
        ss.precision( precision::stream() );
        ss.flags( std::ios::fixed );
        ss << __from;
        ss >> __to;
      }
    };
    
    
  public:
    /**
     * @brief Generic casting operator.
     * @param from The variable to convert from
     * @param to The variable to convert to
     * @note  Conversion takes place throught a stringstream, so conversion
     *        will happen if both types have compatible string representations.
     *        This function is inefficient, so if conversion between types is
     *        speed critical then a native function should be prefered.
     */
    void operator( )( param_type __from, _To& __to ){
      using boost::is_same;
      using boost::mpl::if_;
      
      /* The point of the stream cast functor is that it should be able to
       * cast between types. However, if the types are the same it is a really
       * bad idea to push the conversion through a stream as that is incredibly
       * inefficient. This bit of template magic ensures that stream conversion
       * is only carried out when the types don't match, and this should be
       * optimised at compile time.
       */
      typedef typename if_<is_same<from_type,to_type>,s_type,d_type >::type op;
      op()( __from, __to );
    }
};


#define __SG_BEGIN_DEF_( FROM, TO )                                           \
    template <>                                                               \
        struct stream_cast <FROM,TO>                                          \
  :public std::binary_function<boost::call_traits<FROM>::param_type,          \
                               TO&, void >                                    \
  {                                                                           \
    typedef stream_cast<FROM,TO> self_type;                                   \
    typedef FROM from_type;                                                   \
    typedef TO to_type;                                                       \
    private:                                                                  \
      typedef boost::call_traits<FROM>::param_type param_type;                \
                                                                              \
    public:                                                                   \
      void operator( )( param_type __from, to_type& __to ) const {

#define __SG_END_DEF_ }};

/**
 * @brief Type to \c std::ostream casting
 */
template <class _From>
    struct stream_cast<_From,std::ostream>
  :public std::binary_function<typename boost::call_traits<_From>::param_type,
                               std::ostream&, void >
{
  //! Object type
  typedef stream_cast<_From,std::ostream> self_type;
  //! The type to convert from
  typedef _From from_type;
  //! The type to convert to
  typedef std::ostream to_type;
  
  private:
    typedef typename boost::call_traits<_From>::param_type param_type;
    
  public:
    void operator( )( param_type __from, to_type& __to ) const {
      __to.precision( precision::stream() );
      __to.flags( std::ios::fixed );
      __to << __from;
    }
};


/**
 * @brief \c std::istream to type casting
 */
template <class _To>
    struct stream_cast<std::istream, _To>
  :public std::binary_function< std::istream&, _To&, void >
{
  //! Object type
  typedef stream_cast<std::istream,_To> self_type;
  //! The type to convert from
  typedef std::istream from_type;
  //! The type to convert to
  typedef _To to_type;
  
  void operator( )( std::istream& __from, to_type& __to ) const {
    __from.precision( precision::stream() );
    __from.flags( std::ios::fixed );
    __from >> __to;
  }
};


// Inbuilt
/**
 * @brief \c double to \c long casting
 */
__SG_BEGIN_DEF_( double, long )
  __to = static_cast<long>( __from );
__SG_END_DEF_

/**
 * @brief \c long to \c double casting
 */
__SG_BEGIN_DEF_( long, double )
  __to = static_cast<double>( __from );
__SG_END_DEF_


#ifdef __GMP_PLUSPLUS__
/**
 * @brief Type to \c mpq_class casting
 */
template <class _From>
    struct stream_cast<_From,mpq_class>
  :public std::binary_function<typename boost::call_traits<_From>::param_type,
                               mpq_class&, void >
{
  //! Object type
  typedef stream_cast<_From,mpq_class> self_type;
  //! The type to convert from
  typedef _From from_type;
  //! The type to convert to
  typedef mpq_class to_type;
  
  private:
    typedef typename boost::call_traits<_From>::param_type param_type;
    
  public:
    void operator( )( from_type __from, to_type& __to ) const {
      mpf_class f;
      stream_cast<from_type,mpf_class>( )( __from, f );
      __to = f;
    }
};

/**
 * @brief \c mpt_class to type casting
 */
template <class _To>
    struct stream_cast<mpq_class, _To>
  :public std::binary_function<const mpq_class&, _To&, void >
{
  //! Object type
  typedef stream_cast<mpq_class,_To> self_type;
  //! The type to convert from
  typedef mpq_class from_type;
  //! The type to convert to
  typedef _To to_type;
  
  void operator( )( const from_type& __from, to_type& __to ) const {
    mpf_class f( __from );
    stream_cast<mpf_class,to_type>( )( f, __to );
  }
};

// Note to avoid the compiler getting confused we need to supply same type
// casting.
__SG_BEGIN_DEF_( mpq_class, mpq_class )
  __to = __from;
__SG_END_DEF_


// Double
/**
 * @brief \c mpz_class to \c double casting
 */
__SG_BEGIN_DEF_( mpz_class, double )
  __to = mpz_get_d( __from.get_mpz_t() );
__SG_END_DEF_

/**
 * @brief \c double to \c mpz_class casting

 */
__SG_BEGIN_DEF_( double, mpz_class )
  mpz_set_d( __to.get_mpz_t(), __from );
__SG_END_DEF_

/**
 * @brief \c mpf_class to \c double casting
 */
__SG_BEGIN_DEF_( mpf_class, double )
  __to = mpf_get_d( __from.get_mpf_t() );
__SG_END_DEF_

/**
 * @brief \c double to \c mpf_class casting
 */
__SG_BEGIN_DEF_( double, mpf_class )
  mpf_set_d( __to.get_mpf_t(), __from );
__SG_END_DEF_

/**
 * @brief \c mpq_class to \c double casting
 */
__SG_BEGIN_DEF_( mpq_class, double )
  __to = mpq_get_d( __from.get_mpq_t() );
__SG_END_DEF_

/**
 * @brief \c double to \c mpq_class casting
 */
__SG_BEGIN_DEF_( double, mpq_class )
  mpq_set_d( __to.get_mpq_t(), __from );
__SG_END_DEF_


// Long
/**
 * @brief \c long to \c mpz_class casting
 */
__SG_BEGIN_DEF_( long, mpz_class )
    mpz_set_si( __to.get_mpz_t(), __from );
__SG_END_DEF_

/**
 * @brief \c mpz_class to \c long casting
 */
__SG_BEGIN_DEF_( mpz_class, long )
  __to = mpz_get_si( __from.get_mpz_t() );
__SG_END_DEF_

/**
 * @brief \c long to \c mpf_class casting
 */
__SG_BEGIN_DEF_( long, mpf_class )
    mpf_set_si( __to.get_mpf_t(), __from );
__SG_END_DEF_

/**
 * @brief \c mpf_class to \c long casting
 */
__SG_BEGIN_DEF_( mpf_class, long )
  __to = mpf_get_si( __from.get_mpf_t() );
__SG_END_DEF_

/**
 * @brief \c long to \c mpq_class casting
 */
__SG_BEGIN_DEF_( long, mpq_class )
    mpq_set_si( __to.get_mpq_t(), __from, 1 );
__SG_END_DEF_

/**
 * @brief \c mpq_class to \c long casting
 */
__SG_BEGIN_DEF_( mpq_class, long )
  __to = static_cast<long>( mpq_get_d( __from.get_mpq_t() ) );
__SG_END_DEF_


// MPZ
/**
 * @brief \c mpf_class to \c mpz_class casting
 */
__SG_BEGIN_DEF_( mpf_class, mpz_class )
    mpz_set_f( __to.get_mpz_t(), __from.get_mpf_t() );
__SG_END_DEF_

/**
 * @brief \c mpz_class to \c mpf_class casting
 */
__SG_BEGIN_DEF_( mpz_class, mpf_class )
    mpf_set_z( __to.get_mpf_t(), __from.get_mpz_t() );
__SG_END_DEF_


/**
 * @brief \c mpq_class to \c mpz_class casting
 */
__SG_BEGIN_DEF_( mpq_class, mpz_class )
    mpz_set_q( __to.get_mpz_t(), __from.get_mpq_t() );
__SG_END_DEF_

/**
 * @brief \c mpz_class to \c mpf_class casting
 */
__SG_BEGIN_DEF_( mpz_class, mpq_class )
    mpq_set_z( __to.get_mpq_t(), __from.get_mpz_t() );
__SG_END_DEF_


// MPF
/**
 * @brief \c mpf_class to \c mpq_class casting
 */
__SG_BEGIN_DEF_( mpf_class, mpq_class )
    mpq_set_f( __to.get_mpq_t(), __from.get_mpf_t() );
__SG_END_DEF_

/**
 * @brief \c mpq_class to \c mpf_class casting
 */
__SG_BEGIN_DEF_( mpq_class, mpf_class )
    mpf_set_q( __to.get_mpf_t(), __from.get_mpq_t() );
__SG_END_DEF_
#endif // __GMP_PLUSPLUS__

#undef __SG_BEGIN_DEF_
#undef __SG_END_DEF_

}// functors
}// utility
}// sg
#endif
