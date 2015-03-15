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
 *  @file     static_tree.hpp
 *  @brief    An inline header file for the \c static_tree class.
 *  @note     A static binary tree of arrays of a given type.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-05-19
 */


#ifndef _SG_STATIC_TREE_H
#define _SG_STATIC_TREE_H 1

// Global includes
#include <cstddef>

#include <boost/array.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/comparison.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/assert.hpp>

// Local includes
#include "utility/math/pow.hpp"


namespace sg
{
namespace structure
{
namespace geometric
{
namespace hyperbolic
{
namespace detail
{
/**
 *  @class   static_tree
 *  @brief   A binary tree of fixed sized arrays.
 *  @param   Type The type of object to store in the tree.
 *  @param   Depth The depth of the tree.
 *  @author  Brian Tyler
 *  @version 1.0
 *  @date    2008-05-19
 *  @note    The structure may be best thought of as a pyramid since there is
 *           no direct connection between levels of the tree. And despite the
 *           amount of code, most of it is mpl stuff, so all of the logic is
 *           moved to compile-time, rather than runtime.
 */
template <class _Type, std::size_t _Depth> class static_tree
{
  public:
    //! The type of the object
    typedef static_tree<_Type,_Depth> self_type;
    //! The value type.
    typedef _Type value_type;
    //! The size type
    typedef std::size_t size_type;
    //! The size of an array at this level.
    static const size_type depth = _Depth;
    //! The size of an array at this level.
    static const size_type branches = 2;
    //! The size of an array at this level.
    static const size_type static_size
        = utility::math::powi_ct<branches,_Depth>::value;
    //! The type of the fixed size boost array at this level.
    typedef boost::array<value_type,static_size> array_type;
    
  private:
    typedef boost::mpl::int_<_Depth> _depth;
    typedef boost::mpl::int_<0> _zero;
    
    BOOST_MPL_ASSERT(( boost::mpl::less_equal< _zero, _depth > ));
    
    typedef boost::mpl::equal_to<_depth, _zero> _eq;
    typedef boost::mpl::if_< _eq, void*, static_tree<_Type,_Depth - 1> > _if1;
    
  public:
    //! Type of a child object.
    typedef typename _if1::type child_type;
    
    
    /**
     *  @brief  Get a constant reference root object.
     *  @return A constant reference to the root object.
     */
    const value_type& root( void ) const {return at<0>(0);}
    
    /**
     *  @brief  Get a reference root object.
     *  @return A reference to the root object.
     */
     value_type& root( void ) {return at<0>(0);}
    
    
    /**
     *  @brief  Get a constant reference to the top level array.
     *  @return A constant reference to the top level array.
     */
    const array_type& get( void ) const { return array_; }
    
    /**
     *  @brief  Get a reference to the top level array.
     *  @return A reference to the top level array.
     */
    array_type& get( void ) { return array_; }
    
    
    /**
     *  @brief  Get a constant reference top level object at \c loc.
     *  @return A constant reference to the top level object at \c loc
     */
    const value_type& at( size_type __loc ) const {return array_.at(__loc);}
    
    /**
     *  @brief  Get a reference top level object at \c loc.
     *  @return A reference to the top level object at \c loc
     */
    value_type& at( size_type __loc ) {return array_.at(__loc);}
    
    
    /**
     *  @brief  Get the size of the array at level \c Level
     *  @return The size of the array at level \c Level \$ = 2^{Level} \$
     */
    template <size_type _Level> static const size_type size( void ) {
      using namespace boost::mpl;
      
      // Define mpl integeger types from the template arguments
      typedef int_<_Level> level_int;
      
      // Require 0 <= _Level <= _Depth
      BOOST_MPL_ASSERT(( less_equal< level_int, _depth > ));
      BOOST_MPL_ASSERT(( less_equal< _zero, level_int > ));
      
      return static_tree<void*,_Level>::static_size;
    }
    
    
    /**
     *  @brief  Get a constant reference to the array at level \c Level
     *  @param  Level The level of the array to return
     *  @return A constant reference to the array at level \c Level
     *  @note   It would not be possible to do this using standard template
     *          partial specialization because member functions of template
     *          classes cannot be partially specialized.
     */
    template <size_type _Level>
        const typename static_tree<value_type,_Level>::array_type&
        get( void ) const
    {
      using namespace boost::mpl;
      
      // Define mpl integeger types from the template arguments
      typedef int_<_Level> level_int;
      
      // Require 0 <= _Level <= _Depth
      BOOST_MPL_ASSERT(( less_equal< level_int, _depth > ));
      BOOST_MPL_ASSERT(( less_equal< _zero, level_int > ));
      
      // If _Level == _Depth then return the array, otherwise call get on the
      // child. The functions are accessed through function objects and are
      // completely optimised out at compile time so this is about as efficient
      // a structure as it is possible to get.
      typedef equal_to< level_int, _depth > _eq;
      typedef typename if_< _eq, arrayfunc, childfunc<_Level> >::type _fn;
      
      return _fn()( *this );
    }
    
    
    /**
     *  @brief  Get a reference to the array at level \c Level
     *  @param  Level The level of the array to return
     *  @return A reference to the array at level \c Level
     */
    template <size_type _Level>
        typename static_tree<value_type,_Level>::array_type&
        get( void )
    {
      typedef typename static_tree<value_type,_Level>::array_type& _return_type;
      return const_cast<_return_type>(
               static_cast<const self_type&>( *this ).template get<_Level>()
             );
    }
    
    
    /**
     *  @brief  Get a constant reference to the object at \c loc at given Level
     *  @param  Level The level of the tree to get the object from
     *  @param  loc The location of the object in the array
     *  @return A constant reference to the object at \c <Level>.at(loc)
     */
    template <size_type _Level>
        const value_type& at( size_type __loc ) const
    { return get<_Level>().at( __loc ); }
    
    /**
     *  @brief  Get a reference to the object at \c loc at given Level
     *  @param  Level The level of the tree to get the object from
     *  @param  loc The location of the object in the array
     *  @return A reference to the object at \c <Level>.at(loc)
     */
    template <size_type _Level>
        value_type& at( size_type __loc )
    { return get<_Level>().at( __loc ); }
    
    
  private:
    /**
     *  @class   arrayfunc
     *  @brief   A function object helper that accesses the array member.
     *  @param   t The static_tree to get the array from.
     *  @return  The array of the input object.
     */
    struct arrayfunc {
      const array_type& operator()( const self_type& __t ) const
      { return __t.array_; }
    };
    
    /**
     *  @class   childfunc
     *  @brief   A function object helper to access the array member of child.
     *  @param   Level Determines the level of the array.
     *  @return  The array of the child (possibly recersively) at \c Level
     */
    template <size_type _Level>
    struct childfunc {
      const typename static_tree<value_type,_Level>::array_type&
          operator()( const self_type& __t ) const
      { return __t.child_.template get<_Level>() ; }
    };
    
    
    //! The array at this level in the tree
    array_type array_;
    //! The next layer down in the tree
    child_type child_;
};

}// detail
}// hyperbolic
}// geometric
}// structure
}// sg
#endif
