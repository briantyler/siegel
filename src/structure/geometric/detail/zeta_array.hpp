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
 *  @file     zeta_array.hpp
 *  @brief    An inline header file for the \c zeta_array class.
 *  @note     The \c zeta_array class is an array wrapper for zeta coordinates
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-23
 */


#ifndef _SG_ZETA_ARRAY_H
#define _SG_ZETA_ARRAY_H 1

// Global includes
#include <algorithm>

#include <boost/array.hpp>
#include <boost/operators.hpp>
#include <boost/iterator/indirect_iterator.hpp>

// Local includes
#include "structure/geometric/detail/hyperbolic_base.hpp"
#include "structure/geometric/detail/zeta_coordinate.hpp"
#include "structure/geometric/detail/iterators/zeta_ref_iterator.hpp"

#include "utility/io/container_to_string.hpp"
#include "utility/io/string_to_array.hpp"


namespace sg
{
namespace structure
{
namespace geometric
{
namespace detail
{
/**
 *  @class    zeta_array
 *  @brief    A common interface for all zeta_arrays
 *  @param    N The hyperbolic dimension of the zeta array
 *  @param    Accessor The accessor object type
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-02-23
 */
  template < std::size_t N, class _Accessor >
    class zeta_array
  : public hyperbolic_base<N>,
    private boost::equality_comparable< zeta_array<N,_Accessor > >
{
  public:
    //! Object's type
    typedef zeta_array<N,_Accessor> self_type;
    //! The zeta array type
    typedef self_type zeta_array_type;
    //! Dimension type
    typedef std::size_t size_type;
     //! Zeta coordinate type
    typedef zeta_coordinate<_Accessor> zeta_coordinate_type;
    //! The raw type of a coordinate of zeta
    typedef typename zeta_coordinate_type::zeta_type zeta_type;
    //! Type of the real component of zeta
    typedef typename zeta_coordinate_type::real_type real_type;
    //! Type of the imaginary component of zeta
    typedef typename zeta_coordinate_type::imag_type imag_type;
    //! The common return type; a common type that real and imag cast to.
    typedef typename zeta_coordinate_type::return_type return_type;
    
    
  private:
    //! Zeta array type
    typedef boost::array<zeta_coordinate_type,self_type::zeta_size> array_type;
    //! Pointer array type
    // Here is something really important to realise; This array is created
    // with the knowledge that the larges number of common elements is
    // hyperbolic_size. Having a static array is more efficient than having a
    // vector, so this is a fair assumption. If this were not the case a
    // different approach may want to be taken.
    
  protected:
    typedef boost::array<return_type *, self_type::hyperbolic_size>
            ptr_array_type;
    //! Pointer iterator type
    typedef return_type* const * ptr_iterator;
    //! Const Pointer iterator type
    typedef const return_type* const * const_ptr_iterator;
    
    
  public:
    //! Foward zeta iterator
    typedef typename array_type::iterator zeta_iterator;
    //! Constant forward zeta iterator
    typedef typename array_type::const_iterator const_zeta_iterator;
    //! Reverse zeta iterator
    typedef typename array_type::reverse_iterator reverse_zeta_iterator;
    //! Constant reverse zeta iterator
    typedef typename array_type::const_reverse_iterator
        const_reverse_zeta_iterator;
    
    //! Forward raw zeta iterator
    typedef iterators::zeta_ref_iterator<zeta_iterator> zeta_ref_iterator;
    //! Constant forward raw zeta iterator
    typedef iterators::zeta_ref_iterator<const_zeta_iterator>
        const_zeta_ref_iterator;
    //! Reverse raw zeta iterator
    typedef iterators::zeta_ref_iterator<reverse_zeta_iterator>
        reverse_zeta_ref_iterator;
     //! Constant reverse raw zeta iterator
    typedef iterators::zeta_ref_iterator<const_reverse_zeta_iterator>
        const_reverse_zeta_ref_iterator;
    
    //! Real and imaginary iterator
    typedef boost::indirect_iterator<ptr_iterator> iterator;
    //! Constant real and imaginary iterator
    typedef boost::indirect_iterator<const_ptr_iterator> const_iterator;
    // Note there is no reverse iterator supplied, but it would be straight
    // forward to make one, I just didn't need one, so I didn't write one.
    
    
  public:
    /**
     *  @brief  Default constructor.
     */
    zeta_array( ) : zetaArray_(), pointers_()
    { sys_construct_pointers(); }
    
    
    /**
     *  @brief  Copy constructor.
     *  @param  that The zeta array to copy
     */
    zeta_array( const self_type& __that )
  : zetaArray_( __that.zetaArray_ ), pointers_()
    { sys_construct_pointers(); }
    
    
    /**
     *  @brief  Assignment operator
     *  @param  that The zeta array to assign from
     *  @return \c *this
     */
    self_type& operator=( const self_type& __that ) {
      zetaArray_ = __that.zetaArray_;
      return *this;
    }
    
    
    /**
     *  @brief   Get a constant reference to the zeta coordinate at \c loc
     *  @param   loc The index of the zeta coordinate to return
     *  @return  A constant reference to the zeta coordinate at \c loc
     */
    const zeta_coordinate_type& zeta_at( size_type __loc ) const {
      return zetaArray_.at(__loc);
    }
    
    /**
     *  @brief   Get a reference to the zeta coordinate at \c loc
     *  @param   loc The index of the zeta coordinate to return
     *  @return  A reference to the zeta coordinate at \c loc
     */
    zeta_coordinate_type& zeta_at( size_type __loc ) {
      return zetaArray_.at(__loc);
    }
    
    
    /**
     *  @brief   Get a constant reference to the first zeta coordinate.
     *  @return  A constant reference to the first zeta coordinate.
     */
    const zeta_coordinate_type& front( void ) const {
      return zetaArray_.front();
    }
    
    /**
     *  @brief   Get a reference to the first zeta coordinate.
     *  @return  A reference to the first zeta coordinate.
     */
    zeta_coordinate_type& front( void ) {
      return zetaArray_.front();
    }
    
    
    /**
     *  @brief   Get a constant reference to the last zeta coordinate.
     *  @return  A constant reference to the last zeta coordinate.
     */
    const zeta_coordinate_type& back( void ) const {
      return zetaArray_.back();
    }
    
    /**
     *  @brief   Get a reference to the last zeta coordinate.
     *  @return  A reference to the last zeta coordinate.
     */
    zeta_coordinate_type& back( void ) {
      return zetaArray_.back();
    }
    
    
    /**
     *  @brief   Get a constant reference to the raw real object at \c loc
     *  @param   loc The index of the real object to return
     *  @return  A constant reference to the raw real object at \c loc
     */
    const real_type& real_ref_at( size_type __loc ) const {
      return zetaArray_.at(__loc).real_ref( );
    }
    
    /**
     *  @brief   Get a reference to the raw real object at \c loc
     *  @param   loc The index of the real object to return
     *  @return  A reference to the raw real object at \c loc
     */
    real_type& real_ref_at( size_type __loc ) {
      return zetaArray_.at(__loc).real_ref( );
    }
    
    
    /**
     *  @brief   Get a constant reference to the raw object value at \c loc
     *  @param   loc The index of the object value to return
     *  @return  A constant reference to the raw object value at \c loc
     */
    const imag_type& imag_ref_at( size_type __loc ) const {
      return zetaArray_.at(__loc).imag_ref( );
    }
    
    /**
     *  @brief   Get a reference to the raw imag object at \c loc
     *  @param   loc The index of the imag object to return
     *  @return  A reference to the raw imag object at \c loc
     */
    imag_type& imag_ref_at( size_type __loc ) {
      return zetaArray_.at(__loc).imag_ref( );
    }
    
    
    /**
     *  @brief   Get a const reference to the common typed  object at \c loc
     *  @param   loc The index of the object to return
     *  @return  A const reference to the common typed object at \c loc
     *  @note    Index \in [0,2*zeta_size - 1]
     */
    const return_type& common_at( size_type __loc ) const {
      return *( pointers_.at( __loc ) );
    }
    
    
    /**
     *  @brief   Get a reference to the common typed  object at \c loc
     *  @param   loc The index of the object to return
     *  @return  A reference to the common typed object at \c loc
     *  @note    Index \in [0,2*zeta_size - 1]
     */
    return_type& common_at( size_type __loc ) {
      return const_cast<real_type&>
                (static_cast<const self_type&>(*this).common_at(__loc) );
    }
    
    
    /**
     *  @brief   Get a const reference to the common typed real object at \c loc
     *  @param   loc The index of the real object to return
     *  @return  A const reference to the common typed real object at \c loc
     *  @note    Index \in [0,zeta_size - 1]
     */
    const return_type& real_common_at( size_type __loc ) const {
      return zetaArray_.at(__loc).real_common( );
    }
    
    /**
     *  @brief   Get a reference to the common typed real object at \c loc
     *  @param   loc The index of the real object to return
     *  @return  A reference to the common typed real object at \c loc
     *  @note    Index \in [0,zeta_size - 1]
     */
    return_type& real_common_at( size_type __loc ) {
      return zetaArray_.at(__loc).real_common( );
    }
    
    
    /**
     *  @brief   Get a const reference to the common typed imag object at \c loc
     *  @param   loc The index of the object to return
     *  @return  A const reference to the common typed imag object at \c loc
     *  @note    Index \in [0,zeta_size - 1]
     */
    const return_type& imag_common_at( size_type __loc ) const {
      return zetaArray_.at(__loc).imag_common( );
    }
    
    /**
     *  @brief   Get a reference to the common typed imag object at \c loc
     *  @param   loc The index of the imag object to return
     *  @return  A reference to the common typed imag object at \c loc
     *  @note    Index \in [0,zeta_size - 1]
     */
    return_type& imag_common_at( size_type __loc ) {
      return zetaArray_.at(__loc).imag_common( );
    }
    
    
    /**
     *  @def     __SG_ZETA_ITERATOR_DEF__(ITERATOR,DIRECTION,CONST)
     *  @brief   Defines the iterator access functions to zeta coordinates
     *  @note    Use this when a common interface is needed.
     *  @note    Defines:
     *           <ul>
     *              <li>\c zeta_begin
     *              <li>\c zeta_end
     *              <li>\c zeta_rbegin
     *              <li>\c zeta_rend
     *           </ul>
     *           Including \c const versions.
     */
#define __SG_ZETA_ITERATOR_DEF__(ITERATOR,DIRECTION,CONST)                    \
        ITERATOR zeta_##DIRECTION( void ) CONST {                             \
          return zetaArray_.DIRECTION( );                                     \
        }
    __SG_ZETA_ITERATOR_DEF__(zeta_iterator,begin,);
    __SG_ZETA_ITERATOR_DEF__(zeta_iterator,end,);
    __SG_ZETA_ITERATOR_DEF__(const_zeta_iterator,begin,const);
    __SG_ZETA_ITERATOR_DEF__(const_zeta_iterator,end,const);
    __SG_ZETA_ITERATOR_DEF__(reverse_zeta_iterator,rbegin,);
    __SG_ZETA_ITERATOR_DEF__(reverse_zeta_iterator,rend,);
    __SG_ZETA_ITERATOR_DEF__(const_reverse_zeta_iterator,rbegin,const);
    __SG_ZETA_ITERATOR_DEF__(const_reverse_zeta_iterator,rend,const);
#undef __SG_ZETA_ITERATOR_DEF__
    
    
    /**
     *  @def     __SG_REF_ITERATOR_DEF__(ITERATOR,DIRECTION,CONST)
     *  @brief   Defines the iterator access functions to ref zeta values.
     *  @note    This iterates through the underlying zeta values NOT the
     *           coordinates themselves. Use this when the actual zeta values
     *           are needed.
     *  @note    Defines:
     *           <ul>
     *              <li>\c zeta_ref_begin
     *              <li>\c zeta_ref_end
     *              <li>\c zeta_ref_rbegin
     *              <li>\c zeta_ref_rend
     *           </ul>
     *           Including \c const versions.
     */
#define __SG_REF_ITERATOR_DEF__(ITERATOR,DIRECTION,CONST)                     \
        ITERATOR zeta_ref_##DIRECTION( void ) CONST {                         \
          return ITERATOR( zetaArray_.DIRECTION( ) );                         \
        }
    __SG_REF_ITERATOR_DEF__(zeta_ref_iterator,begin,);
    __SG_REF_ITERATOR_DEF__(zeta_ref_iterator,end,);
    __SG_REF_ITERATOR_DEF__(const_zeta_ref_iterator,begin,const);
    __SG_REF_ITERATOR_DEF__(const_zeta_ref_iterator,end,const);
    __SG_REF_ITERATOR_DEF__(reverse_zeta_ref_iterator,rbegin,);
    __SG_REF_ITERATOR_DEF__(reverse_zeta_ref_iterator,rend,);
    __SG_REF_ITERATOR_DEF__(const_reverse_zeta_ref_iterator,rbegin,const);
    __SG_REF_ITERATOR_DEF__(const_reverse_zeta_ref_iterator,rend,const);
#undef __SG_REF_ITERATOR_DEF__
    
    
    /**
     *  @brief  Provide iterator access to the common structure.
     *  @return An iterator pointing the real part of the first zeta component
     */
    iterator begin( void ) {
      return iterator( pointers_.begin() );
    }
    
    /**
     *  @brief  Provide constant iterator access to the common structure.
     *  @return A constant iterator pointing the real part of the first zeta
     *          component.
     */
    const_iterator begin( void ) const {
      return const_iterator( pointers_.begin() );
    }
    
    /**
     *  @brief  Provide a sentinal iterator for the end of the common structure.
     *  @return An iterator pointing one past the end of the common structure.
     */
    iterator end( void ) { return zeta_reim_end(); }
    
    /**
     *  @brief  Provide a sentinal iterator for the end of the common structure.
     *  @return A const iterator pointing one past the end of the common
     *          structure.
     */
    const_iterator end( void ) const { return zeta_reim_end(); }
    
    /**
     *  @brief  Provide a sentinal iterator for the end of the common zeta
     *          structure.
     *  @return An iterator pointing one past the end of the common zeta
     *          structure.
     */
    iterator zeta_reim_end( void ) {
      return iterator( pointers_.begin() + self_type::zeta_real_size );
    }
    
    /**
     *  @brief  Provide a constant sentinal iterator for the end of the common
     *          zeta structure.
     *  @return A constant iterator pointing one past the end of the common
     *          zeta structure.
     */
    const_iterator zeta_reim_end( void ) const {
      return const_iterator( pointers_.begin() + self_type::zeta_real_size );
    }
    
    
    /**
     *  @brief  Compares two \c zeta_array objects for equality.
     *  @param  rhs A \c zeta_array
     */
    bool operator==( const self_type& __rhs ) const {
      return bool( std::equal( zeta_begin(), zeta_end(), __rhs.zeta_begin()) );
    }
      
      
    /**
     *  @brief  Defines the output stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::ostream& operator<<( std::ostream& __os, const self_type& __t ){
      typedef utility::io::container_to_string<const_zeta_iterator> _cts;
      
      __os << _cts( )( __t.zeta_begin(), __t.zeta_end() );
      return __os;
    }
    
    /**
     *  @brief  Defines the input stream operator
     *  @param  __os An ostream.
     *  @param  __t The object to stream out.
     */
    friend std::istream& operator>>( std::istream& __is, self_type& __t ){
      typedef  utility::io::string_to_array<zeta_iterator> _sta;
      
      std::string str;
      __is >> str;
      _sta( )( str, __t.zeta_begin(), __t.zeta_end() );
      return __is;
    }
    
  private:
    /**
     *  @brief  Initialize the pointers to the common objects.
     *  @note   Should only be called at construction
     */
    void sys_construct_pointers( void ) {
      typename array_type::iterator ait = zetaArray_.begin();
      typename ptr_array_type::iterator pit
          = const_cast<ptr_array_type* const>( &pointers_ )->begin();
      for( ; ait != zetaArray_.end(); ++ait, ++pit ) {
        *pit = &( ait->real_common() );
        ++pit;
        *pit = &( ait->imag_common() );
      }
    }
    
    //! The underlying array of zeta values
    array_type zetaArray_;
    
  protected:
    /**
     *  @brief  Add one of the initialization pointers
     *  @param  loc The location to add at
     *  @param  ptr The pointer to add
     */
    void add_ptr_at( size_type __loc, return_type* __ptr ) {
      (*(const_cast<ptr_array_type* const>( &pointers_ )))[__loc] = __ptr;
    }
    
    //! The array of pointers to the common objects
    const ptr_array_type pointers_;
    
};

}// detail
}// geometric
}// structure
}// sg
#endif
