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

//
//  Copyright (c) 2000-2002
//  Joerg Walter, Mathias Koch
//
//  Permission to use, copy, modify, distribute and sell this software
//  and its documentation for any purpose is hereby granted without fee,
//  provided that the above copyright notice appear in all copies and
//  that both that copyright notice and this permission notice appear
//  in supporting documentation.  The authors make no representations
//  about the suitability of this software for any purpose.
//  It is provided "as is" without express or implied warranty.
//
//  The authors gratefully acknowledge the support of
//  GeNeSys mbH & Co. KG in producing this work.
//


// NOTE: There is a small bug in the boost version which prevernts reverse
//       iteeffectn through negative ranges during debugging, this fixes that
//       problem.

#ifndef SG_RANGE_FIX_H
#define SG_RANGE_FIX_H

#include <algorithm>
#ifdef BOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR
#include <boost/shared_array.hpp>
#endif

#include <boost/numeric/ublas/exception.hpp>
#include <boost/numeric/ublas/detail/iterator.hpp>

namespace boost { namespace numeric { namespace ublas {
    // Range class
    template <class Z, class D>
    class basic_range_fix {
        typedef basic_range_fix<Z, D> self_type;
    public:
        typedef Z size_type;
        typedef D difference_type;
        typedef size_type value_type;
        typedef value_type const_reference;
        typedef const_reference reference;
        typedef const value_type *const_pointer;
        typedef value_type *pointer;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        basic_range_fix ():
            start_ (0), size_ (0) {}
        BOOST_UBLAS_INLINE
        basic_range_fix (size_type start, size_type stop):
            start_ (start), size_ (stop - start) {
            BOOST_UBLAS_CHECK (start_ <= stop, bad_index ());
        }

        BOOST_UBLAS_INLINE
        size_type start () const {
            return start_;
        }
        BOOST_UBLAS_INLINE
        size_type size () const {
            return size_;
        }

        // Random Access Container
        BOOST_UBLAS_INLINE
        size_type max_size () const {
            return size_;
        }
        
        BOOST_UBLAS_INLINE
        bool empty () const {
            return size_ == 0;
        }
            
        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator () (size_type i) const {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
            return start_ + i;
        }

        // Composition
        BOOST_UBLAS_INLINE
        basic_range_fix compose (const basic_range_fix &r) const {
            return basic_range_fix (start_ + r.start_, start_ + r.start_ + r.size_);
        }

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator == (const basic_range_fix &r) const {
            return start_ == r.start_ && size_ == r.size_;
        }
        BOOST_UBLAS_INLINE
        bool operator != (const basic_range_fix &r) const {
            return ! (*this == r);
        }

        // Iterator types
    private:
        // Use and index
        typedef size_type const_subiterator_type;

    public:
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
        typedef indexed_const_iterator<self_type, std::random_access_iterator_tag> const_iterator;
#else
        class const_iterator:
            public container_const_reference<basic_range_fix>,
            public random_access_iterator_base<std::random_access_iterator_tag,
                                               const_iterator, value_type> {
        public:
            typedef typename basic_range_fix::value_type value_type;
            typedef typename basic_range_fix::difference_type difference_type;
            typedef typename basic_range_fix::const_reference reference;
            typedef typename basic_range_fix::const_pointer pointer;

            // Construction and destruction
            BOOST_UBLAS_INLINE
            const_iterator ():
                container_const_reference<basic_range_fix> (), it_ () {}
            BOOST_UBLAS_INLINE
            const_iterator (const basic_range_fix &r, const const_subiterator_type &it):
                container_const_reference<basic_range_fix> (r), it_ (it) {}

            // Arithmetic
            BOOST_UBLAS_INLINE
            const_iterator &operator ++ () {
                ++ it_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator &operator -- () {
                -- it_;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator &operator += (difference_type n) {
                it_ += n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            const_iterator &operator -= (difference_type n) {
                it_ -= n;
                return *this;
            }
            BOOST_UBLAS_INLINE
            difference_type operator - (const const_iterator &it) const {
                return it_ - it.it_;
            }

            // Dereference
            BOOST_UBLAS_INLINE
            const_reference operator * () const {
                BOOST_UBLAS_CHECK ((*this) ().start () <= it_, bad_index ());
                BOOST_UBLAS_CHECK (it_ < (*this) ().start () + (*this) ().size (), bad_index ());
                return it_;
            }

            BOOST_UBLAS_INLINE
            const_reference operator [] (difference_type n) const {
                return *(*this + n);
            }

            // Index
            BOOST_UBLAS_INLINE
            size_type index () const {
                BOOST_UBLAS_CHECK ((*this) ().start () <= it_, bad_index ());
                BOOST_UBLAS_CHECK (it_ < (*this) ().start () + (*this) ().size (), bad_index ());
                return it_ - (*this) ().start ();
            }

            // Assignment
            BOOST_UBLAS_INLINE
            const_iterator &operator = (const const_iterator &it) {
                // Comeau recommends...
                this->assign (&it ());
                it_ = it.it_;
                return *this;
            }

            // Comparison
            BOOST_UBLAS_INLINE
            bool operator == (const const_iterator &it) const {
                BOOST_UBLAS_CHECK ((*this) () == it (), external_logic ());
                return it_ == it.it_;
            }
            BOOST_UBLAS_INLINE
            bool operator < (const const_iterator &it) const {
                BOOST_UBLAS_CHECK ((*this) () == it (), external_logic ());
                return it_ < it.it_;
            }

        private:
            const_subiterator_type it_;
        };
#endif

        BOOST_UBLAS_INLINE
        const_iterator begin () const {
            return const_iterator (*this, start_);
        }
        BOOST_UBLAS_INLINE
        const_iterator end () const {
            return const_iterator (*this, start_ + size_);
        }

        // Reverse iterator
        typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

        BOOST_UBLAS_INLINE
        const_reverse_iterator rbegin () const {
            return const_reverse_iterator (end ());
        }
        BOOST_UBLAS_INLINE
        const_reverse_iterator rend () const {
            return const_reverse_iterator (begin ());
        }

        BOOST_UBLAS_INLINE
        basic_range_fix preprocess (size_type size) const {
            if (this != &all_)
                return *this;
            return basic_range_fix (0, size);
        }
        static
        BOOST_UBLAS_INLINE
        const basic_range_fix &all () {
            return all_;
        }

    private:
        size_type start_;
        size_type size_;
        static const basic_range_fix all_;
    };

    template <class Z, class D>
    const basic_range_fix<Z,D> basic_range_fix<Z,D>::all_  (0, size_type (-1));

}}}

#endif
