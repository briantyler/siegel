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
 *  @file     mem_fun_adaptor.hpp
 *  @brief    An inline header file for the \c mem_fun_adapter classes.
 *  @note     These functor adaptors create functor objects that allow the
 *            member functions to be called by a functor during algorithmic
 *            iteration.
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-25
 *  @note     I wrote these functors when I was into functional programming for
 *            its own sake; considering how long it takes to get the syntax
 *            right to use these functors it is, in general, better to write a
 *            for loop. Either way they are a bit of fun. (Unbilievably I
 *            recently found some guy who thought code should never contain
 *            for and while loops or if statements!!! Some people are crazy)
 */


#ifndef _SG_MEM_FUN_ADAPTOR_H
#define _SG_MEM_FUN_ADAPTOR_H 1


namespace sg
{
namespace utility
{
namespace functors
{
/**
 *  @class    mem_fun_adaptor
 *  @brief    A function object to allow a member function taking no arguments
 *            to be called by a functor during an algorithmic iteration.
 *  @param    _Ret Return type of the member function.
 *  @param    _Tp The type of the object the member function points to.
 *  @param    _Fn The unary functor.
 *  @param    _FnRet The result type of _Fn
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-25
 */
template <class _Ret, class _Tp, class _Fn, class _FnRet>
    class mem_fun_adaptor
  : public std::unary_function<_Tp*, _FnRet>
{
  public:
    //! Object type
    typedef mem_fun_adaptor<_Ret,_Tp,_Fn,_FnRet> self_type;
    
    /**
     *  @brief  Constructor
     *  @param  *__pf Function pointer to member function
     *  @param  fn  The functor to wrap.
     */
    explicit
        mem_fun_adaptor( _Ret (_Tp::*__pf)(), const _Fn& __fn )
  : _M_f(__pf), fn_(__fn) {}
    
    _FnRet operator()(_Tp* __p) const { return fn_( (__p->*_M_f)() ); }
    _FnRet operator()(_Tp* __p) { return fn_( (__p->*_M_f)() ); }
    
  private:
    _Ret (_Tp::*_M_f)();
    _Fn fn_;
};


/**
 *  @class    mem_fun_ref_adaptor
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-28
 *  @note     As \c mem_fun_adaptor but takes reference not pointer
 */
template <class _Ret, class _Tp, class _Fn, class _FnRet>
    class mem_fun_ref_adaptor
  : public std::unary_function<_Tp&, _FnRet>
{
  public:
    //! Object type
    typedef mem_fun_ref_adaptor<_Ret,_Tp,_Fn,_FnRet> self_type;
    
    /**
     *  @brief  Constructor
     *  @param  *__pf Function pointer to member function
     *  @param  fn  The functor to wrap.
     */
    explicit
        mem_fun_ref_adaptor( _Ret (_Tp::*__pf)(), const _Fn& __fn )
  : _M_f(__pf), fn_(__fn) {}
    
    _FnRet operator()(_Tp& __r) const { return fn_( (__r.*_M_f)() ); }
    _FnRet operator()(_Tp& __r) { return fn_( (__r.*_M_f)() ); }
    
  private:
    _Ret (_Tp::*_M_f)();
    _Fn fn_;
};


/**
 *  @class    const_mem_fun_adapter
 *  @brief    A function object to allow a const member function taking no
 *            arguments to be called by a functor during an algorithmic
 *            iteration.
 *  @param    _Ret Return type of the member function.
 *  @param    _Tp The type of the object the member function points to.
 *  @param    _Fn The unary functor.
 *  @param    _FnRet The result type of _Fn
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-25
 */
template <class _Ret, class _Tp, class _Fn, class _FnRet>
    class const_mem_fun_adaptor
  : public std::unary_function<const _Tp*, _FnRet>
{
  public:
    //! Object type
    typedef const_mem_fun_adaptor<_Ret,_Tp,_Fn,_FnRet> self_type;
    
    /**
     *  @brief  Constructor
     *  @param  *__pf Function pointer to member function
     *  @param  fn  The functor to wrap.
     */
    explicit
        const_mem_fun_adaptor(_Ret (_Tp::*__pf)() const, const _Fn& __fn)
  : _M_f(__pf), fn_(__fn) {}
    
    _FnRet operator()(const _Tp* __p) const { return fn_( (__p->*_M_f)() ); }
    _FnRet operator()(const _Tp* __p) { return fn_( (__p->*_M_f)() ); }
    
  private:
    _Ret (_Tp::*_M_f)() const;
    _Fn fn_;
};


/**
 *  @class    const_mem_fun_ref_adaptor
 *  @author   Brian Tyler
 *  @version  1.0
 *  @date     2008-03-28
 *  @note     As \c const_mem_fun_adapter but takes reference not pointer
 */
template <class _Ret, class _Tp, class _Fn, class _FnRet>
    class const_mem_fun_ref_adaptor
  : public std::unary_function<const _Tp&, _FnRet>
{
  public:
    //! Object type
    typedef const_mem_fun_ref_adaptor<_Ret,_Tp,_Fn,_FnRet> self_type;
    
    /**
     *  @brief  Constructor
     *  @param  *__pf Function pointer to member function
     *  @param  fn  The functor to wrap.
     */
    explicit
        const_mem_fun_ref_adaptor(_Ret (_Tp::*__pf)() const, const _Fn& __fn)
  : _M_f(__pf), fn_(__fn) {}
    
    _FnRet operator()(const _Tp& __r) const { return fn_( (__r.*_M_f)() ); }
    _FnRet operator()(const _Tp& __r) { return fn_( (__r.*_M_f)() ); }
    
  private:
    _Ret (_Tp::*_M_f)() const;
    _Fn fn_;
};

}// functors
}// utility
}// sg
#endif
