#pragma once 
#ifndef CONTEXT_HPP
#define CONTEXT_HPP


/***********************************************************************************
 ***********************************************************************************
 **                                                                               **
 **  Copyright (C) 2011 Tom Shorrock <t.h.shorrock@gmail.com> 
 **                                                                               **
 **                                                                               **
 **  This program is free software; you can redistribute it and/or                **
 **  modify it under the terms of the GNU General Public License                  **
 **  as published by the Free Software Foundation; either version 2               **
 **  of the License, or (at your option) any later version.                       **
 **                                                                               **
 **  This program is distributed in the hope that it will be useful,              **
 **  but WITHOUT ANY WARRANTY; without even the implied warranty of               **
 **  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                **
 **  GNU General Public License for more details.                                 **
 **                                                                               **
 **  You should have received a copy of the GNU General Public License            **
 **  along with this program; if not, write to the Free Software                  **
 **  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  **
 **                                                                               **
 ***********************************************************************************
 ***********************************************************************************/
// Include all of Proto
#include <boost/proto/proto.hpp>
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

#include <boost/type_traits.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/placeholders.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/accumulate.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/plus.hpp>
#include <boost/mpl/times.hpp>
#include <boost/bind.hpp>
#include <algorithm>

#include <boost/fusion/include/transform.hpp>
#include <boost/fusion/include/for_each.hpp>
// Create some namespace aliases
namespace mpl = boost::mpl;
namespace fusion = boost::fusion;
namespace proto = boost::proto;


//From this lib
#include "EnsembleLearning/node/Node.hpp"
#include "Placeholder.hpp"

//From boost
#include <boost/call_traits.hpp>

#include <omp.h>
#include <iostream>
#include <map>
namespace ICR{

  namespace EnsembleLearning {
    
    //Forward declaration.
    template<class> class Placeholder;
    
    //Forward declaration.
    namespace detail{
      template<template<class> class ,  class>
      class Deterministic ;
    }

    template<int TerminalExpr>
    class calculator_context
    {
    public:
      void push_back(const std::vector<double>& a)
      {
	std::copy(a.begin(), a.end(), std::back_inserter(m_args));
      }

      calculator_context() { m_args.push_back(0.0); } //m_args(0) is the result placeholder.
  
      template<typename Expr
	       // defaulted template parameters, so we can
	       // specialize on the expressions that need
	       // special handling.
	       , typename Tag = typename proto::tag_of<Expr>::type
	       , typename Arg0 = typename proto::result_of::child_c<Expr, 0>::type
	       > struct eval;
      
    // Handle placeholder terminals here...
    template<typename Expr, int I>
    struct eval<Expr, proto::tag::terminal, placeholder<I> >
    {
      typedef std::pair<double,double> result_type;
      
      result_type operator()(Expr &, calculator_context &ctx) const
      {
	return result_type(ctx.m_args[I],1);
      }
    };

    // Handle other terminals here...
    template<typename Expr, typename Arg0>
    struct eval<Expr, proto::tag::terminal, Arg0>
    {
      typedef std::pair<double,double> result_type;

      result_type operator()(Expr &expr, calculator_context &) const
      {
	return result_type(proto::child(expr),1);
      }
    };

    // Handle addition here...
    template<typename Expr, typename Arg0>
    struct eval<Expr, proto::tag::plus, Arg0>
    {
      typedef std::pair<double,double> result_type;

    
      result_type operator()(Expr &expr, calculator_context &ctx) const
      {
	typedef typename proto::result_of::child_c<Expr, 0>::type ltype_tmp;
	typedef typename proto::result_of::child_c<Expr, 1>::type rtype_tmp;

	typedef typename proto::result_of::child_c<rtype_tmp, 0>::type rtype;
	typedef typename proto::result_of::child_c<ltype_tmp, 0>::type ltype;


	if (boost::is_base_of<mpl::int_<TerminalExpr>,ltype>::value) 
	  {
	    //lhs is equal to the special placeholder.
	    const double rhs1 = proto::eval(proto::right(expr), ctx).first;
	    const double rhs2 = proto::eval(proto::right(expr), ctx).second;
	    return result_type(rhs1, rhs2);
	  }
	else if (boost::is_base_of<mpl::int_<TerminalExpr>,rtype>::value)
	  {
	    //rhs is equal to the special placeholder.
	    const double lhs1 = proto::eval(proto::left(expr), ctx).first;
	    const double lhs2 = proto::eval(proto::left(expr), ctx).second;
	    return result_type(lhs1, lhs2);
	  } 
	else
	  {
	    //no special placeholder
	    const double lhs1 = proto::eval(proto::left(expr), ctx).first;
	    const double lhs2 = proto::eval(proto::left(expr), ctx).second;
	    const double rhs1 = proto::eval(proto::right(expr), ctx).first;
	    const double rhs2 = proto::eval(proto::right(expr), ctx).second;
	    return result_type(lhs1+rhs1, lhs2*rhs2);
	  }
      }
    };

    // Handle multiplication here...
    template<typename Expr, typename Arg0>
    struct eval<Expr, proto::tag::multiplies, Arg0>
    {
      typedef std::pair<double,double> result_type;

      result_type 
      operator()(Expr &expr, calculator_context &ctx) const
      {

    
	typedef typename proto::result_of::child_c<Expr, 0>::type ltype_tmp;
	typedef typename proto::result_of::child_c<Expr, 1>::type rtype_tmp;

	typedef typename proto::result_of::child_c<rtype_tmp, 0>::type rtype;
	typedef typename proto::result_of::child_c<ltype_tmp, 0>::type ltype;


	if (boost::is_same<mpl::int_<TerminalExpr>,mpl::int_<0> >::value) 
	  {
	    //special case for results expression.
	    const double lhs1 = proto::eval(proto::left(expr), ctx).first;
	    const double rhs1 = proto::eval(proto::right(expr), ctx).first;
	    //no special placeholder
	    return result_type(lhs1*rhs1,1);
	  }
	else if (boost::is_base_of<mpl::int_<TerminalExpr>,ltype>::value) 
	  {
	    //lhs is equal to the special placeholder.
	    return result_type(0, proto::eval(proto::right(expr), ctx).first);
	  }
	else if (boost::is_base_of<mpl::int_<TerminalExpr>,rtype>::value)
	  {
	    //rhs is equal to the special placeholder.
	    return result_type(0, proto::eval(proto::left(expr), ctx).first);
	
	  } 
	else if (boost::is_same<proto::tag::terminal,
		 typename proto::tag_of<ltype_tmp>::type>::value)
	  { //lhs is terminal but not special placeholder.
	    const double lhs1 = proto::eval(proto::left(expr), ctx).first;
	    const double rhs1 = proto::eval(proto::right(expr), ctx).first;
	    const double rhs2 = proto::eval(proto::right(expr), ctx).second;
	    return result_type(lhs1*rhs1,lhs1*rhs2);
	  }
	else if (boost::is_same<proto::tag::terminal,
		 typename proto::tag_of<rtype_tmp>::type>::value)
	  { //rhs is terminal but not special placeholder.
	    const double lhs1 = proto::eval(proto::left(expr), ctx).first;
	    const double rhs1 = proto::eval(proto::right(expr), ctx).first;
	    const double lhs2 = proto::eval(proto::left(expr), ctx).second;
	    return result_type(lhs1*rhs1,lhs2*rhs1);
	  }
	else
	  {
	    const double lhs1 = proto::eval(proto::left(expr), ctx).first;
	    const double rhs1 = proto::eval(proto::right(expr), ctx).first;
	    const double lhs2 = proto::eval(proto::left(expr), ctx).second;
	    const double rhs2 = proto::eval(proto::right(expr), ctx).second;
	    //no special placeholder
	    return result_type(lhs1*rhs1,lhs2*rhs2);
	  }
      }
    };


      // Handle the placeholders:
      template<int I>
      double operator()(proto::tag::terminal, placeholder<I>) const
      {
	return this->m_args[I];
      }
private:
  // The values with which we'll replace the placeholders
  std::vector<double> m_args;
    };
    
  // Define the grammar of calculator expressions
  struct calculator_grammar
    : proto::or_<
    proto::plus< calculator_grammar, calculator_grammar >
    , proto::multiplies< calculator_grammar, calculator_grammar >
    , proto::terminal< proto::_ >
    >
  {};

    // Forward-declare an expression wrapper
    template<typename Expr>
    struct calculator;

  // Define a calculator domain. Expression within
  // the calculator domain will be wrapped in the
  // calculator<> expression wrapper.
  struct calculator_domain
    : proto::domain< proto::generator<calculator>, calculator_grammar >//, calculator_grammar 
  {};
  // Define a calculator expression wrapper. It behaves just like
  // the expression it wraps, but with an extra operator() member
  // function that evaluates the expression.    
  template<typename Expr>
struct calculator
  : proto::extends<Expr, calculator<Expr>, calculator_domain>
{
  typedef
  proto::extends<Expr, calculator<Expr>, calculator_domain>
  base_type;

  calculator(Expr const &expr = Expr())
    : base_type(expr)
  {}

  typedef double result_type;

  };




    /** A class that maps the placeholders to the Variables that they represent in an expression.
     *  The same expression can be used for multiple contexts.
     *
     *  @tparam T The type of the stored data.  
     *   This will be either double or float - with the latter used only for memory contrained models.
     *
     *  Example:
     *  @code 
     *
     *  //Create an expression X*Y + Z.
     *  ExpressionFactory<double> Factory;  // Create a Factory.
     *  
     *  Placeholder<double>* X = Factory.placeholder();
     *  Placeholder<double>* Y = Factory.placeholder();
     *  Placeholder<double>* Z = Factory.placeholder();
     *  Expression<double>* XY = Factory.Multiply(X,Y);
     *  Expression<double>* Expr = Factory.Add(XY,Z);
     *  
     *  //Create some variables
     *  Builder<double> builder;
     *  Builder<double>::GaussianNode x = builder.gaussian(0.01,0.01);
     *  Builder<double>::GaussianNode y = builder.gaussian(0.01,0.01);
     *  Builder<double>::GaussianNode z = builder.gaussian(0.01,0.01);
     *  
     *  //Assign this context to the expression.
     *  Context<double> context;
     *  context.Assign(X,x);
     *  context.Assign(Y,y);
     *  context.Assign(Z,z);
     *  
     *  @endcode
     *
     */
    template<class T>
    class Context{
      typedef std::vector<VariableNode<T>* const  > DataContainer;
      typedef VariableNode<T>* const  Datum;
    public:
      /** @name Useful typdefs for types that are exposed to the user.
       */
      ///@{
      
      
      typedef typename boost::call_traits< VariableNode<T>* const>::param_type
      variable_parameter;
      typedef typename boost::call_traits< VariableNode<T>* const>::value_type
      variable_t;

      typedef typename boost::call_traits<size_t>::param_type
      size_parameter;
      
      ///@}
      
      Context(const size_t size) : m_map(size) {}
      
      /** Assign a placeholder a VariableNode.
       *  @param i The Placeholder index.
       *  @param V The VariableNode 
       */
      void 
      Assign(size_t i, variable_parameter  V );


      variable_t
      operator[](size_parameter i) const
      {
	return m_map[i];
      }
      
      size_t
      size() const
      { return m_map.size();}
      
      /** Output the Context to a stream. 
       *  @param c The context.
       *  @param out The output stream.
       *  @return A reference to the output stream.
       *  This function is thread safe.
       */
      template<class U>
      friend
      std::ostream&
      operator<<(std::ostream& out, const Context<U>&  c);

    private:

      template<template<class> class Model,  class U>
      friend class detail::Deterministic;

      mutable DataContainer m_map;
    };


    /*********************************************************
    **********************************************************
    **STREAM OPERATORS IMPLEMENTATION ************************
    **********************************************************
    **********************************************************/

    template<class T>
    std::ostream&
    operator<<(std::ostream& out, 
	       const SubContext<T>& c)
    {
      //lock
#pragma omp critical
      {
	
	typename std::map<const Placeholder<T>*,  T>::const_iterator it;
	for(it = c.m_map.begin();
	    it!= c.m_map.end();
	    ++it)
	  {
	    out<<*it<<" ";
	  }
      }
      return out;
    }


      
    template<class T>
    std::ostream&
    operator<<(std::ostream& out, 
	       const Context<T>&  c)
    { 

#pragma omp critical
      {
	typedef std::map<VariableNode<T>* const ,const Placeholder<T>* > DataContainer;
	typename DataContainer::const_iterator it;
	for(it = c.m_map.begin();
	    it!= c.m_map.end();
	    ++it)
	  {
	    out<<it->first<<" ";
	  }
      }
      return out;
    }


  }
}


/*********************************************************
**********************************************************
***************** IMPLEMENTATION *************************
**********************************************************
**********************************************************/

template<class T>
inline
typename ICR::EnsembleLearning::Context<T>::placeholder_t
ICR::EnsembleLearning::Context<T>::Lookup(variable_parameter V) const
{
  return
    m_map.find(V)->second;
}


template<class T>
void
ICR::EnsembleLearning::Context<T>::Assign(size_t i, 
			     variable_parameter  V )
{
  //make sure not trying to assign to things at once.
#pragma omp critical
  {
    if (m_map.size()<=i) 
      m_map.resize(i+1);
    m_map[i] = V;
  }
}

#endif
