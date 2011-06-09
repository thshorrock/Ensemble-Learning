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


//From this lib
#include "EnsembleLearning/node/Node.hpp"

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


    /** A class that maps the placeholder in an expression
     *    to a set of Moments values of a Varible Nodes of every element in the expression.
     *
     *  @tparam T The type of the stored data.  
     *   This will be either double or float - with the latter used only for memory contrained models.
     *
     *  Example:
     *  @code 
     *  //The context provides the Moments of every element in expression.
     *  Context<double> C;  
     *  SubContext<T> C0 = C[0]; //All the first moments  (the <x>'s of every element in expr) 
     *  SubContext<T> C1 = C[1];  //The second moment (the <x^2> of every element of expression)
     *  @endcode
     *  
     *  @warning SubContext is not thread safe
     *    in that it does not protect its internal data from read/write operations.
     *    It is assumed it will be created and destroyed within a re-entrant member and not stored -
     *    see the example.
     *    
     */
    template<class T>
    class SubContext{
      
      
    public:
      
      /** @name Useful typdefs for types that are exposed to the user.
       */
      ///@{
      
      typedef typename boost::call_traits<T>::param_type
      data_parameter;
      
      typedef typename boost::call_traits<T>::value_type
      data_t;
            
      typedef typename boost::call_traits<T>::const_reference 
      data_const_reference;
      
      typedef typename boost::call_traits<const Placeholder<T>*>::param_type
      placeholder_parameter;

      typedef typename boost::call_traits<const Placeholder<T>*>::value_type
      placeholder_t;
      
      typedef typename boost::call_traits<SubContext<T> >::const_reference
      const_reference;

      typedef typename boost::call_traits<SubContext<T> >::reference
      reference;
      
      typedef typename boost::call_traits<SubContext<T> >::param_type
      parameter;

      typedef typename boost::call_traits<SubContext<T> >::value_type
      type;

      ///@}

      /** Obtain the value associated with for the placeholder.
       *   @param P The placeholder.
       *   @return A value taken from the moment of the Variable node represented by P.
       */
      data_const_reference
      Lookup(placeholder_parameter P) const;
      
      /** Assign a placeholder a value.
       *  @param P The placeholder.
       *  @param V The value 
       */
      void 
      Assign(placeholder_parameter P, 
	     data_parameter V);
      
      /** Multiply to another subcontext.
       *    This is an element-wise multiplication. 
       *    In effect, the variable associated with every placeholder is squared.
       *  @param C The other subcontext.
       *  @return A reference to the product.
       */
      reference
      operator*=(parameter C);

      /** Output the SubContext to a stream. 
       *  @param c The subcontext.
       *  @param out The output stream.
       *  @return A reference to the output stream.
       *  This function is thread safe.
       */
      template<class U>
      friend
      std::ostream&
      operator<<(std::ostream& out, const SubContext<U>&  c);

      /** Multiply two subcontexts together.
       *  @param A The first subcontext.
       *  @param B The second subcontext.
       *  @return The product.
       *
       *  @see operator*=(parameter)
       *  For what it means to multiply subcontexts, 
       */
      friend 
      type
      operator*(parameter A, parameter B)
      {
	SubContext<T> tmp = A;
	return tmp*=B;
      }
    
    private:
      //The store the context
      std::vector<data_t> m_context_data;
      
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
      typedef std::map<VariableNode<T>* const ,const Placeholder<T>* > DataContainer;
      typedef std::pair<  VariableNode<T>*const,const Placeholder<T>*> Datum;
    public:
      /** @name Useful typdefs for types that are exposed to the user.
       */
      ///@{
      
      typedef typename boost::call_traits<SubContext<T> >::value_type
      subcontext_t;
      
      typedef typename boost::call_traits<const Placeholder<T>*>::param_type
      placeholder_parameter;

      typedef typename boost::call_traits<const Placeholder<T>*>::value_type
      placeholder_t;
      
      typedef typename boost::call_traits< VariableNode<T>* const>::param_type
      variable_parameter;
      typedef typename boost::call_traits< VariableNode<T>* const>::value_type
      variable_t;

      typedef typename boost::call_traits<size_t>::param_type
      size_parameter;
      
      ///@}
      
      /** Obtain the placeholder associated  with a particular VariableNoder.
       *   @param V The VariableNode.
       *   @return The placeholder associated with the variable.
       */
      placeholder_t
      Lookup(variable_parameter V) const;
      
      /** Assign a placeholder a VariableNode.
       *  @param P The Placeholder.
       *  @param V The VariableNode 
       */
      void 
      Assign(placeholder_parameter P, 
	     variable_parameter  V );


      /** Splice a Context into Subcontext.
       *  To evaluate an expression the average means of a variable, \f$<x>\f$,
       *   or the average of the square of the means \f$<x^2>\f$ need to be obtained 
       *   from the Varibale Moments that contain all such information.
       *   Subcontext's are obtained from this function.
       *  @param i The index of the subcontext requested
       *  @return The resulting subcontext.
       *
       *  Example: 
       *  @code
       *  //The context provides the Moments of every element in expression.
       *  Context<double>& M
       *  SubContext<double> M0 = M[0];  //All the first moments  (the <x>'s of every element in expr)
       *  SubContext<double> M1 = M[1];  //The second moment (the <x^2> of every element of expression)
       *  //Precision is 1.0/ (<expr(x^2)> - <expr(x)>^2)
       *  double precision = 1.0/(Expr->Evaluate(M1) - Expr->Evaluate(M0*M0) );
       *  @endcode
       */
      subcontext_t
      operator[](size_parameter i) const
      {
	//Defer thread safety to the variable
	// (m_map not being altered here)
	SubContext<T> c;
	for(typename DataContainer::const_iterator it = m_map.begin();
	    it != m_map.end();
	    ++it)
	  {
	    variable_t V = it->first;
	    c.Assign( it->second,V->GetMoments()[i]);
	  }
	return c;
      }

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

      template<template<class> class Model,  class U>
      void
      AddChildFactor(detail::Deterministic<Model,U>* factor ) const
      {
	//Defer thread safety to the variable
	// (m_map not being altered here)
	typename DataContainer::const_iterator it;
	for(  it = m_map.begin();
	      it!= m_map.end();
	      ++it)
	  {
	    it->first->AddChildFactor(factor);
	  }
      }
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
typename ICR::EnsembleLearning::SubContext<T>::data_const_reference
ICR::EnsembleLearning::SubContext<T>::Lookup(placeholder_parameter P) const
{
  return  m_context_data[P->id()];
};


template<class T>
void
ICR::EnsembleLearning::SubContext<T>::Assign(placeholder_parameter P, 
				data_parameter V)
{
  typedef std::map<placeholder_t,data_t> DataContainer;
  typedef std::pair<placeholder_t,data_t> Datum;
      
  typename DataContainer::iterator it;
  std::pair<typename DataContainer::iterator,bool> ret;

  //make sure not trying to assign to things at once.
#pragma omp critical
  {
    if (P->id()<m_context_data.size())
      m_context_data[P->id()] = V;
    else
      {
	m_context_data.resize(P->id()+1);
	m_context_data[P->id()] = V;
      }
    // ret = m_map.insert( Datum(P, V) );
    // if (ret.second == false) //already exists
    //   {
    // 	ret.first->second = V;
    //   }
  }
}


template<class T>
inline
typename ICR::EnsembleLearning::SubContext<T>::reference
ICR::EnsembleLearning::SubContext<T>::operator*=(parameter C)
{
  typedef std::map<placeholder_t,data_t> DataContainer;
  typedef std::pair<placeholder_t,data_t> Datum;
  typename DataContainer::iterator it;
  
  BOOST_ASSERT(m_context_data.size() == C.m_context_data.size());
  for(size_t i=0;i<m_context_data.size();++i){
    m_context_data[i]*=C.m_context_data[i];
  }

  return *this;
}

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
ICR::EnsembleLearning::Context<T>::Assign(placeholder_parameter P, 
			     variable_parameter  V )
{
  std::pair<typename DataContainer::iterator,bool> ret;

  //make sure not trying to assign to things at once.
#pragma omp critical
  {
    ret = m_map.insert( Datum ( V,P) );
    if (ret.second == false) //already exists
      {
	ret.first->second = P;
      }
  }
}

#endif
