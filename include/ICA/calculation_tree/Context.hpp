#pragma once

#include <map>
#include <string>
#include "ICA/node/Node.hpp"

namespace ICR{

  namespace ICA {
    //Forward declaration.
    template<class> class Placeholder;
    
    // template <class>  class FactorNode;  //
    // template<class U>
    // class FactorNode<U>; //ChildFactorWrapper;
    // template<template<class> class Model,  class U>
    // class CalcGaussianFactor;
    // template<class> class GaussianModel;
      
    //Forward declaration.
    namespace Details{
      template<template<class> class ,  class>
      class CalcGaussianFactor ;
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
      //Do the actual subcontext mutliplcation
      struct multiply_subcontext
      {
	typedef std::pair<placeholder_t const,data_t> Datum;
	typedef typename boost::call_traits<Datum&>::param_type data_parameter;
	multiply_subcontext(parameter C) : m_subcontext(C) {}
	void
	operator()(data_parameter d)
	{
	  d.second //the data
	    *= m_subcontext.Lookup( d.first ); // the variable
	}
      private:
	type m_subcontext;
      };
      //The map between the placeholders and the data
      std::map<placeholder_t, data_t> m_map;
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
     *  Builder builder;
     *  Builder::GaussianNode x = Builder.Gaussian(0.01,0.01);
     *  Builder::GaussianNode y = Builder.Gaussian(0.01,0.01);
     *  Builder::GaussianNode z = Builder.Gaussian(0.01,0.01);
     *  
     *  //Assign this context to the expression.
     *  Context context;
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
      
      typedef typename  DataContainer::const_iterator 
      const_iterator;

      placeholder_t
      Lookup(variable_parameter V) const
      {
	return
	  m_map.find(V)->second;
      };
      
      size_t 
      size() const {return m_map.size();}

      // void Assign(boost::shared_ptr<VariableNodeExp>&, const VariableNode<T>*  V );
      void 
      Assign(placeholder_parameter P, 
	     variable_parameter  V )
      {
	std::pair<typename DataContainer::iterator,bool> ret;

	ret = m_map.insert( Datum ( V,P) );
	if (ret.second == false) //already exists
	  {
	    ret.first->second = P;
	  }
      }

      subcontext_t
      operator[](size_parameter i) const
      {
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

      
      /** Obtain an iterator to the beginning of the map;
       */
      const_iterator
      begin() const {return m_map.begin();}
      
      
      /** Obtain an iterator to the end of the map;
       */
      const_iterator
      end() const {return m_map.end();}


      friend
      std::ostream&
      operator<<(std::ostream& out, const Context&  c)
      {
	typename DataContainer::const_iterator it;
	for(it = c.m_map.begin();
	    it!= c.m_map.end();
	    ++it)
	  {
	    out<<it->first<<" ";
	  }
	return out;
      }
    private:

      template<template<class> class Model,  class U>
      friend class Details::CalcGaussianFactor;

      template<template<class> class Model,  class U>
      void
      AddChildFactor(Details::CalcGaussianFactor<Model,U>* factor ) const
      {

	for( const_iterator it = begin();
	     it!= end();
	     ++it)
	  {
	    it->first->AddChildFactor(factor);
	  }
      }
      mutable DataContainer m_map;
    };

    template<class T>
    std::ostream&
    operator<<(std::ostream& out, 
	       const SubContext<T>& c)
    {
      //lock
      boost::mutex mutex;
      boost::lock_guard<boost::mutex> lock(mutex); 
  
      typename std::map<const Placeholder<T>*,  T>::const_iterator it;
      for(it = c.m_map.begin();
	  it!= c.m_map.end();
	  ++it)
	{
	  out<<it->second<<" ";
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
typename ICR::ICA::SubContext<T>::data_const_reference
ICR::ICA::SubContext<T>::Lookup(typename SubContext<T>::placeholder_parameter P) const
{
  return  m_map.find(P)->second;
};


template<class T>
void
ICR::ICA::SubContext<T>::Assign(typename SubContext<T>::placeholder_parameter P, 
				typename SubContext<T>::data_parameter V)
{
  typedef std::map<placeholder_t,data_t> DataContainer;
  typedef std::pair<placeholder_t,data_t> Datum;
      
  typename DataContainer::iterator it;
  std::pair<typename DataContainer::iterator,bool> ret;

  ret = m_map.insert( Datum(P, V) );
  if (ret.second == false) //already exists
    {
      ret.first->second = V;
    }
}


template<class T>
typename ICR::ICA::SubContext<T>::reference
ICR::ICA::SubContext<T>::operator*=(parameter C)
{
  typedef std::map<placeholder_t,data_t> DataContainer;
  typedef std::pair<placeholder_t,data_t> Datum;
  typename DataContainer::iterator it;

  PARALLEL_FOREACH(m_map.begin(),
  		   m_map.end(),
  		   multiply_subcontext(C));

  return *this;
}
