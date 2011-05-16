#pragma once

#include "ICA/node/Node.hpp"

#include <map>
#include <string>

namespace ICR{

  namespace ICA {
    //Forward declaration.
    template<class> class Placeholder;
    
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
      
      typedef typename boost::call_traits<T>::reference  
      data_reference;
      
      typedef typename boost::call_traits<T>::const_reference 
      data_const_reference;
      
      typedef typename boost::call_traits<const Placeholder<T>*>::param_type
      placeholder_parameter;
      
      /** Obtain the value associated with for the placeholder.
       *   @param P The placeholder.
       *   @return A value taken from the moment of the Variable node represented by P.
       */
      data_const_reference
      Lookup(placeholder_parameter P) const;
      
      void Assign(const Placeholder<T>* Expr, const T V )
      {

	typename std::map<const Placeholder<T>*,  T>::iterator it;
	std::pair<typename std::map<const Placeholder<T>*, T>::iterator,bool> ret;

	ret = m_map.insert( std::pair<const Placeholder<T>*, T> (Expr, V) );
	if (ret.second == false) //already exists
	  {
	    ret.first->second = V;
	  }

      }
      

      SubContext&
      operator*=(const SubContext& C)
      {
	
	typename std::map<const Placeholder<T>*,  T>::iterator it;
	for(it = m_map.begin();
	    it!= m_map.end();
	    ++it)
	  {
	    it->second*=C.Lookup(it->first);
	  }
	return *this;
	
      }

      
      friend
      std::ostream&
      operator<<(std::ostream& out, const SubContext&  c)
      {
	typename std::map<const Placeholder<T>*,  T>::const_iterator it;
	for(it = c.m_map.begin();
	    it!= c.m_map.end();
	    ++it)
	  {
	    out<<it->second<<" ";
	  }
	return out;
      }

      friend 
      SubContext<T>
      operator*(const SubContext<T>& A, const SubContext<T>& B)
      {
	SubContext<T> tmp = A;
	return tmp*=B;
      }
    
    private:
      std::map<const Placeholder<T>*, T> m_map;
    };

    template<class T>
    class Context{
      typedef std::map<const VariableNode<T>*,const Placeholder<T>*> DataContainer;
      typedef std::pair< const VariableNode<T>*,const Placeholder<T>*> Datum;
      
    public:
      
      typedef typename  std::map<const VariableNode<T>*,const Placeholder<T>*>::const_iterator const_iterator;
       const_iterator
      begin() const {return m_map.begin();}
      
       const_iterator
      end() const {return m_map.end();}

      const Placeholder<T>* Lookup(const VariableNode<T>* V) const
      {
	return
	  m_map.find(V)->second;
      };
      
      size_t size() const {return m_map.size();}

      // void Assign(boost::shared_ptr<VariableNodeExp>&, const VariableNode<T>*  V );
      void Assign(const Placeholder<T>* Expr, const VariableNode<T>*  V )
      {
	std::pair<typename DataContainer::iterator,bool> ret;

	ret = m_map.insert( Datum ( V,Expr) );
	if (ret.second == false) //already exists
	  {
	    ret.first->second = Expr;
	  }
      }

      SubContext<T>
      operator[](const size_t& i) const
      {
	SubContext<T> c;
	for(typename DataContainer::const_iterator it = m_map.begin();
	    it != m_map.end();
	    ++it)
	  {
	    const VariableNode<T>* V = it->first;
	    c.Assign( it->second,V->GetMoments()[i]);
	  }
	      

	return c;
      }

      
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
      DataContainer m_map;
    };

    

  }
}


template<class T>
typename ICR::ICA::SubContext<T>::data_const_reference
ICR::ICA::SubContext<T>::Lookup(typename SubContext<T>::placeholder_parameter P) const
{
  return  m_map.find(P)->second;
};
