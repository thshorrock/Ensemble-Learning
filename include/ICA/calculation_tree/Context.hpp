#pragma once

#include "ICA/node/Node.hpp"

#include <map>
#include <string>

namespace ICR{

  namespace ICA {

    template<class> class Placeholder;

    template<class T>
    class SubContext{
      
    public:
      
      T 
      Lookup(const Placeholder<T>* P) const
      {
	return
	  m_map.find(P)->second;
      };
      
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
