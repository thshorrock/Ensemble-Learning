#pragma once

#include "VariableExp.hpp"

#include <map>
#include <string>

namespace ICR{

  namespace ICA {

    class Placeholder;

    template<class T>
    class SubContext{
      
    public:
      
      T 
      lookup(const Placeholder* P) const
      {
	return
	  m_map.find(P)->second;
      };
      
      void assign(const Placeholder* Expr, const T V )
      {

	std::map<const Placeholder*,  T>::iterator it;
	std::pair<std::map<const Placeholder*, T>::iterator,bool> ret;

	ret = m_map.insert( std::pair<const Placeholder*, T> (Expr, V) );
	if (ret.second == false) //already exists
	  {
	    ret.first->second = V;
	  }

      }
      
      SubContext&
      operator*=(const SubContext& C)
      {
	
	std::map<const Placeholder*,  T>::iterator it;
	for(it = m_map.begin();
	    it!= m_map.end();
	    ++it)
	  {
	    it->second*=C.lookup(it->first);
	  }

	
      }

      friend 
      SubContext 
      operator*(const SubContext& A, const SubContext& B)
      {
	SubContext tmp = A;
	return A*=B;
      }
    
    private:
      std::map<const Placeholder*, T> m_map;
    };

    template<class T>
    class Context{
      
    public:
      //T lookup(const std::size_t& id) const;
      // void assign(boost::shared_ptr<VariableExp>&, const Variable<T>*  V );
      void assign(const Placeholder* Expr, const Variable<T>*  V )
      {
	std::map<const Placeholder*, const Variable<T>*>::iterator it;
	std::pair<std::map<const Placeholder*, const Variable<T>*>::iterator,bool> ret;

	ret = m_map.insert( std::pair<const Placeholder*, const Variable<T>*> (Expr, V) );
	if (ret.second == false) //already exists
	  {
	    ret.first->second = V;
	  }
      }

      SubContext<T>
      operator[](const size_t& i)
      {
	SubContext c;
	for_each(m_map.begin(), m_map.end(), 
		  boost::bind(&SubContext::assign, boost::ref(c),
			      _1.first, (_1.second)->GetMoments[i])
		 );
	return c;
      }

    private:
      std::map<const Placeholder*,const Variable<T>* > m_map;
    };

  }
}
