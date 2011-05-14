#pragma once
//#include "NaturalParameters.hpp"
#include "MomentsIterator.hpp"
#include "ICA/node/Node.hpp"
#include "ICA/detail/parallel_algorithms.hpp"
#include<boost/assign/list_of.hpp>
#include <boost/assign/std/vector.hpp>
#include<iostream>

namespace ICR{
  namespace ICA{
    
    template<class> class Moments;
    
    template<class T=double>
    class 
    Moments
    {
    public:
      typedef  VariableNode<T>* Variable;
      typedef  VariableNode<T> const * ConstVariable;
      
      typedef typename boost::call_traits<std::vector<T> >::param_type 
      vector_parameter;
      
      typedef typename boost::call_traits<T>::param_type
      data_parameter;
      
      typedef typename boost::call_traits<T>::reference  
      data_reference;
      
      typedef typename boost::call_traits<T>::const_reference 
      data_const_reference;
      
      typedef typename boost::call_traits<T>::value_type
      data_type;
      
      
      typedef typename boost::call_traits<Moments<T> >::param_type 
      parameter;
      
      typedef typename boost::call_traits<Moments<T> >::reference  
      reference;
      
      typedef typename boost::call_traits<Moments<T> >::const_reference  
      const_reference;
      
      typedef typename boost::call_traits<Moments<T> >::value_type  
      type;
      
      typedef typename boost::call_traits<size_t>::param_type  
      size_parameter;
      
      typedef typename boost::call_traits<size_t>::value_type  
      size_type;
      
      typedef typename boost::call_traits<size_t>::reference  
      size_reference;
      
      typedef typename boost::call_traits<size_t>::const_reference  
      size_const_reference;
      
      //typedef typename std::vector<T>::iterator 
      typedef  MomentsIterator<type, data_type>  
      iterator;
      
      //typedef typename std::vector<T>::const_iterator 
      typedef  MomentsIterator<const type, const data_type>  
      const_iterator;
      
      // Moments();
      Moments(size_parameter size = 0);
      Moments(vector_parameter data);
      Moments(data_parameter  d1, data_parameter d2);
      Moments(parameter other);
      
      void operator=(parameter other);
      
      iterator
      begin();

      const_iterator
      begin() const;

      iterator
      end();

      const_iterator
      end() const;

      /** The number of Moments stored.*/
      size_type
      size() const {return m_data.size();}
      
      
      /** Access the moments data.*/
      data_const_reference 
      operator[](size_parameter index) const;
      
      /** Access the moments data.*/
      data_reference 
      operator[](size_parameter index);


      /** Add another Moments to this.
       * @param other The Moments to add.
       * @return A reference to the current Moments.
       */
      reference
      operator+=(parameter other);
      
      // Moments<T>& 
      // operator-=(const Moments<T>& other);

      reference
      operator*=(const data_parameter other);
      
      reference
      operator*=(const parameter other);



    private:
      
      
      struct plus
      {
	typename Moments<T>::data_type
	operator()(data_parameter i, 
		   data_parameter j) 
	{
	  return i+j;
	}
      };

      struct times
      {
      	typename Moments<T>::data_type
      	operator()(data_parameter i, 
      		   data_parameter j) 
      	{
      	  return i*j;
      	}
      };

      struct times_by{
	times_by(data_parameter t) : m_t(t) {};
  
	data_type
	operator()(data_parameter d)
	{ 
	  return d*m_t; 
	}
      private:
	T m_t;
      };


      std::vector<T> m_data;
      mutable boost::mutex m_mutex;

    };



    template <class T>
    inline
    std::ostream&
    operator<<( std::ostream& out, const Moments<T>& m)
    {
      for(size_t i=0;i<m.size();++i)
    	{
    	  out<<m[i]<<" ";
    	}
      return out;
    }

    template<class T>
    inline
    const Moments<T>
    operator+(const Moments<T>& a, const Moments<T>& b)
    {
      Moments<T> tmp = a;
      return tmp+=b;
    }  
    
    template<class T>  
    inline 
    const Moments<T>
    operator*(const Moments<T>& m,  
	      const T d)
    {
      Moments<T> tmp = m;
      return tmp*=d;
    }
      

    template<class T>  
    inline 
    const Moments<T>
    operator*(const T d, 
	      const Moments<T>& m)
    {
      Moments<T> tmp = m;
      return tmp*=d;
    }
    //typename Moments<T>::data_parameter d, 
    //typename Moments<T>::parameter m)
    
    // typedef  MomentsIterator< Moments<double>, double>       DoubleMomentsIterator ;
    // // typedef  typename MomentsIterator<const Moments<T>, T> const_iterator ;
    
  }
}


// template<class T>
// inline
// ICR::ICA::Moments<T>::Moments()
//   : m_data()
// {}




template<class T>
inline
ICR::ICA::Moments<T>::Moments( typename ICR::ICA::Moments<T>::vector_parameter data)
  : m_data(data),
    m_mutex()
{}

template<class T> 
inline   
ICR::ICA::Moments<T>::Moments( typename ICR::ICA::Moments<T>::size_parameter size)
  : m_data(std::vector<T>(size)),
    m_mutex() //non-copiable
{}

template<class T> 
inline   
ICR::ICA::Moments<T>::Moments(typename ICR::ICA::Moments<T>::data_parameter d1,
			      typename ICR::ICA::Moments<T>::data_parameter d2)
  : m_data(),
    m_mutex() //non-copiable
{
  using namespace boost::assign;
  m_data += d1,d2; // insert values at the end of the container
}
      
template<class T> 
inline   
ICR::ICA::Moments<T>::Moments(typename ICR::ICA::Moments<T>::parameter other)
  : m_data(other.m_data),
    m_mutex() //non-copiable
{}
      
template<class T> 
inline   
void 
ICR::ICA::Moments<T>::operator=(typename ICR::ICA::Moments<T>::parameter other)
{    
  if (this!= &other) {
    m_data = (other.m_data);
    // m_mutex non-copiable  
  }
}

template<class T>  
inline  
typename ICR::ICA::Moments<T>::iterator
ICR::ICA::Moments<T>::begin()
{
  return iterator(this,0);
}

template<class T>
inline
typename ICR::ICA::Moments<T>::const_iterator
ICR::ICA::Moments<T>::begin() const
{
  return const_iterator (this,0);
}

template<class T>
inline
typename ICR::ICA::Moments<T>::iterator
ICR::ICA::Moments<T>::end()
{
  return iterator(this,m_data.size());
}


template<class T>
inline
typename ICR::ICA::Moments<T>::const_iterator
ICR::ICA::Moments<T>::end() const
{
  return const_iterator(this,m_data.size());
}


template<class T>
inline
typename ICR::ICA::Moments<T>::data_const_reference
ICR::ICA::Moments<T>::operator[](typename ICR::ICA::Moments<T>::size_parameter index) const
{
  boost::lock_guard<boost::mutex> lock(m_mutex); 
  return m_data[index];
}
      
template<class T>
inline
typename ICR::ICA::Moments<T>::data_reference
ICR::ICA::Moments<T>::operator[](typename ICR::ICA::Moments<T>::size_parameter index)
{
  boost::lock_guard<boost::mutex> lock(m_mutex); 
  return m_data[index];
}
      

template<class T>
typename ICR::ICA::Moments<T>::reference  
ICR::ICA::Moments<T>::operator+=(typename Moments<T>::parameter other)
{
  //This could be called by different threads, so lock
  boost::lock_guard<boost::mutex> lock(m_mutex);  
  PARALLEL_TRANSFORM(m_data.begin(), m_data.end(), other.begin(), m_data.begin(), plus());
  return *this;
}
      
template<class T>  
inline 
typename ICR::ICA::Moments<T>::reference
ICR::ICA::Moments<T>::operator*=(typename Moments<T>::parameter other)
{
  //This could be called by different threads, so lock
  boost::lock_guard<boost::mutex> lock(m_mutex); 
  PARALLEL_TRANSFORM(m_data.begin(), m_data.end(), other.begin(), m_data.begin(), times() );
  return *this;
}

template<class T>  
inline 
typename ICR::ICA::Moments<T>::reference
ICR::ICA::Moments<T>::operator*=(typename Moments<T>::data_parameter d)
{
  //This could be called by different threads, so lock
  boost::lock_guard<boost::mutex> lock(m_mutex); 
  PARALLEL_TRANSFORM(m_data.begin(), m_data.end(),  m_data.begin(), times_by(d));
  return *this;
}
      
   


