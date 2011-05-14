#pragma once
#include <iostream>

#include "Moments.hpp"

#include<vector>


#include "ICA/detail/parallel_algorithms.hpp"

namespace ICR{
  namespace ICA{

    template<class T>
    class NaturalParameters
    {
    public:

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

      typedef typename boost::call_traits<NaturalParameters<T> >::param_type 
      parameter;
      
      typedef typename boost::call_traits<NaturalParameters<T> >::reference  
      reference;
      
      typedef typename boost::call_traits<NaturalParameters<T> >::const_reference  
      const_reference;
      
      typedef typename boost::call_traits<NaturalParameters<T> >::value_type  
      type;
      
      typedef typename boost::call_traits<size_t>::param_type  
      size_parameter;
      
      typedef typename boost::call_traits<size_t>::value_type  
      size_type;
      
   
      typedef typename std::vector<T>::iterator
      iterator;
      
      typedef typename std::vector<T>::const_iterator
      const_iterator;

      /** Constructor */
      NaturalParameters(size_parameter size = 0);
      
      NaturalParameters(vector_parameter NP);
      // NaturalParameters(const Moments<T>& m);

      NaturalParameters(data_parameter d1, 
			data_parameter d2);

      data_const_reference 
      operator[](size_parameter) const;
      
      data_reference
      operator[](size_parameter) ;

      iterator
      begin();

      const_iterator
      begin() const;

      iterator
      end();

      const_iterator
      end() const;

      size_type
      size() const {return m_data.size();}
      
      reference
      operator+=(parameter other);
      
      reference
      operator-=(parameter other);

      reference
      operator*=(data_parameter other);
      
      template<class U>
      friend
      U
      operator*(const NaturalParameters<U>&  a, 
		const Moments<U>& b);
    private:

      struct plus
      {
	data_type
	operator()(data_parameter i, 
		   data_parameter j) 
	{
	  return i+j;
	}
      };

      struct minus
      {
	data_type
	operator()(data_parameter i, 
		   data_parameter j) 
	{
	  return i-j;
	}
      };

      struct times
      {
      	data_type
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
	data_type m_t;
      };

      std::vector<data_type> m_data;
      

    };


    template<class T>
    inline
    NaturalParameters<T>
    operator+(const NaturalParameters<T>& a,
	      const NaturalParameters<T>& b)
    {
      NaturalParameters<T> tmp = a;
      return tmp+=b;
    }  
    
    template<class T>
    inline
    NaturalParameters<T>
    operator-(const NaturalParameters<T>& a,
	      const NaturalParameters<T>& b)
    {
      ICR::ICA::NaturalParameters<T> tmp = a;
      return tmp-=b;
    }  

    
 
    template<class T>
    inline
    T
    operator*(const NaturalParameters<T>&  a, 
	      const Moments<T>& b)
    { 
      
      std::vector<T> prod(a.size());
      PARALLEL_TRANSFORM(a.begin(), a.end(), b.begin(), prod.begin(), 
			 typename NaturalParameters<T>::times());
      T sum = 0.0;
      return PARALLEL_ACCUMULATE(prod.begin(), prod.end(), sum);
    }  
      


    template<class T>
    inline
    T
    operator*(const Moments<T>& m,const NaturalParameters<T>&  p)
    {
      return p*m;
    }  


    template<class T>
    inline
    const NaturalParameters<T>
    operator*(const NaturalParameters<T>& n,const T  d)
    {
      NaturalParameters<T> tmp(n);
      return tmp*=d;
    }  


    // template<class T>
    // inline
    // typename NaturalParameters<T>::data_type
    // operator*(typename NaturalParameters<T>::data_parameter a,
    // 	      typename NaturalParameters<T>::data_parameter b)
    // {
    //   NaturalParameters<T> tmp = a;
    //   return tmp*=b;
    // }  

  }
}

template<class T>
inline
ICR::ICA::NaturalParameters<T>::NaturalParameters(typename NaturalParameters<T>::size_parameter size)
  : //m_mutex_ptr(new boost::mutex), 
  m_data(size)
{
}

template<class T>
inline
ICR::ICA::NaturalParameters<T>::NaturalParameters(typename NaturalParameters<T>::vector_parameter NP)
  : //m_mutex_ptr(new boost::mutex), 
    m_data(NP)
{
}

// template<class T>
// inline
// ICR::ICA::NaturalParameters<T>::NaturalParameters(const ICR::ICA::Moments<T>& m)
// // : m_mutex_ptr(new boost::mutex)
//   : m_data()
// {
//   m_data.resize(m.size());
//   //boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  
//   PARALLEL_COPY(m.begin(),m.end(),m_data.begin());
// }
 
template<class T>
inline
ICR::ICA::NaturalParameters<T>::NaturalParameters(typename NaturalParameters<T>::data_parameter d1,
						  typename NaturalParameters<T>::data_parameter d2 )
//  : 
  //  m_mutex_ptr(new boost::mutex)
  : m_data(2)
{
  //boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  
  m_data[0] = d1;
  m_data[1] = d2;
}

template<class T>
typename ICR::ICA::NaturalParameters<T>::iterator
ICR::ICA::NaturalParameters<T>::begin()
{
  //  boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  
  return m_data.begin();
}


template<class T>
typename ICR::ICA::NaturalParameters<T>::const_iterator
ICR::ICA::NaturalParameters<T>::begin() const
{
  //  boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  
  return m_data.begin();
}


template<class T>
typename ICR::ICA::NaturalParameters<T>::iterator
ICR::ICA::NaturalParameters<T>::end()
{
  // boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  
  return m_data.end();
}


template<class T>
typename ICR::ICA::NaturalParameters<T>::const_iterator
ICR::ICA::NaturalParameters<T>::end() const
{
  // boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  
  return m_data.end();
}


 

template<class T>
inline
typename ICR::ICA::NaturalParameters<T>::data_const_reference 
ICR::ICA::NaturalParameters<T>::operator[](typename NaturalParameters<T>::size_parameter i) const
{
  //boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  
  return m_data[i];
}
   
template<class T>
inline
typename ICR::ICA::NaturalParameters<T>::data_reference
ICR::ICA::NaturalParameters<T>::operator[](typename NaturalParameters<T>::size_parameter i) 
{
  //boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  
  return m_data[i];
}
   


template<class T>
inline
typename ICR::ICA::NaturalParameters<T>::reference
ICR::ICA::NaturalParameters<T>::operator+=(typename ICR::ICA::NaturalParameters<T>::parameter other)
{
  //This could be called by different threads, so lock
  //  boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  //DO KEEP THIS ONE!
  //  std::cout<<"size = "<<m_data.size()<<"other = "<<m_data.size<<std::endl;

  PARALLEL_TRANSFORM(m_data.begin(), m_data.end(), other.begin(), m_data.begin(), plus());
  // for(size_t i=0;i<other.size();++i){
  //   m_data[i]+=other[i];
  // }

  return *this;
}
   

template<class T>  
inline 
typename ICR::ICA::NaturalParameters<T>::reference
ICR::ICA::NaturalParameters<T>::operator-=(typename ICR::ICA::NaturalParameters<T>::parameter other)
{
  //This could be called by different threads, so lock
  //  boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  //DO KEEP THIS ONE!
  PARALLEL_TRANSFORM(m_data.begin(), m_data.end(), other.begin(), m_data.begin(), minus());
  return *this;
}
  
template<class T>
inline
typename ICR::ICA::NaturalParameters<T>::reference
ICR::ICA::NaturalParameters<T>::operator*=(typename ICR::ICA::NaturalParameters<T>::data_parameter other)
{
  
  PARALLEL_TRANSFORM(m_data.begin(), m_data.end(), m_data.begin(), times_by(other));
  return *this;
}    


// template<class T>
// inline
// ICR::ICA::NaturalParameters<T>
// operator+(const ICR::ICA::NaturalParameters<T>& a, const ICR::ICA::NaturalParameters<T>& b)
// {
//   ICR::ICA::NaturalParameters<T> tmp = a;
//   return tmp+=b;
// }  




template<class T>
inline
std::ostream&
operator<<(std::ostream& out, typename ICR::ICA::NaturalParameters<T>::parameter NP)
{
  for(size_t i=0;i<NP.size();++i)
    {
      out<<NP[i]<<" ";
    }

  return out;
}
    
