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
      typedef typename std::vector<T>::iterator iterator;

      NaturalParameters(const size_t size = 2);
      NaturalParameters(const std::vector<T>& NP);
      NaturalParameters(const Moments<T>& m);

      NaturalParameters(const T& d1, const T& d2 );
      ~NaturalParameters()
      {
	//std::cout<<"~NP"<<std::endl;

      }

      T operator[](const short& i) const;
      T& operator[](const short& i) ;

      typename std::vector<T>::iterator
      begin();

      typename std::vector<T>::const_iterator
      begin() const;

      typename std::vector<T>::iterator
      end();

      typename std::vector<T>::const_iterator
      end() const;

      size_t
      size() const {return m_data.size();}
      
      // bool
      // operator<(const NaturalParameters& other) const;
      NaturalParameters<T>& 
      operator+=(const NaturalParameters<T>& other);
      
      NaturalParameters<T>& 
      operator-=(const NaturalParameters<T>& other);

      NaturalParameters<T>& 
      operator*=(const T& other);
      
    private:
      //      mutable boost::shared_ptr<boost::mutex> m_mutex_ptr;
      std::vector<T> m_data;
      

    };


    template<class T>
    inline
    ::ICR::ICA::NaturalParameters<T>
    operator+(const ::ICR::ICA::NaturalParameters<T>& a, const ::ICR::ICA::NaturalParameters<T>& b)
    {
      ::ICR::ICA::NaturalParameters<T> tmp = a;
      return tmp+=b;
    }  
    
    template<class T>
    inline
    T
    operator*(const ICR::ICA::NaturalParameters<T>& a, const ICR::ICA::Moments<T>& b)
    {
      std::vector<T> prod(a.size());
      PARALLEL_TRANSFORM(a.begin(), a.end(), b.begin(), prod.begin(), times<T>);
      T sum = 0;
      return PARALLEL_ACCUMULATE(prod.begin(), prod.end(), sum);
    }  



    template<class T>
    inline
    T
    operator*(const ICR::ICA::Moments<T>& a, const ICR::ICA::NaturalParameters<T>& b)
    {
      return b*a;
    }  


    template<class T>
    inline
    ::ICR::ICA::NaturalParameters<T>
    operator*(const ::ICR::ICA::NaturalParameters<T>& a,const T& b)
    {
      ::ICR::ICA::NaturalParameters<T> tmp = a;
      return tmp*=b;
    }  

  }
}

template<class T>
inline
ICR::ICA::NaturalParameters<T>::NaturalParameters(const size_t size)
  : //m_mutex_ptr(new boost::mutex), 
    m_data(size)
{
  for(size_t i=0;i<size;++i){
    m_data[i] = 0;
  }

}

template<class T>
inline
ICR::ICA::NaturalParameters<T>::NaturalParameters(const std::vector<T>& NP)
  : //m_mutex_ptr(new boost::mutex), 
    m_data(NP)
{
}

template<class T>
inline
ICR::ICA::NaturalParameters<T>::NaturalParameters(const ICR::ICA::Moments<T>& m)
// : m_mutex_ptr(new boost::mutex)
  : m_data()
{
  m_data.resize(m.size());
  //boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  
  PARALLEL_COPY(m.begin(),m.end(),m_data.begin());
}
 
template<class T>
inline
ICR::ICA::NaturalParameters<T>::NaturalParameters(const T& d1 , const T& d2 )
//  : 
  //  m_mutex_ptr(new boost::mutex)
  : m_data()
{
  //boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  
  m_data.resize(2);
  m_data[0] = d1;
  m_data[1] = d2;
}

template<class T>
typename std::vector<T>::iterator
ICR::ICA::NaturalParameters<T>::begin()
{
  //  boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  
  return m_data.begin();
}


template<class T>
typename std::vector<T>::const_iterator
ICR::ICA::NaturalParameters<T>::begin() const
{
  //  boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  
  return m_data.begin();
}


template<class T>
typename std::vector<T>::iterator
ICR::ICA::NaturalParameters<T>::end()
{
  // boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  
  return m_data.end();
}


template<class T>
typename std::vector<T>::const_iterator
ICR::ICA::NaturalParameters<T>::end() const
{
  // boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  
  return m_data.end();
}


 

template<class T>
inline
T 
ICR::ICA::NaturalParameters<T>::operator[](const short& i) const
{
  //boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  
  return m_data[i];
}
   
template<class T>
inline
T& 
ICR::ICA::NaturalParameters<T>::operator[](const short& i) 
{
  //boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  
  return m_data[i];
}
   
// namespace {

// template<class T>
//   T plus(T i, T j) 
//   { 
//     return i+j; 
//   }

// template<class T>
//   T minus(T i, T j) 
//   { 
//     return i-j; 
//   }

// template<class T>
//   T times(T i, T j) 
//   { 
//     return i*j; 
//   }
// }


template<class T>
ICR::ICA::NaturalParameters<T>& 
ICR::ICA::NaturalParameters<T>::operator+=(const NaturalParameters<T>& other)
{
  //This could be called by different threads, so lock
  //  boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  //DO KEEP THIS ONE!
  //  std::cout<<"size = "<<m_data.size()<<"other = "<<m_data.size<<std::endl;

 PARALLEL_TRANSFORM(m_data.begin(), m_data.end(), other.begin(), m_data.begin(), plus<T>);
  // for(size_t i=0;i<other.size();++i){
  //   m_data[i]+=other[i];
  // }

  return *this;
}
   

template<class T>  
inline 
ICR::ICA::NaturalParameters<T>& 
ICR::ICA::NaturalParameters<T>::operator-=(const NaturalParameters<T>& other)
{
  //This could be called by different threads, so lock
  //  boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  //DO KEEP THIS ONE!
  PARALLEL_TRANSFORM(m_data.begin(), m_data.end(), other.begin(), m_data.begin(), minus<T>);
  return *this;
}
  
template<class T>
ICR::ICA::NaturalParameters<T>& 
ICR::ICA::NaturalParameters<T>::operator*=(const T& other)
{
  PARALLEL_TRANSFORM(m_data.begin(), m_data.end(), m_data.begin(), times_by<T>(other));
 
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
ICR::ICA::NaturalParameters<T>
operator-(const ICR::ICA::NaturalParameters<T>& a, const ICR::ICA::NaturalParameters<T>& b)
{
  ICR::ICA::NaturalParameters<T> tmp = a;
  return tmp-=b;
}  



template<class T>
inline
std::ostream&
operator<<(std::ostream& out, const ICR::ICA::NaturalParameters<T>& NP)
{
  for(size_t i=0;i<NP.size();++i)
    {
      out<<NP[i]<<" ";
    }

  return out;
}
    
