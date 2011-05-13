#pragma once
//#include "NaturalParameters.hpp"
#include "ICA/node/Node.hpp"
#include "ICA/detail/parallel_algorithms.hpp"

#include<iostream>

namespace ICR{
  namespace ICA{
    
    template<class T=double>
    class 
    Moments
    {
    public:
      typedef  VariableNode<T>* Variable;
      typedef  VariableNode<T> const * ConstVariable;
      
      // Moments();
      Moments(const size_t size = 0);
      Moments(const std::vector<T> data);

      Moments(const T& d1, const T& d2);
      
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
      
      
      T operator[](const short& index) const;
      T& operator[](const short& index);

      // ConstVariable
      // source() const;
      
      Moments<T>& 
      operator+=(const Moments<T>& other);
      
      Moments<T>& 
      operator-=(const Moments<T>& other);

      Moments<T>& 
      operator*=(const double other);

    private:
      std::vector<T> m_data;
      //      mutable boost::shared_ptr<boost::mutex> m_mutex_ptr;
      // ConstVariable m_source;
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
    ::ICR::ICA::Moments<T>
    operator+(const ::ICR::ICA::Moments<T>& a, const ::ICR::ICA::Moments<T>& b)
    {
      ::ICR::ICA::Moments<T> tmp = a;
      return tmp+=b;
    }  


    
  }
}


// template<class T>
// inline
// ICR::ICA::Moments<T>::Moments()
//   : m_data()
// {}

template<class T>
inline
ICR::ICA::Moments<T>::Moments(const std::vector<T> data)
  : 
  //m_mutex_ptr(new boost::mutex), 
  m_data(data)
{
  
}

template<class T> 
inline   
ICR::ICA::Moments<T>::Moments(const size_t size)
//  : m_mutex_ptr(new boost::mutex)
{
  m_data.resize(size);
  
  for(size_t i=0;i<m_data.size();++i){
    m_data[i] = 1;
  }
}
template<class T> 
inline   
ICR::ICA::Moments<T>::Moments(const T& d1, const T& d2)
//  : m_mutex_ptr(new boost::mutex)
{
  m_data.resize(2);
  m_data[0] = d1;
  m_data[1] = d2;
}
      
template<class T>  
inline  
typename std::vector<T>::iterator
ICR::ICA::Moments<T>::begin()
{
  return m_data.begin();
}

template<class T>
inline
typename 
std::vector<T>::const_iterator
ICR::ICA::Moments<T>::begin() const
{
  return m_data.begin();
}

template<class T>
inline
typename 
std::vector<T>::iterator
ICR::ICA::Moments<T>::end()
{
  return m_data.end();
}


template<class T>
inline
typename std::vector<T>::const_iterator
ICR::ICA::Moments<T>::end() const
{
  return m_data.end();
}


template<class T>
inline
T 
ICR::ICA::Moments<T>::operator[](const short& index) const
{
  return m_data[index];
}
      
template<class T>
inline
T& 
ICR::ICA::Moments<T>::operator[](const short& index)
{
  return m_data[index];
}
      

namespace {

  template<class T>
  inline
  T plus(T i, T j) 
  { 
    return i+j; 
  }

  template<class T>
  inline
  T minus(T i, T j) 
  { 
    return i-j; 
  }

  template<class T>
  inline
  T times(T i, T j) 
  { 
    return i*j; 
  }

  template<class T>
  struct times_by{
    times_by(T t) : m_t(t) {};
  
    double 
    operator()(const double d)
    { 
      return d*m_t; 
    }
    T m_t;

  };
}


template<class T>
ICR::ICA::Moments<T>& 
ICR::ICA::Moments<T>::operator+=(const Moments<T>& other)
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
ICR::ICA::Moments<T>& 
ICR::ICA::Moments<T>::operator-=(const Moments<T>& other)
{
  //This could be called by different threads, so lock
  //  boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  //DO KEEP THIS ONE!
  PARALLEL_TRANSFORM(m_data.begin(), m_data.end(), other.begin(), m_data.begin(), minus<T>);
  return *this;
}
      


template<class T>  
inline 
ICR::ICA::Moments<T>& 
ICR::ICA::Moments<T>::operator*=(const double d)
{
  //This could be called by different threads, so lock
  //  boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  //DO KEEP THIS ONE!
  //  std::cout<<"m_data.size = "<<m_data.size()<<std::endl;

  for(size_t i=0;i<m_data.size();++i){
    m_data[i]*=d;
  }

  return *this;
}
      

template<class T>  
inline 
ICR::ICA::Moments<T> 
operator*(const ICR::ICA::Moments<T>& m, const double d)
{
  ICR::ICA::Moments<T> tmp = m;
  return tmp*=d;
}
      

template<class T>  
inline 
ICR::ICA::Moments<T> 
operator*( const double d, const ICR::ICA::Moments<T>& m)
{
  ICR::ICA::Moments<T> tmp = m;
  return tmp*=d;
}
      
