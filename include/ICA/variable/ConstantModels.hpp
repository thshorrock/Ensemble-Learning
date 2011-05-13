#pragma once
#include "ICA/Node.hpp"
#include "ICA/Moments.hpp"

//No variable data - no updating after initialisation these functions are thread safe.

namespace ICR{
  namespace ICA{
    
    /** Normal Constant */
    template <class T = double>
    class
    NormalConstant
    {
    public:
      NormalConstant(const T d)
	: m_Moments(d, 0)
      {}
	
      static
      T
      GetDataPenalty(const T x)
      {
	return 0;
      }
	       
      const Moments<T>&
      GetMoments() const
      {
	return m_Moments;
      }
	
    protected:
      const Moments<T> m_Moments;
    };

    /** Gaussian Constant */
    template <class T = double>
    class
    GaussianConstant
    {
    public:
      static
      T
      GetDataPenalty(const T x)
      {
	return 0;
      }

      GaussianConstant(const T d)
	: m_Moments(d, d*d)
      {}
		       
      const Moments<T>&
      GetMoments() const
      {
	return m_Moments;
      }
	
    protected:
      const Moments<T> m_Moments;
    };
    
    /** Rectified Gaussian Constant */
    template <class T = double>
    class
    RectifiedGaussianConstant
    {
    public:
      RectifiedGaussianConstant(const T d)
	: m_Moments(d, d*d)
      {}
	 static
      T
      GetDataPenalty(const T x)
      {
	if (x<0) return (-1.0/0);
	else return 0;
      }	       
      const Moments<T>&
      GetMoments() const
      {
	return m_Moments;
      }
	
    protected:
      const Moments<T> m_Moments;
    };
    /** Gamma Constant */
    template <class T = double>
    class
    GammaConstant
    {
    public:
      GammaConstant(const T d)
	: m_Moments(d, std::log(d))
      {
	BOOST_ASSERT(d>0);
      }
		       
      const Moments<T>&
      GetMoments() const
      {
	return m_Moments;
      }
	
    protected:
      const Moments<T> m_Moments;
    };


    
    /** Gamma Constant */
    template <class T = double>
    class
    DirichletConstant
    {
    public:
      DirichletConstant(const size_t s, const T d)
	: m_Moments(std::vector<T>(s,d))
      {}
	
      static
      T
      GetDataPenalty(const T x)
      {
	return 0;
      }	       
      const Moments<T>&
      GetMoments() const
      {
	return m_Moments;
      }
	
    protected:
      const Moments<T> m_Moments;
    };




  }
}

