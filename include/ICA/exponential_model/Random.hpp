#pragma once
//#include "rng.hpp"

#include <gsl/gsl_rng.h>     
#include <gsl/gsl_randist.h>
#include <ctime>

namespace ICR{
  namespace ICA{

    
  class rng{
  private:
    unsigned long int m_seed;
    gsl_rng* m_rng;
    rng(const rng& other); //non-copyable
    rng& operator= (const rng & other); //non-copyable
    
  public:
    rng()
      : m_seed(static_cast<unsigned int>(std::time(0))),
	m_rng(gsl_rng_alloc (gsl_rng_taus))
    {
      gsl_rng_set (m_rng,m_seed);
    }
    rng(unsigned long int seed)
      : m_seed(seed),
	m_rng(gsl_rng_alloc (gsl_rng_taus))
    {
       gsl_rng_set (m_rng,m_seed);
    }
    
    unsigned long int 
    get_seed() const {return m_seed;};
    
    ~rng(){gsl_rng_free (m_rng);};

    double uniform() {return gsl_rng_uniform (m_rng);};
    
    double uniform(double a, double b) {return gsl_ran_flat (m_rng, a, b);}
    
    double gaussian(const double sigma = 1, const double mean = 0){
      return gsl_ran_gaussian (m_rng, sigma) + mean;
    }
    double gaussian_tail(const double sigma = 1, const double mean = 0, const double min = 0){
      return gsl_ran_gaussian_tail (m_rng, min-mean, sigma) + mean;
    } 

    double exponential(const double mean = 1){
      return gsl_ran_exponential (m_rng, mean );
    }
    double gamma(const double shape =1, const double scale =1){
      return gsl_ran_gamma (m_rng, shape, scale);
    } 
    void dirichlet(const size_t size, const double* alpha, double* theta){

      gsl_ran_dirichlet(m_rng,size, alpha, theta);
    } 
  };
  


    //  using ICR::maths::rng;

    template<class T>
    class SingletonDestroyer { 
    public: 
      SingletonDestroyer(T* = 0); 
      ~SingletonDestroyer(); 
      void SetSingleton(T* s); 
    private: 
      T* m_singleton;
    };

    class Random{
      
    public:
      static rng* Instance(unsigned int seed)
      {
	if (m_rng == 0)
	  m_rng = new rng(100);
	m_Destroyer.SetSingleton(m_rng); 
	return m_rng;
      }
      static rng* Instance()
      {
	if (m_rng == 0)
	  m_rng = new rng();
	m_Destroyer.SetSingleton(m_rng); 
	return m_rng;
      }
      
      
      static rng* Restart() 
      {
	if (m_rng!=0) {
	  delete m_rng;
	  m_rng = new rng();
	  m_Destroyer.SetSingleton(m_rng); 
	}
	return m_rng;
      }; 
    private:
      Random(); 
      ~Random(){}

      friend class SingletonDestroyer<rng>; 
      static rng *m_rng;
      static SingletonDestroyer<rng> m_Destroyer;
    };


  }
}

//initialise
ICR::ICA::rng* ICR::ICA::Random::m_rng = 0;

ICR::ICA::SingletonDestroyer<ICR::ICA::rng> 
ICR::ICA::Random::m_Destroyer(ICR::ICA::Random::m_rng);




template<class T>
ICR::ICA::SingletonDestroyer<T>::SingletonDestroyer(T* s) 
{ 
  m_singleton = s; 
}
template<class T>
ICR::ICA::SingletonDestroyer<T>::~SingletonDestroyer() 
{
  if (m_singleton!=0) {
    delete m_singleton;
    m_singleton = 0;
  } 
}

template<class T>
void
ICR::ICA::SingletonDestroyer<T>::SetSingleton(T *s) 
{
  m_singleton = s;
}

    
