#pragma once
#include "rng.hpp"


namespace ICR{
  namespace ICA{

    using ICR::maths::rng;

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
      static rng* Instance()
      {
	if (m_rng == 0)
	  m_rng = new rng(100);
	m_Destroyer.SetSingleton(m_rng); 
	return m_rng;
      }
      
      static rng* Restart() 
      {
	if (m_rng!=0) {
	  delete m_rng;
	  m_rng = new rng(100);
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
  delete m_singleton; 
}

template<class T>
void
ICR::ICA::SingletonDestroyer<T>::SetSingleton(T *s) 
{
  m_singleton = s;
}

    
