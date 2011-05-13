#pragma once

#include "Random.hpp"
#include "ExponentialModel.hpp"

#include <gsl/gsl_sf_psi.h> //for Digamma function gsl_sf_psi
#include <gsl/gsl_sf_gamma.h> //for gamma function

#include <vector>
#include "ICA/parallel_algorithms.hpp"

#include <boost/assert.hpp> 

namespace ICR{
  namespace ICA{

    template<class T>
    class DirichletModel // : public ExponentialModel<T>
    {
    public:

      //moments from parents (the u values)
      static
      T
      CalcLogNorm(const Moments<T>& Us) ;

      static
      T
      CalcLogNorm(const NaturalParameters<T>& NP) ;
      
      static
      T
      CalcAvLog(const Moments<T>& values,
		const Moments<T>& Data) ;

      
      struct take_log
      {
	double operator()(const double d){ return std::log(d);}
      };

      static
      Moments<T>
      CalcSample(const VariableNode<T>* V) 
      {
	ICR::maths::rng* random = Random::Instance();
	Moments<T> Mus = V->GetMoments();
	const size_t     size   = Mus.size();

	double u[size];
	double x[size];
	PARALLEL_COPY(Mus.begin(),Mus.end(),u);
	random->dirichlet(size, u, x);
	
	
	std::vector<double> M(size);
	PARALLEL_TRANSFORM(x,x+size,M.begin(),take_log());
	return Moments<T>(M);
      }


      static
      Moments<T>
      CalcMoments(const NaturalParameters<T>& NP);

      static
      NaturalParameters<T>
      CalcNP2Data(const Moments<T>& Us);


    private:
      static
      T
      CalcLogNorm(const std::vector<T>& u) ;
    };

  }
}

namespace {
  
  struct plus_one
  {
    double operator()(const double d) {return d+1.0;}
  };

  struct minus_one
  {
    double operator()(const double d) {return d-1.0;}
  };

  class
  calculate_digamma_minus_digamma
  {
  public:
    calculate_digamma_minus_digamma(const double U)
      : m_subtract_me(gsl_sf_psi(U))
    {}
    double operator()(const double u) 
    {
      return gsl_sf_psi(u) - m_subtract_me;
    }
  private:
    double m_subtract_me;
  };
  
}

template<class T>
inline
T 
ICR::ICA::DirichletModel<T>::CalcLogNorm(const std::vector<T>& u)  
{
  // std::cout<<"DIRICHLET u = "<<Moments<T>(u)<<std::endl;

  const T U = PARALLEL_ACCUMULATE(u.begin(), u.end(), 0.0);
    
  //calculate Gamma(m_value[i]) for each i
  std::vector<T> Gamma_u(u.size());
  PARALLEL_TRANSFORM(u.begin(), u.end(), Gamma_u.begin(), gsl_sf_lngamma);
    
  const T SumLnGamma_u = PARALLEL_ACCUMULATE(Gamma_u.begin(), Gamma_u.end(), 0.0);
    
  return  gsl_sf_lngamma(U) - SumLnGamma_u;
}


template<class T>
inline
T 
ICR::ICA::DirichletModel<T>::CalcLogNorm(const Moments<T>& Us)  
{
  // std::cout<<"DIRICHLET Mu = "<<Us<<std::endl;
  std::vector<T> u(Us.size());
  PARALLEL_COPY(Us.begin(),Us.end(), u.begin());
  return CalcLogNorm(u);
}



template<class T>
inline
T 
ICR::ICA::DirichletModel<T>::CalcLogNorm(const NaturalParameters<T>& NP)  
{
  // std::cout<<"DIRICHLET NPu = "<<NP<<std::endl;
  std::vector<T> u(NP.size());
  PARALLEL_TRANSFORM(NP.begin(),NP.end(), u.begin(), plus_one());
  return CalcLogNorm(u);
}




template<class T>
inline
ICR::ICA::Moments<T>
ICR::ICA::DirichletModel<T>::CalcMoments(const NaturalParameters<T>& NP)
{
  std::vector<T> u(NP.size());
  PARALLEL_TRANSFORM( NP.begin(), NP.end(), u.begin(), plus_one());

  const T U = PARALLEL_ACCUMULATE(u.begin(),u.end(), 0.0); //maybe save this value need to profile to see if worthwhile
  std::vector<T> the_moments(u.size());
  PARALLEL_TRANSFORM(u.begin(), u.end(), the_moments.begin(), calculate_digamma_minus_digamma(U) );
  // std::cout<<"Dirichlet Moments = "<<Moments<T>(the_moments)<<std::endl;


  return Moments<T>(the_moments);

}


template<class T>
inline
ICR::ICA::NaturalParameters<T>
ICR::ICA::DirichletModel<T>::CalcNP2Data(const Moments<T>& Us)
{
  // std::cout<<"NP.size() = "<< NP.size()<<std::endl;
  // std::cout<<"values.size() = "<< m_values.size()<<std::endl;
  //BOOST_ASSERT(NP.size() == m_values.size());
  ICR::ICA::NaturalParameters<T> NP = NaturalParameters<T>(Us.size());
  PARALLEL_TRANSFORM( Us.begin(), Us.end(), NP.begin(), minus_one());
  //std::cout<<"DIRICHLET NP2DISCRETE = "<<NP<<std::endl;
  return NP;

}
