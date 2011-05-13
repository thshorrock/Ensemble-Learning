#pragma once

#include "ExponentialModel.hpp"

#include <gsl/gsl_sf_psi.h> //for Digamma function gsl_sf_psi
#include <gsl/gsl_sf_gamma.h> //for gamma function

#include "ICA/parallel_algorithms.hpp"
#include <vector>
#include <boost/assert.hpp> 
namespace ICR{
  namespace ICA{

    template<class T=double>
    class DiscreteModel //: public ExponentialModel<T>
    {
    public:

      static
      T
      CalcLogNorm(const Moments<T>& );
      
      static
      T
      CalcLogNorm(const NaturalParameters<T>& NP);
      
      struct exponentiate
      {
	double operator()(const double d) {return std::exp(d);}
      };
      
      static
      Moments<T>
      CalcSample(const VariableNode<T>* prior) 
      {
	Moments<T> PM = prior->GetMoments();
	std::vector<double> M(PM.size());
	std::transform(PM.begin(),PM.end(),M.begin(),exponentiate());
	return Moments<T>(M);
      }

      static
      Moments<T>
      CalcMoments(const NaturalParameters<T>& NP);

      static
      NaturalParameters<T>
      CalcNP2Prior(const Moments<T>& Discrete);
      
      static
      NaturalParameters<T>
      CalcNP2Data(const Moments<T>& Dirichlet);


    private:

      static
      T
      CalcLogNorm(const std::vector<T>& LogProbs) ;
      // void SetLogNorm()
      // {
	
      // 	T norm = 0;
      // 	for(size_t i=0;i<m_LogPi.size();++i){
      // 	  T lp = (m_LogPi[i]);
      // 	  T e = std::exp(m_LogPi[i]);
      // 	  std::cout<<"LP["<<i<<"]="<<lp<<std::endl;
      // 	  std::cout<<"P["<<i<<"]="<<e<<std::endl;
    
      // 	  norm+=e;
      // 	}
      // 	m_LogNorm  = -std::log(norm);
      // 	std::cout<<"log norm = "<<m_LogNorm<<std::endl;

      // }


      // std::vector<T> m_LogQ;
      // std::vector<T> m_data;
    };

  }
}

namespace {
  struct exponentiate
  {
    // exponentiate(double norm) :m_norm(norm) {}
    double operator()(const double d) {return std::exp(d);}
    // double m_norm;
  };
  struct take_log
  {
    double operator()(const double d){ return std::log(d);}
  };
  struct divide_by
  {
    divide_by(const double d) : m_d(d) {};
    double operator()(const double n){ return n/m_d;}
    double m_d;
  };
 
  struct subtract
  {
    subtract(const double d) : m_d(d) {};
    double operator()(const double m){ return m-m_d;}
    double m_d;
  };
 
 
}

template<class T>
inline
T 
ICR::ICA::DiscreteModel<T>::CalcLogNorm(const std::vector<T>& unLogProbs) 
{

  // std::cout<<"DISCRETE "<<std::endl;

  //un normalised
  
  //If all the probs are very small then can easily get servere numerical errors,
  // eg. norm = 0.
  // To solve this we subtract most significant log before exponetating (and add it again after).
  T LogMax = *PARALLEL_MAX(unLogProbs.begin(),unLogProbs.end());

  std::vector<T> unLogProbsTmp(unLogProbs.size());  
  PARALLEL_TRANSFORM( unLogProbs.begin(),unLogProbs.end(), unLogProbsTmp.begin(), subtract(LogMax));



  //exponentiate log (prob/max)
  std::vector<T> unProbs(unLogProbsTmp.size());
  PARALLEL_TRANSFORM( unLogProbsTmp.begin(),unLogProbsTmp.end(), unProbs.begin(), exponentiate());
  T norm = PARALLEL_ACCUMULATE(unProbs.begin(), unProbs.end(),0.0) ;
  
  return -std::log(norm)- LogMax;

}

template<class T>
inline
T 
ICR::ICA::DiscreteModel<T>::CalcLogNorm(const Moments<T>& Dirichlet) 
{
  //the log probs are provided by Dirichlet, need to pass them on
  std::vector<T> unLogProbs(Dirichlet.size());
  PARALLEL_COPY( Dirichlet.begin(), Dirichlet.end(), unLogProbs.begin());
  
  return  CalcLogNorm(unLogProbs);
}

template<class T>
inline
T 
ICR::ICA::DiscreteModel<T>::CalcLogNorm(const NaturalParameters<T>& NP) 
{
  //The NP are the log probs
  std::vector<T> unLogProbs(NP.size());
  PARALLEL_COPY( NP.begin(), NP.end(), unLogProbs.begin());
  
  return  CalcLogNorm(unLogProbs);
}

template<class T>
inline
ICR::ICA::Moments<T>
ICR::ICA::DiscreteModel<T>::CalcMoments(const NaturalParameters<T>& NP)
{
  //NPs are unnormalised log probabilities.

  std::vector<T> unLogProbs(NP.size());
  std::vector<T> LogProbs(NP.size());
  std::vector<T> Probs(NP.size());
  PARALLEL_COPY( NP.begin(), NP.end(), unLogProbs.begin());
  
  T LogNorm = -CalcLogNorm(unLogProbs);
  PARALLEL_TRANSFORM(  unLogProbs.begin(),  unLogProbs.end(), LogProbs.begin(), subtract(LogNorm));
  
  PARALLEL_TRANSFORM(  LogProbs.begin(),  LogProbs.end(), Probs.begin(), exponentiate());

  return Moments<T>(Probs);
}

  //std::cout<<"NP size = "<<NP.size()<<std::endl;
  // std::vector<T> probs_star(NP.size());
  // std::vector<T> probs(NP.size()); //normalised
  //std::cout<<"m_LogPi.size() = "<<m_LogPi.size()<<std::endl;
  //BOOST_ASSERT(NP.size() == m_LogPi.size());

  //update values
  //std::cout<<"NP = "<<NP<<std::endl;
  //PARALLEL_COPY( NP.begin(), NP.end(), m_LogPi.begin());
  //std::cout<<"m_LogPi.size = "<<m_LogPi.size()<<std::endl;

  // PARALLEL_TRANSFORM(  probs_star.begin(),  probs_star.end(), probs, divide_by(norm));
  
  // std::cout<<"log norm = "<<std::log(norm)<<std::endl;
 
  // SetLogNorm();

      // for(size_t i=0;i<m_NP2Weights.size();++i){
      // 	m_NP2Weights[i] -=lognorm;
      // }

  // double norm2 = 0;
  // T lognorm = 1;
  //  for(size_t i=0;i<m_LogPi.size();++i){
    
  //    lognorm*=m_LogPi[i];
  //    //std::cout<<"m_val = "<<m_LogPi[i]<<std::endl;
  //  }
  // std::cout<<"log norm = "<<lognorm<<std::endl;
  // std::cout<<"norm2 = "<<norm2<<std::endl;

  // for(size_t i=0;i<m_LogPi.size();++i){
    
  //   m_LogPi[i]/=norm2;
  //   //    norm2+=m_LogPi[i];
  //   //std::cout<<"m_val = "<<m_LogPi[i]<<std::endl;
  // }
  //No change for discrete the_moments = m_LogPi

  // Moments<T> m =  Moments<T>(m_LogPi);
  // for(size_t i=0;i<m.size();++i){
  //   std::cout<<"m[i] + m_LogNorm = "<<m[i] - m_LogNorm<<std::endl;

  //   m[i] = std::exp(m[i]+m_LogNorm); //
  // }


  // PARALLEL_TRANSFORM( m.begin(), m.end(), m.begin(), exponentiate());
  

  // std::cout<<"DISCRETE MOMENTS = "<<Moments<T>(m)<<std::endl;


// template<class T>
// inline 
// void
// ICR::ICA::DiscreteModel<T>::SetData(const Moments<T>& child)
// {
//   PARALLEL_COPY(child.begin(), child.end(), m_data.begin());
// }

template<class T>  
inline 
ICR::ICA::NaturalParameters<T>
ICR::ICA::DiscreteModel<T>::CalcNP2Prior(const Moments<T>& Discrete)
{
  //The probabilities are provided by Discrete and need to be passsed onto Prior
  //std::cout<<"CalcNP2Prior = "<<Discrete<<std::endl;

  return NaturalParameters<T>(Discrete);
}
  
   
template<class T>   
inline 
ICR::ICA::NaturalParameters<T>
ICR::ICA::DiscreteModel<T>::CalcNP2Data(const Moments<T>& Dirichlet)
{
  //the log probs are provided by Dirichlet, need to pass them on
  
  return  NaturalParameters<T>(Dirichlet);
  //  BOOST_ASSERT(NP.size() == m_LogPi.size());
  //PARALLEL_COPY( m_LogPi.begin(), m_LogPi.end(), NP.begin(), take_log() );
  //PARALLEL_TRANSFORM( m_LogPi.begin(), m_LogPi.end(), NP.begin(), take_log() );
  //std::cout<<"  DISCRETE NP2DATA = "<<NP<<std::endl;

}
   
