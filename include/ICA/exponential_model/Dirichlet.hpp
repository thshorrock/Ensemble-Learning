#pragma once

#include "Random.hpp"
// #include "ExponentialModel.hpp"

#include <gsl/gsl_sf_psi.h> //for Digamma function gsl_sf_psi
#include <gsl/gsl_sf_gamma.h> //for gamma function

#include <vector>
#include "ICA/detail/parallel_algorithms.hpp"

#include <boost/assert.hpp> 

namespace ICR{
  namespace ICA{

    /**  Defines the properties of a Dirichlet distribution.
     *   This includes how to calculate the Moments and natural parameters
     *   between the different factors and variables.
     */
    template<class T>
    class Dirichlet 
    {
    public:

      typedef typename boost::call_traits< Moments<T> >::param_type
      moments_parameter;
      typedef typename boost::call_traits< Moments<T> >::value_type
      moments_t;

      typedef typename boost::call_traits< NaturalParameters<T> >::param_type
      NP_parameter;
      typedef typename boost::call_traits< NaturalParameters<T> >::value_type
      NP_t;

      typedef typename boost::call_traits<T>::value_type
      data_t;
      typedef typename boost::call_traits<T>::param_type
      data_parameter;

      typedef typename boost::call_traits<std::vector<T> >::param_type
      vector_data_parameter;
      typedef typename boost::call_traits<std::vector<T> >::param_type
      vector_data_t;
      
      /** Calculate the natural log of the normalisation factor.
       *  The normalisation is the inverse of the partition factor.
       *  @param Us The moments from the Dirichlet constant node.
       *  @return The log or the normalisation.
        */
      //moments from parents (the u values)
      static
      data_t
      CalcLogNorm(moments_parameter Us) ;

      /** Calculate the natural log of the normalisation factor.
       *  The normalisation is the inverse of the partition factor.
       *  @param NP The natural parameters from which to obtain the mean and precision.
       */
      static
      data_t
      CalcLogNorm(NP_parameter NP) ;
      
      /** Calculate a random sample for a Dirichlet distribution.
       *  This is used for initialisation.
       *  @param V  The Dirichlet Constant VariableNode.
       *  @return A random sample for the moments based on the Dirichlet weights.
       */
      static
      moments_t
      CalcSample(const VariableNode<T>* V);

      /** Calculate the Moments from the Natural Paramters.
       *  @param NP The NaturalParameters from which to calcualate the moments.
       *  @return The Calculated Moments.  
       * 
       *  This function is caled from HiddenMoments or CalculationMoments.
       */
      static
      moments_t
      CalcMoments(NP_parameter NP);

      /** Calculate the Natural Parameters to go the Mixture.
       *  @param Us The moments from the Dirichlet constant Node.
       *  @return The calculated NaturalParameters.  
       */
      static
      NP_t
      CalcNP2Data(moments_parameter Us);


    private:

      struct take_log
      {
	data_t operator()(data_parameter d){ return std::log(d);}
      };

      struct plus_one
      {
	data_t operator()(data_parameter d) {return d+1.0;}
      };

      struct minus_one
      {
	data_t operator()(data_parameter d) {return d-1.0;}
      };

      class
      calculate_digamma_minus_digamma
      {
      public:
	calculate_digamma_minus_digamma(data_parameter U)
	  : m_subtract_me(gsl_sf_psi(U))
	{}
	data_t operator()(data_parameter u) 
	{
	  return gsl_sf_psi(u) - m_subtract_me;
	}
      private:
	data_t m_subtract_me;
      };
      
      static
      data_t
      CalcLogNorm(vector_data_parameter u) ;
    };

  }
}


/**********************************************************
 **********************************************************
 *************** IMPLEMENTATION ***************************
 **********************************************************
 **********************************************************/



template<class T>
inline
typename ICR::ICA::Dirichlet<T>::data_t 
ICR::ICA::Dirichlet<T>::CalcLogNorm(vector_data_parameter u)  
{
  const data_t U = PARALLEL_ACCUMULATE(u.begin(), u.end(), 0.0);
    
  //calculate Gamma(m_value[i]) for each i
  std::vector<T> Gamma_u(u.size());
  PARALLEL_TRANSFORM(u.begin(), u.end(), Gamma_u.begin(), gsl_sf_lngamma);
    
  const data_t SumLnGamma_u = PARALLEL_ACCUMULATE(Gamma_u.begin(), Gamma_u.end(), 0.0);
    
  return  gsl_sf_lngamma(U) - SumLnGamma_u;
}


template<class T>
inline
typename ICR::ICA::Dirichlet<T>::data_t 
ICR::ICA::Dirichlet<T>::CalcLogNorm(moments_parameter Us)  
{
  std::vector<T> u(Us.size());
  PARALLEL_COPY(Us.begin(),Us.end(), u.begin());
  return CalcLogNorm(u);
}



template<class T>
inline
typename ICR::ICA::Dirichlet<T>::data_t 
ICR::ICA::Dirichlet<T>::CalcLogNorm(NP_parameter NP)  
{
  std::vector<data_t> u(NP.size());
  PARALLEL_TRANSFORM(NP.begin(),NP.end(), u.begin(), plus_one());
  return CalcLogNorm(u);
}



template<class T>
inline
typename ICR::ICA::Dirichlet<T>::moments_t
ICR::ICA::Dirichlet<T>::CalcSample(const VariableNode<T>* V) 
{
  rng* random = Random::Instance();
  const moments_t  Mus = V->GetMoments();
  const size_t     size   = Mus.size();

  double u[size]; //double to go to rng
  double x[size];
  PARALLEL_COPY(Mus.begin(),Mus.end(),u);
  random->dirichlet(size, u, x);
	
  std::vector<data_t> M(size);
  PARALLEL_TRANSFORM(x,x+size,M.begin(),take_log());
  return moments_t(M);
}


template<class T>
inline
typename ICR::ICA::Dirichlet<T>::moments_t
ICR::ICA::Dirichlet<T>::CalcMoments(NP_parameter NP)
{
  std::vector<T> u(NP.size());
  PARALLEL_TRANSFORM( NP.begin(), NP.end(), u.begin(), plus_one());
  //maybe save this value need to profile to see if worthwhile?
  const data_t U = PARALLEL_ACCUMULATE(u.begin(),u.end(), 0.0); 
  std::vector<T> the_moments(u.size());
  PARALLEL_TRANSFORM(u.begin(), u.end(), the_moments.begin(), calculate_digamma_minus_digamma(U) );

  return moments_t(the_moments);
}


template<class T>
inline
typename ICR::ICA::Dirichlet<T>::NP_t
ICR::ICA::Dirichlet<T>::CalcNP2Data(moments_parameter Us)
{
  NP_t NP = NP_t(Us.size());
  PARALLEL_TRANSFORM( Us.begin(), Us.end(), NP.begin(), minus_one());
  return NP;
}
