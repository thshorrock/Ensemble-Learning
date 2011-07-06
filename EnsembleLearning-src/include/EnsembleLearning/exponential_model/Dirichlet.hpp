#pragma once
#ifndef DIRICHLET_HPP
#define DIRICHLET_HPP


/***********************************************************************************
 ***********************************************************************************
 **                                                                               **
 **  Copyright (C) 2011 Tom Shorrock <t.h.shorrock@gmail.com> 
 **                                                                               **
 **                                                                               **
 **  This program is free software; you can redistribute it and/or                **
 **  modify it under the terms of the GNU General Public License                  **
 **  as published by the Free Software Foundation; either version 2               **
 **  of the License, or (at your option) any later version.                       **
 **                                                                               **
 **  This program is distributed in the hope that it will be useful,              **
 **  but WITHOUT ANY WARRANTY; without even the implied warranty of               **
 **  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                **
 **  GNU General Public License for more details.                                 **
 **                                                                               **
 **  You should have received a copy of the GNU General Public License            **
 **  along with this program; if not, write to the Free Software                  **
 **  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  **
 **                                                                               **
 ***********************************************************************************
 ***********************************************************************************/




#include "Random.hpp"
#include "EnsembleLearning/detail/parallel_algorithms.hpp"
#include "EnsembleLearning/detail/MACRO_defaults.hpp"//for ENSEMBLE_LEARNING_COMPONENTS
#include "EnsembleLearning/message/Moments.hpp"
#include "EnsembleLearning/message/NaturalParameters.hpp"
#include "EnsembleLearning/node/Node.hpp"

#include <gsl/gsl_sf_psi.h> //for Digamma function gsl_sf_psi
#include <gsl/gsl_sf_gamma.h> //for gamma function

#include <boost/assert.hpp> 
#include <boost/call_traits.hpp> 

#include <vector>
#include <cmath>


namespace ICR{
  namespace EnsembleLearning{

    /**  Defines the properties of a Dirichlet distribution.
     *   This includes how to calculate the Moments and natural parameters
     *   between the different factors and variables.
     */
    template<class T>
    class Dirichlet 
    {
    public:

      /** @name Useful typdefs for types that are exposed to the user.
       */
      ///@{
      typedef typename boost::call_traits< const Moments<T,ENSEMBLE_LEARNING_COMPONENTS>* >::param_type
      moments_parameter;
      typedef typename boost::call_traits<  Moments<T,ENSEMBLE_LEARNING_COMPONENTS> >::value_type
      moments_t;

      typedef typename boost::call_traits< NaturalParameters<T,ENSEMBLE_LEARNING_COMPONENTS> >::param_type
      NP_parameter;
      typedef typename boost::call_traits< NaturalParameters<T,ENSEMBLE_LEARNING_COMPONENTS> >::value_type
      NP_t;

      typedef typename boost::call_traits<T>::value_type
      data_t;
      typedef typename boost::call_traits<T>::param_type
      data_parameter;

      typedef typename boost::call_traits<std::vector<T> >::param_type
      vector_data_parameter;
      typedef typename boost::call_traits<std::vector<T> >::value_type
      vector_data_t;

      typedef typename boost::call_traits<const VariableNode<T,ENSEMBLE_LEARNING_COMPONENTS>*>::param_type
      Variable_parameter;

      
      ///@}

      //moments from parents (the u values)
      /** Calculate the natural log of the normalisation factor.
       *  The normalisation is the inverse of the partition factor.
       *  @param Us The moments from the Dirichlet constant node.
       *  @return The log or the normalisation.
       */
      static
      data_t
      CalcLogNorm(moments_parameter Us) ;

      /** Calculate the natural log of the normalisation factor.
       *  The normalisation is the inverse of the partition factor.
       *  @param NP The natural parameters from which to obtain the mean and precision.
       *  @return The log or the normalisation.
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
      CalcSample(moments_parameter V);

      /** Calculate the Mean from the Natural Paramters.
       *  @param NP The NaturalParameters from which to calcualate the moments.
       *  @return The mean of the distribution.  
       * 
       */
      static
      vector_data_t
      CalcMean(NP_parameter NP)  ;
      
      /** Calculate the standard deviation from the Natural Paramters.
       *  @param NP The NaturalParameters from which to calcualate the moments.
       *  @return The mean of the distribution.  
       * 
       */
      static
      vector_data_t
      CalcPrecision(NP_parameter NP)  ;

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


      static
      NaturalParameters<T,ENSEMBLE_LEARNING_COMPONENTS>
      CalcNP2Prior(moments_parameter Us) {return NP_t();};


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
      
      class
      divide_by
      {
      public:
	divide_by(data_parameter U)
	  : m_divide_me(U)
	{}
	data_t operator()(data_parameter u) 
	{
	  return u/m_divide_me;
	}
      private:
	data_t m_divide_me;
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
typename ICR::EnsembleLearning::Dirichlet<T>::data_t 
ICR::EnsembleLearning::Dirichlet<T>::CalcLogNorm(vector_data_parameter u)  
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
typename ICR::EnsembleLearning::Dirichlet<T>::data_t 
ICR::EnsembleLearning::Dirichlet<T>::CalcLogNorm(moments_parameter Us)  
{
  std::vector<T> u(Us->size());
  PARALLEL_COPY(Us->begin(),Us->end(), u.begin());
  return CalcLogNorm(u);
}



template<class T>
inline
typename ICR::EnsembleLearning::Dirichlet<T>::data_t 
ICR::EnsembleLearning::Dirichlet<T>::CalcLogNorm(NP_parameter NP)  
{
  std::vector<data_t> u(NP.size());
  PARALLEL_TRANSFORM(NP.begin(),NP.end(), u.begin(), plus_one());
  return CalcLogNorm(u);
}



template<class T>
inline
typename ICR::EnsembleLearning::Dirichlet<T>::moments_t
ICR::EnsembleLearning::Dirichlet<T>::CalcSample(moments_parameter Mus) 
{
  rng* random = Random::Instance();
  const size_t     size   = Mus->size();

  double u[size]; //double to go to rng
  double x[size];
  PARALLEL_COPY(Mus->begin(),Mus->end(),u);
  random->dirichlet(size, u, x);
	
  std::vector<data_t> M(size);

  PARALLEL_TRANSFORM(x,x+size,M.begin(),take_log());

  
  return moments_t(M);
}

template<class T>
inline
typename ICR::EnsembleLearning::Dirichlet<T>::vector_data_t
ICR::EnsembleLearning::Dirichlet<T>::CalcMean(NP_parameter NP)
{
  std::vector<T> u(NP.size());
  PARALLEL_TRANSFORM( NP.begin(), NP.end(), u.begin(), plus_one());
  //maybe save this value need to profile to see if worthwhile?
  const data_t U = PARALLEL_ACCUMULATE(u.begin(),u.end(), 0.0); 
  std::vector<T> means(u.size());
  PARALLEL_TRANSFORM(u.begin(), u.end(), means.begin(), divide_by(U));
  return means;
}

template<class T>
inline
typename ICR::EnsembleLearning::Dirichlet<T>::vector_data_t
ICR::EnsembleLearning::Dirichlet<T>::CalcPrecision(NP_parameter NP)
{
  std::vector<T> u(NP.size());
  PARALLEL_TRANSFORM( NP.begin(), NP.end(), u.begin(), plus_one());
  //maybe save this value need to profile to see if worthwhile?
  const data_t U = PARALLEL_ACCUMULATE(u.begin(),u.end(), 0.0); 
  const data_t UplusOne = U+1;
  std::vector<T> means = CalcMean(NP);
  std::vector<T> Prec(means.size());
  for(size_t i=0;i<Prec.size();++i){
    Prec[i] = UplusOne/(means[i]*(1-means[i]));
  }
  
  return Prec;
}


template<class T>
inline
typename ICR::EnsembleLearning::Dirichlet<T>::moments_t
ICR::EnsembleLearning::Dirichlet<T>::CalcMoments(NP_parameter NP)
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
typename ICR::EnsembleLearning::Dirichlet<T>::NP_t
ICR::EnsembleLearning::Dirichlet<T>::CalcNP2Data(moments_parameter Us)
{
  NP_t NP ;
  PARALLEL_TRANSFORM( Us->begin(), Us->end(), NP.begin(), minus_one());
  return NP;
}

#endif  // guard for DIRICHLET_HPP
