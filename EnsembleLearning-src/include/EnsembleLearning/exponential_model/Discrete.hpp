#pragma once
#ifndef DISCRETE_HPP
#define DISCRETE_HPP


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
#include "EnsembleLearning/message/Moments.hpp"
#include "EnsembleLearning/message/NaturalParameters.hpp"
#include "EnsembleLearning/node/Node.hpp"

#include <gsl/gsl_sf_psi.h> //for Digamma function gsl_sf_psi
#include <gsl/gsl_sf_gamma.h> //for gamma function

#include <vector>
#include <boost/assert.hpp> 
#include <boost/call_traits.hpp> 


namespace ICR{
  namespace EnsembleLearning{

    /**  Defines the properties of a Discrete distribution.
     *   This includes how to calculate the Moments and natural parameters
     *   between the different factors and variables.
     */
    template<class T=double>
    class Discrete 
    {
    public:

      /** @name Useful typdefs for types that are exposed to the user.
       */
      ///@{
      typedef typename boost::call_traits< VariableNode<T>* const>::param_type
      variable_parameter;
      
      typedef typename boost::call_traits< Moments<T> >::param_type
      moments_parameter;
      typedef typename boost::call_traits< Moments<T> >::const_reference
      moments_const_reference;
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
      typedef typename boost::call_traits<std::vector<T> >::value_type
      vector_data_t;
      
      ///@}

      /** Calculate the natural log of the normalisation factor.
       *  The normalisation is the inverse of the partition factor.
       *  @param Dirichlet The moments from the Dirichlet
       *     model provide the log probabilities.
       *  @return The log or the normalisation.
        */
      static
      data_t
      CalcLogNorm(moments_parameter Dirichlet);
      
      /** Calculate the natural log of the normalisation factor.
       *  The normalisation is the inverse of the partition factor.
       *  @param NP The natural parameters from which to obtain the mean and precision.
       *  @return The log or the normalisation.
       */
      static
      data_t
      CalcLogNorm(NP_parameter NP);
      
      /** The discrete variable holds the probabilities and is not sampled.
       *  This function returns the probabilities.
       *  This is used for initialisation.
       *  @param prior The log probabilites from the Dirichlet distribution.
       *  @return The probabilities.
       */
      static
      moments_t
      CalcSample(variable_parameter prior);
      
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
 
      /** Calculate the Natural Parameters to go the Dirichlet Prior.
       *  @param Discrete The moments from the Discrete factor
       *    (connected to the Mixture Model).
       *  @return The calculated NaturalParameters.  
       */
      static
      NP_t
      CalcNP2Prior(moments_parameter Discrete);
      
      /** Calculate the Natural Parameters to go the Mixture.
       *  @param Dirichlet The moments from the Dirichlet factor
       *  @return The calculated NaturalParameters.  
       */
      static
      NP_t
      CalcNP2Data(moments_parameter Dirichlet);


    private:

      struct exponentiate
      {
	data_t operator()(data_parameter d) {return std::exp(d);}
      };
      struct take_log
      {
	data_t operator()(data_parameter d){ return std::log(d);}
      };
      struct divide_by
      {
	divide_by(data_parameter d) : m_d(d) {};
	data_t operator()(data_parameter n){ return n/m_d;}
	data_t m_d;
      };
 
      struct subtract
      {
	subtract(data_parameter d) : m_d(d) {};
	data_t operator()(data_parameter m){ return m-m_d;}
	data_t m_d;
      };
 
      static
      data_t
      CalcLogNorm(vector_data_parameter LogProbs) ;
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
typename ICR::EnsembleLearning::Discrete<T>::data_t
ICR::EnsembleLearning::Discrete<T>::CalcLogNorm(vector_data_parameter unLogProbs) 
{
  /* Unnormalised
   *If all the probs are very small then can easily get servere numerical errors,
   * eg. norm = 0.
   * To solve this we subtract most significant log before exponetating 
   *   (and add it again after).
   */
  const data_t LogMax = *PARALLEL_MAX(unLogProbs.begin(),unLogProbs.end());

  std::vector<T> unLogProbsTmp(unLogProbs.size());  
  PARALLEL_TRANSFORM( unLogProbs.begin(),unLogProbs.end(), unLogProbsTmp.begin(), 
		      subtract(LogMax));
  //exponentiate log (prob/max)
  std::vector<T> unProbs(unLogProbsTmp.size());
  PARALLEL_TRANSFORM( unLogProbsTmp.begin(),unLogProbsTmp.end(), unProbs.begin(), 
		      exponentiate());
  const data_t norm = PARALLEL_ACCUMULATE(unProbs.begin(), unProbs.end(),0.0) ;
  
  return -std::log(norm)- LogMax;

}

template<class T>
inline
typename ICR::EnsembleLearning::Discrete<T>::data_t 
ICR::EnsembleLearning::Discrete<T>::CalcLogNorm(moments_parameter Dirichlet) 
{
  //the log probs are provided by Dirichlet, need to pass them on
  std::vector<T> unLogProbs(Dirichlet.size());
  PARALLEL_COPY( Dirichlet.begin(), Dirichlet.end(), unLogProbs.begin());
  
  return  CalcLogNorm(unLogProbs);
}

template<class T>
inline
typename ICR::EnsembleLearning::Discrete<T>::data_t 
ICR::EnsembleLearning::Discrete<T>::CalcLogNorm(NP_parameter NP) 
{
  //The NP are the log probs
  std::vector<T> unLogProbs(NP.size());
  PARALLEL_COPY( NP.begin(), NP.end(), unLogProbs.begin());
  
  return  CalcLogNorm(unLogProbs);
}

template<class T>
inline
typename ICR::EnsembleLearning::Discrete<T>::moments_t
ICR::EnsembleLearning::Discrete<T>::CalcSample(variable_parameter prior) 
{
  const moments_t PM = prior->GetMoments();
  std::vector<data_t> M(PM.size());
  std::transform(PM.begin(),PM.end(),M.begin(),exponentiate());
  return moments_t(M);
}


template<class T>
inline
typename ICR::EnsembleLearning::Discrete<T>::vector_data_t
ICR::EnsembleLearning::Discrete<T>::CalcMean(NP_parameter NP)
{
  //Get the probabilities from the moments
  moments_t moments= CalcMoments(NP);
  data_t mean = 0;
  for(size_t i=0;i<moments.size();++i){
    mean += i*moments[i];
  }
  std::vector<data_t> vmean(1);
  vmean[0]=mean;
  return vmean;
}

template<class T>
inline
typename ICR::EnsembleLearning::Discrete<T>::vector_data_t
ICR::EnsembleLearning::Discrete<T>::CalcPrecision(NP_parameter NP)
{
  //Get the probabilities from the moments
  moments_t moments= CalcMoments(NP);
  data_t mean = CalcMean(NP)[0];
  data_t variance = 0;
  for(size_t i=0;i<moments.size();++i){
    variance += moments[i]*(i-mean)*(i-mean);
  }
  std::vector<data_t> vprec(1);
  vprec[0]=1.0/variance;
  return vprec;
}

template<class T>
inline
typename ICR::EnsembleLearning::Discrete<T>::moments_t
ICR::EnsembleLearning::Discrete<T>::CalcMoments(NP_parameter NP)
{
  //NPs are unnormalised log probabilities.

  std::vector<data_t> unLogProbs(NP.size());  //unnormalised
  std::vector<data_t> LogProbs(NP.size());    //log normalised
  std::vector<data_t> Probs(NP.size());       //normalised
  PARALLEL_COPY( NP.begin(), NP.end(), unLogProbs.begin());
  //calculate the log partition factor, Z.  The normalisation is 1/Z.
  const data_t LogNorm = -CalcLogNorm(unLogProbs); //Thats why there is a minus here
  PARALLEL_TRANSFORM(  unLogProbs.begin(),  unLogProbs.end(), 
		       LogProbs.begin(), subtract(LogNorm));
  
  PARALLEL_TRANSFORM(  LogProbs.begin(),  LogProbs.end(), 
		       Probs.begin(), exponentiate());

  return moments_t(Probs);
}


template<class T>  
inline 
typename ICR::EnsembleLearning::Discrete<T>::NP_t
ICR::EnsembleLearning::Discrete<T>::CalcNP2Prior(moments_parameter Discrete)
{
  //The probabilities are provided by Discrete and need to be passsed onto Prior
  // Want to simply foward these:
  NP_t NP(Discrete.size());
  PARALLEL_COPY(Discrete.begin(), Discrete.end(), NP.begin());
  return NP;
}
  
   
template<class T>   
inline 
typename ICR::EnsembleLearning::Discrete<T>::NP_t
ICR::EnsembleLearning::Discrete<T>::CalcNP2Data(moments_parameter Dirichlet)
{
  //the log probs are provided by Dirichlet, need to pass them on
  //Want to simply copy these
  NP_t NP(Dirichlet.size());
  PARALLEL_COPY(Dirichlet.begin(), Dirichlet.end(), NP.begin());
  return  NP;
}
   
#endif  // guard for DISCRETE_HPP
