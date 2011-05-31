#pragma once
#include "ICA/exception/NotConjugate.hpp"

#include <gsl/gsl_sf_psi.h> //for Digamma function gsl_sf_psi
#include <gsl/gsl_sf_gamma.h> //for gamma function

#include <boost/assert.hpp> 
namespace ICR{
  namespace ICA{


    /**  Defines the properties of a Gamma distribution.
     *   This includes how to calculate the Moments and natural parameters
     *   between the different factors and variables.
     */
    template<class T>
    class Gamma 
    {
    public:
      
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

      typedef typename boost::call_traits< VariableNode<T>* const>::param_type
      variable_parameter;
      typedef typename boost::call_traits< VariableNode<T>* const>::value_type
      variable_t;


      /** Calculate the natural log of the normalisation factor.
       *  The normalisation is the inverse of the partition factor.
       *  @param Shape The moments from the VariableNode that represents the shape
       *  @param IScale The moments from the VariableNode that represents the inverse scale.
       *  @return The log or the normalisation.
        */
      static
      data_t
      CalcLogNorm(moments_parameter Shape, 
		  moments_parameter IScale);
      

      
      /** Calculate the natural log of the normalisation factor.
       *  The normalisation is the inverse of the partition factor.
       *  @param NP The natural parameters from which to obtain the shape and 
       *  inverse scale.
       *  @return The log or the normalisation.
       */
      static
      data_t
      CalcLogNorm(NP_parameter NP);
      

      /** Calculate the average of the log of the distribution.
       *  @param Shape The moments from the VariableNode that represents the shape
       *  @param Iscale The moments from the VariableNode that represents the
       *    inverse scale.
       *  @param Data The moments from the Child(or Data) Variable.
       *  @return The average of the log (<log Gamma(mean,precision) >)
       */
      static
      data_t
      CalcAvLog(moments_parameter Shape,
		moments_parameter Iscale,
		moments_parameter Data) ;
      
      /** Calculate a random sample from the distribution. 
       *  This is used for initialisation.
       *  @param Shape The variable that stores the shape.
       *  @param IScale The variable that stores the invserse scale.
       *  @return A random sample for the moments based on the shape and inverse scale.
       */
      static
      moments_t
      CalcSample(variable_parameter Shape,
		 variable_parameter IScale);
        
      /** Calculate the Moments from the Natural Paramters.
       *  @param NP The NaturalParameters from which to calcualate the moments.
       *  @return The Calculated Moments.  
       * 
       *  This function is caled from HiddenMoments or CalculationMoments.
       */
      static
      moments_t
      CalcMoments(NP_parameter NP);

      /** Calculate the Natural Parameters to go the Parent1 (The shape).
       *   @warning 
       *     The shape is NOT conjugate to the Gamma distribution
       *     and so this function should never be called.
       *     It is required to complete the interface.
       *     Calling this function will result in error.
       *   
       *  @param IScale The moments from the IScale Variable.
       *  @param Data The moments from the Child(or Data) Variable.
       *  @return The calculated NaturalParameters.  
       */
      //should never be called
      static
      NP_t
      CalcNP2Parent1(moments_parameter IScale,
		     moments_parameter Data);
      
      /** Calculate the Natural Parameters to go the Parent2 (The invere scale).
       *  This is evaluated at a Factor from the other two attached Variables.
       *  @param Shape The moments from the Shape Variable.
       *  @param Data The moments from the Child(or Data) Variable.
       *  @return The calculated NaturalParameters.  
       */
      static
      NP_t
      CalcNP2Parent2(moments_parameter Shape,
		     moments_parameter Data);

      /** Calculate the Natural Parameters to go the Child (or data).
       *  This is evaluated at a Factor from the other two attached Variables.
       *  @param Shape  The moments from the Shape Variable.
       *  @param IScale The moments from the IScale Variable.
       *  @return The calculated NaturalParameters.  
       */
      static
      NP_t
      CalcNP2Data(moments_parameter Shape, 
		  moments_parameter IScale);

    private:
      //The Log norm is actually evaluated here.
      static
      data_t
      CalcLogNorm(data_t shape,
		  data_t iscale) ;
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
typename ICR::ICA::Gamma<T>::data_t
ICR::ICA::Gamma<T>::CalcLogNorm(data_t shape, 
				     data_t iscale)  
{
  //These must be greater than zero
  BOOST_ASSERT(iscale>0);
  BOOST_ASSERT(shape>0);
  return shape * std::log(iscale) - gsl_sf_lngamma(shape) ;
}

template<class T>
inline
typename ICR::ICA::Gamma<T>::data_t
ICR::ICA::Gamma<T>::CalcLogNorm(moments_parameter Shape,
				     moments_parameter IScale)
{
  //These must be greater than zero
  BOOST_ASSERT(IScale[0]>0);
  return CalcLogNorm(Shape[0], IScale[0]);
}
      
template<class T>
inline
typename ICR::ICA::Gamma<T>::data_t
ICR::ICA::Gamma<T>::CalcLogNorm(NP_parameter NP)
{
  const data_t shape  = NP[1]+1;
  const data_t iscale = -NP[0];
  BOOST_ASSERT(iscale>0);
  return CalcLogNorm(shape, iscale);
}

      
template<class T>
inline
typename ICR::ICA::Gamma<T>::data_t
ICR::ICA::Gamma<T>::CalcAvLog(moments_parameter Shape,
				   moments_parameter IScale,
				   moments_parameter Data)
{
  return CalcNP2Data(Shape,IScale)*Data +CalcLogNorm(Shape,IScale);
}


template<class T>
inline
typename ICR::ICA::Gamma<T>::moments_t
ICR::ICA::Gamma<T>::CalcSample(variable_parameter Shape,
				    variable_parameter IScale) 
{
  rng* random = Random::Instance();
  const data_t shape  = Shape->GetMoments()[0];
  const data_t iscale = IScale->GetMoments()[0];

  //cannot be zero so add something small incase of numerical error.
  const data_t x=  random->gamma(shape,1.0/iscale) + 1e-16; 

  return Moments<T>(x,std::log(x));
}
      

template<class T>
inline
typename ICR::ICA::Gamma<T>::moments_t
ICR::ICA::Gamma<T>::CalcMoments(NP_parameter NP)
{
  const data_t shape  = NP[1]+1.0;
  const data_t iscale = -NP[0];
  BOOST_ASSERT(iscale>0);

  return moments_t(shape/iscale,
		   gsl_sf_psi(shape) - std::log(iscale) 
		   );
}


template<class T>
inline 
typename ICR::ICA::Gamma<T>::NP_t
ICR::ICA::Gamma<T>::CalcNP2Parent1(moments_parameter IScale,
					moments_parameter Data)
{
  throw Exception::NotConjugate();
}
template<class T>
inline 
typename ICR::ICA::Gamma<T>::NP_t
ICR::ICA::Gamma<T>::CalcNP2Parent2(moments_parameter Shape,
					moments_parameter Data)
{
  return NP_t(-Data[0], Shape[0]);
}

template<class T>
inline 
typename ICR::ICA::Gamma<T>::NP_t
ICR::ICA::Gamma<T>::CalcNP2Data(moments_parameter Shape, 
				     moments_parameter IScale)
{
  return  NP_t(-IScale[0], Shape[0] - 1);
}
