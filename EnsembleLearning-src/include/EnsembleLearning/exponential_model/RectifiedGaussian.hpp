#pragma once
#include "EnsembleLearning/node/variable/Calculation.hpp"
#include "EnsembleLearning/calculation_tree/Context.hpp"
#include "EnsembleLearning/calculation_tree/Expression.hpp"

#include "Random.hpp"
#include <gsl/gsl_sf_erf.h>

namespace ICR{
  namespace EnsembleLearning{
    


    /**  Defines the properties of a Rectified Gaussian distribution.
     *   This includes how to calculate the Moments and natural parameters
     *   between the different factors and variables.
     */
    template<class T=double>
    class RectifiedGaussian
    {
    public:

      typedef typename boost::call_traits< VariableNode<T>* const>::param_type
      variable_parameter;
      typedef typename boost::call_traits< VariableNode<T>* const>::value_type
      variable_t;
      typedef typename boost::call_traits< std::vector<VariableNode<T>*> >::param_type
      variable_vector_parameter;

      
      typedef typename boost::call_traits< Moments<T> >::param_type
      moments_parameter;
      typedef typename boost::call_traits< Moments<T> >::value_type
      moments_t;

      typedef typename boost::call_traits< DeterministicNode<RectifiedGaussian<T>, T>* >::param_type
      deterministic_parameter;
      typedef typename boost::call_traits< Expression<T>* >::param_type
      expression_parameter;
      typedef typename boost::call_traits<Context<T> >::param_type
      context_parameter;
      
      typedef typename boost::call_traits< NaturalParameters<T> >::param_type
      NP_parameter;
      typedef typename boost::call_traits< NaturalParameters<T> >::value_type
      NP_t;

      typedef typename boost::call_traits<T>::value_type
      data_t;
      typedef typename boost::call_traits<T>::param_type
      data_parameter;
      
      typedef typename boost::call_traits<std::vector<T> >::value_type
      vector_data_t;
      
      
      /** Calculate the natural log of the normalisation factor.
       *  The normalisation is the inverse of the partition factor.
       *  @param Mean The moments from the VariableNode that represents Mean
       *  @param Precision The moments from the VariableNode that represents the precision.
       *  @return The log or the normalisation.
       */
      static
      data_t
      CalcLogNorm(moments_parameter Mean,
		  moments_parameter Precision) ;

     
      
      /** Calculate the natural log of the normalisation factor.
       *  The normalisation is the inverse of the partition factor.
       *  @param NP The natural parameters from which to obtain the mean and precision.
       *  @return The log or the normalisation.
       */
      static
      data_t
      CalcLogNorm(NP_parameter NP) ;
      

      /** Calculate the average of the log of the distribution.
       *  @param Mean The moments from the VariableNode that represents Mean
       *  @param Precision The moments from the VariableNode that represents the precision.
       *  @param Data The moments from the Child(or Data) Variable.
       *  @return The average of the log (<log Rectified Gaussian(mean,precision) >)
       */
      static
      data_t
      CalcAvLog(moments_parameter Mean,
		moments_parameter Precision,
		moments_parameter Data) ;

      

      /** Calculate a random sample from the distribution. 
       *  This is used for initialisation.
       *  @param Mean The variable that acts as the mean.
       *  @param Precision The variable that acts as the precision.
       *  @return A random sample for the moments based on the mean and precision.
       */
      static
      moments_t
      CalcSample(variable_parameter Mean,
		 variable_parameter Precision) ;


      /** Calculate a random sample from a mixture.
       *  This is used for initialisation.
       *  @param mean_nodes The vector of variables that acts as the mean.
       *  @param precision_nodes The vector variable that acts as the precision.
       *  @param weights_node The vector variable holds the mixture weights.
       *  @return A random sample for the moments based on the means and precisions.
       */
      static
      moments_t
      CalcSample(variable_vector_parameter mean_nodes,
		 variable_vector_parameter precision_nodes,
		 variable_parameter weights_node);


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
      CalcMoments(NP_parameter NP)  ;
      
      
      /** Calculate the Natural Parameters to go the Parent1 (The mean).
       *  This is evaluated at a Factor from the other two attached Variables.
       *  @param Precision The moments from the Precision Variable.
       *  @param Data The moments from the Child(or Data) Variable.
       *  @return The calculated NaturalParameters.  
       */
      static
      NP_t
      CalcNP2Parent1(moments_parameter Precision, 
		     moments_parameter Data)  ;



      /** Calculate the Natural Parameters to go the Parent2 (The Precision).
       *  This is evaluated at a Factor from the other two attached Variables.
       *  @param Mean The moments from the Mean Variable.
       *  @param Data The moments from the Child(or Data) Variable.
       *  @return The calculated NaturalParameters.  
       */
      static
      NP_t
      CalcNP2Parent2(moments_parameter Mean,
		     moments_parameter Data)  ;
      

      
      /** Calculate the Natural Parameters to go the Child (or data).
       *  This is evaluated at a Factor from the other two attached Variables.
       *  @param Mean The moments from the Mean Variable.
       *  @param Precision The moments from the Precision Variable.
       *  @return The calculated NaturalParameters.  
       */
      static
      NP_t
      CalcNP2Data(moments_parameter Mean, 
		  moments_parameter Precision)  ;




      //Deterministic to Stock
      /** Calculate the Natural Parameters to go to
       *  each of the Parent Variables in a calculation.
       *  This is evaluated from all the other variables and the data.
       *  @param Parent The Node that the messag is to go to.
       *  @param Data The moments from the Child(or Data) Variable.
       *  @param C The context (the parent variables are obtainable from this)
       *  @return The calculated NaturalParameters.  
       */
      static
      NP_t
      CalcNP2Parent(variable_parameter Parent, 
		    deterministic_parameter Data, 
		    context_parameter C);


      //Deterministic to Stock
      /** Calculate the Natural Parameters to go to the Deterministic Variable
       *   in a calculation.
       *  This is evaluated from all the other variables and the data.
       *  @param Expr The expression to evaluate 
       *  @param C The context (the parent variables are obtainable from this)
       *  @return The calculated NaturalParameters.  
       */
      static
      NP_t
      CalcNP2Deterministic(expression_parameter Expr,
			   context_parameter C);

    private:
      
      struct Erfcx
      {
	data_t
	operator()(data_parameter x)
	{
	  return std::exp(x*x)*gsl_sf_erfc(x);
	}
      }; 

      static
      data_t
      CalcLogNorm(data_parameter mean,
		  data_parameter mean_squared,
		  data_parameter precision);
    };

  }
}

template<class T>
inline
typename ICR::EnsembleLearning::RectifiedGaussian<T>::data_t
ICR::EnsembleLearning::RectifiedGaussian<T>::CalcLogNorm(data_parameter mean, 
					    data_parameter mean_squared,
					    data_parameter precision)  
{ 
  // Typo in Miskin's or Winn's thesis for this formula?
  const data_t LN =  0.5* ( std::log(2.0*precision/(M_PI)) - (precision) * (mean_squared))
    - gsl_sf_log_erfc (-mean*std::sqrt(precision/2.0) ); 
  return LN;
}

template<class T>
inline
typename ICR::EnsembleLearning::RectifiedGaussian<T>::data_t
ICR::EnsembleLearning::RectifiedGaussian<T>::CalcLogNorm(moments_parameter Mean,
					    moments_parameter Precision)  
{
  return  CalcLogNorm(Mean[0], Mean[1], Precision[0]);
}

template<class T>
inline
typename ICR::EnsembleLearning::RectifiedGaussian<T>::data_t
ICR::EnsembleLearning::RectifiedGaussian<T>::CalcLogNorm(NP_parameter NP)  
{
  //Obtain data from natural parameters.
  const data_t precision    = -NP[1]*2.0;
  const data_t mean         =  NP[0]/(precision);
  const data_t mean_squared =  (mean)*(mean);
  return  CalcLogNorm(mean, mean_squared, precision);
}
    
template<class T>
inline
typename ICR::EnsembleLearning::RectifiedGaussian<T>::data_t
ICR::EnsembleLearning::RectifiedGaussian<T>::CalcAvLog(moments_parameter Mean,
					  moments_parameter Precision,
					  moments_parameter Data)
{
  //These moments must have a size of 2
  BOOST_ASSERT(Mean.size() == 2);
  BOOST_ASSERT(Precision.size() == 2);
  const NP_t NPData = CalcNP2Data(Mean, Precision);
  return 
    NPData*Data + CalcLogNorm(Mean,Precision);
}
    

template<class T>
inline
typename ICR::EnsembleLearning::RectifiedGaussian<T>::moments_t
ICR::EnsembleLearning::RectifiedGaussian<T>::CalcSample(variable_parameter Mean,
					   variable_parameter Precision) 
{

  rng* random = Random::Instance();
  const data_t mean = Mean->GetMoments()[0];
  const data_t prec = Precision->GetMoments()[0];
  const data_t x=  random->gaussian_tail(1.0/std::sqrt(prec),mean);
  return moments_t(x, x*x +1.0/prec);
	
} 


template<class T>
typename ICR::EnsembleLearning::RectifiedGaussian<T>::moments_t
ICR::EnsembleLearning::RectifiedGaussian<T>::CalcSample(variable_vector_parameter m_mean_nodes,
					   variable_vector_parameter m_precision_nodes,
					   variable_parameter m_weights_node)
	
{
  const moments_t weights = m_weights_node->GetMoments();
  std::vector<moments_t > mean(weights.size());
  std::vector<moments_t > precision(weights.size());

  PARALLEL_TRANSFORM(m_mean_nodes.begin(), m_mean_nodes.end(),
		     mean.begin(), 
		     boost::bind(&VariableNode<T>::GetMoments,
				 _1)
		     );
	
  PARALLEL_TRANSFORM(m_precision_nodes.begin(), m_precision_nodes.end(),
		     precision.begin(), 
		     boost::bind(&VariableNode<T>::GetMoments,
				 _1)
		     );
	
  moments_t AvMean(0,0);
  AvMean= PARALLEL_INNERPRODUCT(weights.begin(), weights.end(),mean.begin(), AvMean);
  moments_t AvPrec(0,0);
  AvPrec = PARALLEL_INNERPRODUCT(weights.begin(), weights.end(),precision.begin(), AvPrec);

  rng* random = Random::Instance();
  const data_t mean0 = AvMean[0];
  const data_t prec0 = AvPrec[0];
  const data_t x=  random->gaussian_tail(1.0/std::sqrt(prec0),mean0);
  return moments_t(x, x*x +1.0/prec0);
}


template<class T>
inline
typename ICR::EnsembleLearning::RectifiedGaussian<T>::vector_data_t
ICR::EnsembleLearning::RectifiedGaussian<T>::CalcMean(NP_parameter NP)
{
  vector_data_t vmean(1);
  vmean[0] =  NP[0]/(CalcPrecision(NP)[0]); 
  return vmean;
}

template<class T>
inline
typename ICR::EnsembleLearning::RectifiedGaussian<T>::vector_data_t
ICR::EnsembleLearning::RectifiedGaussian<T>::CalcPrecision(NP_parameter NP)
{
  vector_data_t vprec(1);
  vprec[0] = -NP[1]*2.0;
  return vprec;
}


template<class T>
inline
typename ICR::EnsembleLearning::RectifiedGaussian<T>::moments_t
ICR::EnsembleLearning::RectifiedGaussian<T>::CalcMoments(NP_parameter NP)
{
  //NP must have a size of 2
  BOOST_ASSERT(NP.size() == 2);
  const data_t precision    =  -NP[1]*2.0;
  const data_t mean         =  NP[0]/(precision); 
  //Mean squared on HiddenVariable  is local and not averaged.
  const data_t mean_squared =  (mean)*(mean); 
  const data_t precision_squared =  (precision)*(precision); 
  
  //The argument to the erfcx function (scaled complimentary error function).
  const data_t arg = -mean*std::sqrt(precision/2.0);

  moments_t M(0,0);
  /* If the mod(argument) becomes too large then numerical errors can force answer to zero,
   * which breaks the gamma distribution.
   * Therefore approximate the function in these cases.
   */
  
  if (arg>15) 
    //will start getting fatal numerical errors above 20 but approximation is not good until about now.
    {
      /** The second moment in this limit is incorrectly given by Miskin.
       *  ( moments_t(-1.0/(mean*precision), 2.0/(mean_squared*precision_squared));
       *  From the article 'Closed-form approximations to the error and
       *  complementary error functions and their applications in
       *  atmospheric science we get (to one approximation better) 
       */
      return  moments_t(-1.0/(mean*precision) + 2.5/(mean*mean_squared*precision_squared),
			 ( mean_squared + 1.0 /(precision))*(1.0-1.0/sqrt(2))
			 + 2.5/(mean_squared*precision_squared));

    }
  if (arg<-15) 
    {
      //its tends to a gaussia
      return  moments_t(mean, mean_squared + 1.0 /(precision));
    }
  else
    {
      //Calculate it properly
      Erfcx erfcx;
      //The following is expensive to calculate and is used twice.
      const data_t iRpt = 1.0/(std::sqrt(M_PI*precision)*erfcx(arg));
      
      return  moments_t(mean + std::sqrt(2.0)*iRpt,
			 mean_squared + 1.0 /(precision) + mean*iRpt);
    }
}




template<class T> 
inline
typename ICR::EnsembleLearning::RectifiedGaussian<T>::NP_t
ICR::EnsembleLearning::RectifiedGaussian<T>::CalcNP2Parent1(moments_parameter Precision,moments_parameter Data)
{
  BOOST_ASSERT(Precision.size() == 2);
  BOOST_ASSERT(Data.size() == 2);
  const data_t precision = Precision[0];
  const data_t data = Data[0];
  return NP_t(precision*data, -0.5*precision);
}

template<class T>
inline
typename ICR::EnsembleLearning::RectifiedGaussian<T>::NP_t
ICR::EnsembleLearning::RectifiedGaussian<T>::CalcNP2Parent2(moments_parameter Mean,moments_parameter Data)
{
  BOOST_ASSERT(Mean.size() == 2);
  BOOST_ASSERT(Data.size() == 2);
  const data_t mean = Mean[0];
  const data_t mean_squared = Mean[1];
  const data_t data = Data[0];
  const data_t data_squared = Data[1];

  //This is going to the precision - therefore we know that:
  BOOST_ASSERT(data_squared - 2*data*mean + mean_squared>0);
  //(This should be a guareentee (should not require run-time check)
  // - if the above condition fails then it should fail noisily )
  return  NP_t(-0.5*(data_squared - 2*data*mean + mean_squared), 0.5 );
}

template<class T>
inline
typename ICR::EnsembleLearning::RectifiedGaussian<T>::NP_t
ICR::EnsembleLearning::RectifiedGaussian<T>::CalcNP2Data(moments_parameter Mean,moments_parameter Precision)
{
  BOOST_ASSERT(Mean.size() == 2);
  BOOST_ASSERT(Precision.size() == 2);
  const data_t mean = Mean[0];
  const data_t precision = Precision[0];
  return NP_t(mean*precision, -0.5*precision);
}

template<class T> 
inline
typename ICR::EnsembleLearning::RectifiedGaussian<T>::NP_t
ICR::EnsembleLearning::RectifiedGaussian<T>::CalcNP2Parent(variable_parameter ParentA, 
					  deterministic_parameter Data, 
					  context_parameter C)
{
  //The moments forwarded from the Deterministic node.
  const moments_parameter FData = Data->GetForwardedMoments();
  const data_t fprec = -1.0/(FData[0]*FData[0]-FData[1]);
  const data_t fdata = FData[0];
  
  //Need to invert the expression based on the context.
  const SubContext<T> C0 = C[0];// The average means: < expr(x_i) >
  const SubContext<T> C1 = C[1];// The average of squares: < expr(x_i)^2 >

  //Find the placeholder associated with the parent that the message is to be sent to.
  const Placeholder<T>* P= C.Lookup(ParentA);
  
  //Invert the expression around P.
  std::pair<T,T> inv_op_data0 = P->Invert(fdata, C0);
  //The data returns in two components:
  //  The subtraction of the all the sums from the forwarded Data.
  const data_t unsummed0 =inv_op_data0.first;
  //  And the product thereafter
  const data_t factor0   =inv_op_data0.second;

  /* Example:  Three parents: X,Y,Z, and forwarded data D,
   *           Expression XY + Z = D.
   *
   *  If ParentA = Z:
   *      unsummed0 = D - XY; factor0 = 1;
   *  If ParentA = Y
   *      unsummed0 = D - Z;  factor0 = X;
   *
   *  That factor0 is X and not 1/X may seem surprising,
   *     but is in result of Miskin's thesis.
   *
   */

  //Now do the same thing for the average of the squares
  std::pair<T,T> inv_op_data1 = P->Invert(fdata, C1);
  const data_t unsummed1 =inv_op_data1.first;
  const data_t factor1   =inv_op_data1.second;

  //The usual (mean*precision, -0.5*precision) (with a scale factor)
  return NP_t(fprec*unsummed0*factor0, -0.5*factor1*fprec);
}


//Deterministic to Stock
template<class T> 
inline
typename ICR::EnsembleLearning::RectifiedGaussian<T>::NP_t
ICR::EnsembleLearning::RectifiedGaussian<T>::CalcNP2Deterministic(expression_parameter Expr,context_parameter C)
{
  //The context provides the Moments of every element in expression.
  const SubContext<T> C0 = C[0];//All the first moments  (the <x>'s of every element in expr)
  const SubContext<T> C1 = C[1];//The second moment (the <x^2> of every element of expression)
  
  //Precision is 1.0/ (<expr(x^2)> - <expr(x)>^2)
  const data_t prec = 1.0/(Expr->Evaluate(C1) - Expr->Evaluate(C0*C0) );
  
  // NP = [<expr(x)> * prec, -0.5*prec]
  return NP_t(  Expr->Evaluate(C0) *prec  , -0.5*prec);

}
