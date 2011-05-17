#pragma once
// #include "ExponentialModel.hpp"
// #include "ICA/NaturalParameters.hpp"
// #include "ICA/NaturalParameters.hpp"
// #include "ICA/Node.hpp"
#include "ICA/node/variable/Calculation.hpp"
#include "ICA/calculation_tree/Context.hpp"
#include "ICA/calculation_tree/Expression.hpp"

#include "Random.hpp"

namespace ICR{
  namespace ICA{
    


    template<class T=double>
    class Gaussian
    {
    public:
      
      typedef typename boost::call_traits< VariableNode<T>* const>::param_type
      variable_parameter;
      typedef typename boost::call_traits< VariableNode<T>* const>::value_type
      variable_t;

      typedef typename boost::call_traits< Moments<T> >::param_type
      moments_parameter;
      typedef typename boost::call_traits< Moments<T> >::const_reference
      moments_const_reference;
      typedef typename boost::call_traits< Moments<T> >::value_type
      moments_t;
      
      typedef typename boost::call_traits< DeterministicNode<Gaussian<T>, T>* >::param_type
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
      
      typedef typename boost::call_traits<const T>::value_type
      const_data_t;
      
      
      /** Calculate the natural log of the normalisation factor.
       *  The normalisation is the inverse of the partition factor.
       *  @param Mean The moments from the VariableNode that represents Mean
       *  @param Precision The moments from the VariableNode that represents the precision.
       *
       */
      static
      data_t
      CalcLogNorm(moments_parameter Mean,
		  moments_parameter Precision) ;

      
      /** Calculate the natural log of the normalisation factor.
       *  The normalisation is the inverse of the partition factor.
       *  @param NP The natural parameters from which to obtain the mean and precision.
       */
      static
      data_t
      CalcLogNorm(NP_parameter NP) ;
      

      /** Calculate the average of the log of the distribution.
       *  @param Mean The moments from the VariableNode that represents Mean
       *  @param Precision The moments from the VariableNode that represents the precision.
       *  @param Data The moments from the Child(or Data) Variable.
       *  @return The average of the log (<log Gaussian(mean,precision) >)
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
		 variable_parameter Precision);
      
      /** Calculate a random sample from a mixture.
       *  This is used for initialisation.
       *  @param mean_nodes The vector of variables that acts as the mean.
       *  @param precision_nodes The vector variable that acts as the precision.
       *  @param weights_node The vector variable holds the mixture weights.
       *  @return A random sample for the moments based on the means and precisions.
       */
      static
      moments_t
      CalcSample(std::vector<variable_t> mean_nodes,
		 std::vector<variable_t> precision_nodes,
		 variable_t weights_node);
      
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
      CalcNP2Data(moments_parameter Mean, moments_parameter Precision)  ;


      //Deterministic to Stock
      /** Calculate the Natural Parameters to go to
       *  each of the Parent Variables in a calculation.
       *  This is evaluated from all the other variables and the data.
       *  @param Parent The Node that the messag is to go to.
       *  @param Data The moments from the Child(or Data) Variable.
       *  @return The calculated NaturalParameters.  
       */
      static
      NP_t
      CalcNP2Parent( variable_parameter  Parent, 
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
      CalcNP2Deterministic(expression_parameter Expr,context_parameter C);

    private:
      //The Log norm is actually evaluated here.
      static
      T
      CalcLogNorm(const T& mean,const T& mean_squared,const T& precision) ;
    };

  }
}

template<class T>
inline
T
ICR::ICA::Gaussian<T>::CalcLogNorm(const T& mean, 
					const T& mean_squared, 
					const T& precision)  
{
  //The log norm
  return   0.5* ( std::log(precision/(2.0*M_PI)) - (precision) * (mean_squared)); 
}

template<class T>
inline
typename ICR::ICA::Gaussian<T>::data_t
ICR::ICA::Gaussian<T>::CalcLogNorm(moments_parameter Mean,
					moments_parameter Precision)  
{
  return  CalcLogNorm(Mean[0], Mean[1], Precision[0]);
}

template<class T>
inline
typename ICR::ICA::Gaussian<T>::data_t
ICR::ICA::Gaussian<T>::CalcLogNorm(NP_parameter NP)  
{
  //Obtain data from natural parameters.
  const_data_t precision    = -NP[1]*2.0;
  const_data_t mean         =  NP[0]/(precision);
  const_data_t mean_squared =  (mean)*(mean);
  return  CalcLogNorm(mean, mean_squared, precision);
}

    
template<class T>
inline
typename ICR::ICA::Gaussian<T>::data_t
ICR::ICA::Gaussian<T>::CalcAvLog(moments_parameter Mean,
				      moments_parameter Precision,
				      moments_parameter Data)
{
  const NP_t NPData = CalcNP2Data(Mean, Precision);
  return 
    NPData*Data + CalcLogNorm(Mean,Precision);
}
    
template<class T>
inline
typename ICR::ICA::Gaussian<T>::moments_t
ICR::ICA::Gaussian<T>::CalcSample(variable_parameter Mean,
				       variable_parameter Precision) 
{
  ICR::maths::rng* random = Random::Instance();
  const_data_t mean = Mean->GetMoments()[0];
  const_data_t prec = Precision->GetMoments()[0];
  const_data_t x=  random->gaussian(1.0/std::sqrt(prec),mean);
  return moments_t(x, x*x +1.0/prec);
}

template<class T>
typename ICR::ICA::Gaussian<T>::moments_t
ICR::ICA::Gaussian<T>::CalcSample(std::vector<variable_t> mean_nodes,
				       std::vector<variable_t> precision_nodes,
				       variable_t weights_node)
{
  moments_t weights = weights_node->GetMoments();
  std::vector<moments_t > mean(weights.size());
  std::vector<moments_t > precision(weights.size());

  PARALLEL_TRANSFORM(mean_nodes.begin(), mean_nodes.end(),
		     mean.begin(), 
		     boost::bind(&VariableNode<T>::GetMoments,
				 _1)
		     );
	
  PARALLEL_TRANSFORM(precision_nodes.begin(), precision_nodes.end(),
		     precision.begin(), 
		     boost::bind(&VariableNode<T>::GetMoments,
				 _1)
		     );
	
  moments_t AvMean(0,0);
  AvMean= PARALLEL_INNERPRODUCT(weights.begin(), weights.end(),mean.begin(), AvMean);
  moments_t AvPrec(0,0);
  AvPrec = PARALLEL_INNERPRODUCT(weights.begin(), weights.end(),precision.begin(), AvPrec);

  ICR::maths::rng* random = Random::Instance();
  const_data_t mean0 = AvMean[0];
  const_data_t prec0 = AvPrec[0];
  const_data_t x=  random->gaussian(1.0/std::sqrt(prec0),mean0);
  return moments_t(x, x*x +1.0/prec0);
}

template<class T>
inline
typename ICR::ICA::Gaussian<T>::moments_t
ICR::ICA::Gaussian<T>::CalcMoments(NP_parameter NP)
{
  //Get the data from the natural paremetes.
  const_data_t precision    = -NP[1]*2.0;
  const_data_t mean         =  NP[0]/(precision); 
  //Mean squared on HiddenVariable  is local and not averaged.
  const_data_t mean_squared =  (mean)*(mean);
  return  moments_t(mean, mean_squared + 1.0 /precision);
}

template<class T> 
inline
typename ICR::ICA::Gaussian<T>::NP_t
ICR::ICA::Gaussian<T>::CalcNP2Parent1(moments_parameter Precision,
					   moments_parameter Data)
{
  //Get the data from the moments.
  const_data_t precision = Precision[0];
  const_data_t data = Data[0];
  return NP_t(precision*data, -0.5*precision);
}


template<class T>
inline
typename ICR::ICA::Gaussian<T>::NP_t
ICR::ICA::Gaussian<T>::CalcNP2Parent2(moments_parameter Mean,
					   moments_parameter Data)
{
  const_data_t mean = Mean[0];
  const_data_t mean_squared = Mean[1];
  const_data_t data = Data[0];
  const_data_t data_squared = Data[1];
  //This is going to the precision - therefore we know that:
  BOOST_ASSERT(data_squared - 2*data*mean + mean_squared > 0);
  //(This should be a guareentee (should not require run-time check)
  // - if the above condition fails then it should fail noisily )
  return  NP_t(-0.5*(data_squared - 2*data*mean + mean_squared), 0.5 );
}

template<class T>
inline
typename ICR::ICA::Gaussian<T>::NP_t
ICR::ICA::Gaussian<T>::CalcNP2Data(moments_parameter Mean,
					moments_parameter Precision)
{
  //Get the data from the moments.
  const_data_t mean = Mean[0];
  const_data_t precision = Precision[0];
  return NP_t(mean*precision, -0.5*precision);
}

//Deterministic to Stock

template<class T> 
inline
typename ICR::ICA::Gaussian<T>::NP_t
ICR::ICA::Gaussian<T>::CalcNP2Parent(variable_parameter ParentA, 
					  deterministic_parameter Data, 
					  context_parameter C)
{
  //The moments forwarded from the Deterministic node.
  moments_const_reference FData = Data->GetForwardedMoments();
  const_data_t fprec = -1.0/(FData[0]*FData[0]-FData[1]);
  const_data_t fdata = FData[0];
  
  //Need to invert the expression based on the context.
  const SubContext<T> C0 = C[0]; // The average means: < expr(x_i) >
  const SubContext<T> C1 = C[1]; // The average of squares: < expr(x_i)^2 >

  //Find the placeholder associated with the parent that the message is to be sent to.
  const Placeholder<T>*  P= C.Lookup(ParentA);
  
  //Invert the expression around P.
  const std::pair<T,T> inv_op_data0 = P->Invert(fdata, C0);
  //The data returns in two components:
  //  The subtraction of the all the sums from the forwarded Data.
  const_data_t unsummed0 =inv_op_data0.first;  
  //  And the product thereafter
  const_data_t factor0   =inv_op_data0.second;

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
  const std::pair<T,T> inv_op_data1 = P->Invert(fdata, C1);
  const_data_t unsummed1 =inv_op_data1.first;
  const_data_t factor1   =inv_op_data1.second;

  //The usual (mean*precision, -0.5*precision) (with a scale factor)
  return NP_t(fprec*unsummed0*factor0, -0.5*factor1*fprec);
}


//Deterministic to Stock
template<class T> 
inline
typename ICR::ICA::Gaussian<T>::NP_t
ICR::ICA::Gaussian<T>::CalcNP2Deterministic(expression_parameter Expr,
						 context_parameter M)
{
  //The context provides the Moments of every element in expression.
  const SubContext<T> M0 = M[0];  //All the first moments  (the <x>'s of every element in expr)
  const SubContext<T> M1 = M[1];  //The second moment (the <x^2> of every element of expression)
  //Precision is 1.0/ (<expr(x^2)> - <expr(x)>^2)
  const_data_t prec = 1.0/(Expr->Evaluate(M1) - Expr->Evaluate(M0*M0) );
  // NP = [<expr(x)> * prec, -0.5*prec]
  return NP_t(  Expr->Evaluate(M0) *prec  , -0.5*prec);
}

