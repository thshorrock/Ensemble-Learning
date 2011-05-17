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
    class GaussianModel
    {
    public:
      
      typedef typename boost::call_traits< VariableNode<T>* const>::param_type
      variable_parameter;
      typedef typename boost::call_traits< VariableNode<T>* const>::value_type
      variable_t;

      static
      T
      CalcLogNorm(const Moments<T>& Mean,const Moments<T>& Precision) ;

      static
      T
      CalcLogNorm(const NaturalParameters<T>& NP) ;
      
      static
      T
      CalcAvLog(const Moments<T>& Mean,const Moments<T>& Precision,const Moments<T>& Data) ;

      static
      Moments<T>
      CalcSample(variable_parameter Mean,variable_parameter Precision) 
      {
	ICR::maths::rng* random = Random::Instance();
	const T mean = Mean->GetMoments()[0];
	const T prec = Precision->GetMoments()[0];
	const T x=  random->gaussian(1.0/std::sqrt(prec),mean);
	return Moments<T>(x, x*x +1.0/prec);
      }
      
      static
      Moments<T>
      CalcSample(std::vector<variable_t> m_mean_nodes,
		 std::vector<variable_t> m_precision_nodes,
		 variable_t m_weights_node)
	
      {
	Moments<T> weights = m_weights_node->GetMoments();
	std::vector<Moments<T> > mean(weights.size());
	std::vector<Moments<T> > precision(weights.size());

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
	
	Moments<T> AvMean(0,0);
	AvMean= PARALLEL_INNERPRODUCT(weights.begin(), weights.end(),mean.begin(), AvMean);
	Moments<T> AvPrec(0,0);
	AvPrec = PARALLEL_INNERPRODUCT(weights.begin(), weights.end(),precision.begin(), AvPrec);

	ICR::maths::rng* random = Random::Instance();
	const T mean0 = AvMean[0];
	const T prec0 = AvPrec[0];
	const T x=  random->gaussian(1.0/std::sqrt(prec0),mean0);
	return Moments<T>(x, x*x +1.0/prec0);
      }
      
      static
      Moments<T>
      CalcMoments(const NaturalParameters<T>& NP)  ;
      

      static
      NaturalParameters<T>
      CalcNP2Parent1(const Moments<T>& Precision, const Moments<T>& Data)  ;

      static
      NaturalParameters<T>
      CalcNP2Data(const Moments<T>& Mean,const  Moments<T>& Precision)  ;

      static
      NaturalParameters<T>
      CalcNP2Parent2(const Moments<T>& Mean,const  Moments<T>& Data)  ;
      


      //Deterministic to Stock
      
      static
      NaturalParameters<T>
      CalcNP2Parent( variable_parameter  ParentA, 
		     DeterministicNode<GaussianModel<T>, T>* Data, 
		     const Context<T>& C);


      //Deterministic to Stock
      static
      NaturalParameters<T>
      CalcNP2Deterministic(const Expression<T>* Expr,const Context<T>& C);

    private:


      
      static
      T
      CalcLogNorm(const T& mean,const T& mean_squared,const T& precision) ;
    };

  }
}

template<class T>
inline
T
ICR::ICA::GaussianModel<T>::CalcLogNorm(const T& mean, const T& mean_squared, const T& precision)  
{
    // std::cout<<"GAUSSIAN "<<std::endl;
  return   0.5* ( std::log(precision/(2.0*M_PI)) - (precision) * (mean_squared)); 
}

template<class T>
inline
T
ICR::ICA::GaussianModel<T>::CalcLogNorm(const Moments<T>& Mean,const Moments<T>& Precision)  
{
  return  CalcLogNorm(Mean[0], Mean[1], Precision[0]);
}

template<class T>
inline
T
ICR::ICA::GaussianModel<T>::CalcLogNorm(const NaturalParameters<T>& NP)  
{
  T precision    = -NP[1]*2.0;
  T mean         =  NP[0]/(precision);
  T mean_squared =  (mean)*(mean);
  return  CalcLogNorm(mean, mean_squared, precision);
}

template<class T>
inline
ICR::ICA::Moments<T>
ICR::ICA::GaussianModel<T>::CalcMoments(const NaturalParameters<T>& NP)
{
  T precision    = -NP[1]*2.0;
  T mean         =  NP[0]/(precision);
  T mean_squared =  (mean)*(mean); //when calculated on Variable this is local and not averaged.

  // std::cout<<"GAUSSIAN MOMENTS = "<<Moments<T>(mean, mean_squared + 1.0 /(precision))<<std::endl;
  return  Moments<T>(mean, mean_squared + 1.0 /(precision));
}

// template<class T>
// inline
// ICR::ICA::NaturalParameters<T>
// ICR::ICA::GaussianModel<T>::CalcNaturalParameters(const Moments<T>& M)
// {

//   double precision = -1.0/(M[0]*M[0]-M[1]);
//   return  NaturalParameters<T>( M[0]*precision, 
// 				-0.5*precision);
// }


//Deterministic to Stock



template<class T> 
inline
ICR::ICA::NaturalParameters<T>
ICR::ICA::GaussianModel<T>::CalcNP2Parent(variable_parameter ParentA, 
					   DeterministicNode<GaussianModel<T>, T>* Data, 
					  const Context<T>& C)
{
  
  const Moments<T>& FData = Data->GetForwardedMoments();
  T fprec = -1.0/(FData[0]*FData[0]-FData[1]);
  T fdata = FData[0];
  
  //Invert the expression based on the context.
  SubContext<T> C0 = C[0];
  SubContext<T> C1 = C[1];

  // std::cout<<"C02P = "<<C[0]<<std::endl;
  // std::cout<<"C12P = "<<C[1]<<std::endl;
  const Placeholder<T>*  P= C.Lookup(ParentA);
  
  std::pair<T,T> inv_op_data0 = P->Invert(fdata, C0);
  T unsummed0 =inv_op_data0.first;
  T factor0   =inv_op_data0.second;

  std::pair<T,T> inv_op_data1 = P->Invert(fdata, C1);
  T unsummed1 =inv_op_data1.first;
  T factor1   =inv_op_data1.second;

  return NaturalParameters<T>(fprec*unsummed0*factor0, -0.5*factor1*fprec);
}


//Deterministic to Stock
template<class T> 
inline
ICR::ICA::NaturalParameters<T>
ICR::ICA::GaussianModel<T>::CalcNP2Deterministic(const Expression<T>* Expr,const Context<T>& M)
{
  //The context provides the Moments of every element in expression.
  SubContext<T> M0 = M[0];  //All the first moments  (the <x>'s of every element in expr)
  SubContext<T> M1 = M[1];  //The second moment (the <x^2> of every element of expression)
  //Precision is 1.0/ (<expr(x^2)> - <expr(x)>^2)
  T prec = 1.0/(Expr->Evaluate(M1) - Expr->Evaluate(M0*M0) );
  // NP = [<expr(x)> * prec, -0.5*prec]
  return NaturalParameters<T>(  Expr->Evaluate(M0) *prec  , -0.5*prec);
}


template<class T> 
inline
ICR::ICA::NaturalParameters<T>
ICR::ICA::GaussianModel<T>::CalcNP2Parent1(const Moments<T>& Precision,const Moments<T>& Data)
{
  T precision = Precision[0];
  T data = Data[0];
  return NaturalParameters<T>(precision*data, -0.5*precision);
}

template<class T>
inline
ICR::ICA::NaturalParameters<T>
ICR::ICA::GaussianModel<T>::CalcNP2Data(const Moments<T>& Mean,const Moments<T>& Precision)
{
  T mean = Mean[0];
  T precision = Precision[0];
  return NaturalParameters<T>(mean*precision, -0.5*precision);
}

template<class T>
inline
ICR::ICA::NaturalParameters<T>
ICR::ICA::GaussianModel<T>::CalcNP2Parent2(const Moments<T>& Mean,const Moments<T>& Data)
{
  T mean = Mean[0];
  T mean_squared = Mean[1];
  T data = Data[0];
  T data_squared = Data[1];
  if (data_squared - 2*data*mean + mean_squared<0) 
    {
      std::cout<<"MeanG = "<<mean<<std::endl;
      std::cout<<"DataG = "<<Data<<std::endl;
    }
  return  NaturalParameters<T>(-0.5*(data_squared - 2*data*mean + mean_squared), 0.5 );
}
    
template<class T>
inline
T
ICR::ICA::GaussianModel<T>::CalcAvLog(const Moments<T>& Mean,const Moments<T>& Precision,const Moments<T>& Data)
{
  NaturalParameters<T> NPData = CalcNP2Data(Mean, Precision);
  return 
    NPData*Data + CalcLogNorm(Mean,Precision);
}
    
