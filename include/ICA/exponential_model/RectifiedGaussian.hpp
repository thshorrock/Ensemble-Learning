#pragma once
// #include "ExponentialModel.hpp"
// #include "ICA/NaturalParameters.hpp"
// #include "ICA/NaturalParameters.hpp"
// #include "ICA/Node.hpp"
// #include "ICA/variable/DeterministicNode.hpp"
// #include "ICA/Deterministic/Context.hpp"
// #include "ICA/Deterministic/Expression.hpp"
#include "ICA/node/variable/Calculation.hpp"
#include "ICA/calculation_tree/Context.hpp"
#include "ICA/calculation_tree/Expression.hpp"

#include "Random.hpp"
#include <gsl/gsl_sf_erf.h>

namespace ICR{
  namespace ICA{
    


    template<class T=double>
    class RectifiedGaussianModel
    {
    public:
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
      CalcSample(const VariableNode<T>* Mean,const VariableNode<T>* Precision) 
      {

	ICR::maths::rng* random = Random::Instance();
	const T mean = Mean->GetMoments()[0];
	const T prec = Precision->GetMoments()[0];
	const T x=  random->gaussian_tail(1.0/std::sqrt(prec),mean);
	return Moments<T>(x, x*x +1.0/prec);
	
      } 
      static
      Moments<T>
      CalcSample(std::vector<VariableNode<T> *> m_mean_nodes,
		 std::vector<VariableNode<T> *> m_precision_nodes,
		 VariableNode<T> * m_weights_node)
	
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
	const T x=  random->gaussian_tail(1.0/std::sqrt(prec0),mean0);
	return Moments<T>(x, x*x +1.0/prec0);
      }
      
      static
      Moments<T>
      CalcMoments(const NaturalParameters<T>& NP)  ;
      
      // static
      // NaturalParameters<T>
      // CalcNaturalParameters(const Moments<T>& M);

      static
      NaturalParameters<T>
      CalcNP2Mean(const Moments<T>& Precision, const Moments<T>& Data)  ;

      static
      NaturalParameters<T>
      CalcNP2Data(const Moments<T>& Mean,const  Moments<T>& Precision)  ;

      static
      NaturalParameters<T>
      CalcNP2Precision(const Moments<T>& Mean,const  Moments<T>& Data)  ;
      


      //Deterministic to Stock
      
      template<template<class> class Op>
      static
      NaturalParameters<T>
      CalcNP2Parent(const VariableNode<T>* OtherParent,const DeterministicNode<RectifiedGaussianModel<T>, T>* Data);

      //Deterministic to Stock
      template<template<class> class Op>
      static
      NaturalParameters<T>
      CalcNP2Deterministic( VariableNode<T>* ParentA,const VariableNode<T>* ParentB);



      static
      NaturalParameters<T>
      CalcNP2Parent(const VariableNode<T>* ParentA, 
		    const DeterministicNode<RectifiedGaussianModel<T>, T>* Data, 
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
ICR::ICA::RectifiedGaussianModel<T>::CalcLogNorm(const T& mean, const T& mean_squared, const T& precision)  
{
  // std::cout<<"GAUSSIAN "<<std::endl;
  const T LN=    0.5* ( std::log(2.0*precision/(M_PI)) - (precision) * (mean_squared))
    - gsl_sf_log_erfc (-mean*std::sqrt(precision/2.0) ); //-mean*std::sqrt(precision/2.0) Typo in Minskins's thesis (... not the first)
  //std::cout<<"LN = "<<LN<<std::endl;
  if (fabs(LN) > 1e6) 
    {
      std::cout<<"LN = "<<LN<<" GP = "<< 0.5* ( std::log(2.0*precision/(M_PI)) - (precision) * (mean_squared))<<" log erfc = "<<- gsl_sf_log_erfc (-mean*std::sqrt(precision/2.0) )<<std::endl;
      std::cout<<"mean = "<<mean<<" prec = "<<precision<<" mean = "<<mean_squared<<std::endl;


    }
  return LN;
}

template<class T>
inline
T
ICR::ICA::RectifiedGaussianModel<T>::CalcLogNorm(const Moments<T>& Mean,const Moments<T>& Precision)  
{
  return  CalcLogNorm(Mean[0], Mean[1], Precision[0]);
}

template<class T>
inline
T
ICR::ICA::RectifiedGaussianModel<T>::CalcLogNorm(const NaturalParameters<T>& NP)  
{
  T precision    = -NP[1]*2.0;
  T mean         =  NP[0]/(precision);
  T mean_squared =  (mean)*(mean);
  return  CalcLogNorm(mean, mean_squared, precision);
}

template <class T>
struct Erfcx
{
  T
  operator()(T x)
  {
    return std::exp(x*x)*gsl_sf_erfc(x);
  }
};
template<class T>
inline
ICR::ICA::Moments<T>
ICR::ICA::RectifiedGaussianModel<T>::CalcMoments(const NaturalParameters<T>& NP)
{
  BOOST_ASSERT(NP.size() == 2);
  const T precision    = -NP[1]*2.0;
  const T mean         =  NP[0]/(precision);
  const T mean_squared =  (mean)*(mean); //when calculated on Variable this is local and not averaged.
  const T precision_squared =  (precision)*(precision); //when calculated on Variable this is local and not averaged.

  const T arg = -mean*std::sqrt(precision/2.0);
  //std::cout<<"arg = "<<arg<<std::endl;

  Moments<T> M(0,0);
  if (arg>10) //will start getting bad numerical errors
    {
      //rememer arg propto -mean
      //its tends to a gamma
      //This is the answer according to Miskin, but the second is clearly wrong (try plotting it)
      //return  Moments<T>(-1.0/(mean*precision), 2.0/(mean_squared*precision_squared));
      //From Closed-form approximations to the error and
      //complementary error functions and their applications in
      //atmospheric science we get (to one approximation better)
      
      return  Moments<T>(-1.0/(mean*precision) + 2.5/(mean*mean_squared*precision_squared),
			 ( mean_squared + 1.0 /(precision))*(1.0-1.0/sqrt(2))
			 + 2.5/(mean_squared*precision_squared));

    }
  if (arg<-10) //will start getting bad numerical errors
    {
      //its tends to a gaussin
      //std::cout<<"M1 "<<Moments<T>(mean, mean_squared + 1.0 /(precision))<<std::endl;

      return  Moments<T>(mean, mean_squared + 1.0 /(precision));
    }
  else
    {
      Erfcx<T> erfcx;
      const T iRpt = 1.0/(std::sqrt(M_PI*precision)*erfcx(arg));
      // std::cout<<"M3 "<<Moments<T>(mean + std::sqrt(2.0)*iRpt,
      // 				  mean_squared + 1.0 /(precision) + mean*iRpt)
      // 	       <<std::endl;

      return  Moments<T>(mean + std::sqrt(2.0)*iRpt,
			 mean_squared + 1.0 /(precision) + mean*iRpt);
    }
}

// template<class T>
// inline
// ICR::ICA::NaturalParameters<T>
// ICR::ICA::RectifiedGaussianModel<T>::CalcNaturalParameters(const Moments<T>& M)
// {

//   double precision = -1.0/(M[0]*M[0]-M[1]);
//   return  NaturalParameters<T>( M[0]*precision, 
// 				-0.5*precision);
// }

//Deterministic to Stock
template<class T> 
template<template<class> class Op>
inline
ICR::ICA::NaturalParameters<T>
ICR::ICA::RectifiedGaussianModel<T>::CalcNP2Parent(const VariableNode<T>* OtherParent,const DeterministicNode<RectifiedGaussianModel<T>, T>* Data)
{
  const Moments<T>& FData = Data->GetForwardedMoments();
  const Moments<T>& Other = OtherParent->GetMoments();
  

  //std::cout<<"Node Det = "<<Data<<std::endl;

  //calc the forwared precision and the forwared data;
  T fprec = -1.0/(FData[0]*FData[0]-FData[1]);
  T fdata = FData[0];
  T other = Other[0];
  
  T inv_op_data = Op<T>::Inverse(fdata,other); 
  
  return NaturalParameters<T>(fprec*inv_op_data, -0.5*fprec);
}

//Deterministic to Stock
template<class T> 
template<template<class> class Op>
inline
ICR::ICA::NaturalParameters<T>
ICR::ICA::RectifiedGaussianModel<T>::CalcNP2Deterministic(VariableNode<T>* ParentA,const VariableNode<T>* ParentB)
{
  //  std::cout<<"Node A = "<<ParentA<<std::endl;


  const Moments<T>&  a = ParentA->GetMoments();
  const Moments<T>&  b = ParentB->GetMoments();
  
  // std::cout<<"Moments A = "<<a<<std::endl;
  // std::cout<<"Moments B = "<<b<<std::endl;
  //calc the forwared precision and the forwared data;
  T prec = 1.0/( Op<T>::Forward( a[1], b[1]) - Op<T>::Forward(a[0]*a[0], b[0]*b[0]) );
  
  return NaturalParameters<T>(  Op<T>::Forward( a[0], b[0])*prec  , -0.5*prec);

}


template<class T> 
inline
ICR::ICA::NaturalParameters<T>
ICR::ICA::RectifiedGaussianModel<T>::CalcNP2Parent(const VariableNode<T>* ParentA, 
					  const DeterministicNode<RectifiedGaussianModel<T>, T>* Data, 
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
  const Placeholder<T>* P= C.Lookup(ParentA);
  
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
ICR::ICA::RectifiedGaussianModel<T>::CalcNP2Deterministic(const Expression<T>* Expr,const Context<T>& C)
{
  //  std::cout<<"Node A = "<<ParentA<<std::endl;
  
  SubContext<T> C0 = C[0];
  SubContext<T> C1 = C[1];
  // std::cout<<"C02Det = "<<C[0]<<std::endl;
  // std::cout<<"C02Det^2 = "<<C[0]*C[0]<<std::endl;
  // std::cout<<"C12Det = "<<C[1]<<std::endl;
  // std::cout<<"Expr->Evaluate(C1) = "<<Expr->Evaluate(C1)<<std::endl;
  // std::cout<<"Expr->Evaluate(C0) = "<<Expr->Evaluate(C0)<<std::endl;
  // std::cout<<"Expr->Evaluate(C0*C0) = "<<Expr->Evaluate(C0)<<std::endl;

  
  T prec = 1.0/(Expr->Evaluate(C1) - Expr->Evaluate(C0*C0) );
  
  return NaturalParameters<T>(  Expr->Evaluate(C0) *prec  , -0.5*prec);

}


template<class T> 
inline
ICR::ICA::NaturalParameters<T>
ICR::ICA::RectifiedGaussianModel<T>::CalcNP2Mean(const Moments<T>& Precision,const Moments<T>& Data)
{
  BOOST_ASSERT(Precision.size() == 2);
  BOOST_ASSERT(Data.size() == 2);
  T precision = Precision[0];
  T data = Data[0];
  return NaturalParameters<T>(precision*data, -0.5*precision);
}

template<class T>
inline
ICR::ICA::NaturalParameters<T>
ICR::ICA::RectifiedGaussianModel<T>::CalcNP2Data(const Moments<T>& Mean,const Moments<T>& Precision)
{
  BOOST_ASSERT(Mean.size() == 2);
  BOOST_ASSERT(Precision.size() == 2);
  T mean = Mean[0];
  T precision = Precision[0];
  return NaturalParameters<T>(mean*precision, -0.5*precision);
}

template<class T>
inline
ICR::ICA::NaturalParameters<T>
ICR::ICA::RectifiedGaussianModel<T>::CalcNP2Precision(const Moments<T>& Mean,const Moments<T>& Data)
{
  BOOST_ASSERT(Mean.size() == 2);
  BOOST_ASSERT(Data.size() == 2);
  T mean = Mean[0];
  T mean_squared = Mean[1];
  T data = Data[0];
  T data_squared = Data[1];
  if (data_squared - 2*data*mean + mean_squared<0) 
    {
      std::cout<<"MeanRG = "<<mean<<std::endl;
      std::cout<<"DataRG = "<<Data<<std::endl;
    }
  return  NaturalParameters<T>(-0.5*(data_squared - 2*data*mean + mean_squared), 0.5 );
}
    
template<class T>
inline
T
ICR::ICA::RectifiedGaussianModel<T>::CalcAvLog(const Moments<T>& Mean,const Moments<T>& Precision,const Moments<T>& Data)
{
  BOOST_ASSERT(Mean.size() == 2);
  BOOST_ASSERT(Precision.size() == 2);
  NaturalParameters<T> NPData = CalcNP2Data(Mean, Precision);
  return 
    NPData*Data + CalcLogNorm(Mean,Precision);
}
    
