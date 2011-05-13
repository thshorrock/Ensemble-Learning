#pragma once

#include "ExponentialModel.hpp"

#include <gsl/gsl_sf_psi.h> //for Digamma function gsl_sf_psi
#include <gsl/gsl_sf_gamma.h> //for gamma function

#include <boost/assert.hpp> 
namespace ICR{
  namespace ICA{

    template<class T=double>
    class GammaModel //: public ExponentialModel<T>
    {
    public:
      static
      T
      CalcLogNorm(const Moments<T>& Shape, const Moments<T> IScale);
      
      static
      T
      CalcLogNorm(const NaturalParameters<T>& NP);
      
      static
      T
      CalcAvLog(const Moments<T>& Shape,const Moments<T>& Iscale,const Moments<T>& Data) ;

      static
      Moments<T>
      CalcSample(const VariableNode<T>* Shape,const VariableNode<T>* IScale) 
      {
	ICR::maths::rng* random = Random::Instance();
	const T shape  = Shape->GetMoments()[0];
	const T iscale = IScale->GetMoments()[0];

	  
	const T x=  random->gamma(shape,1.0/iscale) + 1e-6; //cannot be zero so add something sma
	
	std::cout<<"x = "<<x<<std::endl;
	 Moments<double> M(x,std::log(x));
	 
	 std::cout<<"Init = "<<M<<std::endl;
	 return M;
      }
      
      static
      Moments<T>
      CalcMoments(const NaturalParameters<T>& NP);

      // static
      // NaturalParameters<T>
      // CalcNaturalParameters(const Moments<T>& M);

      static
      NaturalParameters<T>
      CalcNP2IScale(const Moments<T> Shape,const Moments<T>& Data);

      static
      NaturalParameters<T>
      CalcNP2Data(const Moments<T>& Shape, const Moments<T> IScale);


    private:
      static
      T
      CalcLogNorm(const T& shape, const T& iscale) ;
    };

  }
}

template<class T>
inline
T 
ICR::ICA::GammaModel<T>::CalcLogNorm(const T& shape, const T& iscale)  
{
   // std::cout<<"GAMMA "<<std::endl;
  BOOST_ASSERT(iscale>0);
  BOOST_ASSERT(shape!=0);
  //  std::cout<<"shape = "<<shape<<std::endl;

  return shape * std::log(iscale) - gsl_sf_lngamma(shape) ;
}

template<class T>
inline
T
ICR::ICA::GammaModel<T>::CalcLogNorm(const Moments<T>& Shape, const Moments<T> IScale)
{
  BOOST_ASSERT(IScale[0]>0);
  return CalcLogNorm(Shape[0], IScale[0]);
}
      
template<class T>
inline
T
ICR::ICA::GammaModel<T>::CalcLogNorm(const NaturalParameters<T>& NP)
{
  T shape  = NP[1]+1;
  T iscale = -NP[0];
  BOOST_ASSERT(iscale>0);
  return CalcLogNorm(shape, iscale);
}

template<class T>
inline
ICR::ICA::Moments<T>
ICR::ICA::GammaModel<T>::CalcMoments(const NaturalParameters<T>& NP)
{
  T shape    = NP[1]+1.0;
  T iscale    = -NP[0];
  // std::cout<<"shape = "<<m_shape<<" scale = "<<m_iscale<<std::endl;
  BOOST_ASSERT(iscale>0);
 
  //std::cout<<"GAMMA MOMENTS = "<<m<<std::endl;


  return Moments<T>(shape/iscale,
		    gsl_sf_psi(shape) - std::log(iscale) 
		    );
}

// template<class T>
// inline
// ICR::ICA::NaturalParameters<T>
// ICR::ICA::Model<T>::CalcNaturalParameters(const Moments<T>& M)
// {

//   double precision = -1.0/(M[0]*M[0]-M[1]);
//   return  NaturalParameters<T>( M[0]*precision, 
// 				-0.5*precision);
// }

template<class T>
inline 
ICR::ICA::NaturalParameters<T>
ICR::ICA::GammaModel<T>::CalcNP2IScale(const Moments<T> Shape,
					    const Moments<T>& Data)
{
  return NaturalParameters<T>(-Data[0], Shape[0]);
}

template<class T>
inline 
ICR::ICA::NaturalParameters<T>
ICR::ICA::GammaModel<T>::CalcNP2Data(const Moments<T>& Shape, 
					  const Moments<T> IScale)
{
  return  NaturalParameters<T>(-IScale[0], Shape[0] - 1);
}

      
template<class T>
inline
T
ICR::ICA::GammaModel<T>::CalcAvLog(const Moments<T>& Shape,const Moments<T>& IScale,const Moments<T>& Data)
{
  return CalcNP2Data(Shape,IScale)*Data +CalcLogNorm(Shape,IScale);
}
