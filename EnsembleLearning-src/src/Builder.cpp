
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



//builder
#include "EnsembleLearning/Builder.hpp"
//factors
#include "EnsembleLearning/node/factor/Calculation.hpp"
#include "EnsembleLearning/node/factor/Factor.hpp"
#include "EnsembleLearning/node/factor/Mixture.hpp"
//nodes
#include "EnsembleLearning/node/variable/Hidden.hpp"
#include "EnsembleLearning/node/variable/Observed.hpp"
#include "EnsembleLearning/node/variable/Calculation.hpp"
//algorithms
#include "EnsembleLearning/detail/parallel_algorithms.hpp"
//messages
#include "EnsembleLearning/message/Coster.hpp"
//typelists
#include "EnsembleLearning/detail/TypeList.hpp"
//stream
#include <fstream>

template<class T>
ICR::EnsembleLearning::Builder<T>::Builder(const std::string& cost_file)
  : m_PrevCost(-1.0/0.0), // minus infty
    m_Factors(),
    m_Nodes(),
    m_initialised(false),
    m_data_nodes(0),
    m_cost_file(cost_file)
{
  //clear it
  if (m_cost_file != "") { 
    std::ofstream CostFile(m_cost_file.c_str());
    CostFile<<"#Data";
  }
  // CostFile = new std::ofstream("Cost.txt");
  // CostFile.close();
}  

template<class T>
void 
ICR::EnsembleLearning::Builder<T>::set_cost_file(const std::string& cost_file)
{
  m_cost_file = cost_file;
  if (m_cost_file != "") { 
    std::ofstream CostFile(m_cost_file.c_str());
    CostFile<<"#Data";
  }
  
}

template<class T>
ICR::EnsembleLearning::Builder<T>::~Builder()
{
}

template<class T>
typename ICR::EnsembleLearning::Builder<T>::WeightsNode
ICR::EnsembleLearning::Builder<T>::weights()
{
  typedef DirichletConstType prior_t;
  typedef ICR::EnsembleLearning::NoSecondParent blank;
  typedef DirichletType child_t;
  typedef detail::Factor<Dirichlet,T,prior_t,blank,child_t,ENSEMBLE_LEARNING_COMPONENTS> Factor_t;
  

  boost::shared_ptr<DirichletConstType > nDirichletPrior(new DirichletConstType(1.0));
  boost::shared_ptr<DirichletType >      nDirichlet(new DirichletType());
  
  boost::shared_ptr<Factor_t>    
    DirichletF(new Factor_t(nDirichletPrior.get(), nDirichlet.get()));
  
  m_Nodes.push_back(nDirichletPrior);
  m_Nodes.push_back(nDirichlet);
  
  m_Factors.push_back(DirichletF);

  return nDirichlet.get();

	
}

// template<class T>
// typename ICR::EnsembleLearning::Builder<T>::RectifiedGaussianNode
// ICR::EnsembleLearning::Builder<T>::rectified_gaussian(Variable Mean,Variable  Precision)
// {
	
//   boost::shared_ptr<RectifiedGaussianType > RectifiedGaussian(new RectifiedGaussianType());
//   boost::shared_ptr<RectifiedGaussianFactor > RectifiedGaussianF(new RectifiedGaussianFactor(Mean, Precision, RectifiedGaussian.get()));
	
//   m_Factors.push_back(RectifiedGaussianF);
//   m_Nodes.push_back(RectifiedGaussian);
	
//   return RectifiedGaussian.get();
// }

// template<class T>
// typename ICR::EnsembleLearning::Builder<T>::RectifiedGaussianNode
// ICR::EnsembleLearning::Builder<T>::rectified_gaussian(GaussianNode Mean, const T& precision)
// {
	
//   boost::shared_ptr<GammaConstType > Precision(new GammaConstType(precision));
//   m_Nodes.push_back(Precision);
//   return rectified_gaussian(Mean,Precision.get());
// }
    
// template<class T>  
// typename ICR::EnsembleLearning::Builder<T>::RectifiedGaussianNode
// ICR::EnsembleLearning::Builder<T>::rectified_gaussian(const T& mean, GammaNode Precision)
// {
	
//   boost::shared_ptr<GaussianConstType > Mean(new GaussianConstType(mean));
//   m_Nodes.push_back(Mean);
//   return rectified_gaussian(Mean.get(),Precision);
// }
  
// template<class T>
// typename ICR::EnsembleLearning::Builder<T>::RectifiedGaussianNode
// ICR::EnsembleLearning::Builder<T>::rectified_gaussian(const T& mean, const T& precision)
// {
	
//   boost::shared_ptr<GaussianConstType > Mean(new GaussianConstType(mean));
//   boost::shared_ptr<GammaConstType > Precision(new GammaConstType(precision));
	
//   m_Nodes.push_back(Mean);
//   m_Nodes.push_back(Precision);
	
//   return rectified_gaussian(Mean.get(),Precision.get());
// }
      
	


template<class T>	
typename ICR::EnsembleLearning::Builder<T>::RectifiedGaussianNode
ICR::EnsembleLearning::Builder<T>::rectified_gaussian_mixture(std::vector<Variable>& vMean, std::vector<Variable>& vPrecision, WeightsNode Weights)
{
  const size_t number = Weights->size();
	
  typedef DirichletType prior_t;
  typedef ICR::EnsembleLearning::NoSecondParent blank;
  typedef CatagoryType child_t;
  typedef detail::Factor<Discrete,T,prior_t,blank,child_t,ENSEMBLE_LEARNING_COMPONENTS> Factor_t;

  boost::shared_ptr<CatagoryType>         Catagory(new CatagoryType());
  boost::shared_ptr< Factor_t>      CatagoryF(new Factor_t(Weights, Catagory.get()));
  m_Nodes.push_back(Catagory);
  m_Factors.push_back(CatagoryF);
	

  boost::shared_ptr<RectifiedGaussianType> Child(new RectifiedGaussianType());

  m_Nodes.push_back(Child);
	
  // boost::shared_ptr<RectifiedGaussianMixtureFactor> MixtureF(new RectifiedGaussianMixtureFactor(vMean, vPrecision, Catagory.get() , Child.get()));
	
  // m_Factors.push_back(MixtureF);
  return Child.get();
}

template<class T>
typename ICR::EnsembleLearning::Builder<T>::GammaNode
ICR::EnsembleLearning::Builder<T>::gamma(const T& shape, const T& iscale)
{
  if ((shape<=0) || (iscale <=0)) {
    std::cout<<"Shape and iscale must be greater than zero"<<std::endl;
    throw("EXITING");
  }

  typedef detail::Factor<Gamma,T, NormalConstType,GammaConstType,GammaType> Factor_t;
  BOOST_ASSERT(shape>0);
  BOOST_ASSERT(iscale>0);
	
  boost::shared_ptr<NormalConstType > Shape(new NormalConstType(shape));
  boost::shared_ptr<GammaConstType > IScale(new GammaConstType(iscale));
	
  boost::shared_ptr<GammaType > Gamma(new GammaType());
  boost::shared_ptr<Factor_t > GammaF(new Factor_t(Shape.get(), IScale.get(), Gamma.get()));
	
  m_Factors.push_back(GammaF);
  m_Nodes.push_back(Shape);
  m_Nodes.push_back(IScale);
  m_Nodes.push_back(Gamma);
	
  return Gamma.get();
}


template<class T>
typename ICR::EnsembleLearning::Builder<T>::GaussianConstNode
ICR::EnsembleLearning::Builder<T>::gaussian_const(const T value)
{
  boost::shared_ptr<GaussianConstType > Const(new GaussianConstType(value));
  m_Nodes.push_back(Const);
	
  return Const.get();
}
      
template<class T>
typename ICR::EnsembleLearning::Builder<T>::GammaConstNode
ICR::EnsembleLearning::Builder<T>::gamma_const(const T value)
{
  if (value<=0) {
    std::cout<<"Shape and iscale must be greater than zero"<<std::endl;
    throw("EXITING");
  }

  BOOST_ASSERT(value>0);
	
  boost::shared_ptr<GammaConstType > Const(new GammaConstType(value));
  m_Nodes.push_back(Const);
	
  return Const.get();
}
      
template<class T>
typename ICR::EnsembleLearning::Builder<T>::GaussianDataNode
ICR::EnsembleLearning::Builder<T>::gaussian_data(const T data)
{
	
  boost::shared_ptr<GaussianDataType > Data(new GaussianDataType(data));
  ++m_data_nodes;
  m_Nodes.push_back(Data);
	
  return Data.get();
}
      

template<class T>
typename ICR::EnsembleLearning::Builder<T>::GammaDataNode
ICR::EnsembleLearning::Builder<T>::gamma_data(const T data)
{
	
  boost::shared_ptr<GammaDataType > Data(new GammaDataType(data));
  ++m_data_nodes;
  m_Nodes.push_back(Data);
	
  return Data.get();
}

template<class T>
typename ICR::EnsembleLearning::Builder<T>::GaussianResultNode
ICR::EnsembleLearning::Builder<T>::calc_gaussian(Expression<T>* Expr,  Context<T>& context)
{
  // typedef DeterministicNode<Gaussian<T>, T> child_t;

  // boost::shared_ptr<child_t > Child(new child_t ());
  // boost::shared_ptr<DeterministicFactor> ChildF
  //   (new DeterministicFactor(Expr, context,Child.get()));
	
  // m_Nodes.push_back(Child);
  // m_Factors.push_back(ChildF);
  // return Child.get();
}


template<class T>
void 
ICR::EnsembleLearning::Builder<T>::join(T& shape, GammaNode IScale ,  const T& data)
{

  boost::shared_ptr<GammaDataType > Data(new GammaDataType(data));
  ++m_data_nodes;
  boost::shared_ptr<NormalConstType > Shape(new NormalConstType(shape));
	
  
  typedef detail::Factor<Gamma,T,NormalConstType,GammaType,GammaDataType> Factor_t;
  boost::shared_ptr<Factor_t> GammaF (new Factor_t(Shape.get(), IScale, Data.get()));
	 
  m_Factors.push_back(GammaF);
  m_Nodes.push_back(Data);
  m_Nodes.push_back(Shape);
  
	
}
template<class T>
void
ICR::EnsembleLearning::Builder<T>::join(T& shape, T& iscale, GammaDataNode Data  )
{

  typedef ObservedNode<Gamma, T, detail::TypeList::id >    GammaConstType;
  boost::shared_ptr<NormalConstType > Shape(new NormalConstType(shape));
  boost::shared_ptr<GammaConstType > IScale(new GammaConstType(iscale));
  typedef detail::Factor<Gamma,T,NormalConstType,GammaConstType,GammaDataType> Factor_t;
	
  boost::shared_ptr<Factor_t > GammaF (new Factor_t(Shape.get(), IScale.get(), Data));
	 
  m_Factors.push_back(GammaF);
  m_Nodes.push_back(Shape);
  m_Nodes.push_back(IScale);
	
}

template<class T>
void 
ICR::EnsembleLearning::Builder<T>::join(T& shape, GammaNode IScale, GammaDataNode Data  )
{
  typedef detail::Factor<Gamma,T,NormalConstType,GammaType,GammaDataType> Factor_t;

  boost::shared_ptr<NormalConstType > Shape (new NormalConstType(shape));
  boost::shared_ptr<Factor_t > GammaF (new Factor_t(Shape.get(), IScale, Data));
	 
  m_Factors.push_back(GammaF);
  m_Nodes.push_back(Shape);
}

template<class T>
void 
ICR::EnsembleLearning::Builder<T>::join(T& mean, T& precision, GaussianDataNode Data  )
{

  typedef ObservedNode<Gaussian, T, detail::TypeList::id > GaussianConstType   ;
  boost::shared_ptr<GaussianConstType > Mean(new GaussianConstType(mean));
  boost::shared_ptr<GammaConstType > Precision(new GammaConstType(precision));
	
  typedef detail::Factor<Gaussian,T,GaussianConstType,GammaConstType,GaussianDataType> Factor_t;
  boost::shared_ptr<Factor_t > GaussianF(new  Factor_t(Mean.get(), Precision.get(), Data));
	
  m_Factors.push_back(GaussianF);
  m_Nodes.push_back(Mean);
  m_Nodes.push_back(Precision);
}
   

// template<class T>
// void 
// ICR::EnsembleLearning::Builder<T>::join(Variable Mean, GammaNode& Precision, GammaNode& Child  )
// {
//   boost::shared_ptr<GammaFactor> GammaF(new GammaFactor(Mean,Precision,Child));
//   m_Factors.push_back(GammaF);
// }

   
// template<class T>   
// void 
// ICR::EnsembleLearning::Builder<T>::join(Variable Mean, GammaNode& Precision, GaussianNode& Child  )
// {
//   boost::shared_ptr<GaussianFactor> GaussianF(new GaussianFactor(Mean,Precision,Child));
//   m_Factors.push_back(GaussianF);
// }

      


template<class T>
size_t
ICR::EnsembleLearning::Builder<T>::number_of_nodes() const
{
  return m_Nodes.size();
}

template<class T>
size_t
ICR::EnsembleLearning::Builder<T>::number_of_factors() const
{
  return m_Factors.size();
}


template<class T>
void
ICR::EnsembleLearning::Builder<T>::perturb()
{
  PARALLEL_FOREACH(m_Nodes.begin(), m_Nodes.end(),
		   boost::bind(&VariableNode<T>::InitialiseMoments, _1)
		   );
}


template<class T>
double
ICR::EnsembleLearning::Builder<T>::iterate()
{
  // Initialise();
  if (m_Factors.size() == 0)
    {
      std::cout<<"No graph has been built, cannot iterate"<<std::endl;
      return -1.0/0.0;
    }
	
  Coster Cost;
  for(size_t i=0;i<1;++i){
    {
      // std::cout<<"ITERATE"<<std::endl;
      Cost = 0;

      // PARALLEL_FOREACH(m_Factors.begin(), m_Factors.end(),
      // 		     boost::bind(&FactorNode<T>::Iterate, _1 )
      // 		     );
	    
      // std::cout<<"ITERATE Nodes"<<std::endl;
      PARALLEL_FOREACH(m_Nodes.begin(), m_Nodes.end(),
		       boost::bind(&VariableNode<T>::Iterate, _1, boost::ref(Cost))
		       );

    }
  }
  return Cost;

}

template<class T>
bool
ICR::EnsembleLearning::Builder<T>::run(const double& epsilon, const size_t& max_iterations , size_t skip)
{
  //std::cout<<"initialising"<<std::endl;
  for(size_t i=0;i<skip;++i){
    iterate();
  }
	
  for(size_t i=0;i<max_iterations;++i){
    double Cost = iterate()/m_data_nodes;
	  
    if (HasConverged(Cost, epsilon))
      return true;
  }

  return false;

}
    
     
      
template<class T>
bool
ICR::EnsembleLearning::Builder<T>::HasConverged(const T Cost, const T epsilon)
{
#pragma omp critical
  {
    std::cout<<"COST = "<<Cost<<"\t% difference = "<<100.0*(((Cost-m_PrevCost)/std::fabs(Cost)))<<std::endl;
    if (m_cost_file != "") { 
      std::ofstream CostFile(m_cost_file.c_str(),std::ios_base::app);
      CostFile<<Cost<<"\n";
    }
  }
  if (100.0*std::fabs((Cost - m_PrevCost)/Cost)<epsilon) 
    return true;
	
  m_PrevCost = Cost;
  return false;
}
    


template class ICR::EnsembleLearning::Builder<double>;
template class ICR::EnsembleLearning::Builder<float>;

