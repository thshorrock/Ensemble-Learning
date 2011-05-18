#include "ICA.hpp"

// #include "ICA/node/variable/Constant.hpp"
// #include "ICA/node/variable/Data.hpp"
#include "ICA/node/variable/Hidden.hpp"
// #include "ICA/node/variable/Dirichlet.hpp"
#include "ICA/node/variable/Calculation.hpp"
// #include "ICA/variable/Forwarding.hpp"

#include "ICA/node/factor/Calculation.hpp"
#include "ICA/node/factor/Factor.hpp"
#include "ICA/node/factor/Mixture.hpp"
// #include "ICA/factor/DeterministicFactor.hpp"

#include "ICA/exponential_model/RectifiedGaussian.hpp"

#include "ICA/detail/parallel_algorithms.hpp"

// typedef ICR::ICA::Factor<ICR::ICA::RectifiedGaussianModel<T> >      RectifiedGaussianFactor;
// typedef ICR::ICA::Factor<ICR::ICA::GaussianModel<T> >      GaussianFactor;
// typedef ICR::ICA::Factor<ICR::ICA::GammaModel<T> >         GammaFactor;
// typedef ICR::ICA::Factor<ICR::ICA::DirichletModel<T> >     DirichletFactor;
// typedef ICR::ICA::Factor<ICR::ICA::DiscreteModel<T> >      DiscreteFactor;
// typedef ICR::ICA::Mixture<ICR::ICA::RectifiedGaussianModel<T> >     RectifiedGaussianMixtureFactor;
// typedef ICR::ICA::Mixture<ICR::ICA::GaussianModel<T> >     GaussianMixtureFactor;

// typedef ICR::ICA::Details::CalcGaussianFactor<ICR::ICA::GaussianModel, T >  CalcGaussianFactor;

// // typedef DeterministicFactor<GaussianModel,Det::Add, T >     AddFactor;
// // typedef DeterministicFactor<GaussianModel,Det::Multiply, T >     MultiplyFactor;

// typedef ICR::ICA::HiddenNode<ICR::ICA::GaussianModel<T> >      GaussianType;
// typedef ICR::ICA::HiddenNode<ICR::ICA::RectifiedGaussianModel<T> >      RectifiedGaussianType;
// typedef ICR::ICA::HiddenNode<ICR::ICA::GammaModel<T> >         GammaType;
// typedef ICR::ICA::HiddenNode<ICR::ICA::DirichletModel<T> >     DirichletType;
// typedef ICR::ICA::DataNode<ICR::ICA::GaussianConstant<T> >     GaussianDataType;
// typedef ICR::ICA::DataNode<ICR::ICA::GammaConstant<T> >        GammaDataType;
// typedef ICR::ICA::ConstantNode<ICR::ICA::GaussianConstant<T> > GaussianConstType;
// typedef ICR::ICA::ConstantNode<ICR::ICA::GammaConstant<T> >    GammaConstType;
// typedef ICR::ICA::ConstantNode<ICR::ICA::NormalConstant<T> >   NormalConstType;
// typedef ICR::ICA::ConstantNode<ICR::ICA::DirichletConstant<T> >  DirichletConstType;
// typedef ICR::ICA::DeterministicNode<ICR::ICA::GaussianModel<T>, T>    GaussianResultType;

      
// typedef ICR::ICA::HiddenNode<DirichletModel<T> >     WeightsType;
// typedef ICR::ICA::HiddenNode<DiscreteModel<T> >      CatagoryType;
//      typedef ForwardingNode<GaussianModel<T> >   GaussianMixtureType;

      
template<class T>
ICR::ICA::Builder<T>::Builder()
  : m_PrevCost(-1.0/0.0), // minus infty
    m_Factors(),
    m_Nodes(),
    m_initialised(false)
{
  //clear it
  std::ofstream CostFile("Cost.txt");
  CostFile<<"#Data";
  // CostFile = new std::ofstream("Cost.txt");
  // CostFile.close();
}

template<class T>
ICR::ICA::Builder<T>::~Builder()
{
}

template<class T>
typename ICR::ICA::Builder<T>::WeightsNode
ICR::ICA::Builder<T>::weights(const size_t size)
{
  boost::shared_ptr<DirichletConstType > DirichletPrior(new DirichletConstType(size,1.0));
  boost::shared_ptr<DirichletType >      Dirichlet(new DirichletType(size));
  boost::shared_ptr<DirichletFactor >    DirichletF(new DirichletFactor(DirichletPrior.get(), Dirichlet.get()));
  m_Nodes.push_back(DirichletPrior);
  m_Nodes.push_back(Dirichlet);

  m_Factors.push_back(DirichletF);

  return Dirichlet.get();

	
}

template<class T>
typename ICR::ICA::Builder<T>::GaussianNode
ICR::ICA::Builder<T>::gaussian(Variable Mean, Variable Precision)
{
	
  boost::shared_ptr<GaussianType > Gaussian(new GaussianType());
  boost::shared_ptr<GaussianFactor > GaussianF(new GaussianFactor(Mean, Precision, Gaussian.get()));
	
  m_Factors.push_back(GaussianF);
  m_Nodes.push_back(Gaussian);
	
  return Gaussian.get();
}

template<class T>
typename ICR::ICA::Builder<T>::GaussianNode
ICR::ICA::Builder<T>::gaussian(GaussianNode Mean, const T& precision)
{
	
  boost::shared_ptr<GammaConstType > Precision(new GammaConstType(precision));
  m_Nodes.push_back(Precision);
  return gaussian(Mean,Precision.get());
}
    
template<class T>  
typename ICR::ICA::Builder<T>::GaussianNode
ICR::ICA::Builder<T>::gaussian(const T& mean, GammaNode Precision)
{
	
  boost::shared_ptr<GaussianConstType > Mean(new GaussianConstType(mean));
  m_Nodes.push_back(Mean);
  return gaussian(Mean.get(),Precision);
}
  
template<class T>
typename ICR::ICA::Builder<T>::GaussianNode
ICR::ICA::Builder<T>::gaussian(const T& mean, const T& precision)
{
	
  boost::shared_ptr<GaussianConstType > Mean(new GaussianConstType(mean));
  boost::shared_ptr<GammaConstType > Precision(new GammaConstType(precision));
	
  m_Nodes.push_back(Mean);
  m_Nodes.push_back(Precision);
	
  return gaussian(Mean.get(),Precision.get());
}
      

template<class T>
typename ICR::ICA::Builder<T>::RectifiedGaussianNode
ICR::ICA::Builder<T>::rectified_gaussian(Variable Mean,Variable  Precision)
{
	
  boost::shared_ptr<RectifiedGaussianType > RectifiedGaussian(new RectifiedGaussianType());
  boost::shared_ptr<RectifiedGaussianFactor > RectifiedGaussianF(new RectifiedGaussianFactor(Mean, Precision, RectifiedGaussian.get()));
	
  m_Factors.push_back(RectifiedGaussianF);
  m_Nodes.push_back(RectifiedGaussian);
	
  return RectifiedGaussian.get();
}

template<class T>
typename ICR::ICA::Builder<T>::RectifiedGaussianNode
ICR::ICA::Builder<T>::rectified_gaussian(GaussianNode Mean, const T& precision)
{
	
  boost::shared_ptr<GammaConstType > Precision(new GammaConstType(precision));
  m_Nodes.push_back(Precision);
  return rectified_gaussian(Mean,Precision.get());
}
    
template<class T>  
typename ICR::ICA::Builder<T>::RectifiedGaussianNode
ICR::ICA::Builder<T>::rectified_gaussian(const T& mean, GammaNode Precision)
{
	
  boost::shared_ptr<GaussianConstType > Mean(new GaussianConstType(mean));
  m_Nodes.push_back(Mean);
  return rectified_gaussian(Mean.get(),Precision);
}
  
template<class T>
typename ICR::ICA::Builder<T>::RectifiedGaussianNode
ICR::ICA::Builder<T>::rectified_gaussian(const T& mean, const T& precision)
{
	
  boost::shared_ptr<GaussianConstType > Mean(new GaussianConstType(mean));
  boost::shared_ptr<GammaConstType > Precision(new GammaConstType(precision));
	
  m_Nodes.push_back(Mean);
  m_Nodes.push_back(Precision);
	
  return rectified_gaussian(Mean.get(),Precision.get());
}
      
	

template<class T>
typename ICR::ICA::Builder<T>::GaussianNode
ICR::ICA::Builder<T>::gaussian_mixture(std::vector<Variable>& vMean, std::vector<Variable>& vPrecision, WeightsNode Weights)
{
  const size_t number = Weights->size();
	
  boost::shared_ptr<CatagoryType>         Catagory(new CatagoryType(number));
  boost::shared_ptr<DiscreteFactor >      CatagoryF(new DiscreteFactor(Weights, Catagory.get()));
  m_Nodes.push_back(Catagory);
  m_Factors.push_back(CatagoryF);
	

  boost::shared_ptr<GaussianType> Child(new GaussianType());

  m_Nodes.push_back(Child);
	
  boost::shared_ptr<GaussianMixtureFactor> MixtureF(new GaussianMixtureFactor(vMean, vPrecision, Catagory.get() , Child.get()));
	
  m_Factors.push_back(MixtureF);
  return Child.get();
}

template<class T>	
typename ICR::ICA::Builder<T>::RectifiedGaussianNode
ICR::ICA::Builder<T>::rectified_gaussian_mixture(std::vector<Variable>& vMean, std::vector<Variable>& vPrecision, WeightsNode Weights)
{
  const size_t number = Weights->size();
	
  boost::shared_ptr<CatagoryType>         Catagory(new CatagoryType(number));
  boost::shared_ptr<DiscreteFactor >      CatagoryF(new DiscreteFactor(Weights, Catagory.get()));
  m_Nodes.push_back(Catagory);
  m_Factors.push_back(CatagoryF);
	

  boost::shared_ptr<RectifiedGaussianType> Child(new RectifiedGaussianType());

  m_Nodes.push_back(Child);
	
  boost::shared_ptr<RectifiedGaussianMixtureFactor> MixtureF(new RectifiedGaussianMixtureFactor(vMean, vPrecision, Catagory.get() , Child.get()));
	
  m_Factors.push_back(MixtureF);
  return Child.get();
}

template<class T>
typename ICR::ICA::Builder<T>::GammaNode
ICR::ICA::Builder<T>::gamma(const T& shape, const T& iscale)
{
  if ((shape<=0) || (iscale <=0)) {
    std::cout<<"Shape and iscale must be greater than zero"<<std::endl;
    throw("EXITING");
  }

  BOOST_ASSERT(shape>0);
  BOOST_ASSERT(iscale>0);
	
  boost::shared_ptr<NormalConstType > Shape(new NormalConstType(shape));
  boost::shared_ptr<GammaConstType > IScale(new GammaConstType(iscale));
	
  boost::shared_ptr<GammaType > Gamma(new GammaType());
  boost::shared_ptr<GammaFactor > GammaF(new GammaFactor(Shape.get(), IScale.get(), Gamma.get()));
	
  m_Factors.push_back(GammaF);
  m_Nodes.push_back(Shape);
  m_Nodes.push_back(IScale);
  m_Nodes.push_back(Gamma);
	
  return Gamma.get();
}


template<class T>
typename ICR::ICA::Builder<T>::GaussianConstNode
ICR::ICA::Builder<T>::gaussian_const(const T value)
{
  boost::shared_ptr<GaussianConstType > Const(new GaussianConstType(value));
  m_Nodes.push_back(Const);
	
  return Const.get();
}
      
template<class T>
typename ICR::ICA::Builder<T>::GammaConstNode
ICR::ICA::Builder<T>::gamma_const(const T value)
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
typename ICR::ICA::Builder<T>::GaussianDataNode
ICR::ICA::Builder<T>::gaussian_data(const T data)
{
	
  boost::shared_ptr<GaussianDataType > Data(new GaussianDataType(data));
  m_Nodes.push_back(Data);
	
  return Data.get();
}
      

template<class T>
typename ICR::ICA::Builder<T>::GammaDataNode
ICR::ICA::Builder<T>::gamma_data(const T data)
{
	
  boost::shared_ptr<GammaDataType > Data(new GammaDataType(data));
  m_Nodes.push_back(Data);
	
  return Data.get();
}

template<class T>
typename ICR::ICA::Builder<T>::GaussianResultNode
ICR::ICA::Builder<T>::calc_gaussian(Expression<T>* Expr,  Context<T>& context)
{
  boost::shared_ptr<GaussianResultType > Child(new GaussianResultType());
  boost::shared_ptr<CalcGaussianFactor > ChildF
    (new CalcGaussianFactor(Expr, context,Child.get()));
	
  m_Nodes.push_back(Child);
  m_Factors.push_back(ChildF);
  return Child.get();
}


template<class T>
void 
ICR::ICA::Builder<T>::join(T& shape, GammaNode IScale ,  const T& data)
{

  boost::shared_ptr<GammaDataType > Data(new GammaDataType(data));
  boost::shared_ptr<NormalConstType > Shape(new NormalConstType(shape));
	
  boost::shared_ptr<GammaFactor > GammaF (new GammaFactor(Shape.get(), IScale, Data.get()));
	 
  m_Factors.push_back(GammaF);
  m_Nodes.push_back(Data);
  m_Nodes.push_back(Shape);
	
}
template<class T>
void
ICR::ICA::Builder<T>::join(T& shape, T& iscale, GammaDataNode Data  )
{

  boost::shared_ptr<NormalConstType > Shape(new NormalConstType(shape));
  boost::shared_ptr<GammaConstType > IScale(new GammaConstType(iscale));
	
  boost::shared_ptr<GammaFactor > GammaF (new GammaFactor(Shape.get(), IScale.get(), Data));
	 
  m_Factors.push_back(GammaF);
  m_Nodes.push_back(Shape);
  m_Nodes.push_back(IScale);
	
}

template<class T>
void 
ICR::ICA::Builder<T>::join(T& shape, GammaNode IScale, GammaDataNode Data  )
{
  boost::shared_ptr<NormalConstType > Shape (new NormalConstType(shape));
  boost::shared_ptr<GammaFactor > GammaF (new GammaFactor(Shape.get(), IScale, Data));
	 
  m_Factors.push_back(GammaF);
  m_Nodes.push_back(Shape);
}

template<class T>
void 
ICR::ICA::Builder<T>::join(T& mean, T& precision, GaussianDataNode Data  )
{

  boost::shared_ptr<GaussianConstType > Mean(new GaussianConstType(mean));
  boost::shared_ptr<GammaConstType > Precision(new GammaConstType(precision));
	
  boost::shared_ptr<GaussianFactor > GaussianF(new GaussianFactor (Mean.get(), Precision.get(), Data));
	
  m_Factors.push_back(GaussianF);
  m_Nodes.push_back(Mean);
  m_Nodes.push_back(Precision);
}
   
template<class T>   
void 
ICR::ICA::Builder<T>::join(Variable Mean, T& precision, GaussianDataNode Data  )
{

  boost::shared_ptr<GaussianConstType > Precision(new GaussianConstType(precision));
  boost::shared_ptr<GaussianFactor > GaussianF(new GaussianFactor (Mean, Precision.get(), Data));
	
  m_Factors.push_back(GaussianF);
  m_Nodes.push_back(Precision);
}
  
template<class T>    
void 
ICR::ICA::Builder<T>::join(T& mean, GammaNode Precision, GaussianDataNode Data  )
{

  boost::shared_ptr<GaussianConstType > Mean(new GaussianConstType(mean));
  boost::shared_ptr<GaussianFactor > GaussianF(new GaussianFactor(Mean.get(), Precision, Data));
	
  m_Factors.push_back(GaussianF);
  m_Nodes.push_back(Mean);
}
      
  
template<class T>    
void 
ICR::ICA::Builder<T>::join(Variable Mean, GammaNode& Precision, GaussianDataNode& Data  )
{
  boost::shared_ptr<GaussianFactor> GaussianF(new GaussianFactor(Mean,Precision,Data));
  m_Factors.push_back(GaussianF);
}

template<class T>
void 
ICR::ICA::Builder<T>::join(Variable Mean, GammaNode Precision, const T data )
{
	
  boost::shared_ptr<GaussianDataType > Data(new GaussianDataType(data));
  boost::shared_ptr<GaussianFactor> GaussianF(new GaussianFactor(Mean,Precision,Data.get()));
  m_Factors.push_back(GaussianF);
  m_Nodes.push_back(Data);
}


template<class T>
void 
ICR::ICA::Builder<T>::join(Variable Mean, GammaNode& Precision, GammaNode& Child  )
{
  boost::shared_ptr<GammaFactor> GammaF(new GammaFactor(Mean,Precision,Child));
  m_Factors.push_back(GammaF);
}

   
template<class T>   
void 
ICR::ICA::Builder<T>::join(Variable Mean, GammaNode& Precision, GaussianNode& Child  )
{
  boost::shared_ptr<GaussianFactor> GaussianF(new GaussianFactor(Mean,Precision,Child));
  m_Factors.push_back(GaussianF);
}

      
   
template<class T>   
void
ICR::ICA::Builder<T>::join( std::vector<Variable>& vMean, GammaNode& Precision, WeightsNode Weights,const T data )
{

  boost::shared_ptr<GaussianDataType > Data(new GaussianDataType(data));


  const size_t number = Weights->size();
  std::vector<Variable> vPrecision(number);
	
  boost::shared_ptr<CatagoryType>         Catagory(new CatagoryType(number));
  boost::shared_ptr<DiscreteFactor >      CatagoryF(new DiscreteFactor(Weights, Catagory.get()));
  m_Nodes.push_back(Catagory);
  m_Factors.push_back(CatagoryF);
  //make the vector of precision nodes and CatagoryNodes
  for(size_t i=0;i<number;++i){
    vPrecision[i] = Precision;
  }
	
  m_Nodes.push_back(Data);
	
  boost::shared_ptr<GaussianMixtureFactor> MixtureF(new GaussianMixtureFactor(vMean, vPrecision, Catagory.get() , Data.get()));
	
  m_Factors.push_back(MixtureF);
}

template<class T>	
void
ICR::ICA::Builder<T>::join( std::vector<Variable>& vMean, std::vector<Variable>& vPrecision, WeightsNode Weights,const T data )
{

  boost::shared_ptr<GaussianDataType > Data(new GaussianDataType(data));

  const size_t number = Weights->size();
	
  boost::shared_ptr<CatagoryType>         Catagory(new CatagoryType(number));
  boost::shared_ptr<DiscreteFactor >      CatagoryF(new DiscreteFactor(Weights, Catagory.get()));
  m_Nodes.push_back(Catagory);
  m_Factors.push_back(CatagoryF);
  //make the vector of precision nodes and CatagoryNodes
	
  m_Nodes.push_back(Data);
	
  boost::shared_ptr<GaussianMixtureFactor> MixtureF(new GaussianMixtureFactor(vMean, vPrecision, Catagory.get() , Data.get()));
	
  m_Factors.push_back(MixtureF);
}
	

template<class T>
size_t
ICR::ICA::Builder<T>::number_of_nodes() const
{
  return m_Nodes.size();
}

template<class T>
size_t
ICR::ICA::Builder<T>::number_of_factors() const
{
  return m_Factors.size();
}


template<class T>
double
ICR::ICA::Builder<T>::iterate()
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
ICR::ICA::Builder<T>::run(const double& epsilon, const size_t& max_iterations )
{
  Initialise();
	
  for(size_t i=0;i<max_iterations;++i){
    double Cost = iterate();
	  
    if (HasConverged(Cost, epsilon))
      return true;
  }

  return false;

}
      
      
template<class T>
bool
ICR::ICA::Builder<T>::HasConverged(const T Cost, const T epsilon)
{
#pragma omp critical
  {
    std::cout<<"COST = "<<Cost<<" PREV COST = "<<m_PrevCost<<" DIFF = "<<Cost-m_PrevCost<<std::endl;
    std::ofstream CostFile("Cost.txt",std::ios_base::app);
    // CostFile->open("Cost.txt",std::ios_base::app);
    CostFile<<Cost<<"\n";
  }
  // CostFile->close();
  // if (Cost<m_PrevCost) 
  //   std::cerr<<"WARNING:: BOUND HAS DECREASED - THIS IS A BAD SIGN - try making the priors less commital!\n";
  if (std::fabs(Cost - m_PrevCost)<epsilon) 
    return true;
	
  m_PrevCost = Cost;
  return false;
}
    
template<class T>  
void
ICR::ICA::Builder<T>::Initialise()
{
  if (!m_initialised) {
    PARALLEL_FOREACH(m_Nodes.begin(), m_Nodes.end(),
		     boost::bind(&VariableNode<T>::InitialiseMoments, _1)
		     );
    m_initialised=true;
  }

}


template class ICR::ICA::Builder<double>;
//template class ICR::ICA::Builder<float>;

