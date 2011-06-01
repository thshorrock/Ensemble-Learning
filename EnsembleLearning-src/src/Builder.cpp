#include "EnsembleLearning.hpp"

// #include "EnsembleLearning/node/variable/Constant.hpp"
// #include "EnsembleLearning/node/variable/Data.hpp"
#include "EnsembleLearning/node/variable/Hidden.hpp"
// #include "EnsembleLearning/node/variable/Dirichlet.hpp"
#include "EnsembleLearning/node/variable/Calculation.hpp"
// #include "EnsembleLearning/variable/Forwarding.hpp"

#include "EnsembleLearning/node/factor/Calculation.hpp"
#include "EnsembleLearning/node/factor/Factor.hpp"
#include "EnsembleLearning/node/factor/Mixture.hpp"
// #include "EnsembleLearning/factor/DeterministicFactor.hpp"

#include "EnsembleLearning/exponential_model/RectifiedGaussian.hpp"

#include "EnsembleLearning/detail/parallel_algorithms.hpp"

// typedef ICR::EnsembleLearning::Factor<ICR::EnsembleLearning::RectifiedGaussianModel<T> >      RectifiedGaussianFactor;
// typedef ICR::EnsembleLearning::Factor<ICR::EnsembleLearning::GaussianModel<T> >      GaussianFactor;
// typedef ICR::EnsembleLearning::Factor<ICR::EnsembleLearning::GammaModel<T> >         GammaFactor;
// typedef ICR::EnsembleLearning::Factor<ICR::EnsembleLearning::DirichletModel<T> >     DirichletFactor;
// typedef ICR::EnsembleLearning::Factor<ICR::EnsembleLearning::DiscreteModel<T> >      DiscreteFactor;
// typedef ICR::EnsembleLearning::Mixture<ICR::EnsembleLearning::RectifiedGaussianModel<T> >     RectifiedGaussianMixtureFactor;
// typedef ICR::EnsembleLearning::Mixture<ICR::EnsembleLearning::GaussianModel<T> >     GaussianMixtureFactor;

// typedef ICR::EnsembleLearning::Details::CalcGaussianFactor<ICR::EnsembleLearning::GaussianModel, T >  CalcGaussianFactor;

// // typedef DeterministicFactor<GaussianModel,Det::Add, T >     AddFactor;
// // typedef DeterministicFactor<GaussianModel,Det::Multiply, T >     MultiplyFactor;

// typedef ICR::EnsembleLearning::HiddenNode<ICR::EnsembleLearning::GaussianModel<T> >      GaussianType;
// typedef ICR::EnsembleLearning::HiddenNode<ICR::EnsembleLearning::RectifiedGaussianModel<T> >      RectifiedGaussianType;
// typedef ICR::EnsembleLearning::HiddenNode<ICR::EnsembleLearning::GammaModel<T> >         GammaType;
// typedef ICR::EnsembleLearning::HiddenNode<ICR::EnsembleLearning::DirichletModel<T> >     DirichletType;
// typedef ICR::EnsembleLearning::DataNode<ICR::EnsembleLearning::GaussianConstant<T> >     GaussianDataType;
// typedef ICR::EnsembleLearning::DataNode<ICR::EnsembleLearning::GammaConstant<T> >        GammaDataType;
// typedef ICR::EnsembleLearning::ConstantNode<ICR::EnsembleLearning::GaussianConstant<T> > GaussianConstType;
// typedef ICR::EnsembleLearning::ConstantNode<ICR::EnsembleLearning::GammaConstant<T> >    GammaConstType;
// typedef ICR::EnsembleLearning::ConstantNode<ICR::EnsembleLearning::NormalConstant<T> >   NormalConstType;
// typedef ICR::EnsembleLearning::ConstantNode<ICR::EnsembleLearning::DirichletConstant<T> >  DirichletConstType;
// typedef ICR::EnsembleLearning::DeterministicNode<ICR::EnsembleLearning::GaussianModel<T>, T>    GaussianResultType;

      
// typedef ICR::EnsembleLearning::HiddenNode<DirichletModel<T> >     WeightsType;
// typedef ICR::EnsembleLearning::HiddenNode<DiscreteModel<T> >      CatagoryType;
//      typedef ForwardingNode<GaussianModel<T> >   GaussianMixtureType;

      
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
ICR::EnsembleLearning::Builder<T>::weights(const size_t size)
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
typename ICR::EnsembleLearning::Builder<T>::GaussianNode
ICR::EnsembleLearning::Builder<T>::gaussian(Variable Mean, Variable Precision)
{
	
  boost::shared_ptr<GaussianType > Gaussian(new GaussianType());
  boost::shared_ptr<GaussianFactor > GaussianF(new GaussianFactor(Mean, Precision, Gaussian.get()));
	
  m_Factors.push_back(GaussianF);
  m_Nodes.push_back(Gaussian);
	
  return Gaussian.get();
}

template<class T>
typename ICR::EnsembleLearning::Builder<T>::GaussianNode
ICR::EnsembleLearning::Builder<T>::gaussian(GaussianNode Mean, const T& precision)
{
	
  boost::shared_ptr<GammaConstType > Precision(new GammaConstType(precision));
  m_Nodes.push_back(Precision);
  return gaussian(Mean,Precision.get());
}
    
template<class T>  
typename ICR::EnsembleLearning::Builder<T>::GaussianNode
ICR::EnsembleLearning::Builder<T>::gaussian(const T& mean, GammaNode Precision)
{
	
  boost::shared_ptr<GaussianConstType > Mean(new GaussianConstType(mean));
  m_Nodes.push_back(Mean);
  return gaussian(Mean.get(),Precision);
}
  
template<class T>
typename ICR::EnsembleLearning::Builder<T>::GaussianNode
ICR::EnsembleLearning::Builder<T>::gaussian(const T& mean, const T& precision)
{
	
  boost::shared_ptr<GaussianConstType > Mean(new GaussianConstType(mean));
  boost::shared_ptr<GammaConstType > Precision(new GammaConstType(precision));
	
  m_Nodes.push_back(Mean);
  m_Nodes.push_back(Precision);
	
  return gaussian(Mean.get(),Precision.get());
}
      

template<class T>
typename ICR::EnsembleLearning::Builder<T>::RectifiedGaussianNode
ICR::EnsembleLearning::Builder<T>::rectified_gaussian(Variable Mean,Variable  Precision)
{
	
  boost::shared_ptr<RectifiedGaussianType > RectifiedGaussian(new RectifiedGaussianType());
  boost::shared_ptr<RectifiedGaussianFactor > RectifiedGaussianF(new RectifiedGaussianFactor(Mean, Precision, RectifiedGaussian.get()));
	
  m_Factors.push_back(RectifiedGaussianF);
  m_Nodes.push_back(RectifiedGaussian);
	
  return RectifiedGaussian.get();
}

template<class T>
typename ICR::EnsembleLearning::Builder<T>::RectifiedGaussianNode
ICR::EnsembleLearning::Builder<T>::rectified_gaussian(GaussianNode Mean, const T& precision)
{
	
  boost::shared_ptr<GammaConstType > Precision(new GammaConstType(precision));
  m_Nodes.push_back(Precision);
  return rectified_gaussian(Mean,Precision.get());
}
    
template<class T>  
typename ICR::EnsembleLearning::Builder<T>::RectifiedGaussianNode
ICR::EnsembleLearning::Builder<T>::rectified_gaussian(const T& mean, GammaNode Precision)
{
	
  boost::shared_ptr<GaussianConstType > Mean(new GaussianConstType(mean));
  m_Nodes.push_back(Mean);
  return rectified_gaussian(Mean.get(),Precision);
}
  
template<class T>
typename ICR::EnsembleLearning::Builder<T>::RectifiedGaussianNode
ICR::EnsembleLearning::Builder<T>::rectified_gaussian(const T& mean, const T& precision)
{
	
  boost::shared_ptr<GaussianConstType > Mean(new GaussianConstType(mean));
  boost::shared_ptr<GammaConstType > Precision(new GammaConstType(precision));
	
  m_Nodes.push_back(Mean);
  m_Nodes.push_back(Precision);
	
  return rectified_gaussian(Mean.get(),Precision.get());
}
      
	

template<class T>
typename ICR::EnsembleLearning::Builder<T>::GaussianNode
ICR::EnsembleLearning::Builder<T>::gaussian_mixture(std::vector<Variable>& vMean, std::vector<Variable>& vPrecision, WeightsNode Weights)
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
typename ICR::EnsembleLearning::Builder<T>::RectifiedGaussianNode
ICR::EnsembleLearning::Builder<T>::rectified_gaussian_mixture(std::vector<Variable>& vMean, std::vector<Variable>& vPrecision, WeightsNode Weights)
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
typename ICR::EnsembleLearning::Builder<T>::GammaNode
ICR::EnsembleLearning::Builder<T>::gamma(const T& shape, const T& iscale)
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
  boost::shared_ptr<GaussianResultType > Child(new GaussianResultType());
  boost::shared_ptr<CalcGaussianFactor > ChildF
    (new CalcGaussianFactor(Expr, context,Child.get()));
	
  m_Nodes.push_back(Child);
  m_Factors.push_back(ChildF);
  return Child.get();
}


template<class T>
void 
ICR::EnsembleLearning::Builder<T>::join(T& shape, GammaNode IScale ,  const T& data)
{

  boost::shared_ptr<GammaDataType > Data(new GammaDataType(data));
  ++m_data_nodes;
  boost::shared_ptr<NormalConstType > Shape(new NormalConstType(shape));
	
  boost::shared_ptr<GammaFactor > GammaF (new GammaFactor(Shape.get(), IScale, Data.get()));
	 
  m_Factors.push_back(GammaF);
  m_Nodes.push_back(Data);
  m_Nodes.push_back(Shape);
  
	
}
template<class T>
void
ICR::EnsembleLearning::Builder<T>::join(T& shape, T& iscale, GammaDataNode Data  )
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
ICR::EnsembleLearning::Builder<T>::join(T& shape, GammaNode IScale, GammaDataNode Data  )
{
  boost::shared_ptr<NormalConstType > Shape (new NormalConstType(shape));
  boost::shared_ptr<GammaFactor > GammaF (new GammaFactor(Shape.get(), IScale, Data));
	 
  m_Factors.push_back(GammaF);
  m_Nodes.push_back(Shape);
}

template<class T>
void 
ICR::EnsembleLearning::Builder<T>::join(T& mean, T& precision, GaussianDataNode Data  )
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
ICR::EnsembleLearning::Builder<T>::join(Variable Mean, T& precision, GaussianDataNode Data  )
{

  boost::shared_ptr<GaussianConstType > Precision(new GaussianConstType(precision));
  boost::shared_ptr<GaussianFactor > GaussianF(new GaussianFactor (Mean, Precision.get(), Data));
	
  m_Factors.push_back(GaussianF);
  m_Nodes.push_back(Precision);
}
  
template<class T>    
void 
ICR::EnsembleLearning::Builder<T>::join(T& mean, GammaNode Precision, GaussianDataNode Data  )
{

  boost::shared_ptr<GaussianConstType > Mean(new GaussianConstType(mean));
  boost::shared_ptr<GaussianFactor > GaussianF(new GaussianFactor(Mean.get(), Precision, Data));
	
  m_Factors.push_back(GaussianF);
  m_Nodes.push_back(Mean);
}
      
  
template<class T>    
void 
ICR::EnsembleLearning::Builder<T>::join(Variable Mean, GammaNode& Precision, GaussianDataNode& Data  )
{
  boost::shared_ptr<GaussianFactor> GaussianF(new GaussianFactor(Mean,Precision,Data));
  m_Factors.push_back(GaussianF);
}

template<class T>
void 
ICR::EnsembleLearning::Builder<T>::join(Variable Mean, GammaNode Precision, const T data )
{
	
  boost::shared_ptr<GaussianDataType > Data(new GaussianDataType(data));
  ++m_data_nodes;
  boost::shared_ptr<GaussianFactor> GaussianF(new GaussianFactor(Mean,Precision,Data.get()));
  m_Factors.push_back(GaussianF);
  m_Nodes.push_back(Data);
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
void
ICR::EnsembleLearning::Builder<T>::join( std::vector<Variable>& vMean, GammaNode& Precision, WeightsNode Weights,const T data )
{

  boost::shared_ptr<GaussianDataType > Data(new GaussianDataType(data));
  ++m_data_nodes;


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
ICR::EnsembleLearning::Builder<T>::join( std::vector<Variable>& vMean, std::vector<Variable>& vPrecision, WeightsNode Weights,const T data )
{

  boost::shared_ptr<GaussianDataType > Data(new GaussianDataType(data));
  ++m_data_nodes;

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
ICR::EnsembleLearning::Builder<T>::run(const double& epsilon, const size_t& max_iterations )
{
  std::cout<<"initialising"<<std::endl;

  Initialise();
  
  std::cout<<" ... initialised"<<std::endl;
	
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
    std::cout<<"COST = "<<Cost<<" PREV COST = "<<m_PrevCost<<" pc DIFF = "<<100.0*(((Cost-m_PrevCost)/std::fabs(Cost)))<<std::endl;
  if (m_cost_file != "") { 
    std::ofstream CostFile(m_cost_file.c_str(),std::ios_base::app);
    // CostFile->open("Cost.txt",std::ios_base::app);
    CostFile<<Cost<<"\n";
  }
  }
  // CostFile->close();
  // if (Cost<m_PrevCost) 
  //   std::cerr<<"WARNING:: BOUND HAS DECREASED - THIS IS A BAD SIGN - try making the priors less commital!\n";
  if (100.0*std::fabs((Cost - m_PrevCost)/Cost)<epsilon) 
    return true;
	
  m_PrevCost = Cost;
  return false;
}
    
template<class T>  
void
ICR::EnsembleLearning::Builder<T>::Initialise()
{
  if (!m_initialised) {
    PARALLEL_FOREACH(m_Nodes.begin(), m_Nodes.end(),
		     boost::bind(&VariableNode<T>::InitialiseMoments, _1)
		     );
    m_initialised=true;
  }

}


template class ICR::EnsembleLearning::Builder<double>;
template class ICR::EnsembleLearning::Builder<float>;

