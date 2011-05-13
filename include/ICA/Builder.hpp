#ifndef BUILDER_HPP
#define BUILDER_HPP

#include "ICA/variable/ConstantNode.hpp"
#include "ICA/variable/DataNode.hpp"
#include "ICA/variable/HiddenNode.hpp"
#include "ICA/variable/DirichletNode.hpp"
#include "ICA/variable/DeterministicNode.hpp"
// #include "ICA/variable/ForwardingNode.hpp"

#include "ICA/factor/CalcGaussian.hpp"
#include "ICA/factor/Factor.hpp"
#include "ICA/factor/Mixture.hpp"
// #include "ICA/factor/DeterministicFactor.hpp"

#include "ICA/exponential_models/RectifiedGaussianModel.hpp"

#include <fstream>

#include "ICA/parallel_algorithms.hpp"

namespace ICR{
  namespace ICA{

    template<class T>
    class Builder
    {

      typedef Factor<RectifiedGaussianModel<T> >      RectifiedGaussianFactor;
      typedef Factor<GaussianModel<T> >      GaussianFactor;
      typedef Factor<GammaModel<T> >         GammaFactor;
      typedef Factor<DirichletModel<T> >     DirichletFactor;
      typedef Factor<DiscreteModel<T> >      DiscreteFactor;
      typedef Mixture<RectifiedGaussianModel<T> >     RectifiedGaussianMixtureFactor;
      typedef Mixture<GaussianModel<T> >     GaussianMixtureFactor;

      typedef Details::CalcGaussianFactor<GaussianModel, T >  CalcGaussianFactor;

      // typedef DeterministicFactor<GaussianModel,Det::Add, T >     AddFactor;
      // typedef DeterministicFactor<GaussianModel,Det::Multiply, T >     MultiplyFactor;

      typedef HiddenNode<GaussianModel<T> >      GaussianType;
      typedef HiddenNode<RectifiedGaussianModel<T> >      RectifiedGaussianType;
      typedef HiddenNode<GammaModel<T> >         GammaType;
      typedef HiddenNode<DirichletModel<T> >     DirichletType;
      typedef DataNode<GaussianConstant<T> >     GaussianDataType;
      typedef DataNode<GammaConstant<T> >        GammaDataType;
      typedef ConstantNode<GaussianConstant<T> > GaussianConstType;
      typedef ConstantNode<GammaConstant<T> >    GammaConstType;
      typedef ConstantNode<NormalConstant<T> >   NormalConstType;
      typedef ConstantNode<DirichletConstant<T> >  DirichletConstType;
      typedef DeterministicNode<GaussianModel<T>, T>    GaussianResultType;

      
      typedef HiddenNode<DirichletModel<T> >     WeightsType;
      typedef HiddenNode<DiscreteModel<T> >      CatagoryType;
      //      typedef ForwardingNode<GaussianModel<T> >   GaussianMixtureType;
    public:
      typedef VariableNode<T>*                    Variable;
      typedef HiddenNode<RectifiedGaussianModel<T> >*      RectifiedGaussianNode;
      typedef HiddenNode<GaussianModel<T> >*      GaussianNode;
      typedef HiddenNode<GammaModel<T> >*         GammaNode;
      typedef DataNode<GaussianConstant<T> >*     GaussianDataNode;
      typedef DataNode<GammaConstant<T> >*        GammaDataNode;
      typedef ConstantNode<GaussianConstant<T> >* GaussianConstNode;
      typedef ConstantNode<GammaConstant<T> >*    GammaConstNode;
      typedef DeterministicNode<GaussianModel<T>, T>*    GaussianResultNode;
      
      typedef HiddenNode<DirichletModel<T> >*      WeightsNode;
      typedef HiddenNode<DiscreteModel<T> >*       CatagoryNode;
      // typedef ForwardingNode<GaussianModel<T> >*  GaussianMixtureNode;
      
      
      Builder()
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
      ~Builder()
      {
	
	// delete CostFile;

      }

      WeightsNode
      Weights(const size_t size)
      {
	boost::shared_ptr<DirichletConstType > DirichletPrior(new DirichletConstType(size,1.0));
	boost::shared_ptr<DirichletType >      Dirichlet(new DirichletType(size));
	boost::shared_ptr<DirichletFactor >    DirichletF(new DirichletFactor(DirichletPrior.get(), Dirichlet.get()));
		m_Nodes.push_back(DirichletPrior);
	m_Nodes.push_back(Dirichlet);

	m_Factors.push_back(DirichletF);

	return Dirichlet.get();

	
      }

      GaussianNode
      Gaussian(Variable Mean, Variable Precision)
      {
	
	boost::shared_ptr<GaussianType > Gaussian(new GaussianType());
	boost::shared_ptr<GaussianFactor > GaussianF(new GaussianFactor(Mean, Precision, Gaussian.get()));
	
	m_Factors.push_back(GaussianF);
	m_Nodes.push_back(Gaussian);
	
	return Gaussian.get();
      }

      GaussianNode
      Gaussian(GaussianNode Mean, const T& precision)
      {
	
	boost::shared_ptr<GammaConstType > Precision(new GammaConstType(precision));
	m_Nodes.push_back(Precision);
	return Gaussian(Mean,Precision.get());
      }
      
      GaussianNode
      Gaussian(const T& mean, GammaNode Precision)
      {
	
	boost::shared_ptr<GaussianConstType > Mean(new GaussianConstType(mean));
	m_Nodes.push_back(Mean);
	return Gaussian(Mean.get(),Precision);
      }
  
      GaussianNode
      Gaussian(const T& mean, const T& precision)
      {
	
	boost::shared_ptr<GaussianConstType > Mean(new GaussianConstType(mean));
	boost::shared_ptr<GammaConstType > Precision(new GammaConstType(precision));
	
	m_Nodes.push_back(Mean);
	m_Nodes.push_back(Precision);
	
	return Gaussian(Mean.get(),Precision.get());
      }
      

      RectifiedGaussianNode
      RectifiedGaussian(Variable Mean,Variable  Precision)
      {
	
	boost::shared_ptr<RectifiedGaussianType > RectifiedGaussian(new RectifiedGaussianType());
	boost::shared_ptr<RectifiedGaussianFactor > RectifiedGaussianF(new RectifiedGaussianFactor(Mean, Precision, RectifiedGaussian.get()));
	
	m_Factors.push_back(RectifiedGaussianF);
	m_Nodes.push_back(RectifiedGaussian);
	
	return RectifiedGaussian.get();
      }

      RectifiedGaussianNode
      RectifiedGaussian(GaussianNode Mean, const T& precision)
      {
	
	boost::shared_ptr<GammaConstType > Precision(new GammaConstType(precision));
	m_Nodes.push_back(Precision);
	return RectifiedGaussian(Mean,Precision.get());
      }
      
      RectifiedGaussianNode
      RectifiedGaussian(const T& mean, GammaNode Precision)
      {
	
	boost::shared_ptr<GaussianConstType > Mean(new GaussianConstType(mean));
	m_Nodes.push_back(Mean);
	return RectifiedGaussian(Mean.get(),Precision);
      }
  
      RectifiedGaussianNode
      RectifiedGaussian(const T& mean, const T& precision)
      {
	
	boost::shared_ptr<GaussianConstType > Mean(new GaussianConstType(mean));
	boost::shared_ptr<GammaConstType > Precision(new GammaConstType(precision));
	
	m_Nodes.push_back(Mean);
	m_Nodes.push_back(Precision);
	
	 return RectifiedGaussian(Mean.get(),Precision.get());
      }
      
      
      // GaussianNode
      // GaussianMixture( std::vector<Variable>& vMean, const T& precision, WeightsNode Weights)
      // {
      // 	const size_t number = Weights->size();
      // 	std::vector<Variable> vPrecision(number);
	
      // 	make the precision nodes
      // 	for(size_t i=0;i<number;++i){
      // 	  boost::shared_ptr<GammaConstType > Precision(new GammaConstType(precision));
      // 	  vPrecision[i] = Precision.get();
      // 	  m_Nodes.push_back(Precision);
      // 	}
	
      // 	boost::shared_ptr<GaussianType > Gaussian(new GaussianType());
      // 	boost::shared_ptr<GaussianMixtureFactor> MixtureF(new GaussianMixtureFactor(vMean, vPrecision,  Weights, Gaussian.get()));
	
      // 	m_Factors.push_back(MixtureF);
      // 	m_Nodes.push_back(Gaussian);
      // 	return Gaussian.get();
      // }
	

      GaussianNode
      GaussianMixture(std::vector<Variable>& vMean, std::vector<Variable>& vPrecision, WeightsNode Weights)
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
	
      RectifiedGaussianNode
      RectifiedGaussianMixture(std::vector<Variable>& vMean, std::vector<Variable>& vPrecision, WeightsNode Weights)
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
      
      // GaussianNode
      // GaussianMixture( Variable Mean,  Variable Precision, WeightsNode Weights)
      // {
      // 	const size_t number = Weights->size();
      // 	std::vector<Variable> vMean(number);
      // 	std::vector<Variable> vPrecision(number);
      // 	for(size_t i=0;i<number;++i){
      // 	  vMean[i] = Mean;
      // 	  vPrecision[i] = Precision;
      // 	}
      // 	return GaussianMixture(vMean,vPrecision,Weights);
      // }


      GammaNode
      Gamma(const T& shape, const T& iscale)
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

      GaussianConstNode
      GaussianConst(const T value)
      {
	
	boost::shared_ptr<GaussianConstType > Const(new GaussianConstType(value));
	m_Nodes.push_back(Const);
	
	return Const.get();
      }
      
      GammaConstNode
      GammaConst(const T value)
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
      
      GaussianDataNode
      GaussianData(const T data, const size_t skip = 0)
      {
	
	boost::shared_ptr<GaussianDataType > Data(new GaussianDataType(data));
	m_Nodes.push_back(Data);
	
	return Data.get();
      }
      

      GammaDataNode
      GammaData(const T data, const size_t skip = 0)
      {
	
	boost::shared_ptr<GammaDataType > Data(new GammaDataType(data));
	m_Nodes.push_back(Data);
	
	return Data.get();
      }

      GaussianResultNode
      CalcGaussian(Expression<T>* Expr, const Context<T>& context)
      {
	boost::shared_ptr<GaussianResultType > Child(new GaussianResultType());
	boost::shared_ptr<CalcGaussianFactor > ChildF
	  (new CalcGaussianFactor(Expr, context,Child.get()));
	
	m_Nodes.push_back(Child);
	m_Factors.push_back(ChildF);
	return Child.get();
      }

      // GaussianResultNode
      // Add(GaussianNode ParentA, GaussianNode ParentB)
      // {
      // 	std::cout<<"Builder A = "<<ParentA<<std::endl;

      // 	boost::shared_ptr<GaussianResultType > Child(new GaussianResultType());
      // 	boost::shared_ptr<AddFactor > ChildF(new AddFactor(ParentA, ParentB,Child.get()));
	
      // 	std::cout<<"Builder A2= "<<ParentA<<std::endl;
      // 	m_Nodes.push_back(Child);
      // 	m_Factors.push_back(ChildF);
      // 	return Child.get();
      // }

      // GaussianResultNode
      // Multiply(GaussianNode ParentA, GaussianNode ParentB)
      // {
      // 	boost::shared_ptr<GaussianResultType > Child(new GaussianResultType());
      // 	boost::shared_ptr<MultiplyFactor > ChildF(new MultiplyFactor(ParentA, ParentB,Child.get()));
      // 	m_Nodes.push_back(Child);
      // 	m_Factors.push_back(ChildF);
      // 	return Child.get();
      // }
      

      void 
      Join(T& shape, GammaNode IScale ,  const T& data, const size_t skip = 0  )
       {

	boost::shared_ptr<GammaDataType > Data(new GammaDataType(data));
	boost::shared_ptr<NormalConstType > Shape(new NormalConstType(shape));
	
	boost::shared_ptr<GammaFactor > GammaF (new GammaFactor(Shape.get(), IScale, Data.get()));
	 
	m_Factors.push_back(GammaF);
	m_Nodes.push_back(Data);
	m_Nodes.push_back(Shape);
	
       }
      void 
      Join(T& shape, T& iscale, GammaDataNode Data  )
      {

	boost::shared_ptr<NormalConstType > Shape(new NormalConstType(shape));
	boost::shared_ptr<GammaConstType > IScale(new GammaConstType(iscale));
	
	boost::shared_ptr<GammaFactor > GammaF (new GammaFactor(Shape.get(), IScale.get(), Data));
	 
	m_Factors.push_back(GammaF);
	m_Nodes.push_back(Shape);
	m_Nodes.push_back(IScale);
	
       }

      void 
      Join(T& shape, GammaNode IScale, GammaDataNode Data  )
       {
	boost::shared_ptr<NormalConstType > Shape (new NormalConstType(shape));
	boost::shared_ptr<GammaFactor > GammaF (new GammaFactor(Shape.get(), IScale, Data));
	 
	m_Factors.push_back(GammaF);
	m_Nodes.push_back(Shape);
       }

      void 
      Join(T& mean, T& precision, GaussianDataNode Data  )
       {

	boost::shared_ptr<GaussianConstType > Mean(new GaussianConstType(mean));
	boost::shared_ptr<GammaConstType > Precision(new GammaConstType(precision));
	
	boost::shared_ptr<GaussianFactor > GaussianF(new GaussianFactor (Mean.get(), Precision.get(), Data));
	
	m_Factors.push_back(GaussianF);
	m_Nodes.push_back(Mean);
	m_Nodes.push_back(Precision);
       }
      
      void 
      Join(Variable Mean, T& precision, GaussianDataNode Data  )
       {

	boost::shared_ptr<GaussianConstType > Precision(new GaussianConstType(precision));
	boost::shared_ptr<GaussianFactor > GaussianF(new GaussianFactor (Mean, Precision.get(), Data));
	
	m_Factors.push_back(GaussianF);
	m_Nodes.push_back(Precision);
       }
      
      void 
      Join(T& mean, GammaNode Precision, GaussianDataNode Data  )
       {

	boost::shared_ptr<GaussianConstType > Mean(new GaussianConstType(mean));
	boost::shared_ptr<GaussianFactor > GaussianF(new GaussianFactor(Mean.get(), Precision, Data));
	
	m_Factors.push_back(GaussianF);
	m_Nodes.push_back(Mean);
       }
      
      
      void 
      Join(Variable Mean, GammaNode& Precision, GaussianDataNode& Data  )
      {
	boost::shared_ptr<GaussianFactor> GaussianF(new GaussianFactor(Mean,Precision,Data));
	m_Factors.push_back(GaussianF);
      }

      void 
      Join(Variable Mean, GammaNode Precision, const T data )
      {
	
	boost::shared_ptr<GaussianDataType > Data(new GaussianDataType(data));
	boost::shared_ptr<GaussianFactor> GaussianF(new GaussianFactor(Mean,Precision,Data.get()));
	m_Factors.push_back(GaussianF);
	m_Nodes.push_back(Data);
      }


      void 
      Join(Variable Mean, GammaNode& Precision, GammaNode& Child  )
      {
	boost::shared_ptr<GammaFactor> GammaF(new GammaFactor(Mean,Precision,Child));
	m_Factors.push_back(GammaF);
      }

      
      void 
      Join(Variable Mean, GammaNode& Precision, GaussianNode& Child  )
      {
	boost::shared_ptr<GaussianFactor> GaussianF(new GaussianFactor(Mean,Precision,Child));
	m_Factors.push_back(GaussianF);
      }

      
      
      void
      Join( std::vector<Variable>& vMean, GammaNode& Precision, WeightsNode Weights,const T data )
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
	
      void
      Join( std::vector<Variable>& vMean, std::vector<Variable>& vPrecision, WeightsNode Weights,const T data )
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
	


	
      // GaussianMixtureNode
      // Join(const T& mean, const T& precision, WeightsNode Weights,  GaussianNode& Child)
      // {
      // 	const number = Weights->size();

      // 	boost::shared_ptr<GaussianConstType > Mean(new GaussianConstType(mean));
      // 	boost::shared_ptr<GammaConstType > Precision(new GammaConstType(precision));
	
      // 	std::vector<GaussianType> Mixture(number);
      // 	std::vector<GaussianFactor> MixtureF(number);
      // 	for(size_t i=0;i<number;++i){
      // 	  Mixture[i]  = Gaussian(new GaussianType());
      // 	  MixtureF[i] = GaussianF(new GaussianFactor(Mean.get(), Precision.get(), Mixture[i].get()));
      // 	  m_Factors.push_back(GaussianF);
      // 	  m_Nodes.push_back(Data);
      // 	}
	

	
      // 	boost::shared_ptr<GaussianFactor> GaussianMF(Mixture.get(), Weights.get, Child );
	
	
      // 	m_Factors.push_back(GaussianMF);
	
      // }

      size_t
      NumberOfNodes() const
      {
	return m_Nodes.size();
      }

      size_t
      NumberOfFactors() const
      {
	return m_Factors.size();
      }


      double
      Iterate()
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

      bool
      Run(const double& epsilon = 1e-6, const size_t& max_iterations = 100)
      {
	// //Initialise;
	 // for(size_t i=0;i<10;++i){
	 //    Iterate();
	 //  }

	// Coster Cost;
	//     PARALLEL_FOREACH(m_Factors.begin(), m_Factors.end(),
	// 		     boost::bind(&FactorNode<T>::Iterate, _1 )
	// 		     );
	    
	//     // std::cout<<"ITERATE Nodes"<<std::endl;
	//     PARALLEL_FOREACH(m_Nodes.begin(), m_Nodes.end(),
	//     		     boost::bind(&VariableNode<T>::Iterate, _1, boost::ref(Cost))
	//     		     );

	// Initialise();
	
	for(size_t i=0;i<max_iterations;++i){
	  double Cost = Iterate();
	  
	  if (HasConverged(Cost, epsilon))
	    return true;
	}

	return false;


      }
      
      
    private:
      bool
      HasConverged(const T Cost, const T epsilon)
      {
	std::cout<<"COST = "<<Cost<<" PREV COST = "<<m_PrevCost<<" DIFF = "<<Cost-m_PrevCost<<std::endl;
	std::ofstream CostFile("Cost.txt",std::ios_base::app);
	// CostFile->open("Cost.txt",std::ios_base::app);
	CostFile<<Cost<<"\n";
	// CostFile->close();
	// if (Cost<m_PrevCost) 
	//   std::cerr<<"WARNING:: BOUND HAS DECREASED - THIS IS A BAD SIGN - try making the priors less commital!\n";
	if (std::fabs(Cost - m_PrevCost)<epsilon) 
	  return true;
	
	m_PrevCost = Cost;
	return false;
      }
      
      void
      Initialise()
      {
	if (!m_initialised) {
	  PARALLEL_FOREACH(m_Nodes.begin(), m_Nodes.end(),
			   boost::bind(&VariableNode<T>::InitialiseMoments, _1)
			   );
	  m_initialised=true;
	}

      }

      // std::ofstream* CostFile;
      T m_PrevCost;
      std::vector<boost::shared_ptr<FactorNode<T> > > m_Factors;
      std::vector<boost::shared_ptr<VariableNode<T> > > m_Nodes;
      bool m_initialised;
    };

  }
}

#endif //BUILDER_HPP guard
