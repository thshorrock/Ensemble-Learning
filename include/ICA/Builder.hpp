#ifndef BUILDER_HPP
#define BUILDER_HPP

// #include "ICA/node/variable/Constant.hpp"
// #include "ICA/node/variable/Data.hpp"
#include "ICA/node/variable/Hidden.hpp"
#include "ICA/node/variable/Observed.hpp"
// #include "ICA/node/variable/Dirichlet.hpp"
#include "ICA/node/variable/Calculation.hpp"
// #include "ICA/variable/Forwarding.hpp"

#include "ICA/node/factor/Calculation.hpp"
#include "ICA/node/factor/Factor.hpp"
#include "ICA/node/factor/Mixture.hpp"
// #include "ICA/factor/DeterministicFactor.hpp"

#include "ICA/exponential_model/RectifiedGaussian.hpp"

#include <fstream>

#include "ICA/detail/parallel_algorithms.hpp"

namespace ICR{
  namespace ICA{

    template<class T>
    class Builder
    {

      typedef Factor<RectifiedGaussian<T> >      RectifiedGaussianFactor;
      typedef Factor<Gaussian<T> >      GaussianFactor;
      typedef Factor<Gamma<T> >         GammaFactor;
      typedef Factor<Dirichlet<T> >     DirichletFactor;
      typedef Factor<Discrete<T> >      DiscreteFactor;
      typedef Mixture<RectifiedGaussian<T> >     RectifiedGaussianMixtureFactor;
      typedef Mixture<Gaussian<T> >     GaussianMixtureFactor;

      typedef Details::CalcGaussianFactor<Gaussian, T >  CalcGaussianFactor;

      // typedef DeterministicFactor<Gaussian,Det::Add, T >     AddFactor;
      // typedef DeterministicFactor<Gaussian,Det::Multiply, T >     MultiplyFactor;

      typedef HiddenNode<Gaussian<T> >      GaussianType;
      typedef HiddenNode<RectifiedGaussian<T> >      RectifiedGaussianType;
      typedef HiddenNode<Gamma<T> >         GammaType;
      typedef HiddenNode<Dirichlet<T> >     DirichletType;
      typedef ObservedNode<Gaussian<T> >     GaussianDataType;
      typedef ObservedNode<Gamma<T> >        GammaDataType;
      typedef ObservedNode<Gaussian<T> > GaussianConstType;
      typedef ObservedNode<Gamma<T> >    GammaConstType;
      typedef ObservedNode<Gaussian<T> >   NormalConstType;
      typedef ObservedNode<Dirichlet<T> >  DirichletConstType;
      typedef DeterministicNode<Gaussian<T>, T>    GaussianResultType;

      
      typedef HiddenNode<Dirichlet<T> >     WeightsType;
      typedef HiddenNode<Discrete<T> >      CatagoryType;
      //      typedef ForwardingNode<Gaussian<T> >   GaussianMixtureType;
    public:
      typedef VariableNode<T>*                    Variable;
      typedef HiddenNode<RectifiedGaussian<T> >*      RectifiedGaussianNode;
      typedef HiddenNode<Gaussian<T> >*      GaussianNode;
      typedef HiddenNode<Gamma<T> >*         GammaNode;
      typedef ObservedNode<Gaussian<T> >*     GaussianDataNode;
      typedef ObservedNode<Gamma<T> >*        GammaDataNode;
      typedef ObservedNode<Gaussian<T> >* GaussianConstNode;
      typedef ObservedNode<Gamma<T> >*    GammaConstNode;
      typedef DeterministicNode<Gaussian<T>, T>*    GaussianResultNode;
      
      typedef HiddenNode<Dirichlet<T> >*      WeightsNode;
      typedef HiddenNode<Discrete<T> >*       CatagoryNode;
      // typedef ForwardingNode<Gaussian<T> >*  GaussianMixtureNode;
      
      
      Builder();
      ~Builder();

      WeightsNode
      Weights(const size_t size);

      GaussianNode
      gaussian(Variable Mean, Variable Precision);

      GaussianNode
      gaussian(GaussianNode Mean, const T& precision);
      
      GaussianNode
      gaussian(const T& mean, GammaNode Precision);
  
      GaussianNode
      gaussian(const T& mean, const T& precision);
      

      RectifiedGaussianNode
      rectified_gaussian(Variable Mean,Variable  Precision);

      RectifiedGaussianNode
      rectified_gaussian(GaussianNode Mean, const T& precision);
      
      RectifiedGaussianNode
      rectified_gaussian(const T& mean, GammaNode Precision);
  
      RectifiedGaussianNode
      rectified_gaussian(const T& mean, const T& precision);


      GaussianNode
      GaussianMixture(std::vector<Variable>& vMean, std::vector<Variable>& vPrecision, WeightsNode Weights);
	
      RectifiedGaussianNode
      RectifiedGaussianMixture(std::vector<Variable>& vMean, std::vector<Variable>& vPrecision, WeightsNode Weights);
      
      GammaNode
      gamma(const T& shape, const T& iscale);

      GaussianConstNode
      GaussianConst(const T value);
      
      GammaConstNode
      GammaConst(const T value);
      
      GaussianDataNode
      GaussianData(const T data);

      GammaDataNode
      GammaData(const T data);

      GaussianResultNode
      CalcGaussian(Expression<T>* Expr,  Context<T>& context);

      void 
      Join(T& shape, GammaNode IScale ,  const T& data );

      void 
      Join(T& shape, T& iscale, GammaDataNode Data  );

      void 
      Join(T& shape, GammaNode IScale, GammaDataNode Data  );

      void 
      Join(T& mean, T& precision, GaussianDataNode Data  );
      
      void 
      Join(Variable Mean, T& precision, GaussianDataNode Data  );
      
      void 
      Join(T& mean, GammaNode Precision, GaussianDataNode Data  );
      
      
      void 
      Join(Variable Mean, GammaNode& Precision, GaussianDataNode& Data  );

      void 
      Join(Variable Mean, GammaNode Precision, const T data );


      void 
      Join(Variable Mean, GammaNode& Precision, GammaNode& Child  );
      
      void 
      Join(Variable Mean, GammaNode& Precision, GaussianNode& Child  );
      
      
      void
      Join( std::vector<Variable>& vMean, GammaNode& Precision, WeightsNode Weights,const T data );

      void
      Join( std::vector<Variable>& vMean, std::vector<Variable>& vPrecision, WeightsNode Weights,const T data );
	

      size_t
      NumberOfNodes() const;

      size_t
      NumberOfFactors() const;


      double
      Iterate();

      bool
      Run(const double& epsilon = 1e-6, const size_t& max_iterations = 100);
      
    private:
      bool
      HasConverged(const T Cost, const T epsilon);
      
      void
      Initialise();

      // std::ofstream* CostFile;
      T m_PrevCost;
      std::vector<boost::shared_ptr<FactorNode<T> > > m_Factors;
      std::vector<boost::shared_ptr<VariableNode<T> > > m_Nodes;
      bool m_initialised;
    };

  }
}

#endif //BUILDER_HPP guard
