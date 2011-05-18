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

      typedef Factor<RectifiedGaussian<T>, T >      RectifiedGaussianFactor;
      typedef Factor<Gaussian<T>, T >      GaussianFactor;
      typedef Factor<Gamma<T>, T >         GammaFactor;
      typedef Factor<Dirichlet<T>, T >     DirichletFactor;
      typedef Factor<Discrete<T>, T >      DiscreteFactor;
      typedef Mixture<RectifiedGaussian<T>, T >     RectifiedGaussianMixtureFactor;
      typedef Mixture<Gaussian<T>, T >     GaussianMixtureFactor;

      typedef Details::CalcGaussianFactor<Gaussian, T >  CalcGaussianFactor;

      // typedef DeterministicFactor<Gaussian,Det::Add, T >     AddFactor;
      // typedef DeterministicFactor<Gaussian,Det::Multiply, T >     MultiplyFactor;

      typedef HiddenNode<Gaussian<T>, T >      GaussianType;
      typedef HiddenNode<RectifiedGaussian<T>, T >      RectifiedGaussianType;
      typedef HiddenNode<Gamma<T>, T >         GammaType;
      typedef HiddenNode<Dirichlet<T>, T >     DirichletType;
      typedef ObservedNode<Gaussian<T>, T >     GaussianDataType;
      typedef ObservedNode<Gamma<T>, T >        GammaDataType;
      typedef ObservedNode<Gaussian<T>, T > GaussianConstType;
      typedef ObservedNode<Gamma<T> , T>    GammaConstType;
      typedef ObservedNode<Gaussian<T>, T >   NormalConstType;
      typedef ObservedNode<Dirichlet<T>, T >  DirichletConstType;
      typedef DeterministicNode<Gaussian<T>, T>    GaussianResultType;

      
      typedef HiddenNode<Dirichlet<T>, T >     WeightsType;
      typedef HiddenNode<Discrete<T>, T >      CatagoryType;
      //      typedef ForwardingNode<Gaussian<T> >   GaussianMixtureType;
    public:
      typedef VariableNode<T>*                    Variable;
      typedef HiddenNode<RectifiedGaussian<T>, T >*      RectifiedGaussianNode;
      typedef HiddenNode<Gaussian<T>, T >*      GaussianNode;
      typedef HiddenNode<Gamma<T>, T >*         GammaNode;
      typedef ObservedNode<Gaussian<T>, T >*     GaussianDataNode;
      typedef ObservedNode<Gamma<T>, T >*        GammaDataNode;
      typedef ObservedNode<Gaussian<T>, T >* GaussianConstNode;
      typedef ObservedNode<Gamma<T>, T >*    GammaConstNode;
      typedef DeterministicNode<Gaussian<T>, T>*    GaussianResultNode;
      
      typedef HiddenNode<Dirichlet<T>, T >*      WeightsNode;
      typedef HiddenNode<Discrete<T>, T >*       CatagoryNode;
      // typedef ForwardingNode<Gaussian<T> >*  GaussianMixtureNode;
      
      
      Builder();
      ~Builder();

      WeightsNode
      weights(const size_t size);

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
      gaussian_mixture(std::vector<Variable>& vMean, std::vector<Variable>& vPrecision, WeightsNode Weights);
	
      RectifiedGaussianNode
      rectified_gaussian_mixture(std::vector<Variable>& vMean, std::vector<Variable>& vPrecision, WeightsNode Weights);
      
      GammaNode
      gamma(const T& shape, const T& iscale);

      GaussianConstNode
      gaussian_const(const T value);
      
      GammaConstNode
      gamma_const(const T value);
      
      GaussianDataNode
      gaussian_data(const T data);

      GammaDataNode
      gamma_data(const T data);

      GaussianResultNode
      calc_gaussian(Expression<T>* Expr,  Context<T>& context);

      void 
      join(T& shape, GammaNode IScale ,  const T& data );

      void 
      join(T& shape, T& iscale, GammaDataNode Data  );

      void 
      join(T& shape, GammaNode IScale, GammaDataNode Data  );

      void 
      join(T& mean, T& precision, GaussianDataNode Data  );
      
      void 
      join(Variable Mean, T& precision, GaussianDataNode Data  );
      
      void 
      join(T& mean, GammaNode Precision, GaussianDataNode Data  );
      
      
      void 
      join(Variable Mean, GammaNode& Precision, GaussianDataNode& Data  );

      void 
      join(Variable Mean, GammaNode Precision, const T data );


      void 
      join(Variable Mean, GammaNode& Precision, GammaNode& Child  );
      
      void 
      join(Variable Mean, GammaNode& Precision, GaussianNode& Child  );
      
      
      void
      join( std::vector<Variable>& vMean, GammaNode& Precision, WeightsNode Weights,const T data );

      void
      join( std::vector<Variable>& vMean, std::vector<Variable>& vPrecision, WeightsNode Weights,const T data );
	

      size_t
      number_of_nodes() const;

      size_t
      number_of_factors() const;


      double
      iterate();

      bool
      run(const double& epsilon = 1e-6, const size_t& max_iterations = 100);
      
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
