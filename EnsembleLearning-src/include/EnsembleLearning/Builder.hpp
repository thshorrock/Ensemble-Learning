
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



#ifndef BUILDER_HPP
#define BUILDER_HPP

//node definitions required so inheritence relationship available to user.
#include "EnsembleLearning/node/variable/Hidden.hpp"
#include "EnsembleLearning/node/variable/Observed.hpp"
#include "EnsembleLearning/node/variable/Calculation.hpp"


#include <boost/shared_ptr.hpp>
#include <vector>
#include <string>


namespace ICR{
  namespace EnsembleLearning{

    //Forward declarations of models
    template<class T> class Discrete;
    template<class T> class Dirichlet;
    template<class T> class Gamma;
    template<class T> class Gaussian;
    template<class T> class RectifiedGaussian;
    
    //forward declaration of Interfaces
    template<class T> class VariableNode;
    template<class T> class FactorNode;
    

    //Forward declarations of Factors
    namespace detail{
      template<class Model, class T> class Factor;
      template<class Model, class T> class Mixture;
      template<template<class> class Model, class T> class Deterministic;
    }


    //forward declaration of Expressions
    template<class T> class Expression;
    template<class T> class Context;
    

    /** Build the components for ensemble learning.
     *  @tparam T The datatype used by the model.
     *   This will either be float or double.
     *  
     *  The member functions to the builder class returns pointers
     *  to VariableNodes, from which the inferred moments can be found.
     *  The Factor nodes that link the variable nodes are created and stored internally.
     *  The memory management of the VariableNode and the Factor nodes are entirely
     *  handled by the Builder class.
     *
     *  Example of use:
     *  The following infers the mean and variance from a set of data.
     *  @code
     *  #include "EnsembleLearning.hpp"
     *  #include <vector>
     *  #include <iostream>
     *  using namespace ICR::EnsembleLearning;
     *  
     *  int
     *  main  (int ac, char **av)
     *  {
     *    //Create the data
     *    rng* random = Random::Instance(); //a random number generator.
     *    std::vector<double> data(10);
     *    for(size_t i = 0; i<10; ++i)
     *    {
     *      data[i] = random->Gaussian(2, 3); //mean = 3, standard deviation = 2
     *    }
     *
     *    //Build the model
     *    Builder<double> build;  
     *    
     *    typedef Builder<double>::Variable Variable;  //A convenient alias
     *
     *    //Nodes to hold the learnt mean and precision
     *    Variable mean = build.gaussian(0.01,0.01);
     *    Variable prec = build.gamma(0.01,0.1);
     *
     *    //Model The data as Gaussian distributed, 
     *    //  each data point is modelled indepenantly
     *    for(size_t i; i<data.size(); ++i)
     *    {
     *      build.join(mean,prec,data[i];
     *    }
     *  
     *    //The model is complete, now need to do the inference.
     *    //  Iterate until convergance to within 0.1% or 100 iterations
     *    build.run(0.1,100);
     *
     *    //output the inferred mean and precision
     *    std::cout<<"mean      = "<<Mean(mean)<<" +- "<<StandardDeviation(mean)<<std::endl;
     *    std::cout<<"precision = "<<Mean(prec)<<" +- "<<StandardDeviation(prec)<<std::endl;
     *
     *  }
     *  @endcode
     *
     *  @ingroup UserInterface
     *
     */


    template<class T>
    class Builder
    {

      typedef detail::Factor<RectifiedGaussian<T>, T >      RectifiedGaussianFactor;
      typedef detail::Factor<Gaussian<T>, T >      GaussianFactor;
      typedef detail::Factor<Gamma<T>, T >         GammaFactor;
      typedef detail::Factor<Dirichlet<T>, T >     DirichletFactor;
      typedef detail::Factor<Discrete<T>, T >      DiscreteFactor;
      typedef detail::Mixture<RectifiedGaussian<T>, T >     RectifiedGaussianMixtureFactor;
      typedef detail::Mixture<Gaussian<T>, T >     GaussianMixtureFactor;

      typedef detail::Deterministic<Gaussian, T >  DeterministicFactor;

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
    public:
      
      /** @name Useful typdefs for types that are exposed to the user.
       */
      ///@{
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
      ///@}

      /** A constructor.
       *  @param cost_file The output filename to plot the evidence.
       *   If the filename is blank (the default) 
       *   then the evidence is not plotted.
       */
      Builder(const std::string& cost_file = "");
      

      /**
       *  The destructor takes care of all the memory management of the created 
       *  Nodes.
       *  There is no need to use the delete keyword.
       */
      ~Builder();


      /** @name Create a Gaussian VariableNode.
       */
      ///@{
      
      /** Create a Gaussian VariableNode.
       * @param Mean The node that stores the Momets for the Mean.
       *  The model for the Mean node must be a Gaussian or Rectified Gaussian.
       * @param Precision The node that stores that Moments for the Precicision.
       *  The model for the Precision node must be a Gamma.
       * @return The GaussianNode that holds the inferred Gaussian Moments.
       */
      GaussianNode
      gaussian(Variable Mean, Variable Precision);

      /** Create a Gaussian VariableNode.
       * @param Mean The node that stores the Momets for the Mean.
       *  The model for the Mean node must be a Gaussian or Rectified Gaussian.
       * @param precision The value of the precision.
       * @return The GaussianNode that holds the inferred Gaussian Moments.
       */
      GaussianNode
      gaussian(GaussianNode Mean, const T& precision);

      /** Create a Gaussian VariableNode.
       * @param mean The value of the mean.
       *  The model for the Mean node must be a Gaussian or Rectified Gaussian.
       * @param Precision The node that stores that Moments for the Precicision.
       *  The model for the Precision node must be a Gamma.
       * @return The GaussianNode that holds the inferred Gaussian Moments.
       */
      GaussianNode
      gaussian(const T& mean, GammaNode Precision);

      /** Create a Gaussian VariableNode.
       * @param mean The value of the mean.
       * @param precision The value of the precision.
       * @return The GaussianNode that holds the inferred Gaussian Moments.
       */
      GaussianNode
      gaussian(const T& mean, const T& precision);
      
      ///@}
      
      /** @name Create a RectifiedGaussian VariableNode.
       */
      ///@{ 
      /** Create a RectifiedGaussian VariableNode.
       * @param Mean The node that stores the Momets for the Mean.
       *  The model for the Mean node must be a Gaussian or Rectified Gaussian.
       * @param Precision The node that stores that Moments for the Precicision.
       *  The model for the Precision node must be a Gamma.
       * @return The GaussianNode that holds the inferred Rectified Gaussian Moments.
       */
      RectifiedGaussianNode
      rectified_gaussian(Variable Mean,Variable  Precision);

      /** Create a RectifiedGaussian VariableNode.
       * @param Mean The node that stores the Momets for the Mean.
       *  The model for the Mean node must be a Gaussian or Rectified Gaussian.
       * @param precision The value of the precision.
       * @return The GaussianNode that holds the inferred Rectified Gaussian Moments.
       */
      RectifiedGaussianNode
      rectified_gaussian(GaussianNode Mean, const T& precision);
       /** Create a RectifiedGaussian VariableNode.
       * @param mean The value of the mean.
       *  The model for the Mean node must be a Gaussian or Rectified Gaussian.
       * @param Precision The node that stores that Moments for the Precicision.
       *  The model for the Precision node must be a Gamma.
       * @return The GaussianNode that holds the inferred Rectified Gaussian Moments.
       */
      RectifiedGaussianNode
      rectified_gaussian(const T& mean, GammaNode Precision);
  
      /** Create a RectifiedGaussian VariableNode.
       * @param mean The value of the mean.
       * @param precision The value of the precision.
       * @return The GaussianNode that holds the inferred Rectified Gaussian Moments.
       */
      RectifiedGaussianNode
      rectified_gaussian(const T& mean, const T& precision);

      ///@}
      /** @name Create a Gamma VariableNode.
       */
      ///@{ 
      /** Create a Gamma VariableNode
       *  @param shape The value of the shape of the GammaNode.
       *  @param iscale The value of the inverse scale of the GammaNode.
       *  @return The GammaNode that stores the inferred moments to the Gamma Distribution.
       */
      GammaNode
      gamma(const T& shape, const T& iscale);
      ///@}


      /** @name Create a Mixture Model.
       */
      ///@{
      
      /** Create a weights node.
       *  The weights node is used to weight Mixture models.
       *
       *  The moments of the Weights node are the natural logarithm of the
       *  probabilities of the weights.
       *  @param size The number of components in the mixture.
       * @return The WeightsNode that holds the inferred logarithm of the weights.
       */
      WeightsNode
      weights(const size_t size);

      /** Gaussian Mixture Model.
       *  @param vMean The vector container containing all the  variables representing the means.
       *  @param vPrecision The vector container containing all the  variables representing the precisions.
       *  @param Weights A Weights variable that stores the weights to the means and precisions in vMean and vPrecision.
       *  @attention The size of vMean and vPrecision must be identical, 
       *   and must be the same as the size of the Weights Node.
       * @return The GaussianNode that holds the inferred  Gaussian Moments.
       */
      GaussianNode
      gaussian_mixture( std::vector<Variable>& vMean, 
		        std::vector<Variable>& vPrecision, 
		       WeightsNode Weights);
	

  
      /** Rectified Gaussian Mixture Model.
       *  @param vMean The vector containing all the  variables representing the means.
       *  @param vPrecision The vector containing all the  variables representing the precisions.
       *  @param Weights A Weights variable that stores the weights to the means and precisions in vMean and vPrecision.
       *  @attention The size of vMean and vPrecision must be identical, 
       *   and must be the same as the size of the Weights Node.
       * @return The GaussianNode that holds the inferred Rectified Gaussian Moments.
       */
      RectifiedGaussianNode
      rectified_gaussian_mixture( std::vector<Variable>& vMean, 
				  std::vector<Variable>& vPrecision,
				 WeightsNode Weights);
      
      /** Gaussian Mixture Model.
       *  @param MeanBegin The iterator at the beginning of the the container containing all the  variables representing the means.
       *  @param PrecisionBegin The iterator at the beginning of the the container containing all the  variables representing the precisions.
       *  @param Weights A Weights variable that stores the weights to the means and precisions in vMean and vPrecision.
       *  @attention The size of vMean and vPrecision must be identical, 
       *   and must be the same as the size of the Weights Node.
       * @return The GaussianNode that holds the inferred  Gaussian Moments.
       */
      template<class MeanIterator, class PrecIterator>
      GaussianNode
      gaussian_mixture(const MeanIterator& MeanBegin, 
		       const PrecIterator& PrecisionBegin, 
		       WeightsNode Weights);
	

  
      /** Rectified Gaussian Mixture Model.
       *  @param MeanBegin The iterator at the beginning of the the container containing all the  variables representing the means.
       *  @param PrecisionBegin The iterator at the beginning of the the container containing all the  variables representing the precisions.
       *  @param Weights A Weights variable that stores the weights to the means and precisions in vMean and vPrecision.
       *  @attention The size of vMean and vPrecision must be identical, 
       *   and must be the same as the size of the Weights Node.
       * @return The GaussianNode that holds the inferred Rectified Gaussian Moments.
       */
      template<class MeanIterator, class PrecIterator>
      RectifiedGaussianNode
      rectified_gaussian_mixture(const MeanIterator& MeanBegin, 
				 const PrecIterator& PrecisionBegin,
				 WeightsNode Weights);
      
      ///@}
      

      /** @name Create a constant node.
       */
      ///@{
      /** Create a Gaussian Constant.
       *  This is appropriate input for VariableNode's that expect input in Gaussian form.
       *  For example, the mean of a Gaussian or Rectified Gaussian model
       *  @param value The constant value of the node.
       *  @return The GaussianConstNode that holds the constant value.
       */
      GaussianConstNode
      gaussian_const(const T value);
      
      /** Create a Gamma Constant.
       *  This is appropriate input for VariableNode's that expect input in Gamma form.
       *  For example, the precision of a Gaussian or the inverse scale of a Gamma distribution.
       *  @param value The constant value of the node.
       *  @return The GammaConstNode that holds the constant value.
       */
      GammaConstNode
      gamma_const(const T value);
      
      /** Create a Gaussian Data Constant.
       *  This is appropriate input for data that is modelled according to a Gaussian/RectififedGaussian (Mixture or not) distribution.
       *  @param data The value of the data.
       *  @return The GaussianDataNode that holds the data.
       */
      GaussianDataNode
      gaussian_data(const T data);


      /** Create a Gamma Data Constant.
       *  This is appropriate input for data that is modelled according to a Gamma distribution.
       *  @param data The value of the data.
       *  @return The GammaDataNode that holds the data.
       */
      GammaDataNode
      gamma_data(const T data);

      ///@}

      /** @name Create a Calculation Node
       */
      ///@{
      /** Create a calculation node.
       *  The calculation node does not make an inference but rather 
       *  multiplies or adds other moments.
       *  For example in Independant Component Analysis the inferred sources 
       *  are multiplied by an inferred Mixing Matrix to generate the modelled data.
       *  This node stores such a product for Gaussian Data types.
       *  @param Expr The generic Expression in which the nodes are multiplied.
       *  @param context The explicit Context (the actual nodes involved) in
       *  which the calculation takes place.
       *  @return The GaussianResultsNode that holds the calculated Moments.
       */
      GaussianResultNode
      calc_gaussian(Expression<T>* Expr,  Context<T>& context);
      ///@}
      
      /** @name Join Existing Nodes
       */
      ///@{

      /** Join Gamma Data with a shape and inverse scale.
       *  @param shape The value of the shape.
       *  @param IScale The node that infers the inverse scale.
       *  @param data The value of the data 
       */
      void 
      join(T& shape, GammaNode IScale ,  const T& data );

      /** Join Gamma Data with a shape and inverse scale.
       *  @param shape The value of the shape.
       *  @param iscale The value of the inverse scale
       *  @param Data The Data node that holds the data.
       */
      void 
      join(T& shape, T& iscale, GammaDataNode Data  );

      /** Join Gamma Data with a shape and inverse scale.
       *  @param shape The value of the shape.
       *  @param IScale The node that infers the inverse scale.
       *  @param Data The Data node that holds the data.
       */
      void 
      join(T& shape, GammaNode IScale, GammaDataNode Data  );

      /** Join Gaussian Data with the mean and precision.
       *  @param mean The value of the mean.
       *  @param precision The value of the precision
       *  @param Data The GaussianDataNode that holds the data.
       */
      void 
      join(T& mean, T& precision, GaussianDataNode Data  );
      
      /** Join Gaussian Data with the mean and precision.
       *  @param Mean The VariableNode that models the Mean.
       *  @param precision The value of the precision
       *  @param Data The GaussianDataNode that holds the data.
       */
      void 
      join(Variable Mean, T& precision, GaussianDataNode Data  );
      
      /** Join Gaussian Data with the mean and precision.
       *  @param mean The value of the mean.
       *  @param Precision The VariableNode that models the Precision.
       *  @param Data The GaussianDataNode that holds the data.
       */
      void 
      join(T& mean, GammaNode Precision, GaussianDataNode Data  );
      
      /** Join Gaussian Data with the mean and precision.
       *  @param Mean The VariableNode that models the mean.
       *  @param Precision The VariableNode that models the Precision.
       *  @param Data The GaussianDataNode that holds the data.
       */
      void 
      join(Variable Mean, GammaNode& Precision, GaussianDataNode& Data  );
 
      /** Join Gaussian Data with the mean and precision.
       *  @param Mean The VariableNode that models the mean.
       *  @param Precision The VariableNode that models the Precision.
       *  @param data The value of the data.
       */
      void 
      join(Variable Mean, GammaNode Precision, const T data );


      /** Join Gaussian Modelled data to a mixture model.
       *  @param vMean The vector of VariableNode's that models the mean's of the Gaussian Mixture.
       *  @param Precision The VariableNode that models the common Precision to each Gaussian in the mixture.
       *  @param Weights The Weights node that stores the weights.
       *  @param data The value of the data.
       */
      void
      join( std::vector<Variable>& vMean, 
	    GammaNode& Precision, 
	    WeightsNode Weights,
	    const T data );

      /** Join Gaussian Modelled data to a mixture model.
       *  @param vMean The vector of VariableNode's that models the mean's of the Gaussian Mixture.
       *  @param vPrecision The vector of VariableNode's that models the  Precision to each Gaussian in the mixture.
       *  @param Weights The Weights node that stores the weights.
       *  @param data The value of the data.
       */
      void
      join( std::vector<Variable>& vMean, 
	    std::vector<Variable>& vPrecision, 
	    WeightsNode Weights,
	    const T data );
	

      ///@}

      /** @name Auxillary Member functions
       */
      ///@{
      
      /** Set the output file for the evidence
       *  @param cost_file The output filename
       */
      void
      set_cost_file(const std::string& cost_file);

      /** The number of variable nodes in the model.
       *  @return The number of variable nodes used in the model.
       */
      size_t
      number_of_nodes() const;

      /** The number of factor nodes in the model.
       *  @return The number of factor nodes used in the model.
       */
      size_t
      number_of_factors() const;
      
      ///@}
      
      
      /** @name Run the inference.
       */
      ///@{
      
      /** Run the inference.
       *  @param epsilon The percentage difference in the cost (per data point) for convergence.
       *    The model will stop running once the increase in the cost reduced to this threshold.
       *  @param max_iterations The maximum number of iterations.
       *   The model will stop running once this number of iterations has been surpassed.
       *  @return Whether the convergance criterium was met.  
       *   If false, then the number of iterations exceeded the maximum.
       *  
       */
      bool
      run(const double& epsilon = 1e-6, const size_t& max_iterations = 100);
      
      ///@}
    private:
      double
      iterate();

      bool
      HasConverged(const T Cost, const T epsilon);
      
      void
      Initialise();

      T m_PrevCost;
      std::vector<boost::shared_ptr<FactorNode<T> > > m_Factors;
      std::vector<boost::shared_ptr<VariableNode<T> > > m_Nodes;
      bool m_initialised;
      size_t m_data_nodes;
      std::string m_cost_file;
    };

  }
}


/********************************************************
 ********************************************************
 ** Implementation
 ********************************************************
 ********************************************************/




template<class T>
template<class MeanIterator, class PrecIterator>
typename ICR::EnsembleLearning::Builder<T>::GaussianNode
ICR::EnsembleLearning::Builder<T>::gaussian_mixture(const MeanIterator& MeanBegin,
						    const PrecIterator& PrecBegin,
						    WeightsNode Weights)
{
  const size_t number = Weights->size();
  //put in vector form
  std::vector<Variable> vMean;
  std::vector<Variable> vPrec;
	
  std::copy(MeanBegin,MeanBegin+number, vMean.begin());
  std::copy(PrecBegin,PrecBegin+number, vPrec.begin());
  
  return gaussian_mixture(vMean, vPrec, Weights);
}

template<class T>
template<class MeanIterator, class PrecIterator>	
typename ICR::EnsembleLearning::Builder<T>::RectifiedGaussianNode
ICR::EnsembleLearning::Builder<T>::rectified_gaussian_mixture(const MeanIterator& MeanBegin,
							      const PrecIterator& PrecBegin,
							      WeightsNode Weights)
{
  
  const size_t number = Weights->size();
  //put in vector form
  std::vector<Variable> vMean;
  std::vector<Variable> vPrec;
	
  std::copy(MeanBegin,MeanBegin+number, vMean.begin());
  std::copy(PrecBegin,PrecBegin+number, vPrec.begin());
  
  return rectified_gaussian_mixture(vMean, vPrec, Weights);
}



#endif //BUILDER_HPP guard
