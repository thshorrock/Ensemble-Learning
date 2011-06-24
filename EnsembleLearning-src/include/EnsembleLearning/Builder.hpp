
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

#include "EnsembleLearning/node/factor/Factor.hpp"
#include "EnsembleLearning/node/factor/Mixture.hpp"
#include "EnsembleLearning/detail/MixtureVector.hpp"
#include "EnsembleLearning/detail/TernaryOp.hpp"

#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/arithmetic/sub.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>
#include <boost/preprocessor/iteration/local.hpp>
#  include <boost/preprocessor/control/if.hpp>
#  include <boost/preprocessor/comparison/equal.hpp>
#  include <boost/preprocessor/facilities/empty.hpp>

#include <boost/fusion/algorithm/iteration/for_each.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/container/generation/make_vector.hpp>
#include <boost/fusion/include/make_vector.hpp>
#include <boost/fusion/sequence/intrinsic/at.hpp>
#include <boost/fusion/include/at.hpp>


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
    //template<class T> class FactorNode;
    

    //Forward declarations of Factors
    namespace detail{
      // template<template<class> class Model, class T,
      // 	        class parent1_t, class parent2_t, class child_t> class 
      // Factor;
      
      // template<class Model, class T
      // 	       //,class parent1_t, class parent2_t, class child_t
      // 	       > class Mixture;
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

      // typedef detail::Factor<RectifiedGaussian, T >      RectifiedGaussianFactor;
      // typedef detail::Factor<Gaussian, T >      GaussianFactor;
      // typedef detail::Factor<Gamma, T >         GammaFactor;
      // typedef detail::Factor<Dirichlet, T >     DirichletFactor;
      // typedef detail::Factor<Discrete, T >      DiscreteFactor;
      //typedef detail::Mixture<RectifiedGaussian<T>, T >     RectifiedGaussianMixtureFactor;
      //typedef detail::Mixture<Gaussian<T>, T >     GaussianMixtureFactor;
      
      typedef detail::Deterministic<Gaussian, T >  DeterministicFactor;

      typedef HiddenNode<Gaussian, T >      GaussianType;
      typedef HiddenNode<RectifiedGaussian, T >      RectifiedGaussianType;
      typedef HiddenNode<Gamma, T >         GammaType;
      typedef HiddenNode<Dirichlet, T >     DirichletType;
      typedef ObservedNode<Gaussian, T, detail::TypeList::zeros >     GaussianDataType;
      typedef ObservedNode<Gamma, T, detail::TypeList::zeros >        GammaDataType;
      typedef ObservedNode<Gaussian, T, detail::TypeList::zeros > GaussianConstType;
      typedef ObservedNode<Gamma, T, detail::TypeList::zeros >    GammaConstType;
      typedef ObservedNode<Gaussian, T, detail::TypeList::zeros >   NormalConstType;
      typedef ObservedNode<Dirichlet, T, detail::TypeList::zeros >  DirichletConstType;
      typedef DeterministicNode<Gaussian<T>, T>    GaussianResultType;

      
      typedef HiddenNode<Dirichlet, T >     WeightsType;
      typedef HiddenNode<Discrete, T >      CatagoryType;
    public:
      
      /** @name Useful typdefs for types that are exposed to the user.
       */
      ///@{
      typedef VariableNode<T>*                    Variable;
      typedef HiddenNode<RectifiedGaussian, T >*      RectifiedGaussianNode;
      typedef HiddenNode<Gaussian, T >*      GaussianNode;
      typedef HiddenNode<Gamma, T >*         GammaNode;
      typedef ObservedNode<Gaussian, T, detail::TypeList::zeros >*     GaussianDataNode;
      typedef ObservedNode<Gamma, T, detail::TypeList::zeros >*        GammaDataNode;
      typedef ObservedNode<Gaussian, T, detail::TypeList::zeros >* GaussianConstNode;
      typedef ObservedNode<Gamma, T, detail::TypeList::zeros >*    GammaConstNode;
      typedef DeterministicNode<Gaussian<T>, T>*    GaussianResultNode;
      
      typedef HiddenNode<Dirichlet, T >*      WeightsNode;
      typedef HiddenNode<Discrete, T >*       CatagoryNode;
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
      template<template<template<class> class,class,class,class> class MeanType,
	       class MeanId, class MeanEnabler,
	       template<template<class> class,class,class,class> class PrecType,
	       class PrecId, class PrecEnabler>
      HiddenNode<Gaussian, T, typename detail::TypeList::incr_id<MeanId>::type >*
      gaussian(MeanType<Gaussian,T,MeanId,MeanEnabler>* Mean, 
	       PrecType<Gamma,T,PrecId,PrecEnabler>   * Precision);

    
      template<template<template<class> class,class,class,class> class MeanType,
	       class MeanId, class MeanEnabler,
	       template<template<class> class,class,class,class> class PrecType,
	       class PrecId, class PrecEnabler>
      HiddenNode<Gaussian, T, typename detail::TypeList::incr_id<MeanId>::type >*
      gaussian(MeanType<RectifiedGaussian,T,MeanId,MeanEnabler>* Mean, 
	       PrecType<Gamma,T,PrecId,PrecEnabler>   * Precision);
      /** Create a Gaussian VariableNode.
       * @param Mean The node that stores the Momets for the Mean.
       *  The model for the Mean node must be a Gaussian or Rectified Gaussian.
       * @param precision The value of the precision.
       * @return The GaussianNode that holds the inferred Gaussian Moments.
       */
      
      template<template<template<class> class,class,class,class> class MeanType,
	       class MeanId, class MeanEnabler>
      HiddenNode<Gaussian, T, typename detail::TypeList::incr_id<MeanId>::type >*
      gaussian(MeanType<Gaussian,T,MeanId,MeanEnabler>* Mean, const T& precision);

      
      template<template<template<class> class,class,class,class> class MeanType,
	       class MeanId, class MeanEnabler>
      HiddenNode<Gaussian, T, typename detail::TypeList::incr_id<MeanId>::type >*
      gaussian(MeanType<RectifiedGaussian,T,MeanId,MeanEnabler>* Mean, const T& precision);

      /** Create a Gaussian VariableNode.
       * @param mean The value of the mean.
       *  The model for the Mean node must be a Gaussian or Rectified Gaussian.
       * @param Precision The node that stores that Moments for the Precicision.
       *  The model for the Precision node must be a Gamma.
       * @return The GaussianNode that holds the inferred Gaussian Moments.
       */
    
      template<template<template<class> class,class,class,class> class PrecType,
	       class PrecId, class PrecEnabler>
      HiddenNode<Gaussian,T,typename detail::TypeList::incr_id<detail::TypeList::zeros>::type >*
      gaussian(const T& mean, PrecType<Gamma,T,PrecId,PrecEnabler>* Precision);

      /** Create a Gaussian VariableNode.
       * @param mean The value of the mean.
       * @param precision The value of the precision.
       * @return The GaussianNode that holds the inferred Gaussian Moments.
       */
      HiddenNode<Gaussian,T,typename detail::TypeList::incr_id<detail::TypeList::zeros >::type >*
      gaussian(const T mean, const T precision);
      
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
      

      template<template<template<class> class,class,class,class> class MeanType,
	       class MeanId, class MeanEnabler,
	       template<template<class> class,class,class,class> class PrecType,
	       class PrecId, class PrecEnabler>
      HiddenNode<RectifiedGaussian, T, typename detail::TypeList::incr_id<MeanId>::type >*
      rectified_gaussian(MeanType<Gaussian,T,MeanId,MeanEnabler>* Mean, 
			 PrecType<Gamma,T,PrecId,PrecEnabler>   * Precision);


      template<template<template<class> class,class,class,class> class MeanType,
	       class MeanId, class MeanEnabler,
	       template<template<class> class,class,class,class> class PrecType,
	       class PrecId, class PrecEnabler>
      HiddenNode<RectifiedGaussian, T, typename detail::TypeList::incr_id<MeanId>::type >*
      rectified_gaussian(MeanType<RectifiedGaussian,T,MeanId,MeanEnabler>* Mean, 
			 PrecType<Gamma,T,PrecId,PrecEnabler>   * Precision);

      /** Create a RectifiedGaussian VariableNode.
       * @param Mean The node that stores the Momets for the Mean.
       *  The model for the Mean node must be a Gaussian or Rectified Gaussian.
       * @param precision The value of the precision.
       * @return The GaussianNode that holds the inferred Rectified Gaussian Moments.
       */
      template<template<template<class> class,class,class,class> class MeanType,
	       class MeanId, class MeanEnabler>
      HiddenNode<RectifiedGaussian, T, typename detail::TypeList::incr_id<MeanId>::type >*
      rectified_gaussian(MeanType<Gaussian,T,MeanId,MeanEnabler>* Mean, const T& precision);
      
      template<template<template<class> class,class,class,class> class MeanType,
	       class MeanId, class MeanEnabler>
      HiddenNode<RectifiedGaussian, T, typename detail::TypeList::incr_id<MeanId>::type >*
      rectified_gaussian(MeanType<RectifiedGaussian,T,MeanId,MeanEnabler>* Mean, const T& precision);
       /** Create a RectifiedGaussian VariableNode.
       * @param mean The value of the mean.
       *  The model for the Mean node must be a Gaussian or Rectified Gaussian.
       * @param Precision The node that stores that Moments for the Precicision.
       *  The model for the Precision node must be a Gamma.
       * @return The GaussianNode that holds the inferred Rectified Gaussian Moments.
       */
      
    
      template<template<template<class> class,class,class,class> class PrecType,
	       class PrecId, class PrecEnabler>
      HiddenNode<RectifiedGaussian,T, typename detail::TypeList::incr_id<detail::TypeList::zeros>::type >*
      rectified_gaussian(const T& mean,PrecType<Gamma,T,PrecId,PrecEnabler>*  Precision);
  
    /** Create a RectifiedGaussian VariableNode.
     * @param mean The value of the mean.
     * @param precision The value of the precision.
       * @return The GaussianNode that holds the inferred Rectified Gaussian Moments.
       */
      
      HiddenNode<RectifiedGaussian,T, typename detail::TypeList::incr_id<detail::TypeList::zeros>::type>*
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

      

#       define ENSEMBLE_LEARNING_VEC_print(z, n, data) data

      template<template<class> class Model, class fusion_vector1, class fusion_vector2>
      MixtureVector<Model,T>
      mixture_vector(
		     const fusion_vector1& v0,
		     const fusion_vector2& v1
		     )
      {
	typedef  MixtureVector<Model,T>
      	  MV_t;
	
	MV_t MV;

	typedef fusion_vector1 p0_t;
	typedef fusion_vector2 p1_t;
	BOOST_MPL_ASSERT_RELATION(boost::fusion::result_of::size<p0_t>::type::value, ==,boost::fusion::result_of::size<p1_t>::type::value);
	const size_t size = boost::fusion::result_of::size<p0_t>::type::value;
	
	// detail::TernaryOp(v0,v1,MV.data(), MakeMixtureModel<Model>(m_Factors, m_Nodes));
	// //boost::fusion::transform( v0,v1, MakeMixtureModel() );

	//create a new set of nodes.
	//(Iterate with the preprocessor)
#       include <boost/preprocessor/iteration/iterate.hpp>
#       define BOOST_PP_ITERATION_LIMITS (0, ENSEMBLE_LEARNING_COMPONENTS - 1)
#       define BOOST_PP_FILENAME_1       "EnsembleLearning/detail/MixtureVector_create.hpp"
#       include BOOST_PP_ITERATE()
	
	return MV;
      }


      //Gaussian or Rectified Gaussian case.
      template<template<class> class Model>
      MixtureVector<Model,T> 
      mixture_vector(const std::vector<T>& v1, 
      		     const std::vector<T>& v2)
      {

	typedef ObservedNode<Gaussian, T, detail::TypeList::zeros> p0_t;
	typedef ObservedNode<Gamma,    T, detail::TypeList::zeros> p1_t;
	const size_t size = ENSEMBLE_LEARNING_COMPONENTS;
	
	// std::vector<p0_t*> p_0(size);
	// std::vector<p1_t*> p_1(size);
	// for(size_t i=0;i<size;++i){
	//   boost::shared_ptr<p0_t> tmp0(new p0_t(v1[i]));
	//   boost::shared_ptr<p1_t> tmp1(new p1_t(v2[i]));
	//   //store the reference so it doesn't get deleted till end of execution.
	//   m_Nodes.push_back(tmp0);
	//   m_Nodes.push_back(tmp1);
	  
	//   //store so that references can be put easily into fusion vector
	//   p_0[i] = tmp0.get();
	//   p_1[i] = tmp1.get();
	// }
	//The fusion vectors
	// Vector<ObservedNode,Gaussian,T,size> vt0;
	// Vector<ObservedNode,Gamma,T,size> vt1;
	// //Assign
	// boost::fusion::for_each(vt0.data(), AssignVector<p0_t>(p_0));
	// boost::fusion::for_each(vt1.data(), AssignVector<p1_t>(p_1));
	
	// //pass on the pointers
	// return mixture_vector<Model>(vt0.data(),vt1.data());
				
	
	//Create a set of priors
	//(Iterate with the preprocessor)
#       include <boost/preprocessor/iteration/iterate.hpp>
#       define BOOST_PP_ITERATION_LIMITS (0, ENSEMBLE_LEARNING_COMPONENTS - 1)
#       define BOOST_PP_FILENAME_1       "EnsembleLearning/detail/MixtureVector_create_priors.hpp"
#       include BOOST_PP_ITERATE()
	

	//place into a boost fusion vector.
	//shared version
	boost::fusion::vector<
	BOOST_PP_ENUM(BOOST_PP_SUB(ENSEMBLE_LEARNING_COMPONENTS,0), ENSEMBLE_LEARNING_VEC_print,boost::shared_ptr<p0_t>)
	    > shared_vt0(
		   BOOST_PP_ENUM_PARAMS(ENSEMBLE_LEARNING_COMPONENTS,shared_p0_) 
		 );

	boost::fusion::vector<
	BOOST_PP_ENUM(BOOST_PP_SUB(ENSEMBLE_LEARNING_COMPONENTS,0), ENSEMBLE_LEARNING_VEC_print,boost::shared_ptr<p1_t>)
	    > shared_vt1(
		   BOOST_PP_ENUM_PARAMS(ENSEMBLE_LEARNING_COMPONENTS,shared_p1_) 
		 );
	//pointers to shared vectors
	boost::fusion::vector<
	BOOST_PP_ENUM(BOOST_PP_SUB(ENSEMBLE_LEARNING_COMPONENTS,0), ENSEMBLE_LEARNING_VEC_print,p0_t*)
	    > vt0(
		   BOOST_PP_ENUM_PARAMS(ENSEMBLE_LEARNING_COMPONENTS,p0_) 
		 );

	boost::fusion::vector<
	BOOST_PP_ENUM(BOOST_PP_SUB(ENSEMBLE_LEARNING_COMPONENTS,0), ENSEMBLE_LEARNING_VEC_print,p1_t*)
	    > vt1(
		   BOOST_PP_ENUM_PARAMS(ENSEMBLE_LEARNING_COMPONENTS,p1_) 
		 );

	//copy the shared pointers to container so they don't get deleted.
	boost::fusion::for_each(shared_vt0, boost::bind(&std::vector<boost::shared_ptr<VariableNode<T> > >::push_back, 
						 boost::ref(m_Nodes), 
						 _1));

	boost::fusion::for_each(shared_vt1, boost::bind(&std::vector<boost::shared_ptr<VariableNode<T> > >::push_back, 
						 boost::ref(m_Nodes), 
						 _1));

	//pass on the pointers
	return mixture_vector<Model>(vt0,vt1);

      }

      template<template<class> class Model>
      MixtureVector<Model,T>
      mixture_vector(const T d1, 
      		     const T d2)
      {
	
	std::vector<T> v1(ENSEMBLE_LEARNING_COMPONENTS,d1);
	std::vector<T> v2(ENSEMBLE_LEARNING_COMPONENTS,d2);
	
	return mixture_vector<Model>(v1,v2);

      }


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
      template<class vMean_t, class vPrec_t>
      GaussianNode
      gaussian_mixture( vMean_t& vMean, 
		        vPrec_t& vPrecision, 
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
      // template<class MeanIterator, class PrecIterator>
      // GaussianNode
      // gaussian_mixture(const MeanIterator& MeanBegin, 
      // 		       const PrecIterator& PrecisionBegin, 
      // 		       WeightsNode Weights);
	

  
      /** Rectified Gaussian Mixture Model.
       *  @param MeanBegin The iterator at the beginning of the the container containing all the  variables representing the means.
       *  @param PrecisionBegin The iterator at the beginning of the the container containing all the  variables representing the precisions.
       *  @param Weights A Weights variable that stores the weights to the means and precisions in vMean and vPrecision.
       *  @attention The size of vMean and vPrecision must be identical, 
       *   and must be the same as the size of the Weights Node.
       * @return The GaussianNode that holds the inferred Rectified Gaussian Moments.
       */
      // template<class MeanIterator, class PrecIterator>
      // RectifiedGaussianNode
      // rectified_gaussian_mixture(const MeanIterator& MeanBegin, 
      // 				 const PrecIterator& PrecisionBegin,
      // 				 WeightsNode Weights);
      
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
         
      template<template<template<class> class,class,class,class> class MeanType,
	       class MeanId,
	       class Enabler>
      void 
      join(MeanType<Gaussian,T,MeanId,Enabler>* Mean, T& precision, GaussianDataNode Data  );
      
    template<template<template<class> class,class,class,class> class MeanType,
	     class MeanId,
	     class Enabler>
    void 
    join(MeanType<RectifiedGaussian,T,MeanId,Enabler>* Mean, T& precision, GaussianDataNode Data  );

  /** Join Gaussian Data with the mean and precision.
   *  @param mean The value of the mean.
   *  @param Precision The VariableNode that models the Precision.
       *  @param Data The GaussianDataNode that holds the data.
       */
      template<template<template<class> class,class,class,class> class PrecType,
	       class PrecId,
	       class PrecEnabler>
      void 
      join(T& mean,  PrecType<Gamma,T,PrecId,PrecEnabler>* Precision, GaussianDataNode Data  );
      
      /** Join Gaussian Data with the mean and precision.
       *  @param Mean The VariableNode that models the mean.
       *  @param Precision The VariableNode that models the Precision.
       *  @param Data The GaussianDataNode that holds the data.
       */
      // template<template<template<class> class,class,class,class> class MeanType,
      // 	       template<class> class Model1,
      // 	       class MeanId,
      // 	       class MeanEnabler,
      // 	       template<template<class> class,class,class,class> class PrecType,
      // 	       template<class> class  Model2,
      // 	       class PrecId, class PrecEnabler> 
      // void 
      // join(MeanType<Model1,T,MeanId,MeanEnabler>* Mean, PrecType<Model2,T,PrecId,PrecEnabler>* Precision, GaussianDataNode& Data  );  //needs to be specialised

      template<template<template<class> class,class,class,class> class MeanType,
	       class MeanId,
	       class MeanEnabler,
	       template<template<class> class,class,class,class> class PrecType,
	       class PrecId, class PrecEnabler> 
      void 
      join(MeanType<Gaussian,T,MeanId,MeanEnabler>* Mean, PrecType<Gamma,T,PrecId,PrecEnabler>* Precision, GaussianDataNode& Data  );
 
      template<template<template<class> class,class,class,class> class MeanType,
	       class MeanId,class MeanEnabler,
	       template<template<class> class,class,class,class> class PrecType,
	       class PrecId, class PrecEnabler> 
      void 
      join(MeanType<RectifiedGaussian,T,MeanId,MeanEnabler>* Mean, PrecType<Gamma,T,PrecId,PrecEnabler>* Precision, GaussianDataNode& Data  );
      /** Join Gaussian Data with the mean and precision.
       *  @param Mean The VariableNode that models the mean.
       *  @param Precision The VariableNode that models the Precision.
       *  @param data The value of the data.
       */
	   template<template<template<class> class,class,class,class> class MeanType,
	   class MeanId,class MeanEnabler,
	   template<template<class> class,class,class,class> class PrecType,
	   class PrecId, class PrecEnabler> 
	   void 
      join(MeanType<Gaussian,T,MeanId,MeanEnabler>* Mean, PrecType<Gamma,T,PrecId,PrecEnabler>* Precision, const T data );

      // template<template<template<class> class,class,class,class> class MeanType,
      // 	       class MeanId,class MeanEnabler,
      // 	       template<template<class> class,class,class,class> class PrecType,
      // 	       class PrecId, class PrecEnabler> 
      // void 
      // 	   join(MeanType<Gaussian,T,MeanId,MeanEnabler>* Mean, 
      // 		PrecType<Gamma,T,PrecId,PrecEnabler>* Precision, 
      // 		const T data );
	   
      template<template<template<class> class,class,class,class> class MeanType,
	       class MeanId,class MeanEnabler,
	       template<template<class> class,class,class,class> class PrecType,
	       class PrecId, class PrecEnabler> 
      void 
      join(MeanType<RectifiedGaussian,T,MeanId,MeanEnabler>* Mean, PrecType<Gamma,T,PrecId,PrecEnabler>* Precision, const T data );


      /** Join Gaussian Modelled data to a mixture model.
       *  @param vMean The vector of VariableNode's that models the mean's of the Gaussian Mixture.
       *  @param Precision The VariableNode that models the common Precision to each Gaussian in the mixture.
       *  @param Weights The Weights node that stores the weights.
       *  @param data The value of the data.
       */
      // template<class vMean_t>
      // void
      // join( vMean_t& vMean, 
      // 	    GammaNode& Precision, 
      // 	    WeightsNode Weights,
      // 	    const T data );

      /** Join Gaussian Modelled data to a mixture model.
       *  @param vMean The vector of VariableNode's that models the mean's of the Gaussian Mixture.
       *  @param vPrecision The vector of VariableNode's that models the  Precision to each Gaussian in the mixture.
       *  @param Weights The Weights node that stores the weights.
       *  @param data The value of the data.
       */
	   template<class vMean_t, class vPrec_t>
	   void
	   join( vMean_t& vMean, 
		 vPrec_t& vPrecision, 
		 WeightsNode Weights,
		 const T data );
	
	   /** Join Gaussian Modelled data to a mixture model.
	    *  @param MeanBegin The iterator at the beginning of the the container containing all the  variables representing the means.
	    *  @param PrecBegin The iterator at the beginning of the the container containing all the  variables representing the precisions.
	    *  @param Weights The Weights node that stores the weights.
	    *  @param data The value of the data.
	    */
	   // template<class MeanIterator, class PrecIterator>
	   // void
	   // join( MeanIterator MeanBegin, 
	   // 	    PrecIterator PrecBegin, 
	   // 	    WeightsNode Weights,
	   // 	    const T data );


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
	    *  @param skip Skip n iterations at the beginning while the initial moments are getting flushed.  This value should be the span of the graph.
	    *  @return Whether the convergance criterium was met.  
	    *   If false, then the number of iterations exceeded the maximum.
	    *  
	    */
	   bool
	   run(const double& epsilon = 1e-6, const size_t& max_iterations = 100, size_t skip = 1);

	   /** Reset all the moments based on their parents current variables.
	    *  @attention This is an experimental feature,
	    *   it is not recommended that you actually do perturb your variables.
	    */
	   void
	   perturb();
      
	   ///@}
      private:
      
      template<class U>
      struct AssignVector
      {
	AssignVector(const std::vector<U*>& vp) : m_index(0), m_vp(vp) {}
	template <class Node>
	void operator()(Node& node) const
	{
	  node = m_vp[m_index];
	  ++m_index;
	}
	mutable size_t m_index;
	const std::vector<U*>& m_vp;
      };
      
      template<template<class> class Model>
      struct MakeMixtureModel
      {
	MakeMixtureModel(std::vector<boost::shared_ptr<FactorNode_basic> >& F,
			 std::vector<boost::shared_ptr<VariableNode<T> > >& N)
	  : m_F(F), m_N(N)
	{}

	template <class p0_t, class p1_t, class pc_t>
	void operator()(p0_t& p0, p1_t& p1, pc_t& ch) const
	//it0, Itr1 it1, ItrChild itC) const
	{
	  typedef typename boost::remove_pointer<typename boost::remove_reference<p0_t>::type>::type n0_t;
	  typedef typename boost::remove_pointer<typename boost::remove_reference<p1_t>::type>::type n1_t;
	  typedef typename boost::remove_pointer<typename boost::remove_reference<pc_t>::type>::type child_t;
	  ;
	  
	  boost::shared_ptr<child_t> child( new child_t() );
	  
	  typedef detail::Factor<Model,T,n0_t,n1_t,child_t> Factor_t;
	  boost::shared_ptr<Factor_t > F(new Factor_t(p0, 
						      p1, 
						      child.get()));

	  m_F.push_back(F);
	  m_N.push_back(child);
	  ch = child.get();
	}
	mutable std::vector<boost::shared_ptr<FactorNode_basic> >& m_F;
	mutable std::vector<boost::shared_ptr<VariableNode<T> > >& m_N;
      };
	// MV_t& m_MV;
		//};
      
	   double
	   iterate();

	   bool
	   HasConverged(const T Cost, const T epsilon);
      

	   T m_PrevCost;
	   std::vector<boost::shared_ptr<FactorNode_basic> > m_Factors;
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
	   template<template<template<class> class,class,class,class> class MeanType,
	   class MeanId, class MeanEnabler,
	   template<template<class> class,class,class,class> class PrecType,
	   class PrecId, class PrecEnabler>
	   ICR::EnsembleLearning::HiddenNode<ICR::EnsembleLearning::Gaussian,
	   T, typename ICR::EnsembleLearning::detail::TypeList::incr_id<MeanId>::type
	   >*
	   ICR::EnsembleLearning::Builder<T>::gaussian(MeanType<Gaussian,T,MeanId,MeanEnabler>* Mean, 
						       PrecType<Gamma,T,PrecId,PrecEnabler>* Precision)
						       {
  typedef HiddenNode<Gaussian, T, typename detail::TypeList::incr_id<MeanId>::type > c_t;
  typedef MeanType<Gaussian,T,MeanId,MeanEnabler> p0_t;
  typedef PrecType<Gamma,T,PrecId,PrecEnabler>    p1_t;
  typedef detail::Factor<Gaussian,T,p0_t,p1_t,c_t> Factor_t;
  
  boost::shared_ptr<c_t > GaussianN(new c_t());
  boost::shared_ptr<Factor_t > GaussianF(new Factor_t(Mean, Precision, GaussianN.get()));
	
  m_Factors.push_back(GaussianF);
  m_Nodes.push_back(GaussianN);
	
  return GaussianN.get();
}

template<class T>
template<template<template<class> class,class,class,class> class MeanType,
	 class MeanId, class MeanEnabler,
	 template<template<class> class,class,class,class> class PrecType,
	 class PrecId, class PrecEnabler>
ICR::EnsembleLearning::HiddenNode<ICR::EnsembleLearning::Gaussian,
				  T,typename  ICR::EnsembleLearning:: detail::TypeList::incr_id<MeanId>::type
				  >*
ICR::EnsembleLearning::Builder<T>::gaussian(MeanType<RectifiedGaussian,T,MeanId,MeanEnabler>* Mean, 
					    PrecType<Gamma,T,PrecId,PrecEnabler>* Precision)
{
  typedef HiddenNode<Gaussian, T, typename detail::TypeList::incr_id<MeanId>::type > c_t;
  typedef MeanType<RectifiedGaussian,T,MeanId,MeanEnabler> p0_t;
  typedef PrecType<Gamma,T,PrecId,PrecEnabler>    p1_t;
  typedef detail::Factor<Gaussian,T,p0_t,p1_t,c_t> Factor_t;
  
  boost::shared_ptr<c_t > GaussianN(new c_t());
  boost::shared_ptr<Factor_t > GaussianF(new Factor_t(Mean, Precision, GaussianN.get()));
	
  m_Factors.push_back(GaussianF);
  m_Nodes.push_back(GaussianN);
	
  return GaussianN.get();
}


template<class T>
template<template<template<class> class,class,class,class> class MeanType,
	 class MeanId, class MeanEnabler>
ICR::EnsembleLearning::HiddenNode<ICR::EnsembleLearning::Gaussian,
				  T, 
				  typename ICR::EnsembleLearning::detail::TypeList::incr_id<MeanId>::type >*
ICR::EnsembleLearning::Builder<T>::gaussian(MeanType<Gaussian,T,MeanId,MeanEnabler>* Mean, 
					    const T& precision)
{
  typedef ObservedNode<Gamma, T, ICR::EnsembleLearning::detail::TypeList::zeros > p1_t;
  boost::shared_ptr<p1_t > Precision(new p1_t(precision));
  m_Nodes.push_back(Precision);
  return gaussian(Mean,Precision.get());
}

template<class T>
template<template<template<class> class,class,class,class> class MeanType,
	 class MeanId, class MeanEnabler>
ICR::EnsembleLearning::HiddenNode<ICR::EnsembleLearning::Gaussian,
				  T, 
				  typename ICR::EnsembleLearning::detail::TypeList::incr_id<MeanId>::type>*
ICR::EnsembleLearning::Builder<T>::gaussian(MeanType<RectifiedGaussian,T,MeanId,MeanEnabler>* Mean, 
					    const T& precision)
{
  typedef ObservedNode<Gamma, T, ICR::EnsembleLearning::detail::TypeList::zeros > p1_t;
  boost::shared_ptr<p1_t > Precision(new p1_t(precision));
  m_Nodes.push_back(Precision);
  return gaussian(Mean,Precision.get());
}

    
template<class T>  
template<template<template<class> class,class,class,class> class PrecType,
	 class PrecId, class PrecEnabler>
ICR::EnsembleLearning::HiddenNode<ICR::EnsembleLearning::Gaussian,
				  T,
				   typename ICR::EnsembleLearning::detail::TypeList::incr_id<ICR::EnsembleLearning::detail::TypeList::zeros>::type
				  >*
ICR::EnsembleLearning::Builder<T>::gaussian(const T& mean, 
					    PrecType<Gamma,T,PrecId,PrecEnabler>* Precision)
{
  
  typedef ObservedNode<Gaussian, T, ICR::EnsembleLearning::detail::TypeList::zeros > p0_t;

  boost::shared_ptr<p0_t> Mean(new p0_t(mean));
  m_Nodes.push_back(Mean);
  return gaussian(Mean.get(),Precision);
}
  
template<class T>
inline
ICR::EnsembleLearning::HiddenNode<ICR::EnsembleLearning::Gaussian,
				  T, 
				  typename ICR::EnsembleLearning::detail::TypeList::incr_id<ICR::EnsembleLearning::detail::TypeList::zeros>::type
				  >*
ICR::EnsembleLearning::Builder<T>::gaussian(const T mean, const T precision)
{
	
  typedef ObservedNode<Gaussian, T, detail::TypeList::zeros > p0_t;
  typedef ObservedNode<Gamma, T, detail::TypeList::zeros >    p1_t;
  boost::shared_ptr<p0_t> Mean(new p0_t(mean));
  boost::shared_ptr<p1_t> Precision(new p1_t(precision));
	
  m_Nodes.push_back(Mean);
  m_Nodes.push_back(Precision);
	
  return gaussian(Mean.get(),Precision.get());
}
      

//Rectified Gaussian



template<class T>
template<template<template<class> class,class,class,class> class MeanType,
					    class MeanId, class MeanEnabler,
	 template<template<class> class,class,class,class> class PrecType,
	 class PrecId, class PrecEnabler>
ICR::EnsembleLearning::HiddenNode<ICR::EnsembleLearning::RectifiedGaussian,
				  T, typename ICR::EnsembleLearning::detail::TypeList::incr_id<MeanId>::type
				  >*
					    ICR::EnsembleLearning::Builder<T>::rectified_gaussian(MeanType<Gaussian,T,MeanId,MeanEnabler>* Mean, 
						      PrecType<Gamma,T,PrecId,PrecEnabler>* Precision)
{
  typedef HiddenNode<RectifiedGaussian, T, typename detail::TypeList::incr_id<MeanId>::type > c_t;
  typedef MeanType<Gaussian,T,MeanId,MeanEnabler> p0_t;
  typedef PrecType<Gamma,T,PrecId,PrecEnabler>    p1_t;
  typedef detail::Factor<RectifiedGaussian,T,p0_t,p1_t,c_t> Factor_t;
  
  boost::shared_ptr<c_t > RectifiedGaussianN(new c_t());
  boost::shared_ptr<Factor_t > RectifiedGaussianF(new Factor_t(Mean, Precision, RectifiedGaussianN.get()));
	
  m_Factors.push_back(RectifiedGaussianF);
  m_Nodes.push_back(RectifiedGaussianN);
	
  return RectifiedGaussianN.get();
}

template<class T>
template<template<template<class> class,class,class,class> class MeanType,
	 class MeanId,class MeanEnabler,
	 template<template<class> class,class,class,class> class PrecType,
	 class PrecId, class PrecEnabler>
ICR::EnsembleLearning::HiddenNode<ICR::EnsembleLearning::RectifiedGaussian,
				  T, typename ICR::EnsembleLearning::detail::TypeList::incr_id<MeanId>::type
				  >*
ICR::EnsembleLearning::Builder<T>::rectified_gaussian(MeanType<RectifiedGaussian,T,MeanId,MeanEnabler>* Mean, 
						      PrecType<Gamma,T,PrecId,PrecEnabler>* Precision)
{
  typedef HiddenNode<RectifiedGaussian, T, typename detail::TypeList::incr_id<MeanId>::type > c_t;
  typedef MeanType<RectifiedGaussian,T,MeanId,MeanEnabler> p0_t;
  typedef PrecType<Gamma,T,PrecId,PrecEnabler>    p1_t;
  typedef detail::Factor<RectifiedGaussian,T,p0_t,p1_t,c_t> Factor_t;
  
  boost::shared_ptr<c_t > RectifiedGaussianN(new c_t());
  boost::shared_ptr<Factor_t > RectifiedGaussianF(new Factor_t(Mean, Precision, RectifiedGaussianN.get()));
	
  m_Factors.push_back(RectifiedGaussianF);
  m_Nodes.push_back(RectifiedGaussianN);
	
  return RectifiedGaussianN.get();
}


template<class T>
template<template<template<class> class,class,class,class> class MeanType,
	 class MeanId, class MeanEnabler>
ICR::EnsembleLearning::HiddenNode<ICR::EnsembleLearning::RectifiedGaussian,
				  T, 
				  typename ICR::EnsembleLearning::detail::TypeList::incr_id<MeanId>::type >*
ICR::EnsembleLearning::Builder<T>::rectified_gaussian(MeanType<Gaussian,T,MeanId,MeanEnabler>* Mean, 
						      const T& precision)
{
  typedef ObservedNode<Gamma, T, ICR::EnsembleLearning::detail::TypeList::zeros > p1_t;
  boost::shared_ptr<p1_t > Precision(new p1_t(precision));
  m_Nodes.push_back(Precision);
  return rectified_gaussian(Mean,Precision.get());
}

template<class T>
template<template<template<class> class,class,class,class> class MeanType,
	 class MeanId, class MeanEnabler>
ICR::EnsembleLearning::HiddenNode<ICR::EnsembleLearning::RectifiedGaussian,
				  T, 
				  typename ICR::EnsembleLearning::detail::TypeList::incr_id<MeanId>::type >*
ICR::EnsembleLearning::Builder<T>::rectified_gaussian(MeanType<RectifiedGaussian,T,MeanId,MeanEnabler>* Mean, 
						      const T& precision)
{
  typedef ObservedNode<Gamma, T, ICR::EnsembleLearning::detail::TypeList::zeros > p1_t;
  boost::shared_ptr<p1_t > Precision(new p1_t(precision));
  m_Nodes.push_back(Precision);
  return rectified_gaussian(Mean,Precision.get());
}

    
template<class T>  
template<template<template<class> class,class,class,class> class PrecType,
	 class PrecId, class PrecEnabler>
ICR::EnsembleLearning::HiddenNode<ICR::EnsembleLearning::RectifiedGaussian,
				  T, 
				   typename ICR::EnsembleLearning::detail::TypeList::incr_id<ICR::EnsembleLearning::detail::TypeList::zeros>::type
				  >*
ICR::EnsembleLearning::Builder<T>::rectified_gaussian(const T& mean, 
						      PrecType<Gamma,T,PrecId,PrecEnabler>* Precision)
{
  
  typedef ObservedNode<RectifiedGaussian, T, ICR::EnsembleLearning::detail::TypeList::zeros > p0_t;

  boost::shared_ptr<p0_t> Mean(new p0_t(mean));
  m_Nodes.push_back(Mean);
  return rectified_gaussian(Mean.get(),Precision);
}
  
template<class T>
ICR::EnsembleLearning::HiddenNode<ICR::EnsembleLearning::RectifiedGaussian,
				  T, 
				   typename ICR::EnsembleLearning::detail::TypeList::incr_id<ICR::EnsembleLearning::detail::TypeList::zeros>::type
				  >*
ICR::EnsembleLearning::Builder<T>::rectified_gaussian(const T& mean, const T& precision)
{
	
  typedef ObservedNode<RectifiedGaussian, T, ICR::EnsembleLearning::detail::TypeList::zeros > p0_t;
  typedef ObservedNode<Gamma, T, ICR::EnsembleLearning::detail::TypeList::zeros >    p1_t;
  boost::shared_ptr<p0_t> Mean(new p0_t(mean));
  boost::shared_ptr<p1_t> Precision(new p1_t(precision));
	
  m_Nodes.push_back(Mean);
  m_Nodes.push_back(Precision);
	
  return rectified_gaussian(mean,Precision.get());
}
     














//Join Gaussian


template<class T>   
template<template<template<class> class,class,class,class> class MeanType,
	 class MeanId, class MeanEnabler>
void 
ICR::EnsembleLearning::Builder<T>::join(MeanType<Gaussian,T,MeanId,MeanEnabler>* Mean, T& precision, GaussianDataNode Data  )
{
  typedef MeanType<RectifiedGaussian,T,MeanId,MeanEnabler> p0_t;
  typedef detail::Factor<Gaussian,T,p0_t,GammaConstType,GaussianDataType> Factor_t;
  
  boost::shared_ptr<GammaConstType > Precision(new GammaConstType(precision));
  boost::shared_ptr<Factor_t > GaussianF(new Factor_t (Mean, Precision.get(), Data));
	
  m_Factors.push_back(GaussianF);
  m_Nodes.push_back(Precision);
}
template<class T>   
template<template<template<class> class,class,class,class> class MeanType,
	 class MeanId, class MeanEnabler>
void 
ICR::EnsembleLearning::Builder<T>::join(MeanType<RectifiedGaussian,T,MeanId,MeanEnabler>* Mean, T& precision, GaussianDataNode Data  )
{
  typedef MeanType<RectifiedGaussian,T,MeanId,MeanEnabler> p0_t;
  typedef detail::Factor<RectifiedGaussian,T,p0_t,GammaConstType,GaussianDataType> Factor_t;
  
  boost::shared_ptr<GammaConstType > Precision(new GammaConstType(precision));
  boost::shared_ptr<Factor_t > RGaussianF(new Factor_t (Mean, Precision.get(), Data));
	
  m_Factors.push_back(RGaussianF);
  m_Nodes.push_back(Precision);
}
  
  
template<class T> 
template<template<template<class> class,class,class,class> class PrecType,
	 class PrecId, class PrecEnabler>
void 
ICR::EnsembleLearning::Builder<T>::join(T& mean, PrecType<Gamma,T,PrecId,PrecEnabler>* Precision, GaussianDataNode Data  )
{
  typedef PrecType<Gamma,T,PrecId,PrecEnabler> p1_t;
  typedef detail::Factor<Gaussian,T,GaussianConstType,p1_t,GaussianDataType> Factor_t;

  boost::shared_ptr<GaussianConstType > Mean(new GaussianConstType(mean));
  boost::shared_ptr<Factor_t > GaussianF(new Factor_t(Mean.get(), Precision, Data));
	
  m_Factors.push_back(GaussianF);
  m_Nodes.push_back(Mean);
}
      
					
					template<class T>   
					template<template<template<class> class,class,class,class> class MeanType,
					class MeanId, class MeanEnabler,
					template<template<class> class,class,class,class> class PrecType,
					class PrecId, class PrecEnabler>
					void
					ICR::EnsembleLearning::Builder<T>::join(MeanType<Gaussian,T,MeanId,MeanEnabler>* Mean,PrecType<Gamma,T,PrecId,PrecEnabler>* Precision, GaussianDataNode& Data  )
{
  typedef MeanType<Gaussian,T,MeanId,MeanEnabler> p0_t;
  typedef PrecType<Gamma,T,PrecId,PrecEnabler> p1_t;
  typedef detail::Factor<Gaussian,T,p0_t,p1_t,GaussianDataType> Factor_t;

  boost::shared_ptr<Factor_t> GaussianF(new Factor_t(Mean,Precision,Data));
  m_Factors.push_back(GaussianF);
}

template<class T>   
      template<template<template<class> class,class,class,class> class MeanType,
	       class MeanId, class MeanEnabler,
	       template<template<class> class,class,class,class> class PrecType,
	       class PrecId, class PrecEnabler>
void 
ICR::EnsembleLearning::Builder<T>::join(MeanType<RectifiedGaussian,T,MeanId,MeanEnabler>* Mean,PrecType<Gamma,T,PrecId,PrecEnabler>* Precision, GaussianDataNode& Data  )
{
  typedef MeanType<RectifiedGaussian,T,MeanId,MeanEnabler> p0_t;
  typedef PrecType<Gamma,T,PrecId,PrecEnabler> p1_t;
  typedef detail::Factor<RectifiedGaussian,T,p0_t,p1_t,GaussianDataType> Factor_t;

  boost::shared_ptr<Factor_t> RGaussianF(new Factor_t(Mean,Precision,Data));
  m_Factors.push_back(RGaussianF);
}


template<class T>
      template<template<template<class> class,class,class,class> class MeanType,
	       class MeanId, class MeanEnabler,
	       template<template<class> class,class,class,class> class PrecType,
	       class PrecId, class PrecEnabler>  
void 
ICR::EnsembleLearning::Builder<T>::join(MeanType<Gaussian,T,MeanId,MeanEnabler>* Mean,
					PrecType<Gamma,T,PrecId,PrecEnabler>* Precision, 
					const T data )
{
  typedef MeanType<Gaussian,T,MeanId,MeanEnabler> p0_t;
  typedef PrecType<Gamma,T,PrecId,PrecEnabler> p1_t;
	
  typedef detail::Factor<Gaussian,T,p0_t,p1_t,GaussianDataType> Factor_t;

  boost::shared_ptr<GaussianDataType > Data(new GaussianDataType(data));
  ++m_data_nodes;
  boost::shared_ptr<Factor_t> GaussianF(new Factor_t(Mean,Precision,Data.get()));
  m_Factors.push_back(GaussianF);
  m_Nodes.push_back(Data);
}

template<class T>
      template<template<template<class> class,class,class,class> class MeanType,
	       class MeanId, class MeanEnabler,
	       template<template<class> class,class,class,class> class PrecType,
	       class PrecId, class PrecEnabler>  
void 
ICR::EnsembleLearning::Builder<T>::join(MeanType<RectifiedGaussian,T,MeanId,MeanEnabler>* Mean,
					PrecType<Gamma,T,PrecId,PrecEnabler>* Precision, 
					const T data )
{
  typedef MeanType<RectifiedGaussian,T,MeanId,MeanEnabler> p0_t;
  typedef PrecType<Gamma,T,PrecId,PrecEnabler> p1_t;
	
  typedef detail::Factor<RectifiedGaussian,T,p0_t,p1_t,GaussianDataType> Factor_t;

  boost::shared_ptr<GaussianDataType > Data(new GaussianDataType(data));
  ++m_data_nodes;
  boost::shared_ptr<Factor_t> RGaussianF(new Factor_t(Mean,Precision,Data.get()));
  m_Factors.push_back(RGaussianF);
  m_Nodes.push_back(Data);
}

//Mixture

template<class T>
template<class vMean_t, class vPrec_t>
typename ICR::EnsembleLearning::Builder<T>::GaussianNode
ICR::EnsembleLearning::Builder<T>::gaussian_mixture(
						    vMean_t& vMean, 
						    vPrec_t& vPrecision,
						    WeightsNode Weights)
{
  
  typedef HiddenNode<Dirichlet, T > prior_t;
  typedef ICR::EnsembleLearning::NoSecondParent blank;
  typedef CatagoryType weights_t;
  typedef detail::Factor<Discrete,T,prior_t,blank,weights_t> Factor_t;

  const size_t number = Weights->size();
	
  boost::shared_ptr<CatagoryType>         Catagory(new CatagoryType(number));
  boost::shared_ptr<Factor_t >      CatagoryF(new Factor_t(Weights, Catagory.get()));
  m_Nodes.push_back(Catagory);
  m_Factors.push_back(CatagoryF);
	

  boost::shared_ptr<GaussianType> Child(new GaussianType());

  m_Nodes.push_back(Child);
	
  typedef detail::Mixture<Gaussian<T>,
    T,
    typename vMean_t::type, 
    typename vPrec_t::type,
    GaussianType >     GaussianMixtureFactor;
  
  boost::shared_ptr<GaussianMixtureFactor> MixtureF(new GaussianMixtureFactor(vMean, vPrecision, Catagory.get() , Child.get()));
	
   m_Factors.push_back(MixtureF);
  return Child.get();
}
   
template<class T>	
template<class vMean_t, class vPrec_t>
void
ICR::EnsembleLearning::Builder<T>::join( vMean_t& vMean, vPrec_t& vPrecision, WeightsNode Weights,const T data )
{
  typedef HiddenNode<Dirichlet, T > prior_t;
  typedef ICR::EnsembleLearning::NoSecondParent blank;
  typedef DirichletType p0_t;
  typedef CatagoryType weights_t;
  typedef detail::Factor<Discrete,T,prior_t,blank,weights_t> Factor_t;

  boost::shared_ptr<GaussianDataType > Data(new GaussianDataType(data));
  ++m_data_nodes;

  const size_t number = Weights->size();
	
  boost::shared_ptr<CatagoryType>         Catagory(new CatagoryType(number));
  boost::shared_ptr< Factor_t>      CatagoryF(new Factor_t(Weights, Catagory.get()));
  m_Nodes.push_back(Catagory);
  m_Factors.push_back(CatagoryF);
  //make the vector of precision nodes and CatagoryNodes
	
  m_Nodes.push_back(Data);
	
  typedef detail::Mixture<Gaussian<T>, T,vMean_t,  vPrec_t, GaussianDataType >     GaussianMixtureFactor;
  
   boost::shared_ptr<GaussianMixtureFactor> MixtureF(new GaussianMixtureFactor(vMean, vPrecision, Catagory.get() , Data.get()));
	
   m_Factors.push_back(MixtureF);
}
	

// template<class T>   
// void
// ICR::EnsembleLearning::Builder<T>::join( std::vector<Variable>& vMean, GammaNode& Precision, WeightsNode Weights,const T data )
// {
//   BOOST_AUTO( vprec, build.mixture_vector<Gamma>   (1.00,0.001) ); 
// }


// template<class T>
// template<class MeanIterator, class PrecIterator>
// typename ICR::EnsembleLearning::Builder<T>::GaussianNode
// ICR::EnsembleLearning::Builder<T>::gaussian_mixture(const MeanIterator& MeanBegin,
// 						    const PrecIterator& PrecBegin,
// 						    WeightsNode Weights)
// {
//   const size_t number = Weights->size();
//   //put in vector form
//   std::vector<Variable> vMean(number);
//   std::vector<Variable> vPrec(number);
	
//   std::copy(MeanBegin,MeanBegin+number, vMean.begin());
//   std::copy(PrecBegin,PrecBegin+number, vPrec.begin());
  
//   return gaussian_mixture(vMean, vPrec, Weights);
// }








// template<class T>
// template<class MeanIterator, class PrecIterator>	
// typename ICR::EnsembleLearning::Builder<T>::RectifiedGaussianNode
// ICR::EnsembleLearning::Builder<T>::rectified_gaussian_mixture(const MeanIterator& MeanBegin,
// 							      const PrecIterator& PrecBegin,
// 							      WeightsNode Weights)
// {
  
//   const size_t number = Weights->size();
//   //put in vector form
//   std::vector<Variable> vMean(number);
//   std::vector<Variable> vPrec(number);
	
//   std::copy(MeanBegin,MeanBegin+number, vMean.begin());
//   std::copy(PrecBegin,PrecBegin+number, vPrec.begin());
  
//   return rectified_gaussian_mixture(vMean, vPrec, Weights);
// }



// template<class T>
// template<class MeanIterator, class PrecIterator>
// void
// ICR::EnsembleLearning::Builder<T>::join( MeanIterator MeanBegin, 
// 					 PrecIterator PrecBegin, 
// 					 WeightsNode Weights,
// 					 const T data )
// {
//   const size_t number = Weights->size();

//   //put in vector form
//   std::vector<Variable> vMean(number);
//   std::vector<Variable> vPrec(number);
	
//   std::copy(MeanBegin,MeanBegin+number, vMean.begin());
//   std::copy(PrecBegin,PrecBegin+number, vPrec.begin());
  
//   join(vMean, vPrec, Weights,data);
// }

#endif //BUILDER_HPP guard
