

#  define ENSEMBLE_LEARNING_PLACEHOLDERS 2
#include "EnsembleLearning.hpp"
#include <vector>
#include <iostream>
#include <boost/typeof/typeof.hpp>
using namespace ICR::EnsembleLearning;
int
main  (int ac, char **av)
{
 
  //Create the data
  rng* random = Random::Instance(10); //a random number generator with seed 10.
  const size_t data_points = 1000; //make 20 data points.
  std::vector<double> data(data_points);  //the true source data
  for(size_t i = 0; i<data_points; ++i)
    { 
      data[i] = random->gaussian(std::sqrt(1.0/10.0), 3);  //mean = 3, precision = 10;
    }
 
  //Build the model with the Builder (using double precision)
  Builder<double> build;  

  //Some convenient aliases
  // typedef Builder<double>::GaussianNode       GaussianNode;   
  // typedef Builder<double>::RectifiedGaussianNode       RectifiedGaussianNode;   
   typedef Builder<double>::GammaNode          GammaNode; 
  typedef Builder<double>::GaussianResultNode ResultNode;
  typedef Builder<double>::Variable Variable;
  typedef Builder<double>::WeightsNode        WeightsNode;
 

  
  const size_t components = ENSEMBLE_LEARNING_PLACEHOLDERS; //model a mixture with 2 components..

  typedef  CalculationVector<Gaussian,double> CV_t;
  //CV_t vmean =  build.calculation_vector<RectifiedGaussian>(0.0,0.001);
  BOOST_AUTO(vmean, build.calculation_vector<Gaussian>(0.0,0.001));
  // std::vector<RectifiedGaussianNode> vmean(components);
  
  // for(size_t i=0;i<components;++i){
  //   vmean[i] = build.rectified_gaussian(0.0,0.001);
  // }
  
  
  GammaNode precision = build.gamma(1.0,0.01); //Broad priors on the gamma, should be approx 0.32 from propagation of errors calculation.

  
  
  //Build the expression to infer.
  //ExpressionFactory<double> factory;
  //some convenient aliases
  // typedef ExpressionFactory<double>::placeholder_t placeholder_t;   
  // typedef ExpressionFactory<double>::expression_t  expression_t;   

  PlaceholderFactory::make_c<1> g1;   //placeholder representing the 1st gaussian.
  PlaceholderFactory::make_c<2> g2;   //placeholder representing the 2nd gaussian.
  //typedef typename Expression::add<g1,g2>::type expr;
  
  //expression_t expr = factory.Add(g1, g2);  //scale*data+offset;

  
  //first give a particular context to the expression 
  // (i.e. give a value to the placeholder)
  
  Context<double, CV_t> context(vmean);
  // context.Assign(g1,vmean.get<0>());
  // context.Assign(g2,vmean.get<1>());

  //Model the data (each data point is modelled indepenantly)
  for(size_t i=0; i<data.size(); ++i)
    {
      ResultNode result = build.calc_gaussian(g1+g2,context);  
      build.join(result,precision, data[i]);  
    }
  
  //The model is complete, now need to do the inference.
  //  Iterate until convergance of the evidence bound (per data point)
  //  to within 0.01% or  at most 100 iterations.
  build.run(0.0001,1000);
 
  //output the inferred mean and precision

  std::vector<Variable> means = to_std_vector(vmean);
  // std::cout<<"mean0 = "<<means[0]<<std::endl;
  // std::cout<<"mean1 = "<<means[1]<<std::endl;

  //output the inferred weights, means and precisions
  std::cout<<"\nInferred results:\n"<<std::endl;

  for(size_t i=0;i<2;++i){
    std::cout<<"mean["<<i<<"]      = "<<Mean(means[i])<<"\t +- "<<StandardDeviation(means[i])<<std::endl;
  }
}
