
#include "EnsembleLearning.hpp"
#include <vector>
#include <iostream>
using namespace ICR::EnsembleLearning;

int
main  (int ac, char **av)
{
 
  //Create the data
  rng* random = Random::Instance(10); //a random number generator with seed 10.
  const size_t data_points = 20; //make 20 data points.
  std::vector<double> data(data_points);  //the true source data
  for(size_t i = 0; i<data_points; ++i)
    { 
      data[i] = random->gaussian(std::sqrt(1.0/10.0), 3);  //mean = 3, precision = 10;
    }
 
  //Build the model with the Builder (using double precision)
  Builder<double> build;  

  //Some convenient aliases
  typedef Builder<double>::GaussianNode       GaussianNode;   
  typedef Builder<double>::RectifiedGaussianNode       RectifiedGaussianNode;   
  typedef Builder<double>::GammaNode          GammaNode; 
  typedef Builder<double>::GaussianResultNode ResultNode;
  typedef Builder<double>::WeightsNode        WeightsNode;
 

  
  size_t components = 2; //model a mixture with 2 components.
  
  //make 5 means and precisions.
  std::vector<RectifiedGaussianNode> vmean(components);
  
  for(size_t i=0;i<components;++i){
    vmean[i] = build.rectified_gaussian(0.0,0.001);
  }
  
  
  GammaNode precision = build.gamma(1.0,0.01); //Broad priors on the gamma, should be approx 0.32 from propagation of errors calculation.

  
  
  //Build the expression to infer.
  ExpressionFactory<double> factory;
  //some convenient aliases
  typedef ExpressionFactory<double>::placeholder_t placeholder_t;   
  typedef ExpressionFactory<double>::expression_t  expression_t;   

  placeholder_t g1 = factory.placeholder();   //placeholder representing the 1st gaussian.
  placeholder_t g2 = factory.placeholder();   //placeholder representing the 2nd gaussian.
  expression_t expr = factory.Add(g1, g2);  //scale*data+offset;

  
  //first give a particular context to the expression 
  // (i.e. give a value to the placeholder)
  Context<double> context;
  context.Assign(g1, vmean[0]);
  context.Assign(g2, vmean[1]);

  //Model the data (each data point is modelled indepenantly)
  for(size_t i=0; i<data.size(); ++i)
    {
      ResultNode result = build.calc_gaussian(expr,context);  
      build.join(result,precision, data[i]);  
    }
  
  //The model is complete, now need to do the inference.
  //  Iterate until convergance of the evidence bound (per data point)
  //  to within 0.01% or  at most 100 iterations.
  build.run(0.0001,100);
 
  //output the inferred mean and precision

  
  //output the inferred weights, means and precisions
  std::cout<<"\nInferred results:\n"<<std::endl;

  for(size_t i=0;i<components;++i){
    std::cout<<"mean["<<i<<"]      = "<<Mean(vmean[i])<<"\t +- "<<StandardDeviation(vmean[i])<<std::endl;
  }
}
