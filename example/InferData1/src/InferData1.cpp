
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
  std::vector<double> data(data_points);
  for(size_t i = 0; i<data_points; ++i)
    { 
      //The random number generator accepts the standard deviation as an argument
      //  We work with the precision, which is 1.0/variance.
      data[i] = random->gaussian(std::sqrt(1.0/10.0), 3);  //mean = 3, precision = 10;
    }
 
  //Build the model with the Builder (using double precision)
  Builder<double> build;  

  //Some convenient aliases
  typedef Builder<double>::GaussianNode GaussianNode ;   
  typedef Builder<double>::GammaNode GammaNode; 
 
  //Nodes to hold the learnt mean and precision
  GaussianNode mean = build.gaussian(0.00,0.001); //Broad priors on the Gaussian.
  GammaNode    prec = build.gamma(1.0,0.01);     //Broad priors on the Gamma model.
 
  //Model The data as Gaussian distributed with the mean an precision determined from the above nodes.
  //  Each data point is modelled indepenantly.
  for(size_t i=0; i<data.size(); ++i)
    {
      build.join(mean,prec,data[i]);  
    }
  
  //The model is complete, now need to do the inference.
  //  Iterate until convergance of the evidence bound (per data point)
  //  to within 0.01% or  at most 100 iterations.
  build.run(0.01,10);
 
  //output the inferred mean and precision
  std::cout<<"mean      = "<<Mean(mean)<<"\t +- "<<StandardDeviation(mean)<<std::endl;
  std::cout<<"precision = "<<Mean(prec)<<"\t +- "<<StandardDeviation(prec)<<std::endl;

}
