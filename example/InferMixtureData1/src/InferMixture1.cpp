
#include "EnsembleLearning.hpp"
#include <vector>
#include <iostream>
using namespace ICR::EnsembleLearning;

int
main  (int ac, char **av)
{
  //Create the data
  rng* random = Random::Instance(10); //a random number generator with seed 10.
  const size_t data_points = 500; //make 40 data points.
  std::vector<double> data(data_points);
  for(size_t i = 0; i<data_points; ++i)
    { 
      //The random number generator accepts the standard deviation as an argument
      //  We work with the precision, which is 1.0/variance.
      double catagory = random->uniform(); //choose a catagory at random
      if (catagory<0.4)       //40% from here
	data[i] = random->gaussian(std::sqrt(1.0/10.0), 3);  //mean = 3, precision = 10;
      else if (catagory<0.7)  //30% from here
	data[i] = random->gaussian(std::sqrt(1.0/6.0), 6);  //mean = 6, precision = 6;
      else                    //30% from here
	data[i] = random->gaussian(std::sqrt(1.0/2.0), 3);  //mean = 3, precision = 2;
      
    }
 
  //Build the model with the Builder (using double precision),
  //  output the cost at each iteration to the file "cost.dat"
  Builder<double> build("cost.dat");  

  //Some convenient aliases
  typedef Builder<double>::GaussianNode GaussianNode ;   
  typedef Builder<double>::GammaNode GammaNode; 
  typedef Builder<double>::WeightsNode WeightsNode;

  size_t components = 5; //model a mixture with 5 components.
  
  
  //Node to hold the learnt weight of each component of the mixture
  WeightsNode weights = build.weights(components);
  
 
  //make 5 means and precisions.
  std::vector<GaussianNode> vmean(components);
  std::vector<GammaNode>    vprec(components);
  
  for(size_t i=0;i<components;++i){
    vmean[i] = build.gaussian(0.0,0.001);
    vprec[i] = build.gamma(1.0, 0.01);
  }

 
  //Model The data as Gaussian distributed with the mean an precision determined from the above nodes.
  //  Each data point is modelled indepenantly.
  for(size_t i=0; i<data.size(); ++i)
    {
      build.join(vmean.begin(),vprec.begin(), weights, data[i]);  
    }
  
  //The model is complete, now need to do the inference.
  //  Iterate until convergance of the evidence bound (per data point)
  //  to within 0.01% or  at most 100 iterations.
  build.run(1e-6,500);
 
  //output the inferred weights, means and precisions
  std::cout<<"\nInferred results:\n"<<std::endl;

  for(size_t i=0;i<components;++i){
    std::cout<<"weight["<<i<<"]  = "<<Mean(weights,i)<<"\t +- "<<StandardDeviation(weights,i)<<std::endl;
    std::cout<<"mean["<<i<<"]      = "<<Mean(vmean[i])<<"\t +- "<<StandardDeviation(vmean[i])<<std::endl;
    std::cout<<"precision["<<i<<"] = "<<Mean(vprec[i])<<"\t +- "<<StandardDeviation(vprec[i])<<std::endl<<std::endl;
  }

}
