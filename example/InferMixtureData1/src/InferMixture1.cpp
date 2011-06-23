
#include "EnsembleLearning.hpp"
#include <vector>
#include <iostream>
#include <boost/typeof/typeof.hpp>
#include "boost/mpl/deref.hpp"
#include <iterator>
using namespace ICR::EnsembleLearning;

//Define the number of components in the mixture model
#define ENSEMBLE_LEARNING_COMPONENTS 5



int
main  (int ac, char **av)
{
  //Create the data
  rng* random = Random::Instance(10); //a random number generator with seed 10.
  const size_t data_points = 10000; //make some data points.
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
  typedef Builder<double>::Variable Variable;
  typedef Builder<double>::WeightsNode WeightsNode;

  
  //Node to hold the learnt weight of each component of the mixture
  WeightsNode weights = build.weights(ENSEMBLE_LEARNING_COMPONENTS);
  
  //vectors that hold the means and precisions
  BOOST_AUTO( vmean, build.mixture_vector<Gaussian>(0.00,0.001) ); 
  BOOST_AUTO( vprec, build.mixture_vector<Gamma>   (1.00,0.01) ); 
  //if using new c++0x can use the auto keyword instead
  // auto vmean = build.mixture_vector<Gaussian>(0.0,0.001 );
  // auto vprec = build.mixture_vector<Gamma>(1.0,0.001 );

  //Model The data as Gaussian distributed with the mean an precision determined from the above nodes.
  //  Each data point is modelled indepenantly.
  for(size_t i=0; i<data.size(); ++i)
    {
      build.join(vmean,vprec, weights, data[i]);  
    }
  
  //The model is complete, now need to do the inference.
  //  Iterate until convergance of the evidence bound (per data point)
  //  to within 0.01% or  at most 100 iterations.
   build.run(1e-6,100);
 
  //output the inferred weights, means and precisions
  std::cout<<"\nInferred results:\n"<<std::endl;
  
  //first collect the learned nodes into std::vectors (which are more convenient to handle).
  std::vector<Variable> means = Mixture2Vector(vmean);
  std::vector<Variable> precs = Mixture2Vector(vprec);
  
  //Then print.
  for(size_t i=0;i<ENSEMBLE_LEARNING_COMPONENTS;++i){
    std::cout<<"weight["<<i<<"]  = "<<Mean(weights,i)<<"\t +- "<<StandardDeviation(weights,i)<<std::endl;
    std::cout<<"mean["<<i<<"]      = "<<Mean(means[i])<<"\t +- "<<StandardDeviation(means[i])<<std::endl;
    std::cout<<"precision["<<i<<"] = "<<Mean(precs[i])<<"\t +- "<<StandardDeviation(precs[i])<<std::endl<<std::endl;
  }

}
