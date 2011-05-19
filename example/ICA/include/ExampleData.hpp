#pragma once
#include "ICA/exponential_model/Random.hpp"


#include <numeric>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

using namespace boost::numeric::ublas;

//Make a vector with the shape of a Gaussian.
vector<double>
make_Gaussian(double mean, double precision, size_t samples , double from = 0, double to = 10)
{
  vector<double> Data(samples);
  double inorm = std::sqrt(2*precision/(2.0*M_PI));
  
  double dx = (to-from)/double(samples);
  for(size_t i=0;i<Data.size();++i){
    double x = dx*i;
    double expen = std::exp(precision*mean*x-0.5*precision*x*x+ -.5*precision*mean*mean);
    Data[i]=inorm*expen;      
  }
  return Data;
}

//Make a vector of noise
vector<double>
make_noise(double precision, size_t samples)
{
  vector<double> Data(samples);
      
  ICR::ICA::rng* random = ICR::ICA::Random::Instance() ;// get the random
  for(size_t i=0;i<samples;++i){
    Data[i] = random->gaussian(1.0/std::sqrt(precision));
  }
  return Data;
}
    

class Sources
{
public:
  Sources(size_t number_of_sources, size_t samples_per_source, bool positive = false)
    : m_sources(number_of_sources,samples_per_source)
  {
    ICR::ICA::rng* random = ICR::ICA::Random::Instance() ;// get the random
    for(size_t i=0;i<number_of_sources;++i){
      matrix_row<matrix<double> > row(m_sources,i);
      double mean = random->uniform(0,10);
      double prec = random->gamma(1,1);
      double amplitude;
      if (positive) 
	amplitude = random->gamma();
      else
	amplitude = random->gaussian();
      row = amplitude * make_Gaussian(mean,prec,samples_per_source);
    }
  }

  operator matrix<double>() const {return m_sources;}

private:
  matrix<double> m_sources;
};

class Mixing
{
  struct times{
    double operator()(double a, double b)
    {
      return a*b;
    }
  };
public:
  Mixing(size_t number_of_sources, //underlying
	 size_t number_of_records, //mixed data
	 bool   positive = false
	 ) 
    : m_mixing(number_of_records, number_of_sources)
  {
	
    ICR::ICA::rng* random = ICR::ICA::Random::Instance() ;// get the random

    for(size_t i=0;i<m_mixing.size1();++i){
      for(size_t j=0;j<m_mixing.size2();++j){
	double element;
	if (positive)
	  element = random->gamma();
	else
	  element = random->gaussian();
	  
	m_mixing(i,j) = element;
      }
    }
    //normalise the mixing matrix
    for(size_t j=0;j<m_mixing.size2();++j){
      matrix_column<matrix<double> > col(m_mixing, j);
      vector<double> col2 = element_prod(col,col);
      double norm = 1.0/std::sqrt(std::accumulate(col2.begin(), col2.end(), 0.0));
      col*=norm;
    }
  }

  operator matrix<double>() const {return m_mixing;}
private:
  matrix<double> m_mixing;
};
    
matrix<double>
prod(const Mixing& M, const Sources& S)
{
  const matrix<double>& m = M;
  const matrix<double>& s = S;
  return prod(m,s);
}

matrix<double> 
DataRecords(size_t number_of_records, 
	    size_t number_of_sources,
	    size_t samples_per_record,
	    bool   positive_sources = false,
	    bool   positive_mixing  = false
	    )
{
  Sources S(number_of_sources,samples_per_record, positive_sources);
  Mixing  M(number_of_sources,number_of_records, positive_mixing);
  
  return  prod(M,S);
}


void
AddNoise(matrix<double> M, double precision)
{
  matrix<double> noise(M.size1(), M.size2());
  for(size_t r=0;r<M.size1();++r){
    matrix_row<matrix<double> > row(M,r);
    row+=make_noise(precision, M.size2());
  }
  M+=noise;
}
    
