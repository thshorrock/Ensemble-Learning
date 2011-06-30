#pragma once
#include <typeinfo>
#include "EnsembleLearning.hpp"

#include <boost/bind.hpp>
#include <boost/typeof/typeof.hpp>
#include "boost/mpl/deref.hpp"
#include <iterator>

#include <deque>
#include <fstream>

#include <boost/mpl/transform.hpp>
#include <boost/mpl/accumulate.hpp>
// template<class data_t>
// struct Model
// {
//   virtual
//   ICR::EnsembleLearning::Builder<data_t>& 
//   get_builder() = 0;
  
//   virtual
//   void 
//   get_normalised_means(matrix<data_t>& A, matrix<data_t>& S) = 0;
  
// };

using namespace ICR::EnsembleLearning;

class size_exception{};


template<class data_t, int assumed_sources, int mixing_components>
class BuildModel //: public Model<data_t>
{
  typedef typename ICR::EnsembleLearning::Builder<data_t>::GaussianDataNode GaussianData;
  typedef typename ICR::EnsembleLearning::Builder<data_t>::GammaDataNode GammaData;
  typedef typename ICR::EnsembleLearning::Builder<data_t>::GaussianNode GaussianNode;
  typedef typename ICR::EnsembleLearning::Builder<data_t>::RectifiedGaussianNode RectifiedGaussianNode;
  typedef typename ICR::EnsembleLearning::Builder<data_t>::WeightsNode WeightsNode;
  typedef typename ICR::EnsembleLearning::Builder<data_t>::GammaNode    GammaNode;
  typedef typename ICR::EnsembleLearning::Builder<data_t>::GaussianDataNode    Datum;
  typedef typename ICR::EnsembleLearning::Builder<data_t>::Variable    Variable;
  typedef typename ICR::EnsembleLearning::Builder<data_t>::GaussianResultNode    ResultNode;
  
  typedef  HiddenNode<Gaussian,data_t,typename ICR::EnsembleLearning::detail::TypeList::incr_id<ICR::EnsembleLearning::detail::TypeList::zeros >::type > Mean_t;
  
  //A set of T calculation vectors
  typedef CalculationVector<Gaussian,data_t,1,assumed_sources> Ag_t;
  typedef CalculationVector<RectifiedGaussian,data_t,1,assumed_sources> Arg_t;
  typedef CalculationVector<Gaussian,data_t,assumed_sources+1,assumed_sources> Sg_t;
  typedef CalculationVector<RectifiedGaussian,data_t,assumed_sources+1,assumed_sources> Srg_t;
  typedef CalculationVector<Gaussian,data_t,2*assumed_sources+1,2> noise_t;

  template<int i>
  struct Int2Type
  {
    enum {value = i};
  };
  
public:
  BuildModel(matrix<double>& data,
	     const bool positive_sources,
	     bool positive_mixing,
	     bool noise_offset,
	     double GaussianPrecision,
	     double GammaPrecision
	     )
    : m_Build(), 
      // m_Factory(),
      m_A(data.size1(),assumed_sources),
      m_S(assumed_sources, data.size2()),
      m_noiseMean(data.size1()),
      //m_noisePrecision(data.size1()), //N
      m_noisePrecision(data.size2()),  //T
      m_positive_sources(positive_sources),
      m_positive_mixing(positive_mixing),
      m_GaussianPrecision(GaussianPrecision),
      m_GammaVar(1.0/GammaPrecision)
  {
 

    const size_t Components = mixing_components;
    const size_t M = assumed_sources; //assumed sources
    const size_t N = data.size1(); //DataSources
    const size_t T  =data.size2(); //DataLength
    
    std::cout<<"C = "<<Components<<std::endl
  	     <<"M = "<<M<<std::endl
  	     <<"N = "<<N<<std::endl
  	     <<"T = "<<T<<std::endl;

    
     
    //m_Build the source approximation

    //build A weights component for each source
    vector<WeightsNode> Weights(M);
    build_vector(Weights, Components);
    
    //std::vector< MixtureVector<Gaussian,data_t,Components> > ShypMean(M);
    std::vector< ObservedMixtureVector<Gaussian,data_t,Components> > ShypMean(M);
    std::vector< MixtureVector<Gamma,data_t,Components>    > ShypPrec(M);
    
    
    std::vector< Sg_t > cv_Sg(T);
    std::vector< Srg_t > cv_Srg(T);
    
    
    std::vector<data_t> vSHypMean(Components,0.0);
    ObservedMixtureVector<Gaussian,data_t,Components> SHypMean_node 
      = m_Build.template observed_mixture_vector<Gaussian>(vSHypMean);
    

    //m_Build hyper mean and hyper precision for each source mixture component.
    for(size_t m=0;m<M;++m)
      { 
    	ShypMean[m] =  SHypMean_node;//m_Build.template mixture_vector<Gaussian>(0.00,m_GaussianPrecision); 
    	ShypPrec[m] =  m_Build.template mixture_vector<Gamma> (data_t(1.00),data_t(m_GammaVar)) ;
      }
    
    std::cout<<"built hyperparams - nodes = "<< m_Build.number_of_nodes() <<std::endl;
    for(size_t t=0;t<T;++t){
      if (positive_sources) 
    	{
    	  cv_Srg[t]=m_Build.template calculation_mixture<RectifiedGaussian,M+1,M>
    	    (ShypMean,
    	     ShypPrec,
    	     Weights );

    	  std::vector<Variable> S_over_m = to_std_vector(cv_Srg[t]);
    	  for(size_t m=0;m<M;++m){
    	    m_S(m,t) = S_over_m[m];
    	  }
    	}
 
      else
    	{
    	  cv_Sg[t]=m_Build.template calculation_mixture<Gaussian,M+1,M>
    	    (ShypMean, 
    	     ShypPrec,
    	     Weights );
    	  std::vector<Variable> S_over_m = to_std_vector(cv_Sg[t]);
      
    	  for(size_t m=0;m<M;++m){
    	    m_S(m,t) = S_over_m[m];
    	  }
    	}
    }



    //m_Build approximation to the mixing matrix
     
    std::vector<Mean_t*> AMean(M);
    std::vector<GammaNode> APrecision(M);
    build_vector<Mean_t>(AMean);
    build_vector(APrecision);
    
    
    //A set of N calculation vectors

    std::vector< Ag_t  > cv_Ag(N);
    std::vector< Arg_t  > cv_Arg(N);
    // std::vector<noise_t> m_noiseMean(N);

    for(size_t n=0;n<N;++n){
      if (positive_mixing) 
  	{
  	  cv_Arg[n] = m_Build.template calculation<RectifiedGaussian,1,M>(AMean, APrecision);
  	  std::vector<Variable> A_over_m = to_std_vector(cv_Arg[n]);
  	  for(size_t m=0;m<M;++m)
  	    m_A(n,m) = A_over_m[m];
  	}
      else
  	{
  	  cv_Ag[n] = m_Build.template calculation<Gaussian,1,M>(AMean, APrecision);
  	  std::vector<Variable> A_over_m = to_std_vector(cv_Ag[n]);
  	  for(size_t m=0;m<M;++m)
  	    m_A(n,m) = A_over_m[m];
  	} 

      // m_noiseMean[n] = m_Build.template calculation_vector<Gaussian,2*M+1,2>(0.0, m_GaussianPrecision);
    }
    
    
    
    //    vector<GammaNode> m_noisePreceision(N);
    build_vector(m_noisePrecision);
    

  typedef PlaceholderFactory::make_vector<1,M> AP;
  typedef PlaceholderFactory::make_vector<M+1,2*M>  SP;
  typedef PlaceholderFactory::make_c<2*M+1> NP;
  
  std::cout<<"size SP = "<<boost::mpl::size<SP>::value<<std::endl;
  std::cout<<"size AP = "<<boost::mpl::size<AP>::value<<std::endl;

  typedef typename boost::mpl::transform<AP,SP,Times>::type the_product;
  typedef typename boost::mpl::accumulate<the_product,zero_t,Plus>::type the_inner_product;
  the_inner_product Expr;
  typedef typename boost::mpl::accumulate<NP,the_inner_product,Plus>::type the_inner_product_with_offset;
  //the_inner_product_with_offset ExprWithNoise;

  matrix<ResultNode> AtimesSplusN(N,T);
  for(size_t n=0;n<N;++n){ 

    // 	for(size_t m=0;m<M;++m){
	  
    // 	  context.push_back(m_A(n,m))
    // 	}
    for(size_t t=0;t<T;++t){ 
      
      /**  Now need to replace the placeholders with given nodes 
       *   I.e. set up a context to the expression
       */
      if (positive_mixing) 
  	{
  	  if (positive_sources) 
  	    {
  	      typedef fuse<Arg_t,Srg_t> AS_t;
	      AS_t AS(cv_Arg[n],cv_Srg[t] );
  	      if (noise_offset) 
  	      	{
		  // typedef fuse<AS_t, noise_t> ASn_t;
  	      	  // ASn_t ASn(AS,m_noiseMean[n]);
		  // Context<data_t,ASn_t> context(ASn);
		  // AtimesSplusN(n,t) = m_Build.calc_gaussian(ExprWithNoise(),context);  
		}
	      else
  	      	{
  	      	   Context<data_t,AS_t> context(AS);
  	      	   AtimesSplusN(n,t) = m_Build.calc_gaussian(Expr,context);  
  	      	}
  	    }
      	  else
      	    {
      	      typedef fuse<Arg_t,Sg_t> AS_t;
      	      AS_t AS(cv_Arg[n],cv_Sg[t] );
      	      if (noise_offset) 
      		{
      		  // typedef fuse<typename AS_t::type, Mean_t> ASn_t;
      		  // ASn_t ASn(AS,m_noiseMean[n]);
      		  // Context<data_t,ASn_t> context(ASn);
      		  // AtimesSplusN(n,t) = m_Build.calc_gaussian(ExprWithNoise(),context);  
      		}
      	      else
      		{
      		  Context<data_t,AS_t> context(AS);
      		  AtimesSplusN(n,t) = m_Build.calc_gaussian(Expr,context);  
      		}
      	    }
      	}
      else
      	{
      	  if (positive_sources) 
      	    {
      	      typedef fuse<Ag_t,Srg_t> AS_t;
      	      AS_t AS(cv_Ag[n],cv_Srg[t] );
      	      if (noise_offset) 
      		{
      		  // typedef fuse<typename AS_t::type, Mean_t> ASn_t;
      		  // ASn_t ASn(AS,m_noiseMean[n]);
      		  // Context<data_t,ASn_t> context(ASn);
      		  // AtimesSplusN(n,t) = m_Build.calc_gaussian(ExprWithNoise(),context); 
		  
      		}
      	      else
      		{
      		  Context<data_t,AS_t> context(AS);
      		  AtimesSplusN(n,t) = m_Build.calc_gaussian(Expr,context);  
      		}
      	    }
      	  else
      	    {
      	      typedef fuse<Ag_t,Sg_t> AS_t;
      	      AS_t AS(cv_Ag[n],cv_Sg[t] );
      	      if (noise_offset) 
      		{
      		  // typedef fuse<typename AS_t::type, Mean_t> ASn_t;
      		  // ASn_t ASn(AS,m_noiseMean[n]);
      		  // Context<data_t,ASn_t> context(ASn);
      		  // AtimesSplusN(n,t) = m_Build.calc_gaussian(ExprWithNoise(),context);  
      		}
      	      else
      		{
      		  Context<data_t,AS_t> context(AS);
		  // Context<data_t,AS_t> c2(3);
      		  AtimesSplusN(n,t) = m_Build.calc_gaussian(Expr,context);  
      		}
      	    }
       	}
    }
  }
  m_AtimesSplusN = AtimesSplusN;


  std::cout<<"built context - nodes = "<< m_Build.number_of_nodes() <<std::endl;
     
  //join to the data
  for(size_t n=0;n<N;++n){
    for(size_t t=0;t<T;++t){
      // std::cout<<"n = "<<n<<" t = "<<t<<std::endl;
      // std::cout<<"RN"<<std::endl;
      // std::cout<<AtimesSplusN(n,t)<<std::endl;
      // std::cout<<"NP"<<std::endl;;
      // std::cout<<m_noisePrecision[n]<<std::endl;
      // std::cout<<"data"<<std::endl;
      // std::cout<<data(n,t)<<std::endl;

      // std::cout<<"join"<<std::endl;

      //m_Build.join(AtimesSplusN(n,t),m_noisePrecision[n], data(n,t));
      m_Build.join(AtimesSplusN(n,t),m_noisePrecision[t], data(n,t));
    }
  }

  std::cout<<"built data - nodes = "<< m_Build.number_of_nodes() <<std::endl;
     


  }
  ICR::EnsembleLearning::Builder<data_t>& get_builder() {return m_Build;}
  matrix<Variable>& get_sources() {return m_S;}
  matrix<Variable>& get_mixing_matrix(){return m_A;}
  vector<Variable>& get_noiseMean(){return m_noiseMean;}
  //vector<Variable>& get_noisePrecision(){return m_noisePrecision;}
  
  void set_means(matrix<data_t>&Means)
  {
    //First check that number of rows and columns less than m_S.

    std::cout<<"setting means"<<std::endl;
    if (
    	(Means.size1()>m_S.size1())
    	||
    	(Means.size2()>m_S.size2())
    	)
      {
    	std::cout<<"Means size = ("<<Means.size1()<<","<<Means.size2()<<")"<<std::endl;
    	std::cout<<"Model size = ("<<m_S.size1()<<","<<m_S.size2()<<")"<<std::endl;
    	throw size_exception();
      }
    
    //for each row 
    for(size_t m=0;m<Means.size1();++m){
      //set m^th column of A matrix to be one.
      for(size_t n=0;n<m_A.size1();++n){
    	if (m_positive_mixing) 
    	  ICR::EnsembleLearning::SetMean(static_cast<RectifiedGaussianNode>(m_A(m,n)), 1.0);
    	else
    	  ICR::EnsembleLearning::SetMean(static_cast<GaussianNode>(m_A(m,n)), 1.0);
	  
      }
      //for each column
      for(size_t t=0;t<Means.size2(); ++t) {
    	//set the mean
    	if (m_positive_sources) 
    	  ICR::EnsembleLearning::SetMean(static_cast<RectifiedGaussianNode>(m_S(m,t)), Means(m,t));
    	else
    	  ICR::EnsembleLearning::SetMean(static_cast<GaussianNode>(m_S(m,t)), Means(m,t));
	
    	m_AtimesSplusN(m,t)->GetMoments(); //update
      }
      
    }
  }

  void set_sigmas(matrix<data_t>&SD)
  {
    std::cout<<"setting Sigmas"<<std::endl;

    //First check that number of rows and columns less than m_S.
    if (
    	(SD.size1()>m_S.size1())
    	||
    	(SD.size2()>m_S.size2())
    	)
      {
    	std::cout<<"SD    size = ("<<SD.size1()<<","<<SD.size2()<<")"<<std::endl;
    	std::cout<<"Model size = ("<<m_S.size1()<<","<<m_S.size2()<<")"<<std::endl;
    	throw size_exception();
      }
    
    //for each row 
    for(size_t m=0;m<SD.size1();++m){
      //set m^th column of A matrix to be the average sd of row.
      matrix_row<matrix<double> > row(SD,m);
      double av_sd = std::accumulate(row.begin(), row.end(), 0.0)/row.size();
      for(size_t n=0;n<m_A.size1();++n){
    	if (m_positive_mixing) 
    	  ICR::EnsembleLearning::SetStandardDeviation(static_cast<RectifiedGaussianNode>(m_A(n,m)), av_sd);
    	else
    	  ICR::EnsembleLearning::SetStandardDeviation(static_cast<GaussianNode>(m_A(n,m)), av_sd);
      }
      //for each column
      for(size_t t=0;t<SD.size2(); ++t) {
    	//set the mean
    	if (m_positive_sources) 
    	  ICR::EnsembleLearning::SetStandardDeviation(static_cast<RectifiedGaussianNode>(m_S(m,t)), SD(m,t));
    	else
    	  ICR::EnsembleLearning::SetStandardDeviation(static_cast<GaussianNode>(m_S(m,t)), SD(m,t));


    	m_AtimesSplusN(m,t)->GetMoments(); //update 
      }
    }
  }

  void set_mixing_mean(matrix<data_t>& M)
  {
    std::cout<<"setting Mixing"<<std::endl;

    //First check that number of rows and columns less than m_A.
    if (
    	(M.size1()>m_A.size1())
    	||
    	(M.size2()>m_A.size2())
    	)
      {
    	std::cout<<"Mixing size = ("<<M.size1()<<","<<M.size2()<<")"<<std::endl;
    	std::cout<<"Model  size = ("<<m_A.size1()<<","<<m_A.size2()<<")"<<std::endl;
    	throw size_exception();
      }
    
    //for each row 
    for(size_t m=0;m<M.size1();++m){
      for(size_t n=0;n<M.size2();++n){
    	if (m_positive_mixing) 
    	  ICR::EnsembleLearning::SetMean(static_cast<RectifiedGaussianNode>(m_A(n,m)), M(n,m));
    	else
    	  ICR::EnsembleLearning::SetMean(static_cast<GaussianNode>(m_A(n,m)), M(n,m));
      }
      //m_AtimesSplusN(m,t)->GetMoments(); //update 
    }
  
    for(size_t m=0;m<M.size1();++m){
      for(size_t t=0;t<M.size2(); ++t) {
    	m_AtimesSplusN(m,t)->GetMoments(); //update 
      }
    }
  }

  void get_normalised_means(matrix<data_t>& A, matrix<data_t>& S)
  {
    //std::cout<<"m_A size = "<<m_A.size1()<<", "<<m_A.size2()<<std::endl;

    for(size_t i=0;i<m_A.size1();++i){
      for(size_t j=0; j<m_A.size2();++j){
    	A(i,j)= Mean(m_A(i,j));
      }
    }

    for(size_t i=0;i<m_S.size1();++i){
      for(size_t j=0; j<m_S.size2();++j){
    	S(i,j)= Mean(m_S(i,j));
      }
    }

    //normalise
    for(size_t j=0;j<A.size2();++j){
      matrix_column<matrix<data_t> > col(A, j);
      const vector<data_t> col2 = element_prod(col,col);
      const data_t norm = 1.0/std::sqrt(std::accumulate(col2.begin(), col2.end(), 0.0));
      col*=norm;
      
      matrix_row<matrix<data_t> > row(S, j);
      row/=norm;
    }
  }
  matrix<data_t> get_results() const
  {
    matrix<data_t> r(m_AtimesSplusN.size1(), m_AtimesSplusN.size2());
    for(size_t i=0;i<m_AtimesSplusN.size1();++i){
      for(size_t j=0; j<m_AtimesSplusN.size2();++j){
    	r(i,j)= Mean(m_AtimesSplusN(i,j));
      }
    }
    return r;
  }

  vector<data_t> get_noise_precision() const
  {

    vector<data_t> r(m_noisePrecision.size());
    std::cout<<"noise precision size = "<<r.size()<<std::endl;
    for(size_t i=0;i<m_noisePrecision.size();++i){
      r(i)= Mean(m_noisePrecision(i));
    }
    return r;
  }

private: 
  template<class node>
  void 
  build_vector(vector<node*>& V)
  {

    for(size_t i=0;i<V.size();++i){
      V[i] = m_Build.gaussian(0.0,m_GaussianPrecision);
    }
  }
  template<class node>
  void 
  build_vector(std::vector<node*>& V)
  {

    for(size_t i=0;i<V.size();++i){
      V[i] = m_Build.gaussian(0.0,m_GaussianPrecision);
    }
  }
  
  void 
  build_vector(vector<GammaNode>& V)
  {
    std::generate(V.begin(), V.end(), boost::bind(&ICR::EnsembleLearning::Builder<data_t>::gamma,
    						  boost::ref(m_Build),
    						  1,m_GammaVar ));
  }
  void 
  build_vector(std::vector<GammaNode>& V)
  {
    std::generate(V.begin(), V.end(), boost::bind(&ICR::EnsembleLearning::Builder<data_t>::gamma,
    						  boost::ref(m_Build),
    						  1,m_GammaVar ));
  }
  void 
  build_vector(vector<WeightsNode>& V, size_t components)
  {
    std::generate(V.begin(), V.end(), boost::bind(&ICR::EnsembleLearning::Builder<data_t>::weights,
    						  boost::ref(m_Build)
    						  ));
  }
  
  

ICR::EnsembleLearning::Builder<data_t> m_Build;
//ICR::EnsembleLearning::ExpressionFactory<data_t> m_Factory;
matrix<Variable> m_A;
matrix<Variable> m_S;
matrix<Variable> m_AtimesSplusN; 
std::vector<noise_t> m_noiseMean;
vector<GammaNode>     m_noisePrecision;
//GammaNode     m_noisePrecision;
bool m_positive_mixing, m_positive_sources;
data_t m_GaussianPrecision;
data_t m_GammaVar;

};

  template<class data_t>
  std::ofstream&
  operator<<(std::ofstream& out, const matrix<data_t>& M)
  {
    for(size_t j=0; j<M.size2();++j){
      for(size_t i=0;i<M.size1();++i){
        out<< M(i,j) << "\t";
      }
      out<<"\n";
    }
    return out;
  }

template<class data_t>
std::ofstream&
operator<<(std::ofstream& out, const vector<data_t>& M)
{
  for(size_t i=0;i<M.size();++i){
    out<< M(i) << "\t";
  }
  out<<"\n";
  return out;
}
