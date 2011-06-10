#pragma once
#include "EnsembleLearning.hpp"

#include <boost/bind.hpp>

#include <deque>
#include <fstream>

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


class size_exception{};


template<class data_t>
class BuildModel //: public Model<data_t>
{
  typedef typename ICR::EnsembleLearning::Builder<data_t>::GaussianNode GaussianNode;
  typedef typename ICR::EnsembleLearning::Builder<data_t>::RectifiedGaussianNode RectifiedGaussianNode;
  typedef typename ICR::EnsembleLearning::Builder<data_t>::WeightsNode WeightsNode;
  typedef typename ICR::EnsembleLearning::Builder<data_t>::GammaNode    GammaNode;
  typedef typename ICR::EnsembleLearning::Builder<data_t>::GaussianDataNode    Datum;
  typedef typename ICR::EnsembleLearning::Builder<data_t>::Variable    Variable;
  typedef typename ICR::EnsembleLearning::Builder<data_t>::GaussianResultNode    ResultNode;
  
  template<int i>
  struct Int2Type
  {
    enum {value = i};
  };
  
public:
  BuildModel(matrix<double>& data,
	     size_t assumed_sources,
	     size_t mixing_components,
	     const bool positive_sources,
	     bool positive_mixing,
	     bool noise_offset,
	     double GaussianPrecision,
	     double GammaPrecision
	     )
    : m_Build(), 
      m_Factory(),
      m_A(data.size1(),assumed_sources),
      m_S(assumed_sources, data.size2()),
      m_noiseMean(data.size1()),
      m_noisePrecision(data.size1()),
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
    if (positive_sources) 
      build_source_matrix(Int2Type<true>(),M,T, Components);
    else
      build_source_matrix(Int2Type<false>(),M,T, Components);
    
    //m_Build approximation to the mixing matrix
    if (positive_mixing) 
      build_mixing_matrix(Int2Type<true>(),N,M);
    else
      build_mixing_matrix(Int2Type<false>(),N,M);
      

    //m_Build the noise mean approximation
    // m_noiseMean(N);
    //  m_noisePrecision(N);
    for(size_t n=0;n<N;++n){
      m_noiseMean[n] = m_Build.gaussian(0.0,m_GaussianPrecision);
      m_noisePrecision[n] = m_Build.gamma(1,m_GammaVar);
    }
    
    //m_noisePrecision = m_Build.gamma(1.0,1.0);
      
    //Deterministic Node.  Need to make the expression.
    //The inner product plus noise
    /** The inner product is summed over M.
     *  Therefore there are 2M + 1 placeholders required
     *  M for Anm, M for Smt and 1 for Noise
     */
    //Make a set of placeholders for all the elements in the expression.
    vector<ICR::EnsembleLearning::Placeholder<data_t>*> SP(M);
    vector<ICR::EnsembleLearning::Placeholder<data_t>*> AP(M);
    ICR::EnsembleLearning::Placeholder<data_t>* NP = m_Factory.placeholder();
    for(size_t m=0;m<M;++m){
      SP[m] = m_Factory.placeholder();  
      AP[m] = m_Factory.placeholder(); 
    }
    //The expression in terms of these placeholders
    ICR::EnsembleLearning::Expression<data_t>* Expr;
    {
      BOOST_ASSERT(M!=0);
      //First do the multiplication
      std::deque<ICR::EnsembleLearning::Expression<data_t>*> prod(M);
      for(size_t m=0;m<M;++m){
    	prod[m] = m_Factory.Multiply(AP[m], SP[m]);
      }
      //then add them up
      Expr = prod[0];
      prod.pop_front();
      while(prod.size()!=0)
    	{
    	  Expr = m_Factory.Add(Expr,prod[0]);
    	  prod.pop_front();
    	}
      if (noise_offset) 
	//Finnaly add the Noise placeholder
	Expr = m_Factory.Add(Expr, NP);
      
    }
    
    
    matrix<ResultNode> AtimesSplusN(N,T);
    for(size_t n=0;n<N;++n){ 
      for(size_t t=0;t<T;++t){ 
    	/**  Now need to replace the placeholders with given nodes 
    	 *   I.e. set up a context to the expression
    	 */
	ICR::EnsembleLearning::Context<data_t> context;
    	for(size_t m=0;m<M;++m){
    	  context.Assign(AP[m], m_A(n,m));
    	  context.Assign(SP[m], m_S(m,t));
    	}
	if (noise_offset) 
	  context.Assign(NP, m_noiseMean[n]);
	
    	AtimesSplusN(n,t) = m_Build.calc_gaussian(Expr,context);  
      }
    }
    m_AtimesSplusN = AtimesSplusN;


    std::cout<<"built context - nodes = "<< m_Build.number_of_nodes() <<std::endl;
     
    //join to the data
    for(size_t n=0;n<N;++n){
      for(size_t t=0;t<T;++t){
    	m_Build.join(AtimesSplusN(n,t),m_noisePrecision(n), data(n,t));
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
  
    for(size_t m=0;m<SD.size1();++m){
      for(size_t t=0;t<SD.size2(); ++t) {
	m_AtimesSplusN(m,t)->GetMoments(); //update 
      }
    }
  }

  void get_normalised_means(matrix<data_t>& A, matrix<data_t>& S)
  {
    
    for(size_t i=0;i<m_A.size1();++i){
      for(size_t j=0; j<m_A.size2();++j){
	A(i,j)= m_A(i,j)->GetMoments()[0];
      }
    }

    for(size_t i=0;i<m_S.size1();++i){
      for(size_t j=0; j<m_S.size2();++j){
	S(i,j)= m_S(i,j)->GetMoments()[0];
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
	r(i,j)= m_AtimesSplusN(i,j)->GetMoments()[0];
      }
    }
    return r;
  }

private: 
  void 
  build_vector(vector<GaussianNode>& V)
  {
    std::generate(V.begin(), V.end(), 
		  boost::bind(static_cast<GaussianNode (ICR::EnsembleLearning::Builder<data_t>::*)
			      (const data_t&, const data_t&)>
			      (&ICR::EnsembleLearning::Builder<data_t>::gaussian),
			      boost::ref(m_Build),
			      0.0, m_GaussianPrecision));
  }
  void 
  build_vector(vector<GammaNode>& V)
  {
    std::generate(V.begin(), V.end(), boost::bind(&ICR::EnsembleLearning::Builder<data_t>::gamma,
						  boost::ref(m_Build),
						  1,m_GammaVar ));
  }
  void 
  build_vector(vector<WeightsNode>& V, size_t components)
  {
    std::generate(V.begin(), V.end(), boost::bind(&ICR::EnsembleLearning::Builder<data_t>::weights,
						  boost::ref(m_Build),
						  components));
  }
  
  void
  build_mixing_matrix(Int2Type<false>, size_t N, size_t M)
  {
    vector<GaussianNode> AMean(M);
    vector<GammaNode> APrecision(M);
    build_vector(AMean);
    build_vector(APrecision);
    
    for(size_t n=0;n<N;++n){
      for(size_t m=0;m<M;++m){ 
	m_A(n,m) = m_Build.gaussian(AMean[m], APrecision[m]);
      }
    }
    std::cout<<"built mixing - nodes = "<< m_Build.number_of_nodes() <<std::endl;
  };

  void
  build_mixing_matrix(Int2Type<true>, size_t N, size_t M)
  {
    vector<GaussianNode> AMean(M);
    vector<GammaNode> APrecision(M);
    build_vector(AMean);
    build_vector(APrecision);
    
    for(size_t n=0;n<N;++n){
      for(size_t m=0;m<M;++m){ 
	m_A(n,m) = m_Build.rectified_gaussian(AMean[m], APrecision[m]);
      }
    }
    std::cout<<"built mixing (RG)- nodes = "<< m_Build.number_of_nodes() <<std::endl;
  };


  
  void
  build_source_matrix(Int2Type<false>, size_t M, size_t T, size_t Components)
  {
    
    //build A weights component for each source
    vector<WeightsNode> Weights(M);
    build_vector(Weights, Components);
    
    //m_Build hyper mean and hyper precision for each source mixture component.
    std::vector< vector<GaussianNode> > ShypMean(M);
    std::vector< vector<GammaNode> > ShypPrec(M);
    for(size_t m=0;m<M;++m){ 
      ShypMean[m].resize(Components);
      ShypPrec[m].resize(Components);
      build_vector(ShypMean[m]);
      build_vector(ShypPrec[m]);
    }
    std::cout<<"built hyperparams - nodes = "<< m_Build.number_of_nodes() <<std::endl;
    
    for(size_t m=0;m<M;++m){ 
      for(size_t t=0;t<T;++t){ 
	m_S(m,t) = m_Build.gaussian_mixture(ShypMean[m].begin(),
					    ShypPrec[m].begin(),
					    Weights[m]);
      }
    }
    
    std::cout<<"built sources - nodes = "<< m_Build.number_of_nodes() <<std::endl;
     
  }
  
  void
  build_source_matrix(Int2Type<true>, size_t M, size_t T, size_t Components)
  {
    
    //build A weights component for each source
    vector<WeightsNode> Weights(M);
    build_vector(Weights, Components);
    
    //m_Build hyper mean and hyper precision for each source mixture component.
    std::vector< vector<GaussianNode> > ShypMean(M);
    std::vector< vector<GammaNode> > ShypPrec(M);
    for(size_t m=0;m<M;++m){ 
      ShypMean[m].resize(Components);
      ShypPrec[m].resize(Components);
      build_vector(ShypMean[m]);
      build_vector(ShypPrec[m]);
    }
    std::cout<<"built hyperparams - nodes = "<< m_Build.number_of_nodes() <<std::endl;
    
    for(size_t m=0;m<M;++m){ 
      for(size_t t=0;t<T;++t){ 
	m_S(m,t) = m_Build.rectified_gaussian_mixture(ShypMean[m].begin(),
						      ShypPrec[m].begin(),
						      Weights[m]);
      }
    }
    
    std::cout<<"built sources  (RG) - nodes = "<< m_Build.number_of_nodes() <<std::endl;
     
  }
  

  ICR::EnsembleLearning::Builder<data_t> m_Build;
  ICR::EnsembleLearning::ExpressionFactory<data_t> m_Factory;
  matrix<Variable> m_A;
  matrix<Variable> m_S;
  matrix<Variable> m_AtimesSplusN; 
  vector<GaussianNode> m_noiseMean;
  vector<GammaNode>     m_noisePrecision;
  //GammaNode     m_noisePrecision;
  bool m_positive_mixing, m_positive_sources;
  double m_GaussianPrecision;
  double m_GammaVar;

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
