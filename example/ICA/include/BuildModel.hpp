#pragma once
#include "ICA.hpp"

struct Model{};

template<class data_t>
class BuildModel : public Model
{
    typedef typename ICR::ICA::Builder<data_t>::GaussianNode GaussianNode;
    typedef typename ICR::ICA::Builder<data_t>::RectifiedGaussianNode RectifiedGaussianNode;
    typedef typename ICR::ICA::Builder<data_t>::WeightsNode WeightsNode;
    typedef typename ICR::ICA::Builder<data_t>::GammaNode    GammaNode;
    typedef typename ICR::ICA::Builder<data_t>::GaussianDataNode    Datum;
    typedef typename ICR::ICA::Builder<data_t>::Variable    Variable;
    typedef typename ICR::ICA::Builder<data_t>::GaussianResultNode    ResultNode;
   
public:
  BuildModel(matrix<double>& data,
	     size_t assumed_sources,
	     size_t mixing_components,
	     bool positive_sources, 
	     bool positive_mixing,
	     bool noise_offset
	     )
    : m_Build(), 
      m_Factory(),
      m_A(data.size1(),assumed_sources),
      m_S(assumed_sources, data.size2()),
      m_noiseMean(data.size1())//,
      //m_noisePrecision(data.size1()),
      // m_noisePrecision()
  {
 

    const size_t Components = mixing_components;
    const size_t M = assumed_sources; //assumed sources
    const size_t N = data.size1(); //DataSources
    const size_t T  =data.size2(); //DataLength
    
    std::cout<<"C = "<<Components<<std::endl
	     <<"M = "<<M<<std::endl
	     <<"N = "<<N<<std::endl
	     <<"T = "<<T<<std::endl;

    
    //build A weights component for each source
    vector<WeightsNode> Weights(M);
    for(size_t m=0;m<M;++m){
      Weights[m] = m_Build.weights(Components);
    }
    std::cout<<"build weights - nodes = "<< m_Build.number_of_nodes() <<std::endl;

    
    //m_Build hyper mean and hyper precision for each source mixture component.
    std::vector< std::vector<Variable> > ShypMean(M);
    std::vector< std::vector<Variable> > ShypPrec(M);
    for(size_t m=0;m<M;++m){ 
      ShypMean[m].resize(Components);
      ShypPrec[m].resize(Components);
      for(size_t c=0;c<Components;++c){
	ShypMean[m][c]=	m_Build.gaussian(0.0,0.1);
	ShypPrec[m][c]= m_Build.gamma(1.0,1.0);
      }
    }
    std::cout<<"built hyperparams - nodes = "<< m_Build.number_of_nodes() <<std::endl;
     
    //m_Build the source approximation
    for(size_t m=0;m<M;++m){ 
      for(size_t t=0;t<T;++t){ 
	if (positive_sources) 
	  m_S(m,t) = m_Build.rectified_gaussian_mixture(ShypMean[m],
						     ShypPrec[m],
						     Weights[m]);
	else
	  m_S(m,t) = m_Build.gaussian_mixture(ShypMean[m],
					   ShypPrec[m],
					   Weights[m]);
      }
    }
    
    std::cout<<"built sources - nodes = "<< m_Build.number_of_nodes() <<std::endl;
     
    //m_Build approximation to the mixing matrix
    vector<GaussianNode> AMean(M);
    vector<GammaNode> APrecision(M);
    for(size_t m=0;m<M;++m){
      AMean[m] = m_Build.gaussian(0.0,0.1);
      APrecision[m] = m_Build.gamma(1.0,1.0);
    }
    for(size_t n=0;n<N;++n){
      for(size_t m=0;m<M;++m){ 
	if (positive_mixing) 
	  m_A(n,m) = m_Build.rectified_gaussian(AMean[m], APrecision[m]);
	else 
	  m_A(n,m) = m_Build.gaussian(AMean[m], APrecision[m]);
      }
    }
    std::cout<<"built mixing - nodes = "<< m_Build.number_of_nodes() <<std::endl;
     

    //m_Build the noise mean approximation
    // m_noiseMean(N);
    //  m_noisePrecision(N);
    for(size_t n=0;n<N;++n){
      m_noiseMean[n] = m_Build.gaussian(0.00,0.1);
      //m_noisePrecision[n] = m_Build.gamma(1.0,1.0);
    }
    
    m_noisePrecision = m_Build.gamma(1.0,1.0);
      
    //Deterministic Node.  Need to make the expression.
    //The inner product plus noise
    /** The inner product is summed over M.
     *  Therefore there are 2M + 1 placeholders required
     *  M for Anm, M for Smt and 1 for Noise
     */
    //Make a set of placeholders for all the elements in the expression.
    vector<ICR::ICA::Placeholder<data_t>*> SP(M);
    vector<ICR::ICA::Placeholder<data_t>*> AP(M);
    ICR::ICA::Placeholder<data_t>* NP = m_Factory.placeholder();
    for(size_t m=0;m<M;++m){
      SP[m] = m_Factory.placeholder();  
      AP[m] = m_Factory.placeholder(); 
    }
    //The expression in terms of these placeholders
    ICR::ICA::Expression<data_t>* Expr;
    {
      BOOST_ASSERT(M!=0);
      //First do the multiplication
      std::deque<ICR::ICA::Expression<data_t>*> prod(M);
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
	ICR::ICA::Context<data_t> context;
    	for(size_t m=0;m<M;++m){
    	  context.Assign(AP[m], m_A(n,m));
    	  context.Assign(SP[m], m_S(m,t));
    	}
	if (noise_offset) 
	  context.Assign(NP, m_noiseMean[n]);
	
    	AtimesSplusN(n,t) = m_Build.calc_gaussian(Expr,context);  
      }
    }
    std::cout<<"built context - nodes = "<< m_Build.number_of_nodes() <<std::endl;
     
    //join to the data
    for(size_t n=0;n<N;++n){
      for(size_t t=0;t<T;++t){
    	m_Build.join(AtimesSplusN(n,t),m_noisePrecision, data(n,t));
      }
    }

    std::cout<<"built data - nodes = "<< m_Build.number_of_nodes() <<std::endl;
     
  }
  ICR::ICA::Builder<data_t>& get_builder() {return m_Build;}
  matrix<Variable>& get_sources() {return m_S;}
  matrix<Variable>& get_mixing_matrix(){return m_A;}
  vector<Variable>& get_noiseMean(){return m_noiseMean;}
  //vector<Variable>& get_noisePrecision(){return m_noisePrecision;}
  
  void get_normalised_means(matrix<data_t>& A, matrix<data_t>& S)
  {
    for(size_t i=0;i<m_A.size1();++i){
      for(size_t j=0; j<m_A.size2();++j){
	A(i,j)= m_A(i,j)->GetMoments()[0];
      }
    }
    std::cout<<"got a"<<std::endl;

    for(size_t i=0;i<m_S.size1();++i){
      for(size_t j=0; j<m_S.size2();++j){
	S(i,j)= m_S(i,j)->GetMoments()[0];
      }
    }
    std::cout<<"got b"<<std::endl;

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

private: 
  ICR::ICA::Builder<data_t> m_Build;
  ICR::ICA::ExpressionFactory<data_t> m_Factory;
  matrix<Variable> m_A;
  matrix<Variable> m_S;
  vector<GaussianNode> m_noiseMean;
  //vector<GammaNode>     m_noisePrecision;
  GammaNode     m_noisePrecision;
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
