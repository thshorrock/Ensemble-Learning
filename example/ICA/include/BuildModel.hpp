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
  typedef typename ICR::EnsembleLearning::Builder<data_t>::GaussianNode GaussianNode;
  typedef typename ICR::EnsembleLearning::Builder<data_t>::RectifiedGaussianNode RectifiedGaussianNode;
  typedef typename ICR::EnsembleLearning::Builder<data_t>::WeightsNode WeightsNode;
  typedef typename ICR::EnsembleLearning::Builder<data_t>::GammaNode    GammaNode;
  typedef typename ICR::EnsembleLearning::Builder<data_t>::GaussianDataNode    Datum;
  typedef typename ICR::EnsembleLearning::Builder<data_t>::Variable    Variable;
  typedef typename ICR::EnsembleLearning::Builder<data_t>::GaussianResultNode    ResultNode;
  
  typedef  HiddenNode<Gaussian,data_t,typename ICR::EnsembleLearning::detail::TypeList::incr_id<ICR::EnsembleLearning::detail::TypeList::zeros >::type > Mean_t;
  
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

    //build A weights component for each source
    vector<WeightsNode> Weights(M);
    build_vector(Weights, Components);
    
    std::vector< MixtureVector<Gaussian,data_t,Components> > ShypMean(M);
    std::vector< MixtureVector<Gamma,data_t,Components>    > ShypPrec(M);
    
    //A set of T calculation vectors
    std::vector< CalculationVector<Gaussian,data_t,M> > cv_Sg(T);
    std::vector< CalculationVector<RectifiedGaussian,data_t,M> > cv_Srg(T);
    
    //m_Build hyper mean and hyper precision for each source mixture component.
    for(size_t m=0;m<M;++m)
      { 
	ShypMean[m] =  m_Build.template mixture_vector<Gaussian>(0.00,0.001); 
	ShypPrec[m] =  m_Build.template mixture_vector<Gamma> (1.00,0.01) ;
      }
    
    std::cout<<"built hyperparams - nodes = "<< m_Build.number_of_nodes() <<std::endl;
    for(size_t t=0;t<T;++t){
      if (positive_sources) 
	{
	  cv_Srg[t]=m_Build.template calculation_mixture<RectifiedGaussian,M>
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
	  cv_Sg[t]=m_Build.template calculation_mixture<Gaussian,M>
	    (ShypMean, 
	     ShypPrec,
	     Weights );
	  std::vector<Variable> S_over_m = to_std_vector(cv_Sg[t]);
      
	  for(size_t m=0;m<M;++m){
	    m_S(m,t) = S_over_m[m];
	  }
	}
    }


    // for(size_t m=0;m<M;++m){
    //   for(size_t t=0;t<T;++t){ 
    // 	 if (positive_sources) 
    // 	   m_S(m,t) = m_Build.rectified_gaussian_mixture(ShypMean[m], ShypPrec[m],Weights[m]);
    // 	 else
    // 	   m_S(m,t) = m_Build.gaussian_mixture(ShypMean[m], ShypPrec[m],Weights[m]);
    //   }
    // }   


    //m_Build approximation to the mixing matrix
     
    std::vector<Mean_t*> AMean(M);
    std::vector<GammaNode> APrecision(M);
    build_vector<Mean_t>(AMean);
    build_vector(APrecision);
    
    
    //A set of T calculation vectors
    std::vector< CalculationVector<Gaussian,data_t,M> > cv_Ag(N);
    std::vector< CalculationVector<RectifiedGaussian,data_t,M> > cv_Arg(N);
    
    for(size_t n=0;n<N;++n){
      if (positive_mixing) 
	{
	  cv_Arg[n] = m_Build.template calculation<RectifiedGaussian,M>(AMean, APrecision);
	  std::vector<Variable> A_over_m = to_std_vector(cv_Arg[n]);
	  for(size_t m=0;m<M;++m)
	    m_A(n,m) = A_over_m[m];
	}
      else
	{
	  cv_Ag[n] = m_Build.template calculation<Gaussian,M>(AMean, APrecision);
	  std::vector<Variable> A_over_m = to_std_vector(cv_Ag[n]);
	  for(size_t m=0;m<M;++m)
	    m_A(n,m) = A_over_m[m];
	} 
    }
    
    //Build noise
    vector<Mean_t*> m_noiseMean(N);
    vector<GammaNode> m_noisePreceision(N);
    build_vector<Mean_t>(m_noiseMean);
    build_vector(APrecision);
    // for(size_t n=0;n<N;++n){
    //   m_noiseMean[n] = m_Build.gaussian(0.0,m_GaussianPrecision);
    //   m_noisePrecision[n] = m_Build.gamma(1,m_GammaVar);
    // }
    //   for(size_t m=0;m<M;++m){
    // 	if (positive_mixing) 
    // 	  m_A(n,m) = m_Build.rectified_gaussian(AMean[m], APrecision[m]);
    // 	else
    // 	  m_A(n,m) = m_Build.gaussian(AMean[m], APrecision[m]);
	
    //   }
    // }
    //m_Build the noise mean approximation
    // m_noiseMean(N);
    //  m_noisePrecision(N);
    // for(size_t n=0;n<N;++n){
    //   m_noiseMean[n] = m_Build.gaussian(0.0,m_GaussianPrecision);
    //   m_noisePrecision[n] = m_Build.gamma(1,m_GammaVar);
    // }
    
    //m_noisePrecision = m_Build.gamma(1.0,1.0);
      

    //Deterministic Node.  Need to make the expression.
    //The inner product plus noise
    /** The inner product is summed over M.
     *  Therefore there are 2M + 1 placeholders required
     *  M for Anm, M for Smt and 1 for Noise
     */
    //Make a set of placeholders for all the elements in the expression.
    

  typedef PlaceholderFactory::make_vector<1,M>  SP;
  typedef PlaceholderFactory::make_vector<M+1,2*M> AP;
  typedef PlaceholderFactory::make_c<2*M+1> NP;
  
  std::cout<<"size SP = "<<boost::mpl::size<SP>::value<<std::endl;
  std::cout<<"size AP = "<<boost::mpl::size<AP>::value<<std::endl;

  typedef typename boost::mpl::transform<AP,SP,Times>::type the_product;
  typedef typename boost::mpl::accumulate<the_product,zero_t,Plus>::type the_dot_product;
  typename boost::mpl::accumulate<NP,the_dot_product,Plus>::type Expr;
  

  // std::vector<double> d(11,3.0);
  
  // Context::vector_t<10> all_ctx;
  // //Context::assign(all_ctx,d);
  // Context::at<Context::vector_t<10>,0> ctx_0;
  // //ctx_0.push_back(d);

  // //Context::make<mpl::integral_c<int,0> > ctx_0 = 
  // fusion::at_c<0>(all_ctx).push_back(d);
  // Context::make<mpl::integral_c<int,0> > ctx_0 =  fusion::at_c<0>(all_ctx);
    
  //  Context::make_c<0> dotprod_ctx;
  //  dotprod_ctx.args = d;

  // std::pair<double,double>  pl= proto::eval(dot_prod,ctx_0);

  // std::cout<<"pl = "<<pl.first<<std::endl;
  // std::cout<<"pl = "<<pl.second<<std::endl;




  //   std::vector< CalculationVector<Gaussian>::type > m_S(T);
  //   std::vector< CalculationVector<RectifiedGaussian>::type > m_Sp(T);
    
    
  //   for(size_t m=0;m<M;++m){
  //     for(size_t t=0;t<T;++t){ 
  // 	if (positive_sources) 
  // 	  boost::fusion::for_each(m_S(t).data(),  
  // 				  MakeMixure<Gaussian>(
  // 				  boost::bind(&Binder<data_t>::rectified_gaussian_mixture, 
  // 					      m_Build, 
  // 					      ShypMean.begin(), 
  // 					      ShypPrec.begin(),
  // 					      Weights.begin()
  // 					      )
  // 	  m_Sp(t) = m_Build.calculation_vector<RectifiedGaussian>
  // 	    (boost::bind(&Binder<data_t>::rectified_gaussian_mixture, 
  // 			 m_Build, 
  // 			 ShypMean.begin(), 
  // 			 ShypPrec.begin(),
  // 			 Weights.begin()
  // 			 )
			 
			 
								 
  // 	  m_Sp(t) = m_Build.rectified_gaussian_mixture(ShypMean[m],
  // 							ShypPrec[m],
  // 							Weights[m]);
  // 	else
    
  // 	  m_S(m,t) = m_Build.gaussian_mixture(ShypMean[m],
  // 					      ShypPrec[m],
  // 					      Weights[m]);
  // 	}
  //     }
    

 

  //   std::cout<<"built mixing - nodes = "<< m_Build.number_of_nodes() <<std::endl;

  //   // if (positive_mixing) 
  //   //   build_mixing_matrix(Int2Type<true>(),N,M);
  //   // else
  //   //   build_mixing_matrix(Int2Type<false>(),N,M);
      

  //   //Deterministic Node.  Need to make the expression.
  //   //The inner product plus noise
  //   /** The inner product is summed over M.
  //    *  Therefore there are 2M + 1 placeholders required
  //    *  M for Anm, M for Smt and 1 for Noise
  //    */
  //   //Make a set of placeholders for all the elements in the expression.
    
    
  // typedef Factory::make_vector<1,M>  v1;


  // typedef Placeholder::make_vector<1,M>  SP;
  // typedef Placeholder::make_vector<M+1,2*M> AP;
  // typedef Placeholder::make<2*M+1> NP;
  
  // std::cout<<"size SP = "<<mpl::size<SP>::value<<std::endl;
  // std::cout<<"size AP = "<<mpl::size<AP>::value<<std::endl;

  // typedef mpl::transform<AP,SP,times_f>::type prod;
  // mpl::accumulate<prod,zero_t,plus_f>::type dot_prod;
  // mpl::accumulate<NP,dot_prod,plus_f>::type Expr;
  

  // // std::vector<double> d(11,3.0);
  
  // // Context::vector_t<10> all_ctx;
  // // //Context::assign(all_ctx,d);
  // // Context::at<Context::vector_t<10>,0> ctx_0;
  // // //ctx_0.push_back(d);

  // // //Context::make<mpl::integral_c<int,0> > ctx_0 = 
  // // fusion::at_c<0>(all_ctx).push_back(d);
  // // Context::make<mpl::integral_c<int,0> > ctx_0 =  fusion::at_c<0>(all_ctx);
    
  // //  Context::make_c<0> dotprod_ctx;
  // //  dotprod_ctx.args = d;

  // // std::pair<double,double>  pl= proto::eval(dot_prod,ctx_0);

  // // std::cout<<"pl = "<<pl.first<<std::endl;
  // // std::cout<<"pl = "<<pl.second<<std::endl;

  
    
  // // vector<ICR::EnsembleLearning::Placeholder<data_t>*> SP(M);
  // // vector<ICR::EnsembleLearning::Placeholder<data_t>*> AP(M);
  // // ICR::EnsembleLearning::Placeholder<data_t>* NP = m_Factory.placeholder();
  // // for(size_t m=0;m<M;++m){
  // //   SP[m] = m_Factory.placeholder();  
  // //   AP[m] = m_Factory.placeholder(); 
  // // }
  // // //The expression in terms of these placeholders
  // // ICR::EnsembleLearning::Expression<data_t>* Expr;
  // // {
  // //   BOOST_ASSERT(M!=0);
  // //   //First do the multiplication
  // //   std::deque<ICR::EnsembleLearning::Expression<data_t>*> prod(M);
  // //   for(size_t m=0;m<M;++m){
  // // 	prod[m] = m_Factory.Multiply(AP[m], SP[m]);
  // //   }
  // //   //then add them up
  // //   Expr = prod[0];
  // //   prod.pop_front();
  // //   while(prod.size()!=0)
  // // 	{
  // // 	  Expr = m_Factory.Add(Expr,prod[0]);
  // // 	  prod.pop_front();
  // // 	}
  // //   if (noise_offset) 
  // // 	//Finnaly add the Noise placeholder
  // // 	Expr = m_Factory.Add(Expr, NP);
      
  // // }
    
    
  // //   //make vectors of size m.
  // //   BOOST_AUTO(AMean,(m_Build.calculation_vector<Gaussian,0>(0.0,0.001)));
  // //   BOOST_AUTO(APrec,(m_Build.calculation_vector<Gamma,0>(1.0,0.001)));
  // //   BOOST_AUTO(As,(m_Build.calculation_vector<Gaussian,0>(function,m_Build, boost:bind(AMean,APrec))));

  // //   contextAs = context(As)
      
  // //   for(size_t t=0;t<T;++t){ 
	
  // // 	WeightsNode Weights = m_Build.weights();
  // // 	BOOST_AUTO(ShypMean,(m_Build.mixture_vector<Gaussian,0>(0.0,0.001)));
  // // 	BOOST_AUTO(ShypPrec,(m_Build.mixture_vector<Gamma,0>(1.0,0.001)));
  // // 	BOOST_AUTO(Ss,(m_Build.calculation_vector<Gaussian,0>(function,m_Build,boost:bind(ShypMean,ShypPrec,weights))));

  // //   }
      
      
  // //   for(size_t i=0;i<10;++i){
  // // 	delete m_vec;
  // // 	m_vec =  new vector(make_vector(*m_vec_old, m_Build.gaussian(0,0.1)));
  // // 	delete m_vec_old;
  // // 	m_vec_old = new vector(m_vec);
	  
  // // 	GaussianNode* A = new GaussianNode(m_Build.gaussian(0,0.1));
	
	
  // //   }


  // matrix<ResultNode> AtimesSplusN(N,T);
  // for(size_t n=0;n<N;++n){ 

  //   // 	for(size_t m=0;m<M;++m){
	  
  //   // 	  context.push_back(m_A(n,m))
  //   // 	}
  //   for(size_t t=0;t<T;++t){ 
      
  //     /**  Now need to replace the placeholders with given nodes 
  //      *   I.e. set up a context to the expression
  //      */
  //     typedef typename boost::mpl::vector< repeat<Gaussian, M> types
  //     cA = m_Builder.calc_vector<boost::mpl::vector<,M,,M>;
  //     cS = m_Builder.calc_vector<Gaussian>;
  //     ICR::EnsembleLearning::Context<data_t> context(AMean,APrec);
  //     for(size_t m=0;m<M;++m){
	  
  // 	context.push_back(m_A(n,m));
  //     }
	
	
  //     for(size_t m=0;m<M;++m){
  // 	context.push_back(m_S(n,m));
  //     }
  //     if (noise_offset) 
  // 	context.push_back(NP, m_noiseMean[n]);
	
  //     AtimesSplusN(n,t) = m_Build.calc_gaussian(Expr,context);  
  //   }
  // }
  // m_AtimesSplusN = AtimesSplusN;


  // std::cout<<"built context - nodes = "<< m_Build.number_of_nodes() <<std::endl;
     
  // //join to the data
  // for(size_t n=0;n<N;++n){
  //   for(size_t t=0;t<T;++t){
  //     m_Build.join(AtimesSplusN(n,t),m_noisePrecision(n), data(n,t));
  //   }
  // }

  // std::cout<<"built data - nodes = "<< m_Build.number_of_nodes() <<std::endl;
     
  }
  ICR::EnsembleLearning::Builder<data_t>& get_builder() {return m_Build;}
  matrix<Variable>& get_sources() {return m_S;}
  matrix<Variable>& get_mixing_matrix(){return m_A;}
  vector<Variable>& get_noiseMean(){return m_noiseMean;}
  //vector<Variable>& get_noisePrecision(){return m_noisePrecision;}
  
  void set_means(matrix<data_t>&Means)
  {
    // //First check that number of rows and columns less than m_S.

    // std::cout<<"setting means"<<std::endl;
    // if (
    // 	(Means.size1()>m_S.size1())
    // 	||
    // 	(Means.size2()>m_S.size2())
    // 	)
    //   {
    // 	std::cout<<"Means size = ("<<Means.size1()<<","<<Means.size2()<<")"<<std::endl;
    // 	std::cout<<"Model size = ("<<m_S.size1()<<","<<m_S.size2()<<")"<<std::endl;
    // 	throw size_exception();
    //   }
    
    // //for each row 
    // for(size_t m=0;m<Means.size1();++m){
    //   //set m^th column of A matrix to be one.
    //   for(size_t n=0;n<m_A.size1();++n){
    // 	if (m_positive_mixing) 
    // 	  ICR::EnsembleLearning::SetMean(static_cast<RectifiedGaussianNode>(m_A(m,n)), 1.0);
    // 	else
    // 	  ICR::EnsembleLearning::SetMean(static_cast<GaussianNode>(m_A(m,n)), 1.0);
	  
    //   }
    //   //for each column
    //   for(size_t t=0;t<Means.size2(); ++t) {
    // 	//set the mean
    // 	if (m_positive_sources) 
    // 	  ICR::EnsembleLearning::SetMean(static_cast<RectifiedGaussianNode>(m_S(m,t)), Means(m,t));
    // 	else
    // 	  ICR::EnsembleLearning::SetMean(static_cast<GaussianNode>(m_S(m,t)), Means(m,t));
	
    // 	m_AtimesSplusN(m,t)->GetMoments(); //update
    //   }
      
    // }
  }

  void set_sigmas(matrix<data_t>&SD)
  {
    // std::cout<<"setting Sigmas"<<std::endl;

    // //First check that number of rows and columns less than m_S.
    // if (
    // 	(SD.size1()>m_S.size1())
    // 	||
    // 	(SD.size2()>m_S.size2())
    // 	)
    //   {
    // 	std::cout<<"SD    size = ("<<SD.size1()<<","<<SD.size2()<<")"<<std::endl;
    // 	std::cout<<"Model size = ("<<m_S.size1()<<","<<m_S.size2()<<")"<<std::endl;
    // 	throw size_exception();
    //   }
    
    // //for each row 
    // for(size_t m=0;m<SD.size1();++m){
    //   //set m^th column of A matrix to be the average sd of row.
    //   matrix_row<matrix<double> > row(SD,m);
    //   double av_sd = std::accumulate(row.begin(), row.end(), 0.0)/row.size();
    //   for(size_t n=0;n<m_A.size1();++n){
    // 	if (m_positive_mixing) 
    // 	  ICR::EnsembleLearning::SetStandardDeviation(static_cast<RectifiedGaussianNode>(m_A(n,m)), av_sd);
    // 	else
    // 	  ICR::EnsembleLearning::SetStandardDeviation(static_cast<GaussianNode>(m_A(n,m)), av_sd);
    //   }
    //   //for each column
    //   for(size_t t=0;t<SD.size2(); ++t) {
    // 	//set the mean
    // 	if (m_positive_sources) 
    // 	  ICR::EnsembleLearning::SetStandardDeviation(static_cast<RectifiedGaussianNode>(m_S(m,t)), SD(m,t));
    // 	else
    // 	  ICR::EnsembleLearning::SetStandardDeviation(static_cast<GaussianNode>(m_S(m,t)), SD(m,t));


    // 	m_AtimesSplusN(m,t)->GetMoments(); //update 
    //   }
    // }
  }

  void set_mixing_mean(matrix<data_t>& M)
  {
    // std::cout<<"setting Mixing"<<std::endl;

    // //First check that number of rows and columns less than m_A.
    // if (
    // 	(M.size1()>m_A.size1())
    // 	||
    // 	(M.size2()>m_A.size2())
    // 	)
    //   {
    // 	std::cout<<"Mixing size = ("<<M.size1()<<","<<M.size2()<<")"<<std::endl;
    // 	std::cout<<"Model  size = ("<<m_A.size1()<<","<<m_A.size2()<<")"<<std::endl;
    // 	throw size_exception();
    //   }
    
    // //for each row 
    // for(size_t m=0;m<M.size1();++m){
    //   for(size_t n=0;n<M.size2();++n){
    // 	if (m_positive_mixing) 
    // 	  ICR::EnsembleLearning::SetMean(static_cast<RectifiedGaussianNode>(m_A(n,m)), M(n,m));
    // 	else
    // 	  ICR::EnsembleLearning::SetMean(static_cast<GaussianNode>(m_A(n,m)), M(n,m));
    //   }
    //   //m_AtimesSplusN(m,t)->GetMoments(); //update 
    // }
  
    // for(size_t m=0;m<M.size1();++m){
    //   for(size_t t=0;t<M.size2(); ++t) {
    // 	m_AtimesSplusN(m,t)->GetMoments(); //update 
    //   }
    // }
  }

  void get_normalised_means(matrix<data_t>& A, matrix<data_t>& S)
  {
    
    // for(size_t i=0;i<m_A.size1();++i){
    //   for(size_t j=0; j<m_A.size2();++j){
    // 	A(i,j)= m_A(i,j)->GetMoments()[0];
    //   }
    // }

    // for(size_t i=0;i<m_S.size1();++i){
    //   for(size_t j=0; j<m_S.size2();++j){
    // 	S(i,j)= m_S(i,j)->GetMoments()[0];
    //   }
    // }

    // //normalise
    // for(size_t j=0;j<A.size2();++j){
    //   matrix_column<matrix<data_t> > col(A, j);
    //   const vector<data_t> col2 = element_prod(col,col);
    //   const data_t norm = 1.0/std::sqrt(std::accumulate(col2.begin(), col2.end(), 0.0));
    //   col*=norm;
      
    //   matrix_row<matrix<data_t> > row(S, j);
    //   row/=norm;
    // }
  }
  matrix<data_t> get_results() const
  {
    // matrix<data_t> r(m_AtimesSplusN.size1(), m_AtimesSplusN.size2());
    // for(size_t i=0;i<m_AtimesSplusN.size1();++i){
    //   for(size_t j=0; j<m_AtimesSplusN.size2();++j){
    // 	r(i,j)= m_AtimesSplusN(i,j)->GetMoments()[0];
    //   }
    // }
    // return r;
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
  
  void
  build_mixing_matrix(Int2Type<false>, size_t N, size_t M)
  {
    // GaussianNode AMean = m_Build.gaussian(0.0,m_GaussianPrecision);
    // GammaNode    APrecision = m_build.gamma(1.0,m_GammaVar);
    // return m_Build.gaussian(AMean,APrec);
    
    
    // // vector<GaussianNode> AMean(M);
    // // vector<GammaNode> APrecision(M);
    // // build_vector(AMean);
    // // build_vector(APrecision);
    
    // // for(size_t n=0;n<N;++n){
    // //   for(size_t m=0;m<M;++m){ 
    // // 	m_A(n,m) = m_Build.gaussian(AMean[m], APrecision[m]);
    // //   }
    // // }
    // // std::cout<<"built mixing - nodes = "<< m_Build.number_of_nodes() <<std::endl;
  };

  void
  build_mixing_matrix(Int2Type<true>, size_t N, size_t M)
  {
    // vector<GaussianNode> AMean(M);
    // vector<GammaNode> APrecision(M);
    // build_vector(AMean);
    // build_vector(APrecision);
    
    // for(size_t n=0;n<N;++n){
    //   for(size_t m=0;m<M;++m){ 
    // 	m_A(n,m) = m_Build.rectified_gaussian(AMean[m], APrecision[m]);
    //   }
    // }
    // std::cout<<"built mixing (RG)- nodes = "<< m_Build.number_of_nodes() <<std::endl;
  };


  
  void
  build_source_matrix(Int2Type<false>, size_t M, size_t T, size_t Components)
  {
    
    // //build A weights component for each source
    // vector<WeightsNode> Weights(M);
    // build_vector(Weights, Components);
    
    // //m_Build hyper mean and hyper precision for each source mixture component.
    // std::vector< vector<GaussianNode> > ShypMean(M);
    // std::vector< vector<GammaNode> > ShypPrec(M);
    // for(size_t m=0;m<M;++m){ 
    //   ShypMean[m].resize(Components);
    //   ShypPrec[m].resize(Components);
    //   build_vector(ShypMean[m]);
    //   build_vector(ShypPrec[m]);
    // }
    // std::cout<<"built hyperparams - nodes = "<< m_Build.number_of_nodes() <<std::endl;
    
    // for(size_t m=0;m<M;++m){ 
    //   for(size_t t=0;t<T;++t){ 
    // 	m_S(m,t) = m_Build.gaussian_mixture(ShypMean[m].begin(),
    // 					    ShypPrec[m].begin(),
    // 					    Weights[m]);
    //   }
    // }
    
    // std::cout<<"built sources - nodes = "<< m_Build.number_of_nodes() <<std::endl;
     
  }
  
  void
  build_source_matrix(Int2Type<true>, size_t M, size_t T, size_t Components)
  {
    
    // //build A weights component for each source
    // vector<WeightsNode> Weights(M);
    // build_vector(Weights, Components);
    
    // //m_Build hyper mean and hyper precision for each source mixture component.
    // std::vector< vector<GaussianNode> > ShypMean(M);
    // std::vector< vector<GammaNode> > ShypPrec(M);
    // for(size_t m=0;m<M;++m){ 
    //   ShypMean[m].resize(Components);
    //   ShypPrec[m].resize(Components);
    //   build_vector(ShypMean[m]);
    //   build_vector(ShypPrec[m]);
    // }
    // std::cout<<"built hyperparams - nodes = "<< m_Build.number_of_nodes() <<std::endl;
    
    // for(size_t m=0;m<M;++m){ 
    //   for(size_t t=0;t<T;++t){ 
    // 	m_S(m,t) = m_Build.rectified_gaussian_mixture(ShypMean[m].begin(),
    // 						      ShypPrec[m].begin(),
    // 			 			      Weights[m]);
    //   }
    // }
    
    // std::cout<<"built sources  (RG) - nodes = "<< m_Build.number_of_nodes() <<std::endl;
     
  }
  

  ICR::EnsembleLearning::Builder<data_t> m_Build;
  //ICR::EnsembleLearning::ExpressionFactory<data_t> m_Factory;
  matrix<Variable> m_A;
  matrix<Variable> m_S;
  matrix<Variable> m_AtimesSplusN; 
  vector<Mean_t*> m_noiseMean;
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
  // for(size_t j=0; j<M.size2();++j){
  //   for(size_t i=0;i<M.size1();++i){
  //     out<< M(i,j) << "\t";
  //   }
  //   out<<"\n";
  // }
  return out;
}
