

#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/preprocessor/arithmetic/mul.hpp>
#  define ENSEMBLE_LEARNING_COMPONENTS 5
#  define ENSEMBLE_LEARNING_SOURCES 10
#  define ENSEMBLE_LEARNING_PLACEHOLDERS BOOST_PP_ADD(BOOST_PP_MUL(ENSEMBLE_LEARNING_SOURCES,2),1)
#  define FUSION_MAX_VECTOR_SIZE 25
#include "ExampleData.hpp"
//#include "EnsembleLearning.hpp"
#include "BuildModel.hpp"

#include <boost/program_options.hpp>
#include "boost/filesystem.hpp"   // includes all needed Boost.Filesystem declarations

#include <iostream>
#include <fstream>
#include <string>

using namespace ICR::EnsembleLearning;
using namespace boost::numeric::ublas;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

int
main  (int ac, char **av)
{
  //start the Random number generator with fixed seed so model reproducable
  ICR::EnsembleLearning::Random::Instance(10);
  
  bool   use_float = false;
  const size_t assumed_sources   = ENSEMBLE_LEARNING_SOURCES;
  const size_t mixing_components = ENSEMBLE_LEARNING_COMPONENTS;

  bool   positive_source = false;
  bool   positive_mixing = false;
  bool   model_noise_offset = false;
  bool   example = false;
  bool   transpose_priors = false;
  bool   transpose_mixing = false;
  double convergence_criterium;
  double GaussianPrecision;
  double GammaPrecision;
  size_t max_iterations;
  
  std::string data_file;
  std::string output_directory;
  std::string mean_file;
  std::string sigma_file;
  std::string mixing_mean_file;
  
  
  //parse the command line
  po::options_description gen("General options");
  gen.add_options()
    ("help", "Produce this help message")
    ("float,f", 
     "Use floating point (rather than double) arithmatic.")
    ("example",  
     "Use example data rather than a data file")
    ;
  
  po::options_description file("File Handling");
  file.add_options()
    ("output-directory,o", po::value<std::string>(&output_directory), 
     "Set the directory for output")
    ("mean", po::value<std::string>(&mean_file)->default_value(""), 
     "filename containing prior knowledge of the means")
    ("sigma", po::value<std::string>(&sigma_file)->default_value(""), 
     "filename containing prior knowledge of the standard-deviations")
    ("mixing-mean", po::value<std::string>(&mixing_mean_file)->default_value(""), 
     "filename containing prior knowledge of the Mixing means")
    ;
  
  // Hidden options, will be allowed both on command line and
  // in config file, but will not be shown to the user.
  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file",  po::value<std::string>(&data_file), 
     "Set the input data file")
    ;    
  
  po::options_description inference("Inference");
  inference.add_options()
    // ("sources,s", po::value<size_t>(&assumed_sources)->default_value(5), 
    //  "The number of sources to infer (a guess)")
    // ("mixing,m",  po::value<size_t>(&mixing_components)->default_value(5), 
    //  "The number of components to use in Mixture Models")
    ("positive-source",  
     "Assume that the sources are strictly positive")
    ("positive-mixing",  
     "Assume that the mixing matrix is strictly positive")
    ("noise-offset", 
     "Allow the noise to have a non-zero mean")
    ("convergence,c", po::value<double>(&convergence_criterium)->default_value(1e-3), 
     "The change in evidence when the model is assumed to have converged")
    ("iterations,i", po::value<size_t>(&max_iterations)->default_value(1000), 
     "The maximum number of interations")
    ("Gaussian-precision", po::value<double>(&GaussianPrecision)->default_value(0.01), 
     "The precision of the Gaussian Priors in the model")
    ("Gamma-precision", po::value<double>(&GammaPrecision)->default_value(0.01), 
     "The precision of the Gamma Priors in the model")
    ("transpose-priors",  
     "Transpose the data in the mean and sigma prior data files")
    ("transpose-mixing",  
     "Transpose the data in the mixing  prior data files")
    ;
  
  po::options_description cmdline_options;
  cmdline_options.add(gen).add(file).add(hidden).add(inference);

  po::options_description visible("Allowed options");
  visible.add(gen).add(file).add(inference);
      
  po::positional_options_description p;
  p.add("input-file", -1);

  po::variables_map vm;
  
  try
    {
      po::store(po::command_line_parser(ac, av).
		options(cmdline_options).positional(p).run(), vm);
      po::notify(vm);    
    }
  catch(boost::exception_detail::clone_impl<boost::exception_detail::error_info_injector<boost::program_options::multiple_occurrences> >& e)
    {
      std::cout << "Problems processing options: " 
		<< e.what() << "\n\n"
		<< visible  << "\n";
      return 1;
    }

  catch(boost::exception_detail::clone_impl<boost::exception_detail::error_info_injector<boost::program_options::ambiguous_option> >& e)
     {
      std::cout << "Problem processing options: " 
		<< e.what() << "\n\n"
		<< visible  << "\n";
      return 1;
    }


  if (vm.count("help")) {
    std::cout << visible  << "\n";
    return 1;
  }

  if (vm.count("float")) 
    use_float = true;
  if (vm.count("example")) 
    example = true;

  if (vm.count("positive-source")) 
    positive_source = true;
  if (vm.count("positive-mixing")) 
    positive_mixing = true;
  if (vm.count("noise-offset")) 
    model_noise_offset = true;
  if (vm.count("transpose-priors")) 
    transpose_priors = true;
  if (vm.count("transpose-mixing")) 
    transpose_mixing = true;


  if ( !fs::exists(output_directory) ) {
    std::cout << "Output directory '"<<output_directory<<"' does not exist!\n\n"
	      << visible  << "\n";
    return 1;
  }

  //check to see if mean file inputted
  matrix<double> Means;
  if ( mean_file != "") {
    if (!fs::exists(mean_file))  { //wrong filename
      std::cout << "The file '"<<mean_file<<"' does not exist!\n\n"
		<< visible  << "\n";
      return 1;
    }
    //load the data.
    std::cout<<"loading means"<<std::endl;
    std::string line;
    std::ifstream myfile (mean_file.c_str());
    std::vector< std::vector<double> > data;
    while ( getline( myfile, line ) )
      {
	std::stringstream ss(line);
	std::vector<double> v;
	std::copy( 
		  std::istream_iterator<double>(ss), 
		  std::istream_iterator<double>( ),
		  std::back_inserter(v)
		   ); // copies all data into bufferdata[i]
	data.push_back(v);
      }
    if (transpose_priors)
      { 
	const size_t rows = data[0].size();
	const size_t cols = data.size(); //assume all the same size.
	Means = matrix<double>(rows,cols);
	for(size_t r=0;r<Means.size1();++r){
	  for(size_t c=0;c<Means.size2();++c){
	    Means(r,c) = data[c][r];
	  }
	}
      }
    else 
      {
	const size_t rows = data.size();
	const size_t cols = data[0].size(); //assume all the same size.
      
	Means = matrix<double>(rows,cols);
	for(size_t r=0;r<Means.size1();++r){
	  for(size_t c=0;c<Means.size2();++c){
	    Means(r,c) = data[r][c];
	  }
	}
      }
  }
  //check to see if sigma file inputted and exists
  matrix<double> Sigmas;
  if ( sigma_file != "" ) { 
    if (!fs::exists(sigma_file) ) {//wrong filename
      std::cout << "The file '"<<sigma_file<<"' does not exist!\n\n"
		<< visible  << "\n";
      return 1;
    }
    //load the data.
    std::cout<<"loading means"<<std::endl;
    std::string line;
    std::ifstream myfile (sigma_file.c_str());
    std::vector< std::vector<double> > data;
    while ( getline( myfile, line ) )
      {
	std::stringstream ss(line);
	std::vector<double> v;
	std::copy( 
		  std::istream_iterator<double>(ss), 
		  std::istream_iterator<double>( ),
		  std::back_inserter(v)
		   ); // copies all data into bufferdata[i]
	data.push_back(v);
      }
    if (transpose_priors)
      { 
	const size_t rows = data[0].size();
	const size_t cols = data.size(); //assume all the same size.
	Sigmas = matrix<double>(rows,cols);
	for(size_t r=0;r<Sigmas.size1();++r){
	  for(size_t c=0;c<Sigmas.size2();++c){
	    Sigmas(r,c) = data[c][r];
	  }
	}
      }
    else 
      {
	const size_t rows = data.size();
	const size_t cols = data[0].size(); //assume all the same size.
	Sigmas = matrix<double>(rows,cols);
	for(size_t r=0;r<Sigmas.size1();++r){
	  for(size_t c=0;c<Sigmas.size2();++c){
	    Sigmas(r,c) = data[r][c];
	  }
	}
      }
  }
    
  //check to see if mixing_mean file inputted and exists
  matrix<double> MixingMean;
  if ( mixing_mean_file != "" ) { 
    if (!fs::exists(mixing_mean_file) ) {//wrong filename
      std::cout << "The file '"<<mixing_mean_file<<"' does not exist!\n\n"
		<< visible  << "\n";
      return 1;
    }
    //load the data.
    std::cout<<"loading mixing"<<std::endl;
    std::string line;
    std::ifstream myfile (mixing_mean_file.c_str());
    std::vector< std::vector<double> > data;
    while ( getline( myfile, line ) )
      {
	std::stringstream ss(line);
	std::vector<double> v;
	std::copy( 
		  std::istream_iterator<double>(ss), 
		  std::istream_iterator<double>( ),
		  std::back_inserter(v)
		   ); // copies all data into bufferdata[i]
	data.push_back(v);
      }
    if (transpose_mixing)
      { 
	const size_t rows = data[0].size();
	const size_t cols = data.size(); //assume all the same size.
	MixingMean = matrix<double>(rows,cols);
	for(size_t r=0;r<MixingMean.size1();++r){
	  for(size_t c=0;c<MixingMean.size2();++c){
	    MixingMean(r,c) = data[c][r];
	  }
	}
      }
    else 
      {
	const size_t rows = data.size();
	const size_t cols = data[0].size(); //assume all the same size.
	MixingMean = matrix<double>(rows,cols);
	for(size_t r=0;r<MixingMean.size1();++r){
	  for(size_t c=0;c<MixingMean.size2();++c){
	    MixingMean(r,c) = data[r][c];
	  }
	}
      }
  }

  matrix<double> Data;
  if ( !fs::exists(data_file) )
    {
      if (example)
	{
	  size_t number_of_records = 5;
	  size_t number_of_sources = 3;
	  size_t samples_per_record = 100;
	  
	  fs::path data_file = fs::path(output_directory)/fs::path("OrigninalData.txt");
	  fs::path source_file = fs::path(output_directory)/fs::path("SourceData.txt");
	  fs::path mixing_file = fs::path(output_directory)/fs::path("MixingData.txt");
	  
	  Sources S(number_of_sources,samples_per_record, positive_source);
	  Mixing  M(number_of_sources,number_of_records, positive_mixing);
	  
	  Data =  prod(M,S);
	  AddNoise(Data,5000);
  
	  std::ofstream data(data_file.string().c_str());
	  std::ofstream source(source_file.string().c_str());
	  std::ofstream mixing(mixing_file.string().c_str());
	  data<<Data;
	  source<<matrix<double>(S);
	  mixing<<matrix<double>(M);
  
	} 
      else
	{
	  std::cout << "Data file  '"<<data_file<<"' does not exist - and the example option is not set\n\n"
		    << visible  << "\n";
	  return 1;
	}
    }
  else
    {
      //load the data.
      std::cout<<"loading data"<<std::endl;

      std::string line;
      std::ifstream myfile (data_file.c_str());
      std::vector< std::vector<double> > data;
      while ( getline( myfile, line ) )
	{
	  std::stringstream ss(line);
	  std::vector<double> v;
	  std::copy( 
	  	    std::istream_iterator<double>(ss), 
	  	    std::istream_iterator<double>( ),
	  	    std::back_inserter(v)
		     ); // copies all data into bufferdata[i]
	  data.push_back(v);
	}

      const size_t rows = data.size();
      const size_t cols = data[0].size(); //assume all the same size.
      
      Data = matrix<double>(rows,cols);
      for(size_t r=0;r<Data.size1();++r){
	for(size_t c=0;c<Data.size2();++c){
	  Data(r,c) = data[r][c];
	}
      }
      std::cout<<"data loaded"<<std::endl;

    }
  

  std::cout<<"output direcotry = "<<output_directory<<std::endl;
  std::cout<<"input file = "<<data_file<<std::endl;


  
  fs::path source_file = fs::path(output_directory)/fs::path("InferedSources.txt");
  fs::path mixing_file = fs::path(output_directory)/fs::path("InferedMixing.txt");
  fs::path cost_file   = fs::path(output_directory)/fs::path("Cost.txt");
  fs::path result_file   = fs::path(output_directory)/fs::path("InferedResult.txt");
  fs::path result_file2   = fs::path(output_directory)/fs::path("InferedResult2.txt");
  fs::path precisions_file   = fs::path(output_directory)/fs::path("InferedNoise.txt");
  

  //   const size_t Components = mixing_components;
  //   const size_t M = assumed_sources; //assumed sources
  //   const size_t N = 5; //DataSources
  //   typedef double data_t;
    
  //   Builder<data_t> m_Build;
    
  // typedef CalculationVector<Gaussian,data_t,1,assumed_sources> Ag_t;
  // typedef CalculationVector<Gaussian,data_t,assumed_sources+1,assumed_sources> Sg_t;
  // typedef  ICR::EnsembleLearning::Builder<data_t>::Variable    Variable;
  
  // typedef PlaceholderFactory::make_vector<1,M> AP;
  // typedef PlaceholderFactory::make_vector<M+1,2*M>  SP;

  // std::cout<<"size SP = "<<boost::mpl::size<SP>::value<<std::endl;
  // std::cout<<"size AP = "<<boost::mpl::size<AP>::value<<std::endl;

  // typedef  boost::mpl::transform<AP,SP,Times>::type the_product;
  // typedef  boost::mpl::accumulate<the_product,zero_t,Plus>::type the_inner_product;
  // the_inner_product Expr;
  
  // std::vector<Sg_t> cv_Sg(1);
  // std::vector<Ag_t> cv_Ag(1);

  // cv_Sg[0]=m_Build.calculation_vector<Gaussian,M+1,M>
  //   (0.0,0.1 );

  // std::vector<Variable> S_over_m = to_std_vector(cv_Sg[0]);
  // for(size_t m=0;m<M;++m){
  //   std::cout<<"S["<<m<<"] = "<<S_over_m[m]<<std::endl;
  //   std::cout<<"S["<<m<<"] = "<<Mean(S_over_m[m])<<std::endl;

  // }

  // cv_Ag[0] = m_Build.calculation_vector<Gaussian,1,M>(0.0,0.1);
  // std::vector<Variable> A_over_m = to_std_vector(cv_Ag[0]);
  // for(size_t m=0;m<M;++m){
  //   std::cout<<"A["<<m<<"] = "<<A_over_m[m]<<std::endl;
  //   std::cout<<"A["<<m<<"] = "<<Mean(A_over_m[m])<<std::endl;
  // }
  
  // typedef fuse<Ag_t,Sg_t> AS_t;
  // AS_t AS(cv_Ag[0],cv_Sg[0] );
  
  // Context<data_t,AS_t> context(AS);
  // Builder<data_t>::GaussianResultNode  AtimesSplusN = m_Build.calc_gaussian(Expr,context); 
  // Builder<data_t>::GammaNode  Prec = m_Build.gamma(1.0,0.1); 
  
  // m_Build.join(AtimesSplusN,Prec, 0.1);
      
  // //std::cout<<"R = "<<AtimesSplusN->GetMoments()->operator[](0)<<std::endl;
  // std::cout<<"R calc = "<<
  //   S_over_m[0]->GetMoments()->operator[](0)*A_over_m[0]->GetMoments()->operator[](0)
  //   +
  //   S_over_m[1]->GetMoments()->operator[](0)*A_over_m[1]->GetMoments()->operator[](0)
  //    +
  //    S_over_m[2]->GetMoments()->operator[](0)*A_over_m[2]->GetMoments()->operator[](0)
  //    +
  //    S_over_m[3]->GetMoments()->operator[](0)*A_over_m[3]->GetMoments()->operator[](0)
  //    +
  //    S_over_m[4]->GetMoments()->operator[](0)*A_over_m[4]->GetMoments()->operator[](0)
  // 	   <<std::endl;

  // std::cout<<"R = "<<Mean(AtimesSplusN)<<std::endl;
  // Coster dummy;
  // S_over_m[0]->Iterate(dummy);
  // S_over_m[1]->Iterate(dummy);
  // S_over_m[2]->Iterate(dummy);
  // S_over_m[3]->Iterate(dummy);
  // S_over_m[4]->Iterate(dummy);
  // std::cout<<S_over_m[0]<<std::endl;

  

  if (use_float) 
    {
      std::cout<<"Building Model"<<std::endl;
      BuildModel<float, assumed_sources, mixing_components> 
	Model(Data, 
	      positive_source,
	      positive_mixing,
	      model_noise_offset,
	      GaussianPrecision,
	      GammaPrecision
	      );
      size_t size1 = Data.size1();
      size_t size2 = Data.size2();
      Data = matrix<float>(); //clear the memory (clear function doesn't do this)

      {
	
  	std::ofstream sources(source_file.string().c_str());
  	std::ofstream mixing_matrix(mixing_file.string().c_str());
  	std::ofstream result_matrix(result_file.string().c_str());
  	std::ofstream result_matrix2(result_file2.string().c_str());
  	matrix<float> A(size1, assumed_sources);
  	matrix<float> S(assumed_sources,size2);
  	Model.get_normalised_means(A,S);
  	sources<<S;
  	mixing_matrix<< matrix<float>(trans(A));
  	result_matrix<<Model.get_results();
  	result_matrix2<< matrix<float>(prod(A,S));
      }
      
      ICR::EnsembleLearning::Builder<float> Build = Model.get_builder();
      Build.set_cost_file(cost_file.string());
      std::cout<<"Running!"<<std::endl;
      bool converged = false;
      size_t count = 0;
      while(!converged && count<500) {
      	converged = Build.run(convergence_criterium,max_iterations);

  	std::ofstream sources(source_file.string().c_str());
  	std::ofstream mixing_matrix(mixing_file.string().c_str());
  	std::ofstream result_matrix(result_file.string().c_str());
  	std::ofstream result_matrix2(result_file2.string().c_str());

  	matrix<float> A(size1, assumed_sources);
  	matrix<float> S(assumed_sources,size2);
      	Model.get_normalised_means(A,S);
      	sources<<S;
      	mixing_matrix<< matrix<float>(trans(A));
      	result_matrix<<Model.get_results();
      	result_matrix2<< matrix<float>(prod(A,S));
      	++count;
      }
    }
  else
    {

  std::cout<<"Building Model"<<std::endl;
  BuildModel<double, assumed_sources, mixing_components> 
    Model(Data, 
  	  positive_source,
  	  positive_mixing,
  	  model_noise_offset,
  	  GaussianPrecision,
  	  GammaPrecision
  	  );
      
  size_t size1 = Data.size1();
  size_t size2 = Data.size2();
  Data = matrix<double>(); //clear the memory (clear function doesn't do this)

  if (mean_file !="") 
    Model.set_means(Means);
  if (sigma_file != "")
    Model.set_sigmas(Sigmas);
  if (mixing_mean_file != "")
    Model.set_mixing_mean(MixingMean);
  {

    std::ofstream sources(source_file.string().c_str());
    std::ofstream mixing_matrix(mixing_file.string().c_str());
    std::ofstream result_matrix(result_file.string().c_str());
    std::ofstream result_matrix2(result_file2.string().c_str());
    std::ofstream precisions(precisions_file.string().c_str());
  
    matrix<double> A(size1, assumed_sources);
    matrix<double> S(assumed_sources,size2);
    Model.get_normalised_means(A,S);
    sources<<S;
    mixing_matrix<< matrix<double>(trans(A));
    result_matrix<<Model.get_results();
    result_matrix2<< matrix<double>(prod(A,S));
    precisions<<Model.get_noise_precision();
  }
      
  ICR::EnsembleLearning::Builder<double>& Build = Model.get_builder();
  
  Build.set_cost_file(cost_file.string());
  std::cout<<"Running!"<<std::endl;
  bool converged = false;
  size_t count = 0;
  while(!converged && count<500) {
    converged = Build.run(convergence_criterium,max_iterations);

    std::ofstream sources(source_file.string().c_str());
    std::ofstream mixing_matrix(mixing_file.string().c_str());
    std::ofstream result_matrix(result_file.string().c_str());
    std::ofstream result_matrix2(result_file2.string().c_str());
    std::ofstream precisions(precisions_file.string().c_str());

    matrix<double> A(size1, assumed_sources);
    matrix<double> S(assumed_sources,size2);
    Model.get_normalised_means(A,S);
    sources<<S;
    mixing_matrix<< matrix<double>(trans(A));
    result_matrix<<Model.get_results();
    result_matrix2<< matrix<double>(prod(A,S));
    precisions<<Model.get_noise_precision();
    ++count;
  }
  }
}
