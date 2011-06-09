#include "ExampleData.hpp"
#include "BuildModel.hpp"

#include <boost/program_options.hpp>
#include "boost/filesystem.hpp"   // includes all needed Boost.Filesystem declarations

#include <iostream>
#include <fstream>
#include <string>

using namespace boost::numeric::ublas;
namespace po = boost::program_options;
namespace fs = boost::filesystem;


int
main  (int ac, char **av)
{
  //start the Random number generator with fixed seed so model reproducable
  ICR::EnsembleLearning::Random::Instance(10);
  
  bool   use_float = false;
  size_t assumed_sources;
  size_t mixing_components;
  bool   positive_source = false;
  bool   positive_mixing = false;
  bool   model_noise_offset = false;
  bool   example = false;
  double convergence_criterium;
  size_t max_iterations;
  
  std::string data_file;
  std::string output_directory;
  
  
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
    ("sources,s", po::value<size_t>(&assumed_sources)->default_value(5), 
     "The number of sources to infer (a guess)")
    ("mixing,m",  po::value<size_t>(&mixing_components)->default_value(5), 
     "The number of components to use in Mixture Models")
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


  if ( !fs::exists(output_directory) ) {
    std::cout << "Output directory '"<<output_directory<<"' does not exist!\n\n"
	      << visible  << "\n";
    return 1;
  }

  matrix<double> Data;
  if ( !fs::exists(data_file) )
    {
      if (example)
	{
	  size_t number_of_records = 5;
	  size_t number_of_sources = 3;
	  size_t samples_per_record = 100;
	  
	  //data examples
	  //Data = DataRecords(5,2,50,positive_source,positive_mixing);
	  
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
	//for(size_t i=0;i<4;++i){


	//getline( myfile, line );
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
  
  
  std::cout<<"sizeof data = "<<sizeof(Data)<<std::endl;


  std::cout<<"output direcotry = "<<output_directory<<std::endl;
  std::cout<<"input file = "<<data_file<<std::endl;


  
  fs::path source_file = fs::path(output_directory)/fs::path("InferedSources.txt");
  fs::path mixing_file = fs::path(output_directory)/fs::path("InferedMixing.txt");
  fs::path cost_file   = fs::path(output_directory)/fs::path("Cost.txt");
  fs::path result_file   = fs::path(output_directory)/fs::path("InferedResult.txt");
  fs::path result_file2   = fs::path(output_directory)/fs::path("InferedResult2.txt");
  
  

  
  if (use_float) 
    {
      std::cout<<"Building Model"<<std::endl;
      BuildModel<float> Model(Data, assumed_sources, mixing_components, 
			      positive_source,
			      positive_mixing,
			      model_noise_offset);
      ICR::EnsembleLearning::Builder<float> Build = Model.get_builder();
      Build.set_cost_file(cost_file.string());
      std::cout<<"Running!"<<std::endl;
      Build.run(convergence_criterium,max_iterations);
      std::ofstream sources(source_file.string().c_str());
      std::ofstream mixing_matrix(mixing_file.string().c_str());
      matrix<float> A(Data.size1(), assumed_sources);
      matrix<float> S(assumed_sources,Data.size2());
      Model.get_normalised_means(A,S);
      sources<<S;
      mixing_matrix<< matrix<float>(trans(A));
    }
  else
    {
      std::cout<<"Building Model"<<std::endl;
      BuildModel<double> Model(Data, assumed_sources, mixing_components, 
			       positive_source,
			       positive_mixing,
			       model_noise_offset);
      ICR::EnsembleLearning::Builder<double> Build = Model.get_builder();
      Build.set_cost_file(cost_file.string());
      std::cout<<"Running!"<<std::endl;
      bool converged = false;
      size_t count = 0;
      while(!converged && count<10) {
	converged = Build.run(convergence_criterium,max_iterations);


	std::ofstream sources(source_file.string().c_str());
	std::ofstream mixing_matrix(mixing_file.string().c_str());
	std::ofstream result_matrix(result_file.string().c_str());
	std::ofstream result_matrix2(result_file2.string().c_str());
	matrix<double> A(Data.size1(), assumed_sources);
	matrix<double> S(assumed_sources,Data.size2());
	Model.get_normalised_means(A,S);
	sources<<S;
	mixing_matrix<< matrix<double>(trans(A));
	result_matrix<<Model.get_results();
	result_matrix2<< matrix<double>(prod(A,S));
	++count;
      }
    }
}
