#pragma once
//#include "rng.hpp"

#include <gsl/gsl_rng.h>     
#include <gsl/gsl_randist.h>
#include <ctime>

namespace ICR{
  namespace ICA{

    /** A Random number generator.
     *  A simple wrapper around various random number generators found in the GNU Scientific Library.  The random samples are used to initialise the Moments.
     */
  class rng{
  private:
    unsigned long int m_seed;
    gsl_rng* m_rng;
    rng(const rng& other); //non-copyable
    rng& operator= (const rng & other); //non-copyable
    
  public:
    /** Constructor.
     *  The seed is based on the current system time.
     *  @attention: This is not a good seed if multiple random number generators are to be constructed, as there is the strong possibility of many generators starting with the same seed.
     */
    rng()
      : m_seed(static_cast<unsigned int>(std::time(0))),
	m_rng(gsl_rng_alloc (gsl_rng_taus))
    {
      gsl_rng_set (m_rng,m_seed);
    }
    /** Constructor.
     *  @param seed The seed for the random number generator.
     */
    rng(unsigned long int seed)
      : m_seed(seed),
	m_rng(gsl_rng_alloc (gsl_rng_taus))
    {
       gsl_rng_set (m_rng,m_seed);
    }
    
    /** Obtain the seed for the current random number generator.
	@return The current seed.
     */
    unsigned long int 
    get_seed() const {return m_seed;};
    
    /** Destructor */
    ~rng(){gsl_rng_free (m_rng);};

    /** Generate a random number in the uniform distribution (between zero and one).
     *@return The random number.
     */
    double uniform() {return gsl_rng_uniform (m_rng);};
    
    /** Generate a random number in the uniform distribution.
     * @param a The lower bound on the random number.
     * @param b The upper bound on the random number	
     * @return The random number.
     */
    double uniform(double a, double b) {return gsl_ran_flat (m_rng, a, b);}
    
    /** Generate a random number in a Gaussian distrubution.
     * @param sigma The standard deviation of the Gaussian distribution.
     * @param mean The mean of the Gaussian distribution.
     * @return The random number.
     */
    double gaussian(const double sigma = 1, const double mean = 0){
      return gsl_ran_gaussian (m_rng, sigma) + mean;
    }
    
    /** Generate a random number in a Gaussian tail distrubution.
     *  This can be useful for generating random numbers for a Rectified Gaussian distribution.
     * @param sigma The standard deviation of the Gaussian distribution.
     * @param mean The mean of the Gaussian distribution.
     * @param min The miniumum of the random number.
     * @return The random number.
     */
    double gaussian_tail(const double sigma = 1, const double mean = 0, const double min = 0){
      return gsl_ran_gaussian_tail (m_rng, min-mean, sigma) + mean;
    } 

    /** Generate a random number in an exponential distrubution.
     * @param mean The mean of the distribution.
     * @return The random number.
     */
    double exponential(const double mean = 1){
      return gsl_ran_exponential (m_rng, mean );
    }
    /** Generate a random number in a Gamma distrubution.
     * @param shape The shape of the distribution.
     * @param scale The scale of the Gamma distribution.
     * @return The random number.
     */
    double gamma(const double shape =1, const double scale =1){
      return gsl_ran_gamma (m_rng, shape, scale);
    } 

    /** Generate a random number in a Dirichlet distrubution.
     * @param size The dimension of the Dirichlet distribution.
     * @param alpha The set of priors for the Dirichlet Distribution.
     * @param theta The set of random numbers.
     */
    void dirichlet(const size_t size, const double* alpha, double* theta){

      gsl_ran_dirichlet(m_rng,size, alpha, theta);
    } 
  };
  
    /** Destroy the Singleton elegently when the program ends.
     * @tparam The class to destroy on the programs exit.
     */
    template<class T>
    class SingletonDestroyer { 
    public: 
      /** A constructor.
       * @param s A pointer to delete on destruction.
       */
      SingletonDestroyer(T* s = 0); 
      /** A destructor.
       * The  pointer stored in the singleton will be deleted.
       */
      ~SingletonDestroyer(); 
      /** Set the reference stored in the singleton.
       * @param s A pointer  to delete on destruction.
       */
      void SetSingleton(T* s); 
    private: 
      T* m_singleton;
    };

    /** A Singleton used to store a random number generator throughout the programs lifetime.
     * The singleton is destroyed at the end of the program by the SingletonDestroyer.
     */
    class Random{
    public:
      /** Obtain an instance of the random number generator. 
       * @param seed The seed for the random number generator (only used on the first call)
       * @return A pointer to the rng.
       */
      static rng* Instance(unsigned int seed)
      {
	if (m_rng == 0)
	  m_rng = new rng(seed);
	m_Destroyer.SetSingleton(m_rng); 
	return m_rng;
      }
      
      /** Obtain an instance of the random number generator. 
       * @return A pointer to the rng.
       */
      static rng* Instance()
      {
	if (m_rng == 0)
	  m_rng = new rng();
	m_Destroyer.SetSingleton(m_rng); 
	return m_rng;
      }
      
      
      /** Restart the random number generator
       * @param seed The seed for the random number generator
       * @return A pointer to the rng.
       */
      static rng* Restart(unsigned int seed) 
      {
	if (m_rng!=0) {
	  delete m_rng;
	  m_rng = new rng(seed);
	  m_Destroyer.SetSingleton(m_rng); 
	}
	return m_rng;
      }; 
    private:
      Random(); 
      ~Random(){}

      friend class SingletonDestroyer<rng>; 
      static rng *m_rng;
      static SingletonDestroyer<rng> m_Destroyer;
    };


  }
}

//initialise
ICR::ICA::rng* ICR::ICA::Random::m_rng = 0;

ICR::ICA::SingletonDestroyer<ICR::ICA::rng> 
ICR::ICA::Random::m_Destroyer(ICR::ICA::Random::m_rng);




template<class T>
ICR::ICA::SingletonDestroyer<T>::SingletonDestroyer(T* s) 
{ 
  m_singleton = s; 
}
template<class T>
ICR::ICA::SingletonDestroyer<T>::~SingletonDestroyer() 
{
  if (m_singleton!=0) {
    delete m_singleton;
    m_singleton = 0;
  } 
}

template<class T>
void
ICR::ICA::SingletonDestroyer<T>::SetSingleton(T *s) 
{
  m_singleton = s;
}

    
