
//builder
#include "EnsembleLearning/exponential_model/Random.hpp"

//initialise
ICR::EnsembleLearning::rng* ICR::EnsembleLearning::Random::m_rng = 0;

ICR::EnsembleLearning::SingletonDestroyer<ICR::EnsembleLearning::rng> 
ICR::EnsembleLearning::Random::m_Destroyer(ICR::EnsembleLearning::Random::m_rng);

