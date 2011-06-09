
#define BOOST_TEST_MODULE ICA_library_tests

#include <gsl/gsl_sf_psi.h> //for Digamma function gsl_sf_psi
#include <gsl/gsl_sf_gamma.h> //for gamma function

#include "EnsembleLearning/message/Coster.hpp"
#include "EnsembleLearning/message/Moments.hpp"
#include "EnsembleLearning/message/NaturalParameters.hpp"
#include "EnsembleLearning/exponential_model/Gaussian.hpp"
#include "EnsembleLearning/exponential_model/RectifiedGaussian.hpp"
#include "EnsembleLearning/exponential_model/Gamma.hpp"
#include "EnsembleLearning/exponential_model/Discrete.hpp"
#include "EnsembleLearning/exponential_model/Dirichlet.hpp"
#include "EnsembleLearning/message/Moments.hpp"
#include "EnsembleLearning/message/NaturalParameters.hpp"
#include "EnsembleLearning/node/factor/Factor.hpp"

#include "EnsembleLearning.hpp"
//#include "rng.hpp"
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <boost/shared_ptr.hpp>

//#include "vec.hpp"
//#include "mat.hpp"

#include <vector>
#include <iostream>
#include <fstream>
#include<boost/assign/list_of.hpp>
#include <boost/assign/std/vector.hpp>
#include <gsl/gsl_sf_erf.h>
using namespace ICR::EnsembleLearning;



#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

using namespace boost::numeric;

// using namespace boost::assign;
//____________________________________________________________________________//

//typedef ICR::maths::vector<boost::units::si::dimensionless> vec;
//typedef ICR::maths::mat mat;


class ScopedRedirect
{
public:
    ScopedRedirect(std::ostream & inOriginal, std::ostream & inRedirect) :
        mOriginal(inOriginal),
        mOldBuffer(inOriginal.rdbuf(inRedirect.rdbuf()))
    { }

    ~ScopedRedirect()
    {
        mOriginal.rdbuf(mOldBuffer);
    }    

private:
    ScopedRedirect(const ScopedRedirect&);
    ScopedRedirect& operator=(const ScopedRedirect&);

    std::ostream & mOriginal;
  std::streambuf * mOldBuffer;
};


/*****************************************************
 *****************************************************
 *****      Random      TEST       *******
 *****************************************************
 *****************************************************/


BOOST_AUTO_TEST_SUITE( Random_test )


BOOST_AUTO_TEST_CASE( Singleton_test  )
{
  Random::Restart(100);
  
  rng* random1 = Random::Instance(100);
  
  
  rng* random2 = Random::Restart(10);
  double d1 = random2->uniform();
  double d2 = random2->uniform(1,2);
  
  BOOST_CHECK(d1<1);
  BOOST_CHECK(d1>=0);
  BOOST_CHECK(d2<2);
  BOOST_CHECK(d2>=1);
  
  
}


BOOST_AUTO_TEST_SUITE_END()


/*****************************************************
 *****************************************************
 *****               COSTER TEST               *******
 *****************************************************
 *****************************************************/

BOOST_AUTO_TEST_SUITE( Coster_test )

BOOST_AUTO_TEST_CASE( constr_test  )
{
  Coster c;
  BOOST_CHECK_EQUAL(double(c), 0.0); 
}

BOOST_AUTO_TEST_CASE( assign_test  )
{
  Coster c = 4;
  BOOST_CHECK_EQUAL(double(c), 4.0);  
}

BOOST_AUTO_TEST_CASE( incr_test  )
{
  Coster c = 4;
  c += 4;
  BOOST_CHECK_EQUAL(double(c), 8.0);  
}

BOOST_AUTO_TEST_SUITE_END()


/*****************************************************
 *****************************************************
 *****               Moments TEST               *******
 *****************************************************
 *****************************************************/

BOOST_AUTO_TEST_SUITE( Moments_test )

BOOST_AUTO_TEST_CASE( constr_test  )
{
  Moments<double> M1;
  Moments<double> M2(2);
  std::vector<double> v = boost::assign::list_of(2.0)(1.5)(4.0);
  Moments<double> M3(v);
  Moments<double> M4(2,2.5);
  Moments<double> M5(M3);
  
  Moments<double> M6;
  M6 = M3;
  
  //also tests size
  BOOST_CHECK_EQUAL(M1.size(), (size_t) 0);  
  BOOST_CHECK_EQUAL(M2.size(), (size_t) 2);  
  BOOST_CHECK_EQUAL(M3.size(), (size_t) 3);  
  BOOST_CHECK_EQUAL(M4.size(), (size_t) 2);  
  BOOST_CHECK_EQUAL(M5.size(), (size_t) 3); 
  BOOST_CHECK_EQUAL(M5.size(), (size_t) 3);  
  
  //also tests operator[] const
  BOOST_CHECK_CLOSE(M2[0], 0.0, 0.0001);
  BOOST_CHECK_CLOSE(M2[1], 0.0, 0.0001);
  BOOST_CHECK_CLOSE(M3[0], 2.0, 0.0001);
  BOOST_CHECK_CLOSE(M3[1], 1.5, 0.0001);
  BOOST_CHECK_CLOSE(M3[2], 4.0, 0.0001);
  BOOST_CHECK_CLOSE(M4[0], 2.0, 0.0001);
  BOOST_CHECK_CLOSE(M4[1], 2.5, 0.0001);
  BOOST_CHECK_CLOSE(M5[0], 2.0, 0.0001);
  BOOST_CHECK_CLOSE(M5[1], 1.5, 0.0001);
  BOOST_CHECK_CLOSE(M5[2], 4.0, 0.0001);
  BOOST_CHECK_CLOSE(M6[0], 2.0, 0.0001);
  BOOST_CHECK_CLOSE(M6[1], 1.5, 0.0001);
  BOOST_CHECK_CLOSE(M6[2], 4.0, 0.0001);
}

BOOST_AUTO_TEST_CASE( op_square_test  )
{
  std::vector<double> v = boost::assign::list_of(2.0)(1.5)(4.0);
  Moments<double> M1(v);
  BOOST_CHECK_CLOSE(M1[0], 2.0, 0.0001);
  BOOST_CHECK_CLOSE(M1[1], 1.5, 0.0001);
  BOOST_CHECK_CLOSE(M1[2], 4.0, 0.0001);
  
  M1[0] = 3.0;
  M1[1] = 3.5;
  M1[2] = 3.7;

  BOOST_CHECK_CLOSE(M1[0], 3.0, 0.0001);
  BOOST_CHECK_CLOSE(M1[1], 3.5, 0.0001);
  BOOST_CHECK_CLOSE(M1[2], 3.7, 0.0001);
  
  //check not pushed to the front or anyting silly
  BOOST_CHECK_EQUAL(M1.size(),  (size_t) 3);  
}

BOOST_AUTO_TEST_CASE( maths_op_test  )
{
  std::vector<double> v1 = boost::assign::list_of(2.0)(1.5)(4.0);
  std::vector<double> v2 = boost::assign::list_of(2.0)(-3.0)(4.0);
  Moments<double> M1(v1);
  BOOST_CHECK_CLOSE(M1[0], 2.0, 0.0001);
  BOOST_CHECK_CLOSE(M1[1], 1.5, 0.0001);
  BOOST_CHECK_CLOSE(M1[2], 4.0, 0.0001);

  Moments<double> M2(v2);

  //Plus
  M1+=M2;
  BOOST_CHECK_CLOSE(M1[0], 4.0, 0.0001);
  BOOST_CHECK_CLOSE(M1[1], -1.5, 0.0001);
  BOOST_CHECK_CLOSE(M1[2], 8.0, 0.0001);
  //check not pushed to the front or anyting silly
  BOOST_CHECK_EQUAL(M1.size(), (size_t) 3);  

  //Times double
  M1 = Moments<double>(v1);  //reset
  M1*=-2.0;
  BOOST_CHECK_CLOSE(M1[0], -4.0, 0.0001);
  BOOST_CHECK_CLOSE(M1[1], -3.0, 0.0001);
  BOOST_CHECK_CLOSE(M1[2], -8.0, 0.0001);
  BOOST_CHECK_EQUAL(M1.size(), (size_t) 3);  
  
  //Times Momenst
  M1 = Moments<double>(v1);  //reset
  Moments<double> M3(v2);
  M1*= M3;
  BOOST_CHECK_CLOSE(M1[0], 4.0, 0.0001);
  BOOST_CHECK_CLOSE(M1[1], -4.5, 0.0001);
  BOOST_CHECK_CLOSE(M1[2], 16, 0.0001);
  BOOST_CHECK_EQUAL(M1.size(), (size_t) 3);  

  //CHECK FREE FUNCTIONS;
  
  //Plus
  Moments<double> M4 = Moments<double>(v1)+Moments<double>(v2);
  BOOST_CHECK_CLOSE(M4[0], 4.0, 0.0001);
  BOOST_CHECK_CLOSE(M4[1], -1.5, 0.0001);
  BOOST_CHECK_CLOSE(M4[2], 8.0, 0.0001);
  //check not pushed to the front or anyting silly
  BOOST_CHECK_EQUAL(M4.size(), (size_t) 3);

  //Times double
  Moments<double> M5 =  Moments<double>(v1)*-2.0;  //reset
  BOOST_CHECK_CLOSE(M5[0], -4.0, 0.0001);
  BOOST_CHECK_CLOSE(M5[1], -3.0, 0.0001);
  BOOST_CHECK_CLOSE(M5[2], -8.0, 0.0001);
  BOOST_CHECK_EQUAL(M5.size(), (size_t) 3); 
  
  Moments<double> M6 = -2.0* Moments<double>(v1);  //reset
  BOOST_CHECK_CLOSE(M6[0], -4.0, 0.0001);
  BOOST_CHECK_CLOSE(M6[1], -3.0, 0.0001);
  BOOST_CHECK_CLOSE(M6[2], -8.0, 0.0001);
  BOOST_CHECK_EQUAL(M6.size(), (size_t) 3);  
  
  //Times Momenst
  // Moments<double> M6 = Moments<double>(v1)*M3;  //reset
  // BOOST_CHECK_CLOSE(M6[0], 4.0, 0.0001);
  // BOOST_CHECK_CLOSE(M6[1], -4.5, 0.0001);
  // BOOST_CHECK_CLOSE(M6[2], 16, 0.0001);
  // BOOST_CHECK_EQUAL(M6.size(), (size_t) 3);  

}

BOOST_AUTO_TEST_CASE( iterator_test  )
{
  typedef Moments<double>::iterator iterator;
  typedef Moments<double>::const_iterator const_iterator;
  
  std::vector<double> v1 = boost::assign::list_of(2.0)(1.5)(4.0);
  Moments<double> M1(v1);
  const Moments<double> cM1(v1);
  

  //check end
  BOOST_CHECK_EQUAL(M1.end() - M1.begin(), 3);

  iterator it1 = M1.begin();
  BOOST_CHECK_CLOSE(*it1, 2.0, 0.0001);
  it1++;
  BOOST_CHECK_CLOSE(*it1, 1.5, 0.0001);
  BOOST_CHECK_CLOSE(*(++it1), 4.0, 0.0001);
  BOOST_CHECK_EQUAL(it1 - M1.begin(), 2);
  it1--;
  BOOST_CHECK_CLOSE(*it1, 1.5, 0.0001);
  BOOST_CHECK_CLOSE(*(--it1), 2.0, 0.0001);

  //check end
  BOOST_CHECK_EQUAL(cM1.end() - cM1.begin(), 3);

  const_iterator cit1 = cM1.begin();
  BOOST_CHECK_CLOSE(*cit1, 2.0, 0.0001);
  cit1++;
  BOOST_CHECK_CLOSE(*cit1, 1.5, 0.0001);
  BOOST_CHECK_CLOSE(*(++cit1), 4.0, 0.0001);
  BOOST_CHECK_EQUAL(cit1 - cM1.begin(), 2);
  cit1--;
  BOOST_CHECK_CLOSE(*cit1, 1.5, 0.0001);
  BOOST_CHECK_CLOSE(*(--cit1), 2.0, 0.0001);

  Moments<double>::const_iterator it2 = M1.begin();
  BOOST_CHECK_CLOSE(*it2, 2.0, 0.0001);
  it2++;
  BOOST_CHECK_CLOSE(*(it2), 1.5, 0.0001);
  BOOST_CHECK_CLOSE(*(++it2), 4.0, 0.0001);
  BOOST_CHECK_EQUAL(it2 - M1.begin(), 2);
  BOOST_CHECK_EQUAL(it2 - cM1.begin(), 2);
  BOOST_CHECK_EQUAL(it2 - it1, 2);
  it2--;
  BOOST_CHECK_CLOSE(*(it2), 1.5, 0.0001);
  BOOST_CHECK_CLOSE(*(--it2), 2.0, 0.0001);

  
  it1 = M1.end();
  it2 = cM1.begin();
  it1--;
  //it1 = cit1; //test const conversion
  it2 = it1; //test const conversion
  BOOST_CHECK_CLOSE(*it1, 4.0, 0.0001);
  BOOST_CHECK_CLOSE(*it2, 4.0, 0.0001);
  
  *it1 = 3.3;
  BOOST_CHECK_CLOSE(M1[2], 3.3, 0.0001);
  

  it1 = M1.end();
  //Equality
  BOOST_CHECK(M1.end() == it1);
  BOOST_CHECK(M1.begin() != it1);
  
  
}


BOOST_AUTO_TEST_SUITE_END()


/*****************************************************
 *****************************************************
 *****      NaturalParameters     TEST         *******
 *****************************************************
 *****************************************************/

BOOST_AUTO_TEST_SUITE( NaturalParameters_test )


BOOST_AUTO_TEST_CASE( constr_test  )
{
  NaturalParameters<double> NP1;
  NaturalParameters<double> NP2(2);
  std::vector<double> v = boost::assign::list_of(2.0)(1.5)(4.0);
  NaturalParameters<double> NP3(v);
  NaturalParameters<double> NP4(2,2.5);
  NaturalParameters<double> NP5(NP3);
  
  NaturalParameters<double> NP6;
  NP6 = NP3;
  
  //also tests size
  BOOST_CHECK_EQUAL(NP1.size(), (size_t) 0);  
  BOOST_CHECK_EQUAL(NP2.size(), (size_t) 2);  
  BOOST_CHECK_EQUAL(NP3.size(), (size_t) 3);  
  BOOST_CHECK_EQUAL(NP4.size(), (size_t) 2);  
  BOOST_CHECK_EQUAL(NP5.size(), (size_t) 3); 
  BOOST_CHECK_EQUAL(NP5.size(), (size_t) 3);  
  
  //also tests operator[] const
  BOOST_CHECK_CLOSE(NP2[0], 0.0, 0.0001);
  BOOST_CHECK_CLOSE(NP2[1], 0.0, 0.0001);
  BOOST_CHECK_CLOSE(NP3[0], 2.0, 0.0001);
  BOOST_CHECK_CLOSE(NP3[1], 1.5, 0.0001);
  BOOST_CHECK_CLOSE(NP3[2], 4.0, 0.0001);
  BOOST_CHECK_CLOSE(NP4[0], 2.0, 0.0001);
  BOOST_CHECK_CLOSE(NP4[1], 2.5, 0.0001);
  BOOST_CHECK_CLOSE(NP5[0], 2.0, 0.0001);
  BOOST_CHECK_CLOSE(NP5[1], 1.5, 0.0001);
  BOOST_CHECK_CLOSE(NP5[2], 4.0, 0.0001);
  BOOST_CHECK_CLOSE(NP6[0], 2.0, 0.0001);
  BOOST_CHECK_CLOSE(NP6[1], 1.5, 0.0001);
  BOOST_CHECK_CLOSE(NP6[2], 4.0, 0.0001);
}

BOOST_AUTO_TEST_CASE( op_square_test  )
{
  std::vector<double> v = boost::assign::list_of(2.0)(1.5)(4.0);
  NaturalParameters<double> NP1(v);
  BOOST_CHECK_CLOSE(NP1[0], 2.0, 0.0001);
  BOOST_CHECK_CLOSE(NP1[1], 1.5, 0.0001);
  BOOST_CHECK_CLOSE(NP1[2], 4.0, 0.0001);
  
  NP1[0] = 3.0;
  NP1[1] = 3.5;
  NP1[2] = 3.7;

  BOOST_CHECK_CLOSE(NP1[0], 3.0, 0.0001);
  BOOST_CHECK_CLOSE(NP1[1], 3.5, 0.0001);
  BOOST_CHECK_CLOSE(NP1[2], 3.7, 0.0001);
  
  //check not pushed to the front or anyting silly
  BOOST_CHECK_EQUAL(NP1.size(),  (size_t) 3);  
}

BOOST_AUTO_TEST_CASE( maths_op_test  )
{
  std::vector<double> v1 = boost::assign::list_of(2.0)(1.5)(4.0);
  std::vector<double> v2 = boost::assign::list_of(2.0)(-3.0)(4.0);
  NaturalParameters<double> NP1(v1);
  BOOST_CHECK_CLOSE(NP1[0], 2.0, 0.0001);
  BOOST_CHECK_CLOSE(NP1[1], 1.5, 0.0001);
  BOOST_CHECK_CLOSE(NP1[2], 4.0, 0.0001);

  NaturalParameters<double> NP2(v2);

  //Plus
  NP1+=NP2;
  BOOST_CHECK_CLOSE(NP1[0], 4.0, 0.0001);
  BOOST_CHECK_CLOSE(NP1[1], -1.5, 0.0001);
  BOOST_CHECK_CLOSE(NP1[2], 8.0, 0.0001);
  //check not pushed to the front or anyting silly
  BOOST_CHECK_EQUAL(NP1.size(), (size_t) 3);  

  //Minus
  NP1-=NP2;
  BOOST_CHECK_CLOSE(NP1[0], 2.0, 0.0001);
  BOOST_CHECK_CLOSE(NP1[1], 1.5, 0.0001);
  BOOST_CHECK_CLOSE(NP1[2], 4.0, 0.0001);
  //check not pushed to the front or anyting silly
  BOOST_CHECK_EQUAL(NP1.size(), (size_t) 3);  

  //Times double
  NP1 = NaturalParameters<double>(v1);  //reset
  NP1*=-2.0;
  BOOST_CHECK_CLOSE(NP1[0], -4.0, 0.0001);
  BOOST_CHECK_CLOSE(NP1[1], -3.0, 0.0001);
  BOOST_CHECK_CLOSE(NP1[2], -8.0, 0.0001);
  BOOST_CHECK_EQUAL(NP1.size(), (size_t) 3);  
  
  
  //FREE FUNCTIONS TEST
  //plus
  NaturalParameters<double> NP3 
    = NaturalParameters<double>(v1)+NaturalParameters<double>(v2);
  BOOST_CHECK_CLOSE(NP3[0], 4.0, 0.0001);
  BOOST_CHECK_CLOSE(NP3[1], -1.5, 0.0001);
  BOOST_CHECK_CLOSE(NP3[2], 8.0, 0.0001);
  //check not pushed to the front or anyting silly
  BOOST_CHECK_EQUAL(NP3.size(), (size_t) 3); 
  
  //Minus
  NaturalParameters<double> NP4 = NP3-NP2;
  BOOST_CHECK_CLOSE(NP4[0], 2.0, 0.0001);
  BOOST_CHECK_CLOSE(NP4[1], 1.5, 0.0001);
  BOOST_CHECK_CLOSE(NP4[2], 4.0, 0.0001);
  //check not pushed to the front or anyting silly
  BOOST_CHECK_EQUAL(NP4.size(), (size_t) 3);  

   //times double
  NaturalParameters<double> NP5 =  NaturalParameters<double>(v1)*-2.0;
  BOOST_CHECK_CLOSE(NP5[0], -4.0, 0.0001);
  BOOST_CHECK_CLOSE(NP5[1], -3.0, 0.0001);
  BOOST_CHECK_CLOSE(NP5[2], -8.0, 0.0001);
  BOOST_CHECK_EQUAL(NP5.size(), (size_t) 3);  
  
  //INNER PRODUCT
  NaturalParameters<double> NPIP = NaturalParameters<double>(v1);
  Moments<double> MIP(v2);
  BOOST_CHECK_CLOSE(NPIP*MIP,4.0-4.5+ 16.0 , 0.0001);
  BOOST_CHECK_CLOSE(MIP*NPIP,4.0-4.5+ 16.0 , 0.0001);
   
  

}

BOOST_AUTO_TEST_CASE( iterator_test  )
{
  typedef NaturalParameters<double>::iterator iterator;
  typedef NaturalParameters<double>::const_iterator const_iterator;
  
  std::vector<double> v1 = boost::assign::list_of(2.0)(1.5)(4.0);
  NaturalParameters<double> M1(v1);
  const NaturalParameters<double> cM1(v1);
  

  //check end
  BOOST_CHECK_EQUAL(M1.end() - M1.begin(), 3);

  iterator it1 = M1.begin();
  BOOST_CHECK_CLOSE(*it1, 2.0, 0.0001);
  it1++;
  BOOST_CHECK_CLOSE(*it1, 1.5, 0.0001);
  BOOST_CHECK_CLOSE(*(++it1), 4.0, 0.0001);
  BOOST_CHECK_EQUAL(it1 - M1.begin(), 2);
  it1--;
  BOOST_CHECK_CLOSE(*it1, 1.5, 0.0001);
  BOOST_CHECK_CLOSE(*(--it1), 2.0, 0.0001);

  //check end
  BOOST_CHECK_EQUAL(cM1.end() - cM1.begin(), 3);

  const_iterator cit1 = cM1.begin();
  BOOST_CHECK_CLOSE(*cit1, 2.0, 0.0001);
  cit1++;
  BOOST_CHECK_CLOSE(*cit1, 1.5, 0.0001);
  BOOST_CHECK_CLOSE(*(++cit1), 4.0, 0.0001);
  BOOST_CHECK_EQUAL(cit1 - cM1.begin(), 2);
  cit1--;
  BOOST_CHECK_CLOSE(*cit1, 1.5, 0.0001);
  BOOST_CHECK_CLOSE(*(--cit1), 2.0, 0.0001);

  NaturalParameters<double>::const_iterator it2 = M1.begin();
  BOOST_CHECK_CLOSE(*it2, 2.0, 0.0001);
  it2++;
  BOOST_CHECK_CLOSE(*(it2), 1.5, 0.0001);
  BOOST_CHECK_CLOSE(*(++it2), 4.0, 0.0001);
  BOOST_CHECK_EQUAL(it2 - M1.begin(), 2);
  BOOST_CHECK_EQUAL(it2 - it1, 2);
  it2--;
  BOOST_CHECK_CLOSE(*(it2), 1.5, 0.0001);
  BOOST_CHECK_CLOSE(*(--it2), 2.0, 0.0001);

  
  it1 = M1.end();
  it2 = cM1.begin();
  it1--;
  //it1 = cit1; //test const conversion
  it2 = it1; //test const conversion
  BOOST_CHECK_CLOSE(*it1, 4.0, 0.0001);
  BOOST_CHECK_CLOSE(*it2, 4.0, 0.0001);
  
  *it1 = 3.3;
  BOOST_CHECK_CLOSE(M1[2], 3.3, 0.0001);
  
  
}




BOOST_AUTO_TEST_SUITE_END()




/*****************************************************
 *****************************************************
 *****      Exponential Models      TEST       *******
 *****************************************************
 *****************************************************/

BOOST_AUTO_TEST_SUITE( ExpModels_test )
  
BOOST_AUTO_TEST_CASE( Gaussian_test  )
{
  //Initialise
  Moments<double> Mean(2,5);
  Moments<double> Precision(3,6);
  Moments<double> Data(4,17);
  NaturalParameters<double> SumNP(6,-1.5);

  //Collect
  NaturalParameters<double> NPMean
    = Gaussian<double>::CalcNP2Parent1(Precision,Data);
  NaturalParameters<double> NPPrec
    = Gaussian<double>::CalcNP2Parent2(Mean,Data);
  NaturalParameters<double> NPData
    = Gaussian<double>::CalcNP2Data(Mean,Precision);
  

  Moments<double> Update =  Gaussian<double>::CalcMoments(SumNP);

  double LogNorm1 = Gaussian<double>::CalcLogNorm(Mean,Precision);
  double LogNorm2 = Gaussian<double>::CalcLogNorm(SumNP);
  double AvLog    = Gaussian<double>::CalcAvLog(Mean,Precision,Data);
  
  

  //Check
  BOOST_CHECK_CLOSE(NPMean[0], 12.0 , 0.0001);  //3*4
  BOOST_CHECK_CLOSE(NPMean[1], -1.5 , 0.0001);  //-3/2
  BOOST_CHECK_CLOSE(NPPrec[0], -3.0 , 0.0001);  // - 0.5*(17-2*4*2 +5)
  BOOST_CHECK_CLOSE(NPPrec[1], 0.5 , 0.0001);  // 0.5
  BOOST_CHECK_CLOSE(NPData[0], 6.0 , 0.0001);  //2*3
  BOOST_CHECK_CLOSE(NPData[1], -1.5 , 0.0001);  //-3/2

  BOOST_CHECK_CLOSE(Update[0],2.0 , 0.0001);  
  BOOST_CHECK_CLOSE(Update[1],4.0+1.0/3.0 , 0.0001);  

  //check mean & precision
  std::vector<double> mean      =  Gaussian<double>::CalcMean(SumNP);
  std::vector<double> precision =  Gaussian<double>::CalcPrecision(SumNP);
  BOOST_CHECK_CLOSE(mean[0],2.0 , 0.0001); 
  BOOST_CHECK_CLOSE(precision[0],3.0  , 0.0001); 
  

  //calc 1000 samples and check mean
  std::vector< Moments<double> > sample(5000);
  double av = 0;
  for(size_t i=0;i<sample.size();++i){
    sample[i] = Gaussian<double>::CalcSample(Mean,Precision);
    av+=sample[i][0];
  }
  av*=(1.0/sample.size());
  BOOST_CHECK_CLOSE(av,  2.0, 5); 
  

  BOOST_CHECK_CLOSE(LogNorm1, 0.5*(std::log(3.) - 3*5 - std::log(2*M_PI)) , 0.0001);  
  BOOST_CHECK_CLOSE(LogNorm2, 0.5*(std::log(3.) - 3*4 - std::log(2*M_PI)), 0.0001);  
  
  BOOST_CHECK_CLOSE(AvLog, NPData[0]*4 + NPData[1]*17 + LogNorm1, 0.0001);  
} 

double argm(double mean, double prec)
{
  return -mean*std::sqrt(prec/2.0);
}
double erfcx(double x)
{
  return std::exp(x*x)* gsl_sf_erfc(x);
  
}

double
mean_0(double mean, double prec)
{
  
  return mean + std::sqrt(2.0/(M_PI*prec))/erfcx(argm(mean,prec));

}

double
mean_1(double mean, double prec)
{  
  return mean*mean + 1.0/prec + std::sqrt(2.0/(2.0*M_PI*prec))*mean/erfcx(argm(mean,prec));
}

BOOST_AUTO_TEST_CASE( RectifiedGaussian_test  )
{
  //Initialise
  Moments<double> Parent1(2,5);
  Moments<double> Parent2(3,6);
  Moments<double> Data(4,17);
  NaturalParameters<double> SumNP(6,-1.5);
  //test large values (numerical errors can occur arg of +-30
  NaturalParameters<double> SumNP2(60,-1);
  NaturalParameters<double> SumNP3(-40,-1);
  //test large values (numerical errors can occur arg of +-15
  NaturalParameters<double> SumNP4(30,-1);
  NaturalParameters<double> SumNP5(-30,-1);


  //Collect
  NaturalParameters<double> NPParent1
    = RectifiedGaussian<double>::CalcNP2Parent1(Parent2,Data);
  NaturalParameters<double> NPPrec
    = RectifiedGaussian<double>::CalcNP2Parent2(Parent1,Data);
  NaturalParameters<double> NPData
    = RectifiedGaussian<double>::CalcNP2Data(Parent1,Parent2);
  

  Moments<double> Update =  RectifiedGaussian<double>::CalcMoments(SumNP);
  Moments<double> Update1 =  RectifiedGaussian<double>::CalcMoments(SumNP2);
  Moments<double> Update2 =  RectifiedGaussian<double>::CalcMoments(SumNP3);
  Moments<double> Update3 =  RectifiedGaussian<double>::CalcMoments(SumNP4);
  Moments<double> Update4 =  RectifiedGaussian<double>::CalcMoments(SumNP5);

  double LogNorm1 = RectifiedGaussian<double>::CalcLogNorm(Parent1,Parent2);
  double LogNorm2 = RectifiedGaussian<double>::CalcLogNorm(SumNP);
  double AvLog    = RectifiedGaussian<double>::CalcAvLog(Parent1,Parent2,Data);
  

  //Check
  BOOST_CHECK_CLOSE(NPParent1[0], 12.0 , 0.0001);  //3*4
  BOOST_CHECK_CLOSE(NPParent1[1], -1.5 , 0.0001);  //-3/2
  BOOST_CHECK_CLOSE(NPPrec[0], -3.0 , 0.0001);  // - 0.5*(17-2*4*2 +5)
  BOOST_CHECK_CLOSE(NPPrec[1], 0.5 , 0.0001);  // 0.5
  BOOST_CHECK_CLOSE(NPData[0], 6.0 , 0.0001);  //2*3
  BOOST_CHECK_CLOSE(NPData[1], -1.5 , 0.0001);  //-3/2

  BOOST_CHECK_CLOSE(Update[0],mean_0(2,3), 0.0001);  
  BOOST_CHECK_CLOSE(Update[1],mean_1(2,3) , 0.0001);  
  BOOST_CHECK_CLOSE(Update1[0],mean_0(30,2), 0.0001);  
  BOOST_CHECK_CLOSE(Update1[1],mean_1(30,2) , 0.0001);  
  BOOST_CHECK_CLOSE(Update2[0],mean_0(-20,2), 0.1);  
  BOOST_CHECK_CLOSE(Update2[1],mean_1(-20,2) , 0.1);  
  BOOST_CHECK_CLOSE(Update3[0],mean_0(15,2), 0.1);  //approximations to within 0.1%
  BOOST_CHECK_CLOSE(Update3[1],mean_1(15,2) , 0.1);  //approximations to within 0.1%
  BOOST_CHECK_CLOSE(Update4[0],mean_0(-15,2), 0.1);   //approximations to within 0.1%
  BOOST_CHECK_CLOSE(Update4[1],mean_1(-15,2) , 0.01);  //approximations to within 0.01%



  //check mean & precision
  std::vector<double> mean      =  RectifiedGaussian<double>::CalcMean(SumNP);
  std::vector<double> precision =  RectifiedGaussian<double>::CalcPrecision(SumNP);
  BOOST_CHECK_CLOSE(mean[0],2.0 , 0.0001); 
  BOOST_CHECK_CLOSE(precision[0],3.0  , 0.0001); 
  

  //calc 1000 samples and check mean
  std::vector< Moments<double> > sample(5000);
  double av = 0;
  for(size_t i=0;i<sample.size();++i){
    sample[i] = RectifiedGaussian<double>::CalcSample(Parent1,Parent2);
    av+=sample[i][0];
  }
  av*=(1.0/sample.size());
  BOOST_CHECK_CLOSE(av,  2.0, 5); 
  

  BOOST_CHECK_CLOSE(LogNorm1, 0.5*(std::log(2.*3) - 3*5 - std::log(M_PI)) - gsl_sf_log_erfc(-2.0*std::sqrt(3.0/2.0)) , 0.0001);  
  BOOST_CHECK_CLOSE(LogNorm2,  0.5*(std::log(2.*3) - 3*4 - std::log(M_PI)) - gsl_sf_log_erfc(-2.0*std::sqrt(3.0/2.0)), 0.0001);  
  
  BOOST_CHECK_CLOSE(AvLog, NPData[0]*4 + NPData[1]*17 + LogNorm1, 0.0001);  
} 

BOOST_AUTO_TEST_CASE( Gamma_test  )
{
  //Initialise
  Moments<double> Shape(2,5);
  Moments<double> IScale(3,0);
  Moments<double> Data(4,7);
  NaturalParameters<double> SumNP(-3,1);

  //Collect
  NaturalParameters<double> NPIScale
    = Gamma<double>::CalcNP2Parent2(Shape,Data);
  NaturalParameters<double> NPData
    = Gamma<double>::CalcNP2Data(Shape,IScale);
  

  Moments<double> Update =  Gamma<double>::CalcMoments(SumNP);

  double LogNorm1 = Gamma<double>::CalcLogNorm(Shape,IScale);
  double LogNorm2 = Gamma<double>::CalcLogNorm(SumNP);
  double AvLog    = Gamma<double>::CalcAvLog(Shape,IScale,Data);
  

  //Check
  BOOST_CHECK_CLOSE(NPIScale[0], -4.0 , 0.0001);  //-4
  BOOST_CHECK_CLOSE(NPIScale[1], 2.0 , 0.0001);  //3-1
  BOOST_CHECK_CLOSE(NPData[0], -3.0 , 0.0001);  //-3
  BOOST_CHECK_CLOSE(NPData[1], 1.0 , 0.0001);  //2-1

  BOOST_CHECK_CLOSE(Update[0], 2.0/3.0, 0.0001);  
  BOOST_CHECK_CLOSE(Update[1], gsl_sf_psi(2.0) - std::log(3.0), 0.0001);  
 
  //check mean & precision
  std::vector<double> mean      =  Gamma<double>::CalcMean(SumNP);
  std::vector<double> precision =  Gamma<double>::CalcPrecision(SumNP);
  BOOST_CHECK_CLOSE(mean[0],2.0/3.0 , 0.0001); 
  BOOST_CHECK_CLOSE(precision[0], 9.0/2.0 , 0.0001); 
  

  //calc 1000 samples and check mean
  std::vector< Moments<double> > sample(5000);
  double av = 0;
  for(size_t i=0;i<sample.size();++i){
    sample[i] = Gamma<double>::CalcSample(Shape,IScale);
    av+=sample[i][0];
  }
  av*=(1.0/sample.size());
  BOOST_CHECK_CLOSE(av,  Shape[0]/IScale[0], 5); 
  

  BOOST_CHECK_CLOSE(LogNorm1, 2.0*std::log(3.0) - gsl_sf_lngamma(2.0) , 0.0001);  
  BOOST_CHECK_CLOSE(LogNorm2, LogNorm1, 0.0001);  
  
  BOOST_CHECK_CLOSE(AvLog, NPData[0]*4.0 + NPData[1]*7.0 + LogNorm1, 0.0001);  
} 

BOOST_AUTO_TEST_CASE( Dirichlet_test  )
{
  

  //Initialise
  std::vector<double> u(3,1.1);
  Moments<double> Us(u);
  NaturalParameters<double> SumNP(std::vector<double>(3,0.1));

  //Collect
  NaturalParameters<double> NPData
    = Dirichlet<double>::CalcNP2Data(Us);
  

  Moments<double> Update =  Dirichlet<double>::CalcMoments(SumNP);

  double LogNorm1 = Dirichlet<double>::CalcLogNorm(Us);
  double LogNorm2 = Dirichlet<double>::CalcLogNorm(SumNP);
  

  //Check
  BOOST_CHECK_CLOSE(NPData[0], 0.1 , 0.0001);  //1.1-1
  BOOST_CHECK_CLOSE(NPData[1], 0.1 , 0.0001);  //1.1-1
  BOOST_CHECK_CLOSE(NPData[2], 0.1 , 0.0001);  //1.1-1

  BOOST_CHECK_CLOSE(Update[0], gsl_sf_psi(1.1)-gsl_sf_psi(3.3), 0.0001); 
  BOOST_CHECK_CLOSE(Update[1], gsl_sf_psi(1.1)-gsl_sf_psi(3.3), 0.0001); 
  BOOST_CHECK_CLOSE(Update[2], gsl_sf_psi(1.1)-gsl_sf_psi(3.3), 0.0001); 
  
  //check mean
  std::vector<double> mean =  Dirichlet<double>::CalcMean(NPData);
  BOOST_CHECK_CLOSE(mean[0], 1.0/3.0, 0.0001); 
  BOOST_CHECK_CLOSE(mean[1], 1.0/3.0, 0.0001); 
  BOOST_CHECK_CLOSE(mean[2], 1.0/3.0, 0.0001); 
  
  //check precision
  std::vector<double> precision =  Dirichlet<double>::CalcPrecision(NPData);
  BOOST_CHECK_CLOSE(precision[0], 1.0/((2.0/9.0)/(4.3)), 0.0001); 
  BOOST_CHECK_CLOSE(precision[1], 1.0/((2.0/9.0)/(4.3)), 0.0001); 
  BOOST_CHECK_CLOSE(precision[2], 1.0/((2.0/9.0)/(4.3)), 0.0001); 
  

  //calc 1000 samples and check mean
  std::vector< Moments<double> > sample(5000);
  std::vector<double> av(3);
  for(size_t i=0;i<sample.size();++i){
    sample[i] = Dirichlet<double>::CalcSample(Us);
    for(size_t j=0;j<3;++j){
      av[j]+=std::exp(sample[i][j]);
      // std::cout<<"av["<<j<<"] = "<<av[j]<<std::endl;

    }
  }
  for(size_t j=0;j<3;++j){
    av[j]*=(1.0/sample.size());
  }
  
  BOOST_CHECK_CLOSE((av[0]), (1.0/3.0), 5.0); 
  BOOST_CHECK_CLOSE((av[1]), (1.0/3.0), 5.0); 
  BOOST_CHECK_CLOSE((av[2]), (1.0/3.0), 5.0); 

  

  //BOOST_CHECK_CLOSE(3*std::exp(Update[0]), 1.0, 0.0001);

  BOOST_CHECK_CLOSE(LogNorm1, gsl_sf_lngamma(3.3) -3.0*gsl_sf_lngamma(1.1) , 0.0001);  
  BOOST_CHECK_CLOSE(LogNorm2, LogNorm1, 0.0001);  
  
} 


BOOST_AUTO_TEST_CASE( Discrete_test  )
{
  //Initialise
  std::vector<double> p(3);
  std::vector<double> logp(3);
  p[0] = 0;
  p[1] = 0.4;
  p[2] = 0.6;
  logp[0] = -500;
  logp[1] = -0.9;
  logp[2] = -0.5;

  Moments<double> dirichlet(logp);
  Moments<double> discrete(p);
  NaturalParameters<double> SumNP(logp);

  //Collect
  NaturalParameters<double> NPData
    = Discrete<double>::CalcNP2Data(dirichlet);
  NaturalParameters<double> NPPrior
    = Discrete<double>::CalcNP2Prior(discrete);
  

  Moments<double> Update =  Discrete<double>::CalcMoments(SumNP);

  double LogNorm1 = Discrete<double>::CalcLogNorm(dirichlet);
  double LogNorm2 = Discrete<double>::CalcLogNorm(SumNP);
  
  //Check
   BOOST_CHECK_CLOSE(NPData[0],  -500 , 0.0001);  //1.1-1
   BOOST_CHECK_CLOSE(NPData[1],   -0.9, 0.0001);  //1.1-1
   BOOST_CHECK_CLOSE(NPData[2],  -0.5 , 0.0001);  //1.1-1
   BOOST_CHECK_CLOSE(NPPrior[0], 0.0 , 0.0001);  //1.1-1
   BOOST_CHECK_CLOSE(NPPrior[1], 0.4 , 0.0001);  //1.1-1
   BOOST_CHECK_CLOSE(NPPrior[2], 0.6 , 0.0001);  //1.1-1

   BOOST_CHECK_CLOSE(Update[0], std::exp(-500+LogNorm1), 0.0001); 
   BOOST_CHECK_CLOSE(Update[1], std::exp(-0.9+LogNorm1), 0.0001); 
   BOOST_CHECK_CLOSE(Update[2], std::exp(-0.5+LogNorm1), 0.0001); 

  //check mean
  std::vector<double> mean =  Discrete<double>::CalcMean(NPData);
  BOOST_CHECK_CLOSE(mean[0], 1.6, 1); 
  
  //check precision
  std::vector<double> precision =  Discrete<double>::CalcPrecision(NPData);
  BOOST_CHECK_CLOSE(precision[0], 1.0/(0.4*0.36+0.6*.16), 1); 
  
  //check calcsample
  Moments<double> sample =  Discrete<double>::CalcSample(dirichlet);
  for(size_t i=0;i<3;++i){
    BOOST_CHECK_CLOSE(discrete[i]+1.0,sample[i]+1.0,  5); 
  }

   
  BOOST_CHECK_CLOSE(Update[0]+Update[1]+Update[2], 1.0, 0.0001);

  BOOST_CHECK_CLOSE(LogNorm1,-0.013015253 , 0.0001);  
  BOOST_CHECK_CLOSE(LogNorm2, LogNorm1, 0.0001);  
  
} 

BOOST_AUTO_TEST_SUITE_END()


/*****************************************************
 *****************************************************
 *****       Observed Nodes TEST               *******
 *****************************************************
 *****************************************************/


BOOST_AUTO_TEST_SUITE( Observed_Node_test )

BOOST_AUTO_TEST_CASE( Observed_Gaussian_test  )
{
  ObservedNode<Gaussian,double> obs2(2.0); //2 elements value 2.0
  
  Moments<double> M2 = obs2.GetMoments();
  BOOST_CHECK_EQUAL(M2.size(), size_t(2));
  BOOST_CHECK_CLOSE(M2[0], 2.0, 0.0001);
  BOOST_CHECK_CLOSE(M2[1], 4.0, 0.0001);


  std::vector<double> mean2 = obs2.GetMean();
  std::vector<double> var2 = obs2.GetVariance();
  
  BOOST_CHECK_CLOSE(mean2[0], 2.0, 0.0001);
  BOOST_CHECK_EQUAL(mean2.size(), size_t(1));
  

  BOOST_CHECK_CLOSE(var2[0], 0.0, 0.0001);
  BOOST_CHECK_EQUAL(var2.size(), size_t(1));
  
  Coster C = 0;
  obs2.Iterate(C);
  BOOST_CHECK_CLOSE(double(C), 0.0, 0.0001);
  
}
BOOST_AUTO_TEST_CASE( Observed_RGaussian_test  )
{
  ObservedNode<RectifiedGaussian,double> obs2(2.0); //2 elements value 2.0
  
  Moments<double> M2 = obs2.GetMoments();
  BOOST_CHECK_EQUAL(M2.size(), size_t(2));
  BOOST_CHECK_CLOSE(M2[0], 2.0, 0.0001);
  BOOST_CHECK_CLOSE(M2[1], 4.0, 0.0001);


  std::vector<double> mean2 = obs2.GetMean();
  std::vector<double> var2 = obs2.GetVariance();
  
  BOOST_CHECK_CLOSE(mean2[0], 2.0, 0.0001);
  BOOST_CHECK_EQUAL(mean2.size(), size_t(1));
  

  BOOST_CHECK_CLOSE(var2[0], 0.0, 0.0001);
  BOOST_CHECK_EQUAL(var2.size(), size_t(1));
  
  Coster C = 0;
  obs2.Iterate(C);
  BOOST_CHECK_CLOSE(double(C), 0.0, 0.0001);
  
}
BOOST_AUTO_TEST_CASE( Observed_Gamma_test  )
{
  ObservedNode<Gamma,double> obs2(2.0); //2 elements value 2.0
  
  Moments<double> M2 = obs2.GetMoments();
  BOOST_CHECK_EQUAL(M2.size(), size_t(2));
  BOOST_CHECK_CLOSE(M2[0], 2.0, 0.0001);
  BOOST_CHECK_CLOSE(M2[1], std::log(2), 0.0001);


  std::vector<double> mean2 = obs2.GetMean();
  std::vector<double> var2  = obs2.GetVariance();
  
  BOOST_CHECK_CLOSE(mean2[0], 2.0, 0.0001);
  BOOST_CHECK_EQUAL(mean2.size(), size_t(1));
  

  BOOST_CHECK_CLOSE(var2[0], 0.0, 0.0001);
  BOOST_CHECK_EQUAL(var2.size(), size_t(1));
  
  Coster C = 0;
  obs2.Iterate(C);
  BOOST_CHECK_CLOSE(double(C), 0.0, 0.0001);
  
}

BOOST_AUTO_TEST_CASE( Observed_Dirichlet_test  )
{
  ObservedNode<Dirichlet,double> obs3(3,2.0); //3 elements value 2.0
  
  Moments<double> M3 = obs3.GetMoments();
  BOOST_CHECK_EQUAL(M3.size(), size_t(3));
  BOOST_CHECK_CLOSE(M3[0], 2.0, 0.0001);
  BOOST_CHECK_CLOSE(M3[1], 2.0, 0.0001);
  BOOST_CHECK_CLOSE(M3[2], 2.0, 0.0001);


  std::vector<double> mean3 = obs3.GetMean();
  std::vector<double> var3  = obs3.GetVariance();
  
  BOOST_CHECK_EQUAL(mean3.size(), size_t(3));
  BOOST_CHECK_CLOSE(mean3[0], 2.0, 0.0001);
  BOOST_CHECK_CLOSE(mean3[1], 2.0, 0.0001);
  BOOST_CHECK_CLOSE(mean3[2], 2.0, 0.0001);
  

  BOOST_CHECK_EQUAL(var3.size(), size_t(3));
  BOOST_CHECK_CLOSE(var3[0], 0.0, 0.0001);
  BOOST_CHECK_CLOSE(var3[1], 0.0, 0.0001);
  BOOST_CHECK_CLOSE(var3[2], 0.0, 0.0001);
  
  Coster C = 0;
  obs3.Iterate(C);
  BOOST_CHECK_CLOSE(double(C), 0.0, 0.0001);
  
}
BOOST_AUTO_TEST_SUITE_END()





/*****************************************************
 *****************************************************
 *****       Factor TEST                       *******
 *****************************************************
 *****************************************************/

BOOST_AUTO_TEST_SUITE( Factor_Node_test )

BOOST_AUTO_TEST_CASE( Factor_Gaussian_test  )
{

  ObservedNode<Gaussian,double> obsGaussian(2.0);
  ObservedNode<Gamma,double>    obsGamma(3.0); 
  ObservedNode<Gaussian,double> obsData(4.0); 
  
  detail::Factor<Gaussian,double> GF(&obsGaussian, &obsGamma, &obsData);
  BOOST_CHECK_CLOSE(GF.CalcLogNorm(), 0.0, 0.0001);
  
  NaturalParameters<double> NPp1 = GF.GetNaturalNot(&obsGaussian);
  NaturalParameters<double> NPp2 = GF.GetNaturalNot(&obsGamma);
  NaturalParameters<double> NPd  = GF.GetNaturalNot(&obsData);
  
  BOOST_CHECK_CLOSE(NPp1[0], 12.0 , 0.0001);
  BOOST_CHECK_CLOSE(NPp1[1], -1.5 ,0.0001 );
  BOOST_CHECK_CLOSE(NPp2[0], -0.5*(16 - 16+ 4),0.0001 );
  BOOST_CHECK_CLOSE(NPp2[1], 0.5 ,0.0001 );
  BOOST_CHECK_CLOSE(NPd[0], 6 ,0.0001 );
  BOOST_CHECK_CLOSE(NPd[1], -1.5 ,0.0001 );
  BOOST_CHECK_CLOSE(GF.CalcLogNorm(), 0.5* ( std::log(3/(2.0*M_PI)) - (3) * (4)), 0.0001);
  
  
  //calc 5000 samples and check mean
  std::vector< Moments<double> > sample(5000);
  double av = 0;
  double av_squared = 0;
  for(size_t i=0;i<sample.size();++i){
    sample[i] = GF.InitialiseMoments();
    const double av_tmp = sample[i][0];
    av+=av_tmp;
    av_squared+=av_tmp*av_tmp;
  }
  av*=(1.0/sample.size());
  double var = av_squared*(1.0/sample.size()) - av*av;
  BOOST_CHECK_CLOSE(av,  2.0, 5);
  BOOST_CHECK_CLOSE(var, 1.0/3.0, 5);
  

}
BOOST_AUTO_TEST_CASE( Factor_Dirichlet_test  )
{

  ObservedNode<Dirichlet,double> obsDirichlet(3,1.1); 
  ObservedNode<Gaussian,double> obsData(4.0); 
  
  detail::Factor<Dirichlet,double> GF(&obsDirichlet,  &obsData);
  BOOST_CHECK_CLOSE(GF.CalcLogNorm(), 0.0, 0.0001);
  
  NaturalParameters<double> NPd = GF.GetNaturalNot(&obsData);
  
  BOOST_CHECK_CLOSE(NPd[0], 0.1 , 0.0001);  //1.1-1
  BOOST_CHECK_CLOSE(NPd[1], 0.1 , 0.0001);  //1.1-1
  BOOST_CHECK_CLOSE(NPd[2], 0.1 , 0.0001);  //1.1-1
  BOOST_CHECK_CLOSE(GF.CalcLogNorm(),gsl_sf_lngamma(3.3) -3.0*gsl_sf_lngamma(1.1) , 0.0001);
  
  
  //calc 5000 samples and check mean
  std::vector< Moments<double> > sample(5000);
  std::vector<double> av(3);
  std::vector<double> av_squared(3);
  for(size_t i=0;i<sample.size();++i){
    sample[i] = GF.InitialiseMoments();
    for(size_t j=0;j<3;++j){
      const double av_tmp = std::exp(sample[i][j]);
      av[j]+=av_tmp;
      av_squared[j]+=av_tmp*av_tmp;
    }
  }
  std::vector<double> var(3);
  for(size_t j=0;j<3;++j){
    av[j]*=(1.0/sample.size());
    var[j] =  av_squared[j]*(1.0/sample.size()) - av[j]*av[j];
    BOOST_CHECK_CLOSE(av[j],  1.0/3.0, 5);
    BOOST_CHECK_CLOSE(var[j], ((2.0/9.0)/(4.3)), 5);
  }
}

BOOST_AUTO_TEST_CASE( Factor_Discrete_test  )
{

  ObservedNode<Dirichlet,double> obsPrior(3,1.0/3.0); 
  ObservedNode<Dirichlet,double> obsData(3,2.0/3.0); 
  
  detail::Factor<Discrete,double> GF(&obsPrior,  &obsData);
  BOOST_CHECK_CLOSE(GF.CalcLogNorm(), 0.0, 0.0001);
  
  NaturalParameters<double> NPp = GF.GetNaturalNot(&obsPrior);
  NaturalParameters<double> NPd = GF.GetNaturalNot(&obsData);
  
  BOOST_CHECK_CLOSE(NPp[0],  2.0/3.0, 0.0001);  //1.1-1
  BOOST_CHECK_CLOSE(NPp[1],  2.0/3.0 , 0.0001);  //1.1-1
  BOOST_CHECK_CLOSE(NPp[2],  2.0/3.0 , 0.0001);  //1.1-1
  BOOST_CHECK_CLOSE(NPd[0],  1.0/3.0, 0.0001);  //1.1-1
  BOOST_CHECK_CLOSE(NPd[1],  1.0/3.0, 0.0001);  //1.1-1
  BOOST_CHECK_CLOSE(NPd[2],  1.0/3.0, 0.0001);  //1.1-1
  BOOST_CHECK_CLOSE(GF.CalcLogNorm(),-1.431945622001443, 0.0001);
  
  
  //calc 5000 samples and check mean
  std::vector< Moments<double> > sample(1);
  double av = 0;
  double av_squared = 0;
  for(size_t i=0;i<sample.size();++i){
    sample[i] = GF.InitialiseMoments();
    const double av_tmp = sample[i][0];
    av+=av_tmp;
    av_squared+=av_tmp*av_tmp;
  }
  av*=(1.0/sample.size());
  double var = av_squared*(1.0/sample.size()) - av*av;
  BOOST_CHECK_CLOSE(av,  1.3956124250860897, 5);
  BOOST_CHECK_CLOSE(var, 0, 5);
  
}

BOOST_AUTO_TEST_SUITE_END()

/*****************************************************
 *****************************************************
 *****       Hidden Nodes TEST                 *******
 *****************************************************
 *****************************************************/


BOOST_AUTO_TEST_SUITE( Hidden_Node_test )

BOOST_AUTO_TEST_CASE( Hidden_Gaussian_test  )
{
  HiddenNode<Gaussian,double> G; //2 elements value 2.0
  Moments<double> M2 = G.GetMoments();
  BOOST_CHECK_EQUAL(M2.size(), size_t(2));
  BOOST_CHECK_CLOSE(M2[0], 0.0, 0.0001);
  BOOST_CHECK_CLOSE(M2[1], 0.0, 0.0001);

  ObservedNode<Gaussian,double> obsGaussian(2.0);
  ObservedNode<Gamma,double>    obsGamma(3.0); 

  Random::Restart(10);
  detail::Factor<Gaussian,double> GF(&obsGaussian, &obsGamma, &G);

  
  Moments<double> M22 = G.GetMoments();
  BOOST_CHECK_EQUAL(M22.size(), size_t(2));
  BOOST_CHECK_CLOSE(M22[0], 1.7337736219042532, 0.0001); //inialised
  BOOST_CHECK_CLOSE(M22[1], 3.3393043053443257, 0.0001); //inialised

  std::vector<double> mean2 = G.GetMean();
  std::vector<double> var2  = G.GetVariance();
  
  BOOST_CHECK_CLOSE(mean2[0], 2.0, 0.0001);
  BOOST_CHECK_EQUAL(mean2.size(), size_t(1));
  

  BOOST_CHECK_CLOSE(var2[0], 1.0/3.0, 0.0001);
  BOOST_CHECK_EQUAL(var2.size(), size_t(1));
  
  Coster C = 0; 
  G.Iterate(C);
  BOOST_CHECK_CLOSE(double(C), 0.0 , 0.0001); //No child so no cost
  
  Moments<double> M23 = G.GetMoments();
  BOOST_CHECK_EQUAL(M23.size(), size_t(2));
  BOOST_CHECK_CLOSE(M23[0], 2.0, 0.0001);
  BOOST_CHECK_CLOSE(M23[1], 13.0/3.0, 0.0001);

  //Next 

  ObservedNode<Gaussian,double> obsData1(4.0); 
  detail::Factor<Gaussian,double> GF2(&G, &obsGamma, &obsData1);
  
  G.Iterate(C);
  double Cost =( -12*3 +  1.5*(9+1.0/6.0) )   -6.36963 +27.0231 ;
  BOOST_CHECK_CLOSE(double(C), Cost , 0.01); //No child so no cost
  
  Moments<double> M33 = G.GetMoments();
  BOOST_CHECK_EQUAL(M33.size(), size_t(2));
  BOOST_CHECK_CLOSE(M33[0], 3.0, 0.0001);
  BOOST_CHECK_CLOSE(M33[1], 9+1.0/6.0, 0.0001);

  
  BOOST_CHECK_EQUAL(G.size(), size_t(2));
}

 BOOST_AUTO_TEST_SUITE_END()




/*****************************************************
 *****************************************************
 *****      Calculation Tree        TEST       *******
 *****************************************************
 *****************************************************/

BOOST_AUTO_TEST_SUITE( Calculation_Tree_test )

BOOST_AUTO_TEST_CASE( ExpressionFactor_test  )
{
  //Create an expression X*Y + Z.
  ExpressionFactory<double> Factory;  // Create a Factory.
       
  Placeholder<double>* X = Factory.placeholder();
  Placeholder<double>* Y = Factory.placeholder();
  Placeholder<double>* Z = Factory.placeholder();
  Expression<double>* XY = Factory.Multiply(X,Y);
  Expression<double>* Expr = Factory.Add(XY,Z);

  
  //Create some variables
  Builder<double> builder;
  Builder<double>::GaussianNode x = builder.gaussian(0.0,0.01);
  Builder<double>::GaussianNode y = builder.gaussian(0.0,0.01);
  Builder<double>::GaussianNode z = builder.gaussian(0.0,0.01);

  Random::Restart(10);
  x->InitialiseMoments();
  Random::Restart(10);
  y->InitialiseMoments();
  Random::Restart(10);
  z->InitialiseMoments();

  Moments<double> Mx = x->GetMoments();
  BOOST_CHECK_EQUAL(Mx.size(), size_t(2));
  BOOST_CHECK_CLOSE(Mx[0], -4.6111761317687545, 0.0001); //inialised
  BOOST_CHECK_CLOSE(Mx[1], 121.26294531819386, 0.0001); //inialised
  Moments<double> My = y->GetMoments();
  Moments<double> Mz = z->GetMoments();

  
  //Assign this context to the expression.
  Context<double> context;
  context.Assign(X,x);
  context.Assign(Y,y);
  context.Assign(Z,z);

  //The context provides the Moments of every element in expression.
  const SubContext<double> M0 = context[0];  //All the first moments  (the <x>'s of every element in expr)
  const SubContext<double> M1 = context[1];  //The second moment (the <x^2> of every element of expression)
  
  BOOST_CHECK_CLOSE(Expr->Evaluate(M0), Mx[0]*My[0]+Mz[0]  , 0.001);
  BOOST_CHECK_CLOSE(Expr->Evaluate(M1), Mx[1]*My[1]+Mz[1]  , 0.001);
  const double prec = 1.0/(Expr->Evaluate(M1) - Expr->Evaluate(M0*M0) );

  BOOST_CHECK_CLOSE(prec,6.9673840417644678e-05, 0.001);
  
  BOOST_CHECK_CLOSE(Gaussian<double>::CalcNP2Deterministic(Expr,context)[0],
		    Expr->Evaluate(M0) *prec, 0.001);
  BOOST_CHECK_CLOSE(Gaussian<double>::CalcNP2Deterministic(Expr,context)[1],
		    -0.5*prec, 0.001);


  // typedef detail::Deterministic<Gaussian, double >  DeterministicFactor;
  // Builder<double>::GaussianResultNode D  = new Builder<double>::GaussianResultNode();
  // DeterministicFactor DF = new DeterministicFactor(Expr,context,D);
  
  
  

  // BOOST_CHECK_CLOSE(Gaussian<double>::CalcNP2Parent(Expr,context)[0],
  // 		    Expr->Evaluate(M0) *prec, 0.001);
  // BOOST_CHECK_CLOSE(Gaussian<double>::CalcNP2Parent(Expr,context)[1],
  // 		    -0.5*prec, 0.001);
}



BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( Add_Test_Suite )



// BOOST_AUTO_TEST_CASE( Add_test  )
// {
//    std::ofstream filestream("Add1.txt");
//    {
//      ScopedRedirect redirect(std::cout, filestream);
//   rng random(10);
//   size_t data_points = 150;//5447
//   typedef Builder<double>::GaussianNode GaussianNode;
//   typedef Builder<double>::GammaNode    GammaNode;
//   typedef Builder<double>::GaussianDataNode    Datum;
//   typedef Builder<double>::GaussianResultNode    ResultNode;
   
//   Builder<double> Build;
//   GammaNode    Precision = Build.Gamma(0.01,0.01);
//   GaussianNode Parent11    = Build.Gaussian(0.0,0.01);
//   GaussianNode Mean2    = Build.Gaussian(0.0,0.01);
  
//   std::cout<<"Mean1 = "<<Mean1<<std::endl;
//   std::cout<<"Mean2 = "<<Mean2<<std::endl;

//   ResultNode SumMeans = Build.Add(Mean1,Mean2);
   
//   for(size_t i=0;i<data_points;++i) 
//     {
//       double splitter = random.uniform();
//       double data;
//       // if (splitter>0.8) {
//       // 	data = random.gaussian(1.0/std::sqrt(50),0);
//       // }
//       // else if(splitter>0.4) {
//       // 	data = random.gaussian(1.0/std::sqrt(1),0);
//       // }
//       // else {
// 	data = random.gaussian(1.0/std::sqrt(0.3),6);
//       // }
// 	//    std::cout<<"data = "<<data<<std::endl;

//       Build.Join(SumMeans, Precision, data);
//     }
  
//   Build.Run(1e-6, 10);
   

//   double mean1= Mean1->GetMoments()[0];
//   double mean2= Mean2->GetMoments()[0];
//   double sum= SumMeans->GetMoments()[0];
//   double prec = Precision->GetMoments()[0];
//   std::cout<<"Gaussian mean1      = "<<mean1<<std::endl;
//   std::cout<<"Gaussian mean2      = "<<mean2<<std::endl;
//   std::cout<<"Gaussian meansum    = "<<sum<<std::endl;
//   std::cout<<"Gaussian precission = "<<prec<<std::endl;
//    }
// }



BOOST_AUTO_TEST_CASE( Add_test_2  )
{

  Random::Restart(10);
   std::ofstream filestream("Add2.txt");
   {
     ScopedRedirect redirect(std::cout, filestream);
  rng random(10);
  size_t data_points = 150;//5447
  typedef Builder<double>::GaussianNode GaussianNode;
  typedef Builder<double>::GammaNode    GammaNode;
  typedef Builder<double>::GaussianDataNode    Datum;
  typedef Builder<double>::GaussianResultNode    ResultNode;
   
  Builder<double> Build;
  GammaNode    Precision = Build.gamma(0.01,0.01);
  GaussianNode Mean1    = Build.gaussian(0.0,0.01);
  GaussianNode Mean2    = Build.gaussian(0.0,0.01);
  
  std::cout<<"Mean1 = "<<Mean(Mean1)<<std::endl;
  std::cout<<"Mean2 = "<<Mean(Mean2)<<std::endl;
  ExpressionFactory<double> Factory;
  Placeholder<double>* P1 = Factory.placeholder();
  Placeholder<double>* P2 = Factory.placeholder();
  Expression<double>* Expr = Factory.Add(P1,P2);
  Context<double> context;
  context.Assign(P1,Mean1);
  context.Assign(P2,Mean2);
  ResultNode SumMeans = Build.calc_gaussian(Expr,context);
   
  for(size_t i=0;i<data_points;++i) 
    {
      double splitter = random.uniform();
      double data;
      // if (splitter>0.8) {
      // 	data = random.gaussian(1.0/std::sqrt(50),0);
      // }
      // else if(splitter>0.4) {
      // 	data = random.gaussian(1.0/std::sqrt(1),0);
      // }
      // else {
	data = random.gaussian(1.0/std::sqrt(0.3),6);
      // }
	// std::cout<<"data = "<<data<<std::endl;

      Build.join(SumMeans, Precision, data);
    }
  
  Build.run(1e-6, 10);
   

  double mean1= Mean1->GetMoments()[0];
  double mean2= Mean2->GetMoments()[0];
  double sum= SumMeans->GetMoments()[0];
  double prec = Precision->GetMoments()[0];
  std::cout<<"Gaussian mean1      = "<<mean1<<std::endl;
  std::cout<<"Gaussian mean2      = "<<mean2<<std::endl;
  std::cout<<"Gaussian meansum    = "<<sum<<std::endl;
  std::cout<<"Gaussian precission = "<<prec<<std::endl;
   }
}
 BOOST_AUTO_TEST_SUITE_END()
 
 BOOST_AUTO_TEST_SUITE( Multiply_Test_Suite )
 

// BOOST_AUTO_TEST_CASE( Multiply_test  )
// {
//   Random::Restart(10);
//    std::ofstream filestream("Multiply1.txt");
//    {
//      ScopedRedirect redirect(std::cout, filestream);
//   rng random(10);
//   size_t data_points = 150;//5447
//   typedef Builder<double>::GaussianNode GaussianNode;
//   typedef Builder<double>::GammaNode    GammaNode;
//   typedef Builder<double>::GaussianDataNode    Datum;
//   typedef Builder<double>::GaussianResultNode    ResultNode;
   
//   Builder<double> Build;
//   GammaNode    Precision = Build.gamma(0.01,1.0);
//   GaussianNode Mean1    = Build.gaussian(0.0,0.01);
//   GaussianNode Mean2    = Build.gaussian(0.0,0.01);
  
//   std::cout<<"Mean1 = "<<Mean1<<std::endl;
//   std::cout<<"Mean2 = "<<Mean2<<std::endl;

//   ResultNode SumMeans = Build.Multiply(Mean1,Mean2);
   
//   for(size_t i=0;i<data_points;++i) 
//     {
//       double splitter = random.uniform();
//       double data;
//       // if (splitter>0.8) {
//       // 	data = random.gaussian(1.0/std::sqrt(50),0);
//       // }
//       // else if(splitter>0.4) {
//       // 	data = random.gaussian(1.0/std::sqrt(1),0);
//       // }
//       // else {
// 	data = random.gaussian(1.0/std::sqrt(0.3),6);
//       // }
// 	//    std::cout<<"data = "<<data<<std::endl;

//       Build.join(SumMeans, Precision, data);
//     }
  
//   Build.Run(1e-6, 10);
   

//   double mean1= Mean1->GetMoments()[0];
//   double mean2= Mean2->GetMoments()[0];
//   double sum= SumMeans->GetMoments()[0];
//   double prec = Precision->GetMoments()[0];
//   std::cout<<"Gaussian mean1      = "<<mean1<<std::endl;
//   std::cout<<"Gaussian mean2      = "<<mean2<<std::endl;
//   std::cout<<"Gaussian meansum    = "<<sum<<std::endl;
//   std::cout<<"Gaussian precission = "<<prec<<std::endl;
//    }
// }



BOOST_AUTO_TEST_CASE( Multiply_test_2  )
{

  Random::Restart(10);
  std::ofstream filestream("Multiply2.txt");
   {
       ScopedRedirect redirect(std::cout, filestream);
  rng random(10);
  size_t data_points = 150;//5447
  typedef Builder<double>::GaussianNode GaussianNode;
  typedef Builder<double>::GammaNode    GammaNode;
  typedef Builder<double>::GaussianDataNode    Datum;
  typedef Builder<double>::GaussianResultNode    ResultNode;
   
  Builder<double> Build;
  GammaNode    Precision = Build.gamma(0.01,0.01);
  GaussianNode Mean1    = Build.gaussian(0.0,0.01);
  GaussianNode Mean2    = Build.gaussian(0.0,0.01);
  
  std::cout<<"Mean1 = "<<Mean1<<std::endl;
  std::cout<<"Mean2 = "<<Mean2<<std::endl;
  ExpressionFactory<double> Factory;
  Placeholder<double>* P1 = Factory.placeholder();
  Placeholder<double>* P2 = Factory.placeholder();
  Expression<double>* Expr = Factory.Multiply(P1,P2);
  Context<double> context;
  context.Assign(P1,Mean1);
  context.Assign(P2,Mean2);
  ResultNode SumMeans = Build.calc_gaussian(Expr,context);
   
  for(size_t i=0;i<data_points;++i) 
    {
      double splitter = random.uniform();
      double data;
      // if (splitter>0.8) {
      // 	data = random.gaussian(1.0/std::sqrt(50),0);
      // }
      // else if(splitter>0.4) {
      // 	data = random.gaussian(1.0/std::sqrt(1),0);
      // }
      // else {
	data = random.gaussian(1.0/std::sqrt(0.3),6);
      // }
	// std::cout<<"data = "<<data<<std::endl;

      Build.join(SumMeans, Precision, data);
    }
  
  Build.run(1e-6, 5);
   

  double mean1= Mean1->GetMoments()[0];
  double mean2= Mean2->GetMoments()[0];
  double sum= SumMeans->GetMoments()[0];
  double prec = Precision->GetMoments()[0];
  std::cout<<"Gaussian mean1      = "<<mean1<<std::endl;
  std::cout<<"Gaussian mean2      = "<<mean2<<std::endl;
  std::cout<<"Gaussian meansum    = "<<sum<<std::endl;
  std::cout<<"Gaussian precission = "<<prec<<std::endl;
   }
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( DotProduct_Test_Suite )
 
BOOST_AUTO_TEST_CASE( DotProduct_test  )
{

  Random::Restart(10);
  std::ofstream filestream("DotProduct.txt");
   {
     //     ScopedRedirect redirect(std::cout, filestream);
  rng random(10);
  size_t data_points = 150;//5447
  typedef Builder<double>::GaussianNode GaussianNode;
  typedef Builder<double>::GammaNode    GammaNode;
  typedef Builder<double>::GaussianDataNode    Datum;
  typedef Builder<double>::GaussianResultNode    ResultNode;
   
  Builder<double> Build;
  GammaNode    Precision = Build.gamma(0.01,0.01);
  GaussianNode S1    = Build.gaussian(0.0,0.01);
  GaussianNode S2    = Build.gaussian(0.0,0.01);
  GaussianNode A1    = Build.gaussian(0.0,0.01);
  GaussianNode A2    = Build.gaussian(0.0,0.01);
  
  ExpressionFactory<double> Factory;
  Placeholder<double>* PA1 = Factory.placeholder();
  Placeholder<double>* PA2 = Factory.placeholder();
  Placeholder<double>* PS1 = Factory.placeholder();
  Placeholder<double>* PS2 = Factory.placeholder();
  Expression<double>* Expr1 = Factory.Multiply(PA1,PS1);
  Expression<double>* Expr2 = Factory.Multiply(PA2,PS2);
  Expression<double>* Expr = Factory.Add(Expr1,Expr2);
  Context<double> context;
  context.Assign(PA1,A1);
  context.Assign(PA2,A2);
  context.Assign(PS1,S1);
  context.Assign(PS2,S2);
  ResultNode SumMeans = Build.calc_gaussian(Expr,context);
   
  for(size_t i=0;i<data_points;++i) 
    {
      double splitter = random.uniform();
      double data;
      if (splitter>0.8) {
       	data = random.gaussian(1.0/std::sqrt(50),0);
      }
      else if(splitter>0.4) {
      	data = random.gaussian(1.0/std::sqrt(1),0);
      }
      else {
	data = random.gaussian(1.0/std::sqrt(0.3),6);
      }
      // std::cout<<"data = "<<data<<std::endl;

      Build.join(SumMeans, Precision, data);
    }
  
  Build.run(1e-6, 10);
   
  std::cout<<"Gaussian S1      = "<< S1->GetMoments()[0] <<std::endl;
  std::cout<<"Gaussian S2      = "<< S2->GetMoments()[0]<<std::endl;
  std::cout<<"Gaussian A1      = "<< A1->GetMoments()[0] <<std::endl;
  std::cout<<"Gaussian A2      = "<< A2->GetMoments()[0]<<std::endl;
  std::cout<<"Gaussian meansum    = "<<SumMeans->GetMoments()[0]<<std::endl;
  std::cout<<"Gaussian precission = "<<Precision->GetMoments()[0]<<std::endl;
   }
}
BOOST_AUTO_TEST_SUITE_END()


//  BOOST_AUTO_TEST_CASE( GaussianConstant_test  )
// {
  
//   typedef Builder::GaussianConstNode Const;
//   Builder Build;
//   Const mu = Build.GaussianConst(2);
//   Moments moments = mu->GetMoments();
//   BOOST_CHECK_CLOSE(moments[0], 2.0, 0.0001);
//   BOOST_CHECK_CLOSE(moments[1], 4.0, 0.0001);
//   BOOST_CHECK(Build.NumberOfNodes()   == 1);
//   BOOST_CHECK(Build.NumberOfFactors() == 0);
  
// }

//  BOOST_AUTO_TEST_CASE( DataConstant_test  )
// {
  
//   typedef Builder::GaussianDataNode Const;
//   Builder Build;
//   Const mu = Build.GaussianData(2);
//   Moments moments = mu->GetMoments();
//   BOOST_CHECK_CLOSE(moments[0], 2.0, 0.0001);
//   BOOST_CHECK_CLOSE(moments[1], 4.0, 0.0001);
//   BOOST_CHECK(Build.NumberOfNodes()   == 1);
//   BOOST_CHECK(Build.NumberOfFactors() == 0);
  

//   std::cout<<"SIZEOF DATA = "<<sizeof( *mu )<<std::endl;


// }

//  BOOST_AUTO_TEST_CASE( GammaConstant_test  )
// {
  
//   typedef Builder::GammaDataNode Const;
//   Builder Build;
//   Const mu = Build.GammaData(2);
//   Moments moments = mu->GetMoments();
//   BOOST_CHECK_CLOSE(moments[0], 2.0, 0.0001);
//   BOOST_CHECK_CLOSE(moments[1], 0.69314718055994529, 0.0001);

//   BOOST_CHECK(Build.NumberOfNodes()   ==  (size_t) 1);
//   BOOST_CHECK(Build.NumberOfFactors() ==  (size_t) 0);
  
// }



//  BOOST_AUTO_TEST_SUITE_END()


//  BOOST_AUTO_TEST_SUITE( FactorConstructionTests_test )

 // BOOST_AUTO_TEST_CASE( GammaFactor_test  )
 // {
 //   Builder<double> Build;

 //   typedef Builder<double>::GaussianNode GaussianNode;
 //   typedef Builder<double>::GammaNode    GammaNode;
   

 //   GammaNode Scale =  Build.Gamma(2,3);

 //   // std::cout<<"scale = "<<Scale<<std::endl;

   
 //   // BOOST_CHECK_EQUAL(Scale->NumberOfNeighbours(), (size_t) 1);
  
 //   BOOST_CHECK(Build.NumberOfNodes()   ==  (size_t) 3);
 //   BOOST_CHECK(Build.NumberOfFactors() ==  (size_t) 1);
 //   //std::cout<<"size of builder = "<<sizeof(Build)<<std::endl;

 // }

// BOOST_AUTO_TEST_CASE( TestSize_test  )
// {
//   boost::shared_ptr<ConstantNode<NormalConstant> > Shape(new ConstantNode<NormalConstant>(2));
//   boost::shared_ptr<ConstantNode<GammaConstant> > IScale(new ConstantNode<GammaConstant>(3));
	
//   boost::shared_ptr<DataNode<GammaConstant> > Data(new DataNode<GammaConstant>(3));
  
//   boost::shared_ptr<Factor<Gamma> > GammaF (new Factor<Gamma>(Shape.get(), IScale.get(), Data.get()));

  
  
//    std::cout<<"size of const  = "<<sizeof(*Shape)<<std::endl;
//    std::cout<<"size of data   = "<<sizeof(*Data)<<std::endl;
//    std::cout<<"size of factor = "<<sizeof(*GammaF)<<std::endl;
//    std::cout<<"total size of each data point = "<<sizeof(*Data)+sizeof(*GammaF)<<std::endl;
   
//    std::cout<<"mutex = "<<sizeof(boost::mutex)<<std::endl;
//    std::cout<<"condition = "<<sizeof(boost::condition)<<std::endl;
//    std::cout<<"thread = "<<sizeof(boost::thread)<<std::endl;
   
   
// }

// BOOST_AUTO_TEST_CASE( GaussianFactor_test  )
// {
//   Builder Build;

//   typedef Builder::GaussianNode GaussianNode;
   

//    GaussianNode Scale =  Build.Gaussian(2,3);

   
 //    // BOOST_CHECK_EQUAL(Scale->NumberOfNeighbours(), (size_t) 1);
  
 //    BOOST_CHECK(Build.NumberOfNodes()   ==  (size_t) 3);
 //    BOOST_CHECK(Build.NumberOfFactors() ==  (size_t) 1);
 //  }


 // BOOST_AUTO_TEST_CASE( Weight_test  )
 // {
 //   ICR::maths::rng random;
 //   size_t data_points = 1;//5447
 //   typedef Builder<double>::GaussianNode GaussianNode;
 //   typedef Builder<double>::WeightsNode WeightsNode;
 //   typedef Builder<double>::GammaNode    GammaNode;
 //   typedef Builder<double>::GaussianDataNode    Datum;
 //   typedef Builder<double>::Variable    Variable;
   
 //   Builder<double> Build;
   
 //   WeightsNode Weights = Build.Weights(4);

 //   BOOST_CHECK(Build.NumberOfNodes()   ==  (size_t) 3);
 //   BOOST_CHECK(Build.NumberOfFactors() ==  (size_t) 2);

   
 //   std::cout<<"Weights   = "<<Weights<<std::endl;
 //   Build.Run(1e-6, 10);
 // }

 //  BOOST_AUTO_TEST_SUITE_END()



// BOOST_AUTO_TEST_CASE( Tree_test  )
// {
//   std::ofstream filestream("Tree.txt");
//   {
//     ScopedRedirect redirect(std::cout, filestream);
//     ICR::maths::rng random(10);
//     size_t data_points = 150;//5447
//     typedef Builder<double>::GaussianNode GaussianNode;
//     typedef Builder<double>::GammaNode    GammaNode;
//     typedef Builder<double>::GaussianDataNode    Datum;
   
//     Builder<double> Build;
//     // GaussianNode HypMean = Build.Gaussian(0.001,0.001);
//     GammaNode    Precision = Build.Gamma(0.01,0.01);
//     GaussianNode Mean    = Build.Gaussian(0.0,0.01);
   
//     for(size_t i=0;i<data_points;++i) 
//       {
// 	double splitter = random.uniform();
// 	double data;
// 	if (splitter>0.8) {
// 	  data = random.gaussian(1.0/std::sqrt(50),0);
// 	}
// 	else if(splitter>0.4) {
// 	  data = random.gaussian(1.0/std::sqrt(1),0);
// 	}
// 	else {
// 	  data = random.gaussian(1.0/std::sqrt(0.3),6);
// 	}
// 	//std::cout<<"data = "<<data<<std::endl;

// 	Build.Join(Mean, Precision, data);
//       }
  
//     //     Build.Join(Mean, Precision, random.gaussian(1.0/(3),2));
  
   
//     // for(size_t i=0;i<data_points;++i){
//     //   Build.Join(Mean, Precision, 1.0*random.gaussian(1.0/(3),2)+0.0*random.gaussian(1.0/(3),5));
//     // }


//     // BOOST_CHECK(Build.NumberOfNodes()   ==  (size_t) 6 + data_points );
//     // BOOST_CHECK(Build.NumberOfFactors() ==  (size_t) 2 + data_points);
   
//     Build.Run(1e-6, 10);
   
//     // for(size_t i=0;i<10000;++i){
//     //   Build.Iterate();

//     // }


//     double mean = Mean->GetMoments()[0];
//     double prec = Precision->GetMoments()[0];
//     std::cout<<"Gaussian mean = "<<mean<<std::endl;
//     std::cout<<"Gaussian precission = "<<prec<<std::endl;
//   }
// }




// BOOST_AUTO_TEST_CASE( MixtureTree_test  )
// {
//   std::ofstream filestream("MixtureTree.txt");
//   {
//     ScopedRedirect redirect(std::cout, filestream);
//     ICR::maths::rng random(10);
//     ICR::maths::rng init(11);
//     size_t data_points = 150;//5447
//     typedef Builder<double>::GaussianNode GaussianNode;
//     typedef Builder<double>::WeightsNode WeightsNode;
//     typedef Builder<double>::GammaNode    GammaNode;
//     typedef Builder<double>::GaussianDataNode    Datum;
//     typedef Builder<double>::Variable    Variable;
   
//     Builder<double> Build;


//     const size_t Components = 5;
//     WeightsNode Weights    = Build.Weights(Components);
//     std::vector<Variable> PrecMixture(Components);
//     std::vector<Variable> MeanMixture(Components);
//     for(size_t i=0;i<Components;++i){
//       PrecMixture[i] = Build.Gamma(0.01,0.01);
//       MeanMixture[i] = Build.Gaussian(0.0,0.01); //
//     }
   

//     //  Build.Run(1e-6, 100);
   
//     std::ofstream DataPoints("DataPnts.txt");
//     for(size_t i=0;i<data_points;++i){ 
//       double splitter = random.uniform();
//       double data;
//       if (splitter>0.8)
// 	data = random.gaussian(1.0/std::sqrt(50),0);
//       else if(splitter>0.4)
// 	data = random.gaussian(1.0/std::sqrt(1),0);
//       else 
// 	data = random.gaussian(1.0/std::sqrt(0.3),6);
//       //	std::cout<<"data = "<<data<<std::endl;
//       DataPoints<<data<<"\n";
//       Build.Join(MeanMixture,PrecMixture,Weights, data);
//       //     Build.Join(Mean, Precision, random.gaussian(1.0/(3),2));
//     }
//     DataPoints.close();

//     // BOOST_CHECK(Build.NumberOfNodes()   ==  (size_t) 6 + data_points );
//     // BOOST_CHECK(Build.NumberOfFactors() ==  (size_t) 2 + data_points);
   
//     Build.Run(1e-10, 100);
   
//     // for(size_t i=0;i<10000;++i){
//     //   Build.Iterate();

//     // }

//     for(size_t i=0;i<Components;++i){
//       std::cout<<"Gaussian Mean["<<i<<"] (weight "<<std::exp(Weights->GetMoments()[i])<<")\t= "<<MeanMixture[i]->GetMoments()[0]<<", Prec\t "<<PrecMixture[i]->GetMoments()[0]<<std::endl;
//     }
//     //double prec = Precision[i]->GetMoments()[0];
//     //std::cout<<"Gaussian precission = "<<prec<<std::endl;
//   }
// }





std::ostream&
operator<<(std::ostream& out, const ublas::matrix<double>& M)
{
  for(size_t j=0; j<M.size2();++j){
    for(size_t i=0;i<M.size1();++i){
      out<< M(i,j) << "\t";
    }
    out<<"\n";
  }
  return out;
}



std::ostream&
operator<<(std::ostream& out, const ublas::matrix_row<ublas::matrix<double> >& M)
{
  for(size_t j=0; j<M.size();++j){
    out<< double(M(j)) << "\t";
  }
  return out;
}




std::ostream&
operator<<(std::ostream& out, const ublas::matrix_column<ublas::matrix<double> >& M)
{
  for(size_t j=0; j<M.size();++j){
      out<< M(j) << "\t";
  }
  return out;
}




ublas::vector<double>
make_Gaussian(double mean, double precision, size_t samples , double from = 0, double to = 10)
{
  ublas::vector<double> Data(samples);

  double inorm = std::sqrt(2*precision/(2.0*M_PI));
  
  double dx = (to-from)/double(samples);
  for(size_t i=0;i<Data.size();++i){
    double x = dx*i;
    double expen = std::exp(precision*mean*x-0.5*precision*x*x+ -.5*precision*mean*mean);
    Data[i]=inorm*expen;      
  }

  return Data;

}


class GaussianComponent
{
public:
  GaussianComponent(double mean, double precision,  size_t samples ,double from = 0, double to = 10)
    : m_component(make_Gaussian(mean,precision,samples,from,to))
  {  }
  
  GaussianComponent(const ublas::vector<double>& C)
    : m_component(C)
  {  }

  const double
  operator[] (const size_t i) const {return m_component[i];}

  operator ublas::vector<double>() const {
    return m_component;
  }

  GaussianComponent& operator+=(const GaussianComponent& C)
  {
    // std::cout<<"C size = "<<size()<<std::endl;
    // std::cout<<"O size = "<<C.size()<<std::endl;
    // for(size_t i=0;i<m_component.size();++i){
    //   m_component[i]+=C[i];
    // } 
    m_component+=ublas::vector<double>(C);

    return *this;
  }

  size_t size() const
  {
    return m_component.size();
  }

  friend 
  GaussianComponent operator+(const GaussianComponent& C1,const GaussianComponent& C2)
  {
    GaussianComponent tmp = C1;
    // std::cout<<"C1 size = "<<C1.size()<<std::endl;
    // std::cout<<"C2 size = "<<C2.size()<<std::endl;

    return tmp+=C2;
  }
  
private:
  ublas::vector<double> m_component;

};

class NoiseComponent
{
public:
  NoiseComponent(size_t seed, double precision, size_t samples)
    : m_component(samples)
  {
    ICR::EnsembleLearning::rng random(seed);
    for(size_t i=0;i<samples;++i){
      m_component[i] = random.gaussian(1.0/std::sqrt(precision));
    }
  }
  operator ublas::vector<double>() const
  {
    return (m_component);
  }
  
  operator GaussianComponent() const
  {
    return GaussianComponent(m_component);
  }
private:
  ublas::vector<double> m_component;
  
};




struct exponentiate
{
  double operator()(const double d) {return std::exp(d);}
};
      

BOOST_AUTO_TEST_CASE( ICA_test  )
{
  std::ofstream filestream("MixtureTree.txt");
  {
    //ScopedRedirect redirect(std::cout, filestream);

    //Make 3 data sources from 5 Gaussians
    size_t samples = 100;
    GaussianComponent G0(0,2,samples);
    GaussianComponent G1(3,5,samples);
    GaussianComponent G2(7,1,samples);
    GaussianComponent G3(3,0.1,samples);
    GaussianComponent G4(7,0.1,samples);

    GaussianComponent s1 = G0+G1+G2;
    GaussianComponent s2 = G3;
    GaussianComponent s3 = G4;
    std::cout<<"size="<<s1.size()<<std::endl;

    // NoiseComponent N1(1,5000);
    // NoiseComponent N2(2,5000);
    // NoiseComponent N3(3,5000);
    // s1+=N1;
    // s2+=N2;
    // s3+=N3;
    {
      std::ofstream source1("ICAsource1.txt");
      std::ofstream source2("ICAsource2.txt");
      std::ofstream source3("ICAsource3.txt");
      std::cout<<s1.size()<<std::endl;

      for(size_t i=0;i<s1.size();++i){
	source1<<s1[i]<<"\n";
	source2<<s2[i]<<"\n";
	source3<<s3[i]<<"\n";
      }
    }

    //Make Mixure Matrix A 5 rows 3 cols
    ublas::matrix<double> ATrue(5,3);
    ICR::EnsembleLearning::rng random(10);
    
    for(size_t i=0;i<ATrue.size1();++i){
      for(size_t j=0;j<ATrue.size2();++j){
	ATrue(i,j) = random.gamma();
      }
    }

    for(size_t j=0;j<ATrue.size2();++j){
      ublas::matrix_column<ublas::matrix<double> > col(ATrue,j);
      std::cout<<"col = "<<col<<std::endl;

      ublas::vector<double> col2 = element_prod(col,col);
      const double norm = 1.0/std::sqrt(std::accumulate(col2.begin(), col2.end(), 0.0));
      std::cout<<"norm "<<norm<<std::endl;
      col*=norm;
      // double norm = 1.0/sqrt(mean(col2));

      // col*=norm;
      // ATrue.set_col(j, col);
    }
    
    //std::cout<<A<<std::endl;
    
    ublas::matrix<double> sources(3,samples); 
    ublas::matrix_row<ublas::matrix<double> >(sources,0) = ublas::vector<double>(s1);
    ublas::matrix_row<ublas::matrix<double> >(sources,1) = ublas::vector<double>(s2);
    ublas::matrix_row<ublas::matrix<double> >(sources,2) = ublas::vector<double>(s3);
    
    // sources.set_row(0,s1);
    // sources.set_row(1,s2);
    // sources.set_row(2,s3);
    
    ublas::matrix<double> Data = prod(ATrue,sources);
    ublas::vector<double> d0 = ublas::matrix_row<ublas::matrix<double> >(Data,0) ;//+  ublas::vector<double>(NoiseComponent(0,5000,samples));
    ublas::vector<double> d1 = ublas::matrix_row<ublas::matrix<double> >(Data,1) ;// +  ublas::vector<double>(NoiseComponent(1,5000,samples));
    ublas::vector<double> d2 = ublas::matrix_row<ublas::matrix<double> >(Data,2) ;//+  ublas::vector<double>(NoiseComponent(2,5000,samples));
    ublas::vector<double> d3 = ublas::matrix_row<ublas::matrix<double> >(Data,3) ;//+  ublas::vector<double>(NoiseComponent(3,5000,samples));
    ublas::vector<double> d4 = ublas::matrix_row<ublas::matrix<double> >(Data,4) ;//+  ublas::vector<double>(NoiseComponent(4,5000,samples));
    
    // std::
    
    {
      std::ofstream Atrue_stream("ICAATrue.txt");
      for(size_t i=0;i<ATrue.size1();++i){
	for(size_t j=0;j<ATrue.size2();++j){
	  Atrue_stream<< ATrue(i,j)<<" ";
	}
	Atrue_stream<<"\n";
      }
    }
    {
      std::ofstream data0("ICAData0.txt");
      std::ofstream data1("ICAData1.txt");
      std::ofstream data2("ICAData2.txt");
      std::ofstream data3("ICAData3.txt");
      std::ofstream data4("ICAData4.txt");
      for(size_t i=0;i<d1.size();++i){
	data0<<d0[i]<<"\n";
	data1<<d1[i]<<"\n";
	data2<<d2[i]<<"\n";
	data3<<d3[i]<<"\n";
	data4<<d4[i]<<"\n";
      }
    }
    //Make 5 Data output
    
    
    

    
    
    // ICR::maths::rng random(10);
    // ICR::maths::rng init(11);
    // size_t data_points = 150;//5447
    typedef Builder<double>::GaussianNode GaussianNode;
    typedef Builder<double>::RectifiedGaussianNode RectifiedGaussianNode;
    typedef Builder<double>::WeightsNode WeightsNode;
    typedef Builder<double>::GammaNode    GammaNode;
    typedef Builder<double>::GaussianDataNode    Datum;
    typedef Builder<double>::Variable    Variable;
    typedef Builder<double>::GaussianResultNode    ResultNode;
   
     ICR::EnsembleLearning::Random::Restart(10);
    Builder<double> Build("Cost.txt");


    const size_t Components = 3;
    const size_t M = 4; //assumed sources
    const size_t N = Data.size1(); //DataSources
    const size_t T  = Data.size2(); //DataLength

    std::cout<<"C = "<<Components<<std::endl
	     <<"M = "<<M<<std::endl
	     <<"N = "<<N<<std::endl
	     <<"T = "<<T<<std::endl;

    std::vector<WeightsNode> Weights(M);
    for(size_t m=0;m<M;++m){
      Weights[m] = Build.weights(Components);
    }
    
    std::vector< std::vector<Variable> > ShypMean(M);
    std::vector< std::vector<Variable> >    ShypPrec(M);
    for(size_t m=0;m<M;++m){ 
      ShypMean[m].resize(Components);
      ShypPrec[m].resize(Components);
      for(size_t c=0;c<Components;++c){
	ShypMean[m][c]=	Build.gaussian(0.0,0.1);
	ShypPrec[m][c]= Build.gamma(0.1,0.1);
      }


    }

    std::cout<<"built hyperparams - nodes = "<< Build.number_of_nodes() <<std::endl;
	
    std::vector< std::vector<RectifiedGaussianNode> > S(M);
    for(size_t m=0;m<M;++m){ 
      S[m].resize(T);
      for(size_t t=0;t<T;++t){ 
	S[m][t] = Build.rectified_gaussian_mixture(ShypMean[m],
						 ShypPrec[m],
						 Weights[m]);
      }
    }

    std::cout<<"built sources - nodes = "<< Build.number_of_nodes() <<std::endl;
    
    std::vector<GaussianNode> AMean(M);
    for(size_t m=0;m<M;++m){
      AMean[m] = Build.gaussian(0.0,0.1);
    }
    std::vector<GammaNode> APrecision(M);
    for(size_t m=0;m<M;++m){
      APrecision[m] = Build.gamma(0.1,0.1);
    }
    
    std::vector< std::vector<RectifiedGaussianNode> > A(N);
    for(size_t n=0;n<N;++n){ 
      A[n].resize(M);
      for(size_t m=0;m<M;++m){ 
	A[n][m] = Build.rectified_gaussian(AMean[m], APrecision[m]);
      }
    }

    
    std::vector<RectifiedGaussianNode> noiseMean(N);
    for(size_t n=0;n<N;++n){
      noiseMean[n] = Build.rectified_gaussian(0.00,0.1);
    }
    
    
    
    std::vector<GammaNode> DPrecision(N);
    for(size_t n=0;n<N;++n){
      DPrecision[n] = Build.gamma(0.1,0.1);
    }
    
    //Deterministic Node.  Need to make the expression.
    //The inner product plus noise
    /** The inner product is summed over M.
     *  Therefore there are 2M + 1 placeholders required
     *  M for Anm, M for Smt and 1 for Noise
     */
    ExpressionFactory<double> Factory;
    //Make a set of placeholders for all the elements in the expression.
    std::vector<Placeholder<double>*> SP(M);
    std::vector<Placeholder<double>*> AP(M);
    for(size_t m=0;m<M;++m){
      SP[m] = Factory.placeholder();  
      AP[m] = Factory.placeholder(); 
    }
    //Placeholder<double>* NP = Factory.placeholder();

    //The expression in terms of these placeholders
    Expression<double>* Expr;
    {
      BOOST_ASSERT(M!=0);
      //First do the multiplication
      std::deque<Expression<double>*> prod(M);
      for(size_t m=0;m<M;++m){
	prod[m] = Factory.Multiply(AP[m], SP[m]);
      }
      //then add them up
      Expr = prod[0];
      prod.pop_front();
      while(prod.size()!=0)
	{
	  Expr = Factory.Add(Expr,prod[0]);
	  prod.pop_front();
	}
      //Finnaly add the Noise placeholder
      //Expr = Factory.Add(Expr, NP);
      
    }
    
    

    std::vector< std::vector<ResultNode> > AtimesSplusN(N);
    for(size_t n=0;n<N;++n){ 
      AtimesSplusN[n].resize(T);
      for(size_t t=0;t<T;++t){ 

	/**  Now need to replace the placeholders with given nodes 
	 *   I.e. set up a context to the expression
	 */
	Context<double> context;
	for(size_t m=0;m<M;++m){
	  context.Assign(AP[m], A[n][m]);
	  context.Assign(SP[m], S[m][t]);
	}
	//context.Assign(NP, noiseMean[n]);
	
	AtimesSplusN[n][t] = Build.calc_gaussian(Expr,context);  
      }
    }

    for(size_t n=0;n<N;++n){
      for(size_t t=0;t<T;++t){
	Build.join(AtimesSplusN[n][t],DPrecision[n], Data(n,t));
      }
    }

    std::cout<<"NODES   = "<<Build.number_of_nodes()<<std::endl;
    std::cout<<"FACTORS = "<<Build.number_of_factors()<<std::endl;

     ICR::EnsembleLearning::Random::Restart(11);
    // // BOOST_CHECK(Build.NumberOfNodes()   ==  (size_t) 6 + data_points );
    // // BOOST_CHECK(Build.NumberOfFactors() ==  (size_t) 2 + data_points);
   
    

    bool not_converged = true;
    size_t max_attempts = 1;

    size_t attempt = 0;
    
    while(not_converged && attempt< max_attempts)
      {

	//check weights to see if they have coallect
	not_converged = false;
	for(size_t m=0;m<M;++m){
	  const Moments<double> W = Weights[m]->GetMoments();
	  std::vector<double> weights(W.size());
	  PARALLEL_TRANSFORM(W.begin(), W.end(), weights.begin(),exponentiate());
	  if (*min_element(weights.begin(), weights.end()) < 0.01) {
	    //Weights[m]->InitialiseMoments();
	    std::cout<<"reinitializing weights "<<m<<std::endl;
	    // for(size_t c=0;c<Components;++c){
	    //   // std::cout<<std::exp(Weights[m]->GetMoments()[c])<<"\t";
	    //   ShypMean[m][c]->InitialiseMoments();
	    //   ShypPrec[m][c]->InitialiseMoments();
	    // }
	    // std::cout<<std::endl;



      
	    // std::ublas::vector<double> hprec1(M);
	    // PARALLEL_TRANSFORM(ShypMean[m].begin(), ShypMean[m].end(), 
	    // 		       hprec1.begin(),
	    // 		       boost::bind(&Variable::GetMoments(), _1)
	    // 		       );
	    // std::ublas::vector<double> hprec2 = hprec1;
	    // PARALLEL_NEXT_PERMUTATION(hprec2.begin(), hprec2.end());
	    
	    // std::ublas::vector<double> hprec_sub;
      
	    // PARALLEL_TRANSFROM(hprec1.begin(), hprec1.end(),
	    // 		       hprec2.begin(), hprec_sub.begin(),
	    // 		       subtract()
	    // 		       );
	    
	    // for(size_t c=0;c<Components;++c){
	    //   if (hprec_sub<0.01)
	    // 	{

	    // 	  ShypMean[m][c]->InitialiseMoments();
	    // 	  ShypPrec[m][c]->InitialiseMoments();
	    // 	}
	    // }

	    // PARALLEL_SORT(weights1.begin(), weights1.end());
	    // weights2 = weights1;
	    // PARALLEL_NEXT_PERMUTATION(weights2.begin(), weights2.end());
	    // PARALLEL_TRANSFROM(weights1.begin(), weights1.end(),
	    // 		       weights2.begin(), weights_sub.begin(),
	    // 		       subtract()
	    // 		       );
	    // if (*min_element(weights_sub.begin(), weights_sub.end() < 0.01) 
			 
	    // 	for(size_t c=0;c<Components;++c){
	    // 	  InferredWeights<<std::exp(Weights[m]->GetMoments()[c])<<"\t";
	    // 	}
	    // 	InfHypMeans <<"\n";
	    // 	InfHypPrec  <<"\n";
	    // 	InferredWeights<<"\n";
	    // 	}
	    


	    // not_converged = true;
	  }
	}
	
	not_converged = !Build.run(0.1,200);

	// std::cout<<"WEIGHTS ["<<attempt<<"] = "<<std::endl;

	// for(size_t m=0;m<M;++m){
	//   for(size_t c=0;c<Components;++c){
	//     std::cout<<std::exp(Weights[m]->GetMoments()[c])<<"\t";
	//   }
	//   std::cout<<std::endl;
	//}

	std::ofstream InferredResult("ICAInferedResult.txt");
	for(size_t t=0;t<T;++t){
	  for(size_t n=0;n<N;++n){
	    InferredResult<<Mean(AtimesSplusN[n][t])<<"\t";
	  }
	  InferredResult<<"\n";

	}
	// mat tAInf(N,M);
	// mat tA2Inf(N,M);
	// vec tAVarianceScale(M);

	// for(size_t n=0;n<N;++n){
	//   for(size_t m=0;m<M;++m){
	//     tAInf(n,m)  = A[n][m]->GetMoments()[0];
	//     tA2Inf(n,m) = A[n][m]->GetMoments()[1];
	//     // AVarianceScale(m,t) = GSL_POW_2(SInf(m,t))/SInf2(m,t);
	//   }
	// }
	// //std::cout<<"A = "<<std::endl;
	// //std::cout<<tAInf<<std::endl;


	// for(size_t m=0;m<M;++m){
	//   const vec a =tAInf.get_col(m);
	//   tAVarianceScale(m) = mean(a*a)/mean(tA2Inf.get_col(m));
      
	//   if (tAVarianceScale(m)<0.9) 
	//     {
	//       std::cout<<"restarting "<<m<<std::endl;
	//       attempt =0;

	//       //restart
	//       for(size_t n=0;n<N;++n){
	// 	A[n][m]->InitialiseMoments();
	//       }
	//       for(size_t t=0;t<T;++t){
	//       	S[m][t]->InitialiseMoments();
	//       }
	//     }
	// }
	// std::cout<<"A variance scale = "<<tAVarianceScale<<"\n";
	++attempt;


	/* Fix weightings on mixing matrix */
	ublas::matrix<double> AInf(N,M);
	// mat A2Inf(N,M);
	// vec AVarianceScale(M);
    
	ublas::matrix<double> SInf(M,T);
	// mat SInf2(M,T);
	// mat SVarianceScale(M,T);
    
	for(size_t n=0;n<N;++n){
	  for(size_t m=0;m<M;++m){
	    AInf(n,m)  = Mean(A[n][m]);
	    // A2Inf(n,m) = A[n][m]->GetMoments()[1];
	    // AVarianceScale(m,t) = GSL_POW_2(SInf(m,t))/SInf2(m,t);
	  }
	}
    
    

	for(size_t m=0;m<M;++m){
	  for(size_t t=0;t<T;++t){
	    SInf(m,t) = Mean(S[m][t]);
	    // SInf2(m,t) = S[m][t]->GetMoments()[1];
	    // SVarianceScale(m,t) = GSL_POW_2(SInf(m,t))/SInf2(m,t);
	  }
	}
    
    
	for(size_t m=0;m<M;++m){
	  ublas::matrix_column<ublas::matrix<double> > col(AInf,m);
	  //ublas::vector<double> col = AInf.get_col(m);
	  ublas::vector<double> col2 = element_prod(col,col);
	  double norm = 1.0/std::sqrt(std::accumulate(col2.begin(), col2.end(), 0.0));
	  //double<double> norm = 1.0/std::sqrt(mean(col2));
	  col*=norm;
	  //AInf.set_col(m, col);
      
	  ublas::matrix_row<ublas::matrix<double> > row(SInf, m);
	  row/=norm;
	  // ublas::vector<double> row = SInf.get_row(m);
	  // row /= norm;
	  // SInf.set_row(m,row);
	}
    

	std::ofstream InferredA("ICAMixingMatrix.txt");
	for(size_t n=0;n<N;++n){
	  //ublas::vector<double> row = AInf.get_row(n);
	  
	  ublas::matrix_row<ublas::matrix<double> > row(AInf, n);
	  InferredA<<row<<"\n";
	}


	std::ofstream InferredSources("ICAInferedSources.txt");
	for(size_t t=0;t<T;++t){
	  InferredSources<< ublas::matrix_column<ublas::matrix<double> >(SInf,t) <<"\n";
	  //InferredSources<< SInf.get_col(t)<<"\n";
	  // InferredSources <<col<<"\n";
	  // InferredSources<<S[m][t]->GetMoments()[0]<<"\t";
	}
	// InferredSources<<"\n";


	std::ofstream InferredNoise("ICInferedNoise.txt");
	for(size_t n=0;n<N;++n){
	  InferredNoise<<Mean(noiseMean[n])<<"\t";
	}

	std::ofstream  InferredWeights("ICAWeights.txt");
	std::ofstream  InferredWeightsSD("ICAWeightsSD.txt");
	std::ofstream  InfHypMeans("ICAHypMean.txt");
	std::ofstream  InfHypPrec("ICAHypPrec.txt");
	for(size_t m=0;m<M;++m){
	  for(size_t c=0;c<Components;++c){
	    InfHypMeans<<Mean(ShypMean[m][c])<<"\t";
	    InfHypPrec <<Mean(ShypPrec[m][c])<<"\t";
	    InferredWeights<<Mean(Weights[m],c)<<"\t";
	    InferredWeightsSD<<StandardDeviation(Weights[m],c)<<"\t";
	  }
	  InfHypMeans <<"\n";
	  InfHypPrec  <<"\n";
	  InferredWeights<<"\n";
	  InferredWeightsSD<<"\n";
	}
    


	
      }
  

    //or if 

    // std::ofstream InferredA("ICAMixingMatrix.txt");
    // for(size_t n=0;n<N;++n){
    //   for(size_t m=0;m<M;++m){

    // 	InferredA<<A[n][m]->GetMoments()[0]<<"\t";
    //   }
    //   InferredA<<"\n";
    // }
    //    std::cout<<<<std::endl;

 

  }
}









// BOOST_AUTO_TEST_CASE( MixtureTree_test  )
// {
//   ICR::maths::rng random;
//   size_t data_points = 1000;//5447
//   typedef Builder<double>::GaussianNode GaussianNode;
//   typedef Builder<double>::WeightsNode WeightsNode;
//   typedef Builder<double>::GammaNode    GammaNode;
//   typedef Builder<double>::GaussianDataNode    Datum;
//   typedef Builder<double>::Variable    Variable;
   
//   Builder<double> Build;
//   const size_t Components = 1;
//   std::vector<Variable> HypMeanMixture(Components);
//   for(size_t i=0;i<Components;++i){
//     HypMeanMixture[i] = Build.Gaussian(0.001,0.001);
//   }

//   //GaussianNode HypMean = Build.Gaussian(0.001,0.001);
   
//   WeightsNode Weights    = Build.Weights(Components);
//   GaussianNode Mean      = Build.GaussianMixture(HypMeanMixture,0.001,Weights);
//   GammaNode    Precision = Build.Gamma(0.001,0.0001);
   
//   std::cout<<"Weights   = "<<Weights<<std::endl;
//   std::cout<<"Mean      = "<<Mean<<std::endl;
//   std::cout<<"Precision = "<<Precision<<std::endl;

   
//   for(size_t i=0;i<data_points;++i){
//     Build.Join(Mean, Precision, random.gaussian(1.0/(3),2));
//   }
//   Build.Iterate();
//   Build.Run(1e-6, 20);
   

//   //   BOOST_CHECK(Build.NumberOfNodes()   ==  (size_t) 6 + data_points );
//   //BOOST_CHECK(Build.NumberOfFactors() ==  (size_t) 2 + data_points);
   
//   // for(size_t i=0;i<10000;++i){
//   //   Build.Iterate();

//   // }


//   // double mean = Mean->GetMoments()[0];
//   // double prec = Precision->GetMoments()[0];
//   // std::cout<<"Gaussian mean = "<<mean<<std::endl;
//   // std::cout<<"Gaussian precission = "<<prec<<std::endl;
// }






// // BOOST_AUTO_TEST_CASE( GaussianFactor_test  )
// {
//   ConstantNode  HypShape(2.0);
//   ConstantNode  HypIScale(3.0);
//   GaussianNode     Scale;
    
//   GaussianFactor TEST(HypShape, HypIScale ,Scale);    
//   BOOST_CHECK_EQUAL(TEST.NumberOfNeighbours(), (size_t) 3);
//   BOOST_CHECK_EQUAL(TEST.IterationStarted(), false);
  
//   std::list<VariableNode*> Neighbours = TEST.GetNeighbours();
  
//   std::list<VariableNode*>::iterator it =  Neighbours.begin();
//   BOOST_CHECK(*it == &HypShape);
//   ++it;
//   BOOST_CHECK(*it == &HypIScale);
//   ++it;
//   BOOST_CHECK(*it == &Scale);
  
//   BOOST_CHECK_EQUAL(HypShape.NumberOfNeighbours(), (size_t)1);
//   BOOST_CHECK_EQUAL(HypIScale.NumberOfNeighbours(),(size_t) 1);
//   BOOST_CHECK_EQUAL(Scale.NumberOfNeighbours(), (size_t) 1);

//   std::list<FactorNode*> NeighboursHypShape = HypShape.GetNeighbours();
//   BOOST_CHECK(*NeighboursHypShape.begin() == &TEST);
// }




// BOOST_AUTO_TEST_SUITE( TestMessages_test )

// BOOST_AUTO_TEST_CASE( Gaussiandoubleest  )
// {
//   ConstantNode  HypShape(2.0);
//   ConstantNode  HypIScale(3.0);
//   GammaNode     Precision;

  
//   ConstantNode  HypMean(3.0);
//   ConstantNode  HypPrec(4.0);
//   GaussianNode  Mean;
  
  
//   DataNode<Model::Gaussian>  Data(5.0);
  
//   GammaFactor PrecisionF(HypShape, HypIScale , Precision);    
//   GaussianFactor MeanF(HypMean,HypPrec,Mean);
  
//   GaussianFactor DataF(Mean, Precision,Data);
    

//     boost::thread_group messages;
//     double Cost;
//     DataF.Iterate(messages);
//     messages.join_all();

//     //Data is aquired in following order
    
//     //Can Get these moments immediately (ConstantNode -> FactorNode)
//     //To PrecisionF
//     BOOST_CHECK_CLOSE(HypShape.GetMoments()[0], 2.0, 0.0001);
//     BOOST_CHECK_CLOSE(HypShape.GetMoments()[1], 0.0, 0.0001);
//     BOOST_CHECK_CLOSE(HypIScale.GetMoments()[0], 3.0, 0.0001);
//     BOOST_CHECK_CLOSE(HypIScale.GetMoments()[1], 0.0, 0.0001);

//     //To MeanF
//     BOOST_CHECK_CLOSE(HypMean.GetMoments()[0], 3.0, 0.0001);
//     BOOST_CHECK_CLOSE(HypMean.GetMoments()[1], 0.0, 0.0001);
//     BOOST_CHECK_CLOSE(HypPrec.GetMoments()[0], 4.0, 0.0001);
//     BOOST_CHECK_CLOSE(HypPrec.GetMoments()[1], 0.0, 0.0001);

//     //To DataF
//     BOOST_CHECK_CLOSE(Data.GetMoments()[0], 5.0, 0.0001);
//     BOOST_CHECK_CLOSE(Data.GetMoments()[1], 25.0, 0.0001);
  


//     //Can now send message from FactorNode -> VariableNode (This is the Natural Parameters)
//     //From Factor PrecisionF to Precision
//     BOOST_CHECK_CLOSE(PrecisionF.GetNaturalNot(Precision)[0], -HypIScale.GetMoments()[0], 0.0001);
//     BOOST_CHECK_CLOSE(PrecisionF.GetNaturalNot(Precision)[1], HypShape.GetMoments()[0]-1.0, 0.0001);
//     // //From Factor MeanF to Mean
//     BOOST_CHECK_CLOSE(MeanF.GetNaturalNot(Mean)[0], HypPrec.GetMoments()[0]*HypMean.GetMoments()[0], 0.0001); //=12
//     BOOST_CHECK_CLOSE(MeanF.GetNaturalNot(Mean)[1], -0.5*HypPrec.GetMoments()[0], 0.0001); // = -2.0

//     // Can now send message from VariableNode -> FactorNode (This is the moments)
//     //From Precision to DataF
//     BOOST_CHECK_CLOSE(Precision.GetMoments()[0], 0.0, 0.0001); //first go so moments not initialised yet!
//     BOOST_CHECK_CLOSE(Precision.GetMoments()[1], 0.0, 0.0001);

//     //From Mean to DataF
//     BOOST_CHECK_CLOSE(Mean.GetMoments()[0],0, 0.0001);//first go so moments not initialised yet!
//     BOOST_CHECK_CLOSE(Mean.GetMoments()[1],0, 0.0001);
    
//     //Can now send message from FactorNode -> VariableNode (This is the Natural Parameters)
//     //From Factor DataF to Precision
//     BOOST_CHECK_CLOSE(DataF.GetNaturalNot(Precision)[0], -0.5*Data.GetMoments()[0]*Data.GetMoments()[0], 0.0001); // -0.5*(x^2 - 2*x<mu> + <m^2>) = -12.5
//     BOOST_CHECK_CLOSE(DataF.GetNaturalNot(Precision)[1], 0.5, 0.0001); //0.5
//     // //From Factor DataF to Mean
//     BOOST_CHECK_CLOSE(DataF.GetNaturalNot(Mean)[0], 0.0, 0.0001); // <gamma>x
//     BOOST_CHECK_CLOSE(DataF.GetNaturalNot(Mean)[1], 0.0, 0.0001); // -<gamma>/2

//     //Now all variable nodes have both sets of messages - Can update and evaluate the cost.
//     DataF.EvaluateCost(Cost);  //this calls update
    
//     //Natual Parameters = 
//     /* Mean 
//      *  [ beta(4) m(3) ]   [ 0 ]    [ 12 ]
//      * =[ -beta(4)/2   ] + [ 0 ]  = [ -2 ]
//      *
//      * It follows that beta' = -2* (-2) = 4, m' = 12/beta' = 3;
//      * and so moments = [m', m'^2 + 1/beta' ] = [3, 9.25]
//      */
//     BOOST_CHECK_CLOSE(Mean.GetMoments()[0],3.0, 0.0001); //Moments are updated!
//     BOOST_CHECK_CLOSE(Mean.GetMoments()[1],9.25, 0.0001);

//     /* Precision 
//      *  [ -b(3)    ]   [ -12.5 ]    [ -15.5 ]
//      * =[ a(2)-1   ] + [ 0.5   ]  = [ 1.5  ]
//      *
//      * It follows that b' = -(-15.5) = 15.5, a' = 1.5+1 = 2.5;
//      * and so moments = [a'/b', digamma(a) - log(b) ] = [ 0.1612903 , -2.037683]
//      */
//     BOOST_CHECK_CLOSE(Precision.GetMoments()[0],0.1612903 , 0.0001); //Moments are updated!
//     BOOST_CHECK_CLOSE(Precision.GetMoments()[1],-2.037683, 0.0001);

//     //COST
//     //Data Cost
//     /*  Lx = <gamma><mu> x - 2*<gamma> x^2 + .5<log gamma> - <gamma><mu^2> - log(2*pi) : <gamma (on F)> = 0, <mu (on F)> = 0
//      *     = -infty
//      *
//      */
    
//     //COST
//     //Data Mean
//     /*  Lmu = (beta m - beta'm')<mu>  + (- 0.5*beta +  0.5*beta')<mu^2> + 0.5(log beta/beta' - beta m^2  + beta'm'^2) : <mu> = 3 <mu^2> = 9.25
//      *      =  4*3    - 4*3                        0                      0.5(log 1 - 0*4*9)  
//      *      = 0
//      */

//     //COST
//     //Data gamma
//     /*  Lmu = (b'-b)<gamma> + (a - a') <log gamma> + a log b -a' log b' - log GAMMA(a)/GAMMA(a') 
//      *      =  -12.5            0.5                         0                      0.5(log 1 - 0*4*9)  
//      *      = 0
//      */


//     // BOOST_CHECK_CLOSE(MeanF.GetNaturalNot(Mean)[0], HypPrec.GetMoments()[0]*HypMean.GetMoments()[0], 0.0001);
//     // BOOST_CHECK_CLOSE(MeanF.GetNaturalNot(Mean)[1], -0.5*HypPrec.GetMoments()[0], 0.0001);
    
//     // BOOST_CHECK_CLOSE(Mean.GetMoments()[0],HypMean.GetMoments()[0], 0.0001);
//     // BOOST_CHECK_CLOSE(Mean.GetMoments()[1],HypMean.GetMoments()[0]*HypMean.GetMoments()[0]+ 1.0/HypPrec.GetMoments()[0], 0.0001);
  



//     //BOOST_CHECK_CLOSE(Precision.GetMoments()[0],HypShape.GetMoments()[0]/HypIScale.GetMoments()[0] , 0.0001);
//     //BOOST_CHECK_CLOSE(Precision.GetMoments()[1],gsl_sf_psi(HypShape.GetMoments()[0]) - std::log10(HypIScale.GetMoments()[0]), 0.0001);




//   std::cout<<"FINISH"<<std::endl;

// }
//   // NaturalParameters NotShape = 
//   //   TEST.GetNaturalNot(HypShape);
//   // NaturalParameters NotHypIScale = 
//   //   PrecisionF.GetNaturalNot(HypIScale);
//   // NaturalParameters NotPrecision = 
//   //   PrecisionF.GetNaturalNot(Precision);
  
//   // Moments MomPrec = Precision.GetMoments();
  

//   // NaturalParameters NotHypMean = 
//   //   MeanF.GetNaturalNot(HypMean);
//   // NaturalParameters NotHypPrec = 
//   //   MeanF.GetNaturalNot(HypPrec);
//   // NaturalParameters NotMean = 
//   //   MeanF.GetNaturalNot(Mean);
   
//   // NaturalParameters NotData = 
//   //   DataF.GetNaturalNot(Data);
//   // NaturalParameters DataNotPrecision = 
//   //   DataF.GetNaturalNot(Precision);
//   // NaturalParameters DataNotMean = 
//   //   DataF.GetNaturalNot(Mean);


//   // //Data is aquired in following order
  
//   // //Can Get these moments immediately (ConstantNode -> FactorNode)
//   // //To PrecisionF
//   // BOOST_CHECK_CLOSE(HypShape.GetMoments()[0], 2.0, 0.0001);
//   // BOOST_CHECK_CLOSE(HypShape.GetMoments()[1], 0.0, 0.0001);
//   // BOOST_CHECK_CLOSE(HypIScale.GetMoments()[0], 3.0, 0.0001);
//   // BOOST_CHECK_CLOSE(HypIScale.GetMoments()[1], 0.0, 0.0001);

//   // //To MeanF
//   // BOOST_CHECK_CLOSE(HypMean.GetMoments()[0], 3.0, 0.0001);
//   // BOOST_CHECK_CLOSE(HypMean.GetMoments()[1], 0.0, 0.0001);
//   // BOOST_CHECK_CLOSE(HypPrec.GetMoments()[0], 4.0, 0.0001);
//   // BOOST_CHECK_CLOSE(HypPrec.GetMoments()[1], 0.0, 0.0001);

//   // //To DataF
//   // BOOST_CHECK_CLOSE(Data.GetMoments()[0], 0.0, 0.0001);
//   // BOOST_CHECK_CLOSE(Data.GetMoments()[1], 0.0, 0.0001);
  
//   // //Can now send message from FactorNode -> VariableNode (This is the Natural Parameters)
//   // //From Factor PrecisionF to Precision
//   // BOOST_CHECK_CLOSE(NotPrecision[0], -HypIScale.GetMoments()[0], 0.0001);
//   // BOOST_CHECK_CLOSE(NotPrecision[1], HypShape.GetMoments()[0]-1.0, 0.0001);
//   // //From Factor MeanF to Mean
//   // BOOST_CHECK_CLOSE(NotMean[0], HypMean.GetMoments()[0]*HypPrec.GetMoments()[0], 0.0001);
//   // BOOST_CHECK_CLOSE(NotMean[1], -0.5* HypPrec.GetMoments()[0], 0.0001);

//   // // Can now send message from VariableNode -> FactorNode (This is the moments)
//   // //From Precision to DataF
//   // //BOOST_CHECK_CLOSE(Precision.GetMoments()[0],HypShape.GetMoments()[0]/HypIScale.GetMoments()[0] , 0.0001);
//   // //BOOST_CHECK_CLOSE(Precision.GetMoments()[1],gsl_sf_psi(HypShape.GetMoments()[0]) - std::log10(HypIScale.GetMoments()[0]), 0.0001);

//   // //From Mean to DataF
//   // BOOST_CHECK_CLOSE(Mean.GetMoments()[0],HypMean.GetMoments()[0], 0.0001);
//   // BOOST_CHECK_CLOSE(Mean.GetMoments()[1],HypMean.GetMoments()[0]*HypMean.GetMoments()[0]+ 1.0/HypPrec.GetMoments()[0], 0.0001);
  


// // BOOST_AUTO_TEST_CASE( Builder_test  )
// // {
// //   ConstantNode  HypShape(2.0);
// //   ConstantNode  HypIScale(3.0);
// //   GammaNode     Precision;

  
// //   ConstantNode  HypMean(3.0);
// //   ConstantNode  HypPrec(4.0);
// //   GaussianNode  Mean;

// //   std::cout<<"MEAN = "<<&Mean<<std::endl;
// //   std::cout<<"PREC = "<<&Precision<<std::endl;


// //   DataNode<Model::Gaussian>  Data(0.0);
  
// //   Builder B;

// //   B.Join(HypShape, HypIScale , Precision);
// //   B.Join(HypMean,HypPrec,Mean);
// //   B.Join(Mean, Precision,Data);
// //   for(size_t i=0;i<10;++i){
// //     std::cout<<"ITERATE"<<std::endl;

// //     B.Iterate();
// //   }
  
  
// //  }


// BOOST_AUTO_TEST_CASE( Builder_two_test  )
// {
//   ConstantNode  HypShape(0.001);
//   ConstantNode  HypIScale(0.001);
//   GammaNode     Precision;

  
//   ConstantNode  HypMean(0.001);
//   ConstantNode  HypPrec(0.001);
//   GaussianNode  Mean;

//   std::cout<<"MEAN = "<<&Mean<<std::endl;
//   std::cout<<"PREC = "<<&Precision<<std::endl;

//   size_t data_points = 1000;
//   typedef boost::shared_ptr<DataNode<Model::Gaussian> > datum;
//   std::vector< datum >  Data(data_points);
//   ICR::maths::rng random;
//   for(size_t i=0;i<data_points;++i){
//     Data[i]   = datum(new DataNode<Model::Gaussian>(random.gaussian(1.0/(3),2),0));
//   }

//   typedef Builder::GaussianNode GaussianNode;
//   typedef Builder::GammaNode    GammaNode;
  
//   Builder Build;
//   GaussianNode HypMean      = Build.Gaussian(0.001,0.001);
//   GammaNode    HypPrecision = Build.Gaussian(0.001,0.001);
//   GaussianNode Mean         = Build.Gaussian(HypMean,HypPrecision);
//   GammaNode    Precision    = Build.Gamma(0.001,0.001);
//   Build.Data<Model::Gaussian>(Mean,Precision, random.gaussian(1.0/(3),2),0);

//   Build.Iterate();

//   DataNode<Model::Gaussian>& datum;
//   std::vector< datum >  Data;
//   std::copy(
//   for_each(Data.begin(), Data.end(), boost::bind(&Builder::Join, boost::ref(B), Mean, Precision, _1) );


//   B.Join(HypShape, HypIScale , Precision);
//   B.Join(HypMean,HypPrec,Mean);
//   for(size_t i=0;i<data_points;++i){
//     B.Join(Mean, Precision,*(Data[i]));
//   }

//   for(size_t i=0;i<50;++i){
//     std::cout<<"ITERATE"<<std::endl;
    
//     B.Iterate();
//   }
  
//   double mean = Mean.GetMoments()[0];
//   double prec = Precision.GetMoments()[0];
//   std::cout<<"Gaussian mean = "<<mean<<std::endl;
//   std::cout<<"Gaussian precission = "<<prec<<std::endl;
//   //std::cout<<"Gaussian mean = "<<Data.variance()<<std::endl;


  
// }




// // BOOST_AUTO_TEST_CASE( Gaussiandoubleest  )
// // {
// //   ConstantNode  HypShape(2.0);
// //   ConstantNode  HypIScale(3.0);
// //   GammaNode     Precision;

  
// //   ConstantNode  HypMean(3.0);
// //   ConstantNode  HypPrec(4.0);
// //   GaussianNode  Mean;
  
  
// //   DataNode<Model::Gaussian>  Data(0.0);
  
// //   GammaFactor PrecisionF(HypShape, HypIScale , Precision);    
// //   GaussianFactor MeanF(HypMean,HypPrec,Mean);
  
// //   GaussianFactor DataF(Mean, Precision,Data);
    
// //   for(size_t i=0;i<2;++i){

// //     boost::thread_group messages;
// //     double Cost;
// //     DataF.Iterate(messages);
// //     messages.join_all();
// //     DataF.EvaluateCost(Cost);
// //     std::cout<<"RESET"<<std::endl;

// //     DataF.Reset();
// //     std::cout<<HypShape.IterationStarted()<<std::endl;
// //     std::cout<<HypIScale.IterationStarted()<<std::endl;
// //     std::cout<<Precision.IterationStarted()<<std::endl;
// //     std::cout<<HypMean.IterationStarted()<<std::endl;
// //     std::cout<<HypPrec.IterationStarted()<<std::endl;
// //     std::cout<<Mean.IterationStarted()<<std::endl;
// //     std::cout<<Data.IterationStarted()<<std::endl;
// //     std::cout<<PrecisionF.IterationStarted()<<std::endl;
// //     std::cout<<MeanF.IterationStarted()<<std::endl;
// //     std::cout<<DataF.IterationStarted()<<std::endl;
// //     std::cout<<"Cost = "<<Cost<<std::endl;

// //   }
// //   std::cout<<"FINISH"<<std::endl;

// //   // NaturalParameters NotShape = 
// //   //   TEST.GetNaturalNot(HypShape);
// //   // NaturalParameters NotHypIScale = 
// //   //   PrecisionF.GetNaturalNot(HypIScale);
// //   // NaturalParameters NotPrecision = 
// //   //   PrecisionF.GetNaturalNot(Precision);
  
// //   // Moments MomPrec = Precision.GetMoments();
  

// //   // NaturalParameters NotHypMean = 
// //   //   MeanF.GetNaturalNot(HypMean);
// //   // NaturalParameters NotHypPrec = 
// //   //   MeanF.GetNaturalNot(HypPrec);
// //   // NaturalParameters NotMean = 
// //   //   MeanF.GetNaturalNot(Mean);
   
// //   // NaturalParameters NotData = 
// //   //   DataF.GetNaturalNot(Data);
// //   // NaturalParameters DataNotPrecision = 
// //   //   DataF.GetNaturalNot(Precision);
// //   // NaturalParameters DataNotMean = 
// //   //   DataF.GetNaturalNot(Mean);


// //   // //Data is aquired in following order
  
// //   // //Can Get these moments immediately (ConstantNode -> FactorNode)
// //   // //To PrecisionF
// //   // BOOST_CHECK_CLOSE(HypShape.GetMoments()[0], 2.0, 0.0001);
// //   // BOOST_CHECK_CLOSE(HypShape.GetMoments()[1], 0.0, 0.0001);
// //   // BOOST_CHECK_CLOSE(HypIScale.GetMoments()[0], 3.0, 0.0001);
// //   // BOOST_CHECK_CLOSE(HypIScale.GetMoments()[1], 0.0, 0.0001);

// //   // //To MeanF
// //   // BOOST_CHECK_CLOSE(HypMean.GetMoments()[0], 3.0, 0.0001);
// //   // BOOST_CHECK_CLOSE(HypMean.GetMoments()[1], 0.0, 0.0001);
// //   // BOOST_CHECK_CLOSE(HypPrec.GetMoments()[0], 4.0, 0.0001);
// //   // BOOST_CHECK_CLOSE(HypPrec.GetMoments()[1], 0.0, 0.0001);

// //   // //To DataF
// //   // BOOST_CHECK_CLOSE(Data.GetMoments()[0], 0.0, 0.0001);
// //   // BOOST_CHECK_CLOSE(Data.GetMoments()[1], 0.0, 0.0001);
  
// //   // //Can now send message from FactorNode -> VariableNode (This is the Natural Parameters)
// //   // //From Factor PrecisionF to Precision
// //   // BOOST_CHECK_CLOSE(NotPrecision[0], -HypIScale.GetMoments()[0], 0.0001);
// //   // BOOST_CHECK_CLOSE(NotPrecision[1], HypShape.GetMoments()[0]-1.0, 0.0001);
// //   // //From Factor MeanF to Mean
// //   // BOOST_CHECK_CLOSE(NotMean[0], HypMean.GetMoments()[0]*HypPrec.GetMoments()[0], 0.0001);
// //   // BOOST_CHECK_CLOSE(NotMean[1], -0.5* HypPrec.GetMoments()[0], 0.0001);

// //   // // Can now send message from VariableNode -> FactorNode (This is the moments)
// //   // //From Precision to DataF
// //   // //BOOST_CHECK_CLOSE(Precision.GetMoments()[0],HypShape.GetMoments()[0]/HypIScale.GetMoments()[0] , 0.0001);
// //   // //BOOST_CHECK_CLOSE(Precision.GetMoments()[1],gsl_sf_psi(HypShape.GetMoments()[0]) - std::log10(HypIScale.GetMoments()[0]), 0.0001);

// //   // //From Mean to DataF
// //   // BOOST_CHECK_CLOSE(Mean.GetMoments()[0],HypMean.GetMoments()[0], 0.0001);
// //   // BOOST_CHECK_CLOSE(Mean.GetMoments()[1],HypMean.GetMoments()[0]*HypMean.GetMoments()[0]+ 1.0/HypPrec.GetMoments()[0], 0.0001);
  

// //   // Can Now Update DataNode



// //   // //Can now send message from FactorNode -> VariableNode
// //   // //From  to DataF to Data;
  
  
// //   // //Can now pass message to Precision
  
// //   // //from HypMean and HypPrecision
// //   // BOOST_CHECK_CLOSE(NotMean[0], , 0.0001);
// //   // BOOST_CHECK_CLOSE(NotMean[1], -0.5*3.0, 0.0001);
// //   // //Can now pass message to Mean

  

// //   // //Can now pass message to Precision
// //   // //from HypShape and HypIScale
// //   // BOOST_CHECK_CLOSE(NotPrecision[0], -3.0, 0.0001);
// //   // BOOST_CHECK_CLOSE(NotPrecision[1], 2.0-1.0, 0.0001);

  
  

// //   // // BOOST_CHECK_CLOSE(NotShape[0], , 0.0001);
// //   // // BOOST_CHECK_CLOSE(NotShape[1], , 0.0001);
// //   // BOOST_CHECK_CLOSE(NotIScale[0], , 0.0001); //- child mean (Scale)
// //   // BOOST_CHECK_CLOSE(NotIScale[1], 2.0, 0.0001);

// //   // BOOST_CHECK_CLOSE(NotHypMean[0],  , 0.0001); 
// //   // BOOST_CHECK_CLOSE(NotHypMean[1], , 0.0001);
// //   // BOOST_CHECK_CLOSE(NotHypPrec[0], , 0.0001);
// //   // BOOST_CHECK_CLOSE(NotHypPrec[1], , 0.0001);
// //   // BOOST_CHECK_CLOSE(NotMean[0], 3.0*5.0, 0.0001);
// //   // BOOST_CHECK_CLOSE(NotMean[1], -0.5*3.0, 0.0001);

  
// //   // BOOST_CHECK_CLOSE(NotData[0], 3.0*5.0, 0.0001);
// //   // BOOST_CHECK_CLOSE(NotData[1], -0.5*3.0, 0.0001);




// //   // BOOST_CHECK_CLOSE(NotData[0], -3  , 0.0001);
// //   // BOOST_CHECK_CLOSE(NotData[1], 1 , 0.0001);


  
// // }

// // BOOST_AUTO_TEST_CASE( GaussianGNN_test  )
// // {
// //   ConstantNode  mean(2.0);
// //   ConstantNode  precision(3.0);
// //   ConstantNode  data(5.0);
    
// //   GaussianFactor TEST(mean, precision ,data);    

// //   boost::thread_group messages;
// //   TEST.Iterate(messages);
// //   messages.join_all();
    
// //    NaturalParameters NotMean = 
// //      TEST.GetNaturalNot(mean);
// //   NaturalParameters NotPrec = 
// //     TEST.GetNaturalNot(precision);
// //   NaturalParameters NotData = 
// //     TEST.GetNaturalNot(data);
  
// //   BOOST_CHECK_CLOSE(NotMean[0], 3.0*5.0, 0.0001);
// //   BOOST_CHECK_CLOSE(NotMean[1], -0.5*3.0, 0.0001);
// //   BOOST_CHECK_CLOSE(NotPrec[0], , 0.0001);
// //   BOOST_CHECK_CLOSE(NotPrec[1], , 0.0001);
// //   BOOST_CHECK_CLOSE(NotData[0], -3  , 0.0001);
// //   BOOST_CHECK_CLOSE(NotData[1], 1 , 0.0001);
  
// // }

// BOOST_AUTO_TEST_SUITE_END()





// // BOOST_AUTO_TEST_SUITE( ConstantNode_test )

// // BOOST_AUTO_TEST_CASE( Constructor_test  )
// // {

// //   ConstantNode mu(2.0);
// //   ConstantNode ivar(3.0);
// //   ConstantNode data(5.0);
  
// // }
// // BOOST_AUTO_TEST_CASE( GetMoments_test  )
// // {

// //   ConstantNode mu(2.0);
// //   Moments moments = mu.GetMoments();
// //   BOOST_CHECK_CLOSE(moments[0], 2.0, 0.0001)
// //   BOOST_CHECK_CLOSE(moments[1], 0.0, 0.0001)
    
// // }

// // // BOOST_AUTO_TEST_CASE( AddFactor_test  )
// // // {

// // //     Variable  HypShape  (new ConstantNode(2.0) );
// // //     Variable  HypIScale (new ConstantNode(3.0) );
// // //     Variable  Scale  (new GammaNode() ); 
    
// // //     GammaFactor TEST(HypShape.get(),HypIScale.get(),Scale.get());
    
    
// // //   BOOST_CHECK_EQUAL(HypShape.NumberOfNeighbours, 1)
// // //   BOOST_CHECK_EQUAL(HypIScale.NumberOfNeighbours, 1)
// // //   BOOST_CHECK_EQUAL(Scale.NumberOfNeighbours, 1)

// // // }




// // BOOST_AUTO_TEST_SUITE_END()


// // BOOST_AUTO_TEST_SUITE( GammaFactor_test )





// // BOOST_AUTO_TEST_CASE( Constructor_test  )
// // {

// //   typedef boost::shared_ptr<VariableNode> Variable;
// //   typedef boost::shared_ptr<FactorNode>   Factor;
  
// //     Variable  HypShape  (new ConstantNode(2.0) );
// //     Variable  HypIScale (new ConstantNode(3.0) );
// //     Variable  Scale  (new GammaNode() ); 
    
// //     GammaFactor TEST(HypShape.get(),HypIScale.get(),Scale.get());
// //     boost::mutex m_mutex;
// //     boost::condition m_cond;
// //     std::cout<<"SIZE OF mutex = "  <<sizeof(m_mutex)<<std::endl;
// //     std::cout<<"SIZE OF condition = "  <<sizeof(m_cond)<<std::endl;

    
// //     Factor hypGF(new GammaFactor(HypShape.get(),HypIScale.get(),Scale.get()) );
    
// //     Variable  Shape  (new ConstantNode(5.0) ); 
// //     Variable  data(new ConstantNode(5.0) );
// //     Factor GF(new GammaFactor(Shape.get(),Scale.get(),data.get()) );

// //     //   for(size_t i=0;i<10;++i){
// //       boost::thread_group messages;
// //       GF->Iterate(messages);
// //       std::cout<<"messages = "<<messages.size()<<std::endl;
// //       messages.join_all();
// //       GF->Reset();
// //       // }
// //     std::cout<<"FINISHED"<<std::endl;

// //     // GF.stop();

// //     //   GaussianFactor GF(mu, ivar);
// //     //   GF.AttachChild(data);
  
// //   //   GF.Init();
  
// // }


// // // BOOST_AUTO_TEST_CASE( Constructor_two_test  )
// // // {

// // ConstantNode HypShape(2.0);
// // ConstantNode HypIScale(3.0);
// // GammaNode    Scale;

// // ConstantNode mean(3);
// // GammaNode    precision;

// // DataNode<Gaussian>  data(4);

// // Builder b;
// // b.join(HypShape,HypIScale, Scale) 
// // b.join(mean, precision, data);
// // b.iterate();


// // //   typedef boost::shared_ptr<VariableNode> Variable;
// // //   typedef boost::shared_ptr<FactorNode>   Factor;
  
// // //     Variable  HypShape  (new ConstantNode(2.0) );
// // //     Variable  HypIScale (new ConstantNode(3.0) );
// // //     Variable  Scale  (new GammaNode() ); 
    
// // //     Factor hypGammaF(new GammaFactor(HypShape.get(),HypIScale.get(),Scale.get()) );
    
// // //     Variable  Shape     (new ConstantNode(5.0) ); 
// // //     Variable  Precision (new GammaNode() ); 
// // //     Factor GammaF(new GammaFactor(Shape.get(),Scale.get(),Precision.get()) );
    
// // //     Variable  Mean(new ConstantNode(3.0) );
// // //     Variable  Data(new ConstantNode(6.5));
// // //     Factor GaussianF(new GaussianFactor(Mean.get(), Precision.get(), Data.get()));
    
    
    
    
// // //     // for(size_t i=0;i<10;++i){
// // //       boost::thread_group messages;
// // //       GaussianF->Iterate(messages);
// // //       std::cout<<"messages = "<<messages.size()<<std::endl;
// // //       messages.join_all();
// // //       GaussianF->Reset();
// // //       //}
// // //     std::cout<<"FINISHED"<<std::endl;

// // //     // GF.stop();

// // //     //   GaussianFactor GF(mu, ivar);
// // //     //   GF.AttachChild(data);
  
// // //   //   GF.Init();
  
// // // }



// // BOOST_AUTO_TEST_SUITE_END()


// // //BOOST_AUTO_TEST_SUITE_END()


// // // BOOST_AUTO_TEST_SUITE( GaussianFactor_test )

// // // BOOST_AUTO_TEST_CASE( Constructor_test  )
// // // {

// // //   typedef boost::shared_ptr<node> Node;
  
// // //   Node mu  (new ConstantNode(2.0) );
// // //   Node ivar(new ConstantNode(3.0) );
// // //   Node data(new ConstantNode(5.0) );
  
// // //   GaussianFactor GF(mu, ivar);
// // //   GF.AttachChild(data);
  
// // //   GF.Init();
  
// // // }



// // //BOOST_AUTO_TEST_SUITE_END()
