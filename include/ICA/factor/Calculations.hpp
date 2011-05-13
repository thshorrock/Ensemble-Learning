#pragma once



#include "ICA/Moments.hpp"


namespace ICR{
  namespace ICA{
    
    template<class Model, class T>
    class Add
    {
    public:
      static
      NaturalParameters<T>
      CalcNP2First(const Moments& second,const Moments& child)
      {
	
	//Undo the add to parent
	return NaturalPamameters<T>(child-second);	
      }
      
      static
      NaturalParameters<T>
      CalcNP2Data(const Moments& first,const Moments& second)
      {
	
	//Undo the add to parent
	return NaturalPamameters<T>(first-second);	
      }
      
    private:

    };

  }
}
