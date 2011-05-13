#pragma once



#include<string>

namespace ICR{

  namespace ICA{
    
    template<class>  class SubContext;

    struct FunctionName{
      enum Value{
	PLUS,
	TIMES,
      };
    };
    
    template<class> class Function;
    
    template<class T>
    class Expression
    {
    public:
      virtual ~Expression(){};
      
      virtual void SetParent(const Function<T>*) = 0;
      virtual const Function<T>* GetParent() const = 0;
      virtual T Evaluate(const SubContext<T>&) const = 0;
    };

    template<class T>
    class Function : public Expression<T>
    {
    public:
      virtual ~Function(){};
      
      virtual
      FunctionName::Value
	GetFunctionType() const = 0;
    };
    

  }
}
