#pragma once


namespace ICR{

  namespace ICA{
    
    template<class>  class SubContext;

    struct FunctionName{
      enum Value{
	PLUS,
	TIMES,
      };
    };
    
    template<class T>
    class Plus;
    
    template<class T>
    class Times;
    
    template<class T>
    class Placeholder;
    
    template<class> class Function;
    
    template<class T>
    class Expression
    {
    public:
      typedef typename boost::call_traits<Function<T>*>::param_type
      function_parameter;
      
      typedef typename boost::call_traits<Function<T>*>::value_type
      function_t;

      typedef typename boost::call_traits<SubContext<T> >::param_type
      subcontext_parameter;

      typedef typename boost::call_traits<T>::value_type
      data_t;
      
      virtual ~Expression(){};
      
      virtual data_t Evaluate(subcontext_parameter) const = 0;

    protected:
      friend class Plus<T>;
      friend class Times<T>;
      friend class Placeholder<T>;
      
      
      virtual void SetParent(function_parameter) = 0;
      virtual function_t GetParent() const = 0;
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
