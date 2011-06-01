#pragma once


namespace ICR{

  namespace ICA{
    
    //Forward declarations
    template<class>  class SubContext;
    template<class>  class FunctionIterator;
    
    template<class T>
    class Plus;
    
    template<class T>
    class Times;
    
    template<class T>
    class Placeholder;
    
    template<class> class Function;
    
    /** A struct containing the Function Names.
     *  These couple the functions together and so will 
     *  hopefully be removed in a future version.
     */
    struct FunctionName{
      /** An enumeration indicating the type of function. */
      enum Value{
	PLUS,  //<!-- Addition function specifier.
	TIMES, //<!-- Mutliplication function specifier.
      };
    };
    
    /** An interface for every expression (and sub expression).
     *  @tparam T The data type: float or double.
     */
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
      
      /** A destructor */
      virtual ~Expression(){};
      
      /** Evaluate the expression for a given context.
       * @param C The context from which to evaluate the expression.
       * @return  The result of the expression.
       */
      virtual data_t Evaluate(subcontext_parameter C) const = 0;

    protected:
      friend class Plus<T>;
      friend class Times<T>;
      friend class Placeholder<T>;
      template<class> friend  class FunctionIterator;
      
      virtual void SetParent(function_parameter) = 0;
      virtual function_t GetParent() const = 0;
    };

    /** An interface to Function expressions.
     * Function expressions are operators - Plus, Times.
     * The interface returns the function type.
     * In a future version this will be done by the compiler directly
     * without the need for the function call.
     */
    template<class T>
    class Function : public Expression<T>
    {
    public:
      /** A destructor */
      virtual ~Function(){};
    
      /** Get the function type.
       *  @return A member of the enumeration FunctionName::Value 
       * that indicates the function type.
       */
      virtual
      FunctionName::Value
      GetFunctionType() const = 0;
    };
    

  }
}
