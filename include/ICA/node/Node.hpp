#pragma once

#include "ICA/message/Coster.hpp"

//#include "ICA/Message.hpp"
#include <boost/shared_ptr.hpp>

#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp>

#include <iostream>
#include <list>
#include <vector>

namespace ICR{
  namespace ICA{
    
    //forward declaration
    template<class T> class Moments;
    template<class T> class NaturalParameters;
    
    // template<class T=double>
    // class LearningNode 
    // {
    // public:
    //   virtual 
    //   ~LearningNode(){};
      
    // };

    // template<class T=double>
    // class CanMakeMixture{
    // public:
    //   virtual 
    //   T
    //   GetAvLog() const = 0;
    // };

    //forward
    template<class T=double>
    class FactorNode;
    
    template<class T=double>
    class VariableNode  // : public Node
    {
    public:
      
      virtual
      const Moments<T>&
      GetMoments() const = 0;

      virtual
      void
      InitialiseMoments()  = 0;

      virtual
      void
      SetParentFactor(FactorNode<T>* f) = 0;
      
      virtual
      void
      AddChildFactor(FactorNode<T>* f) = 0;
      
      virtual 
      void
      Iterate(Coster&) = 0;

      virtual 
      ~VariableNode(){};
    };

    template<class T>
    std::ostream&
    operator<<(std::ostream& out, const VariableNode<T>* v)
    {
      printf("NODE: %p \t", v );
      out<<v->GetMoments();
      return out;
    }

    template<class T>
    class FactorNode 
    {
    public:
      
      typedef typename boost::call_traits< VariableNode<T>* const>::param_type
      variable_parameter;
      typedef typename boost::call_traits< VariableNode<T>* const>::value_type
      variable_t;

      virtual
      NaturalParameters<T>
      GetNaturalNot( variable_parameter ) const = 0;
 
      virtual 
      T
      CalcLogNorm() const = 0;
      
      virtual
      Moments<T>
      InitialiseMoments() const  = 0;
      // virtual 
      // T
      // GetFFunction() const = 0;

      // virtual 
      // void
      // Iterate() = 0;

      virtual 
      ~FactorNode(){};
    };
      
  }
}
