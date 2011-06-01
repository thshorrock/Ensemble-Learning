#ifndef NODE_HPP
#define NODE_HPP

#include "EnsembleLearning/message/Coster.hpp"
#include <boost/call_traits.hpp>
#include <vector>
#include <iostream>

namespace ICR{
  namespace EnsembleLearning{
    
    //forward declaration
    template<class T> class Moments;
    template<class T> class NaturalParameters;
    
    //forward
    template<class T>
    class FactorNode;
    
    /** An interface for the Variable Nodes.
     * Every variable node derives from this interface.
     * @tparam T The data type used by the variable nodes in the model.
     * This is either double or float.
     */
    template<class T>
    class VariableNode  // : public Node
    {
    public:
      
      /** Get the moments stored in this node. 
       *  @return The stored moments */
      virtual
      const Moments<T>&
      GetMoments() const = 0;

      /** Initialise the moments to this node.
       *  The moments are random based upon the prior moments (of the parent nodes)*/
      virtual
      void
      InitialiseMoments()  = 0;

      /** Set the parent factor for a node.
       *  Every variable node has only one parent factor, 
       *  which composes the messages from parent nodes.
       *  For example, the parent factor for Gaussian distributed data will be 
       *  a Gaussian Factor, whose parents in turn are variable nodes that store the mean and precision.
       *  @param f A pointer to the factor node that is the sole parent
       */
      virtual
      void
      SetParentFactor(FactorNode<T>* f) = 0;
      
      /** Add a child factor to a node.
       *  Every variable can have many child factors
       *  to which the moments are passed.
       *  For example, an observed node that stores Gaussian moment with zero mean
       *  might be a parent node to many Gaussian factors, each of which are
       *  added with this function.
       *  @param f A pointer to the factor node that is a child.
       */
      virtual
      void
      AddChildFactor(FactorNode<T>* f) = 0;
      
      /** Collect a vector of means from the node.
       *  @return A vector containing all the means of the node.
       *  Typically this vector will contain only one element,
       *   but Dirichlet models contain the mean for every component of the mixture model.
       */
      virtual
      const std::vector<T>
      GetMean()  = 0;
      
      /** Collect a vector of variances from the node.
       *  @return A vector containing all the variances of the node.
       *  Typically this vector will contain only one element,
       *   but Dirichlet models contain the variance for every component of the mixture model.
       */
      virtual
      const std::vector<T>
      GetVariance()  = 0;

      /** Collect the messages of adjacent factors and update the stored moments.
       *  @param Cost The total cost of the approximation to which this node contributes.
       */
      virtual 
      void
      Iterate(Coster& Cost) = 0;

      /** Destructor. */
      virtual 
      ~VariableNode(){};
    };

    /** Output the moments from a pointer to a variable node to a stream
     * @param out The output stream
     * @param v The pointer to the variable node.
     * @return A reference to the output stream
     */
    template<class T>
    std::ostream&
    operator<<(std::ostream& out, const VariableNode<T>* v)
    {
      printf("NODE: %p \t", v );
      out<<v->GetMoments();
      return out;
    }

    /** The interface to the factor nodes.
     *  Every factor node derives from this.
     * @tparam T The data type used by all the factor nodes, either double or float.
     */
    template<class T>
    class FactorNode 
    {
    public:
      
      typedef typename boost::call_traits< VariableNode<T>* const>::param_type
      variable_parameter;
      typedef typename boost::call_traits< VariableNode<T>* const>::value_type
      variable_t;

      /** Obtain the natural parameter destined for the variable_parameter v.
       * @param v A pointer to the  VariableNode for which the message is destined.
       *  The message is calculated from the moments of every node adjacent to the factor withe exception of v.
       * @return The natural parameter calculated for v.
       */
      virtual
      NaturalParameters<T>
      GetNaturalNot( variable_parameter v) const = 0;
 
      /** Calculate natural logarithm of the Models nomalisation constant required by the child node to evaluate the model's cost.
       * @return The natural logarithm of the Models nomalisation constant.
       */
      virtual 
      T
      CalcLogNorm() const = 0;
      
      /** Calculate the initial moments for the child node.
       *  These are evaluated from Random samples from the model's distribution,
       * based upon the values of the parent nodes.
       *  @return The Inital Moments for the child nodes.
       */
      virtual
      Moments<T>
      InitialiseMoments() const  = 0;

      /** Destructor */
      virtual 
      ~FactorNode(){};
    };
      
  }
}

#endif //NODE_HPP
