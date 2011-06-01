#ifndef OUTPUT_HPP
#define OUTPUT_HPP


#include "ICA/node/Node.hpp"
#include <cmath>  
namespace ICR{
  namespace ICA{
    
    /** Evaluate the mean of a variable node. 
     * @param node A pointer to the VariableNode to evaluate.
     * @param index The index of the mean to return.
     *   Typically there is only one mean in the node (the index defaults to zero),
     *   however, Dirichlet models contain a mean for every component of the mixture.
     * @return The mean for the given node and index.
     * @ingroup UserInterface
     */
    template<class T>
    T
    Mean(VariableNode<T>* node, size_t index = 0)
    {
      return node->GetMean()[index];
    }
    /** Evaluate the standard deviation of a variable node. 
     * @param node A pointer to the VariableNode to evaluate.
     * @param index The index of the mean to return.
     *   Typically there is only one mean in the node (the index defaults to zero),
     *   however, Dirichlet models contain a mean for every component of the mixture.
     * @return The variance for the given node and index.
     * @ingroup UserInterface
     */
    template<class T>
    T
    StandardDeviation(VariableNode<T>* node, size_t index = 0)
    {
      return std::sqrt(node->GetVariance()[index]);
    }
    
  }
}


#endif //OUTPUT_HPP guard
