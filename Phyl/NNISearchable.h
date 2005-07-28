//
// File: NNISearchable.h
// Created by: Julien Dutheil
//             Vincent Ranwez
// Created on: Thu Jul 28 16:32 2005
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004)

This software is a computer program whose purpose is to provide classes
for phylogenetic data analysis.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef _NNISEARCHABLE_H_
#define _NNISEARCHABLE_H_

/**
 * @brief Interface for Nearest Neighbor Interchanged algorithms.
 *
 * This interface defines the methods to work with NNI algorithms.
 * 
 * NNISearchable objects are supposed to work with TreeTemplate objects.
 * NNI are defined as follow:
 * <pre>
 *D      C
 * \    /
 *  \  /
 *   \/
 *   X\    /B
 *     \  /
 *      \/
 *      F\
 *        \
 *         \A
 * 
 * </pre>
 * Where:
 * -F is the focuse (parent) node,
 * -A and B are the son of F
 * -X is the parent of F and so on.
 * We work on rooted trees, whether this root as a biological sens or not.
 * Two NNI's are possible for branch (XF):
 * - swaping B<->C, which is the same as D<->A
 * - swaping A<->C, which is the same as D<->B
 * Because of the rooted representation, we'll consider B<->C and A<->C,
 * which are simpler to perform.
 */
class NNISearchable
{
	public:
		NNISearchable() {}
		~NNISearchable() {}

	public:

		/**
		 * @brief Send the score of a NNI movement, without performing it.
		 *
		 * @param parent The focus node.
		 * @param son    The son node.
		 * @return The score of the NNI.
		 * @throw NodeException If 'son' is not a son of 'parent'.
		 */
		virtual double testNNI(const Node * parent, const Node * son) const throw (NodeException) = 0;

		/**
		 * @brief Perform a NNI movement.
		 *
		 * @param parent The focus node.
		 * @param son    The son node.
		 * @throw NodeException If 'son' is not a son of 'parent'.
		 */
		virtual void doNNI(const Node * parent, const Node * son) throw (NodeException) = 0;
};

#endif //_NNISEARCHABLE_H_

