//
// File: DRTreeLikelihoodTools.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Mon Janv 17 09:56 2005
//

#ifndef _DRTREELIKELIHOODTOOLS_H_
#define _DRTREELIKELIHOODTOOLS_H_

#include "DRHomogeneousTreeLikelihood.h"
#include <Seq/AlignedSequenceContainer.h>

class DRTreeLikelihoodTools {

	public: static vector<int> getMarginalAncestralStates(
							DRHomogeneousTreeLikelihood & drl,
							const Node * node);

	public: static map<const Node *, vector<int> > getAllMarginalAncestralStates(
							DRHomogeneousTreeLikelihood & drl);

	private: static void recursiveMarginalAncestralStates(
							 DRHomogeneousTreeLikelihood & drl,
							 const Node * node,
							 map<const Node *, vector<int> > & ancestors,
							 AlignedSequenceContainer & data);

	public: static VVVdouble getPosteriorProbabilitiesForEachStateForEachRate(
							DRHomogeneousTreeLikelihood & drl,
							const Node * node);
};

#endif //_DRTREELIKELIHOODTOOLS_H_

