//
// File: DRTreeLikelihoodTools.cpp
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Mon Janv 17 09:56 2005
//

#include "DRTreeLikelihoodTools.h"
#include <NumCalc/VectorTools.h>
using namespace VectorFunctions;

//-----------------------------------------------------------------------------------------

vector<int> DRTreeLikelihoodTools::getMarginalAncestralStates(
	DRHomogeneousTreeLikelihood & drl,
	const Node * node)
{
	unsigned int nSites   = drl.getNumberOfDistinctSites();
	unsigned int nClasses = drl.getNumberOfClasses();
	unsigned int nStates  = drl.getNumberOfStates();
	vector<int> ancestors(nSites);
	if(node -> isLeaf()) {
		VVdouble larray = drl.getLeafLikelihoods(node);
		for(unsigned int i = 0; i < nSites; i++) {
			ancestors[i] = (int)posmax(larray[i]);
		}
	} else {
		VVVdouble larray = drl.computeLikelihoodAtNode(node);
		Vdouble freqs = drl.getSubstitutionModel() -> getFrequencies();
		Vdouble rcProbs = drl.getRateDistribution() -> getProbabilities(); 
		for(unsigned int i = 0; i < nSites; i++) {
			Vdouble likelihoods(nStates, 0);
			VVdouble * larray_i = & larray[i];
			for(unsigned int c = 0; c < nClasses; c++) {
				Vdouble * larray_i_c = & (* larray_i)[c];
				double rcp = rcProbs[c];
				for(unsigned int x = 0; x < nStates; x++) {
					likelihoods[x] += (* larray_i_c)[x] * freqs[x] * rcp;
				}
			}
			ancestors[i] = (int)posmax(likelihoods);
		}
	}
	return ancestors;
}

//-----------------------------------------------------------------------------------------

map<const Node *, vector<int> > DRTreeLikelihoodTools::getAllMarginalAncestralStates(
	DRHomogeneousTreeLikelihood & drl)
{
	map<const Node *, vector<int> > ancestors;
	// Clone the data into a AlignedSequenceContainer for more efficiency:
	AlignedSequenceContainer * data = new AlignedSequenceContainer(* drl.getShrunkData());
	recursiveMarginalAncestralStates(drl, drl.getTree() -> getRootNode(), ancestors, *data);
	delete data;
	return ancestors;
}

//-----------------------------------------------------------------------------------------

void DRTreeLikelihoodTools::recursiveMarginalAncestralStates(
	DRHomogeneousTreeLikelihood & drl,
	const Node * node,
	map<const Node *, vector<int> > & ancestors,
	AlignedSequenceContainer & data)
{
	if(node -> isLeaf()) {
		ancestors[node] = data.getSequence(node -> getName()) -> getContent();
	} else {
		ancestors[node] = getMarginalAncestralStates(drl, node);
		for(unsigned int i = 0; i < node -> getNumberOfSons(); i++) {
			recursiveMarginalAncestralStates(drl, node -> getSon(i), ancestors, data);
		}
	}
}

//-----------------------------------------------------------------------------------------

VVVdouble DRTreeLikelihoodTools::getPosteriorProbabilitiesForEachStateForEachRate(
							DRHomogeneousTreeLikelihood & drl,
							const Node * node)
{
	unsigned int nSites   = drl.getNumberOfDistinctSites();
	unsigned int nClasses = drl.getNumberOfClasses();
	unsigned int nStates  = drl.getNumberOfStates();
	VVVdouble postProb(nSites);
	
	const DiscreteDistribution * rDist = drl.getRateDistribution();
	Vdouble rcProbs = rDist -> getProbabilities();
	if(node -> isLeaf()) {
		VVdouble larray = drl.getLeafLikelihoods(node);
		for(unsigned int i = 0; i < nSites; i++) {
			VVdouble * postProb_i = & postProb[i];
			postProb_i -> resize(nClasses);
			Vdouble * larray_i = & larray[i];
			for(unsigned int c = 0; c < nClasses; c++) {
				Vdouble * postProb_i_c = & (* postProb_i)[c];
				postProb_i_c -> resize(nStates);
				double * rcProb = & rcProbs[c];
				for(unsigned int x = 0; x < nStates; x++) {
					(* postProb_i_c)[x] = (* larray_i)[x] * (* rcProb);
				}
			}
		}
	} else {
		VVVdouble larray = drl.computeLikelihoodAtNode(node);
		
		Vdouble likelihoods(nSites, 0);
		Vdouble freqs = drl.getSubstitutionModel() -> getFrequencies();
		Vdouble rcRates = rDist -> getCategories();
		for(unsigned int i = 0; i < nSites; i++) {
			VVdouble * larray_i = & larray[i];
			for(unsigned int c = 0; c < nClasses; c++) {
				Vdouble * larray_i_c = & (* larray_i)[c];
				double rcp = rcProbs[c];
				for(unsigned int s = 0; s < nStates; s++) {
					likelihoods[i] += rcp * freqs[s] * (* larray_i_c)[s];
				}
			}
		}
		
		for(unsigned int i = 0; i < nSites; i++) {
			VVdouble * postProb_i = & postProb[i];
			postProb_i -> resize(nClasses);
			VVdouble * larray_i = & larray[i];
			double likelihood = likelihoods[i];
			for(unsigned int c = 0; c < nClasses; c++) {
				Vdouble * postProb_i_c = & (* postProb_i)[c];
				postProb_i_c -> resize(nStates);
				Vdouble * larray_i_c = & (* larray_i)[c];
				double rcProb = rcProbs[c];
				for(unsigned int x = 0; x < nStates; x++) {
					(* postProb_i_c)[x] = (* larray_i_c)[x] * freqs[x] * rcProb / likelihood;
				}
			}
		}
	}
	return postProb;
}

//-----------------------------------------------------------------------------------------

