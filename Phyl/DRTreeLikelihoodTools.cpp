//
// File: DRTreeLikelihoodTools.cpp
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Mon Janv 17 09:56 2005
//

/*
Copyright ou � ou Copr. CNRS, (16 Novembre 2004) 

Julien.Dutheil@univ-montp2.fr

Ce logiciel est un programme informatique servant � fournir des classes
pour l'analyse de donn�es phylog�n�tiques.

Ce logiciel est r�gi par la licence CeCILL soumise au droit fran�ais et
respectant les principes de diffusion des logiciels libres. Vous pouvez
utiliser, modifier et/ou redistribuer ce programme sous les conditions
de la licence CeCILL telle que diffus�e par le CEA, le CNRS et l'INRIA 
sur le site "http://www.cecill.info".

En contrepartie de l'accessibilit� au code source et des droits de copie,
de modification et de redistribution accord�s par cette licence, il n'est
offert aux utilisateurs qu'une garantie limit�e.  Pour les m�mes raisons,
seule une responsabilit� restreinte p�se sur l'auteur du programme,  le
titulaire des droits patrimoniaux et les conc�dants successifs.

A cet �gard  l'attention de l'utilisateur est attir�e sur les risques
associ�s au chargement,  � l'utilisation,  � la modification et/ou au
d�veloppement et � la reproduction du logiciel par l'utilisateur �tant 
donn� sa sp�cificit� de logiciel libre, qui peut le rendre complexe � 
manipuler et qui le r�serve donc � des d�veloppeurs et des professionnels
avertis poss�dant  des  connaissances  informatiques approfondies.  Les
utilisateurs sont donc invit�s � charger  et  tester  l'ad�quation  du
logiciel � leurs besoins dans des conditions permettant d'assurer la
s�curit� de leurs syst�mes et ou de leurs donn�es et, plus g�n�ralement, 
� l'utiliser et l'exploiter dans les m�mes conditions de s�curit�. 

Le fait que vous puissiez acc�der � cet en-t�te signifie que vous avez 
pris connaissance de la licence CeCILL, et que vous en avez accept� les
termes.
*/

/*
Copyright or � or Copr. CNRS, (November 16, 2004)

Julien.Dutheil@univ-montp2.fr

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

#include "DRTreeLikelihoodTools.h"
#include <NumCalc/VectorTools.h>
using namespace VectorFunctions;

//-----------------------------------------------------------------------------------------

VVVdouble DRTreeLikelihoodTools::getPosteriorProbabilitiesForEachStateForEachRate(
							DRHomogeneousTreeLikelihood & drl,
							const Node * node)
{
	unsigned int nSites   = drl.getLikelihoodData()->getNumberOfDistinctSites();
	unsigned int nClasses = drl.getNumberOfClasses();
	unsigned int nStates  = drl.getNumberOfStates();
	VVVdouble postProb(nSites);
	
	const DiscreteDistribution * rDist = drl.getRateDistribution();
	Vdouble rcProbs = rDist -> getProbabilities();
	if(node -> isLeaf()) {
		VVdouble larray = drl.getLikelihoodData()->getLeafLikelihoods(node);
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

