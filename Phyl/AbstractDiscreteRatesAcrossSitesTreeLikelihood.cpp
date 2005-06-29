//
// File: AbstractDiscreteRateAcrossSitesTreeLikelihood.cpp
// Created by: Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Wue Jun 15 09:42 2005
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

#include "AbstractDiscreteRatesAcrossSitesTreeLikelihood.h"

// From NumCalc:
#include <NumCalc/VectorTools.h>
using namespace VectorFunctions;

// From the STL:
#include <iostream>
using namespace std;

/******************************************************************************/

AbstractDiscreteRatesAcrossSitesTreeLikelihood::AbstractDiscreteRatesAcrossSitesTreeLikelihood(
	DiscreteDistribution * rDist,
	bool verbose
)	throw (Exception):
	AbstractTreeLikelihood(true)
{
	_rateDistribution = rDist;
}

/******************************************************************************/

ParameterList AbstractDiscreteRatesAcrossSitesTreeLikelihood::getRateDistributionParameters() const
{
	return _rateDistribution -> getParameters().getCommonParametersWith(_parameters);
}

/******************************************************************************/

VVdouble AbstractDiscreteRatesAcrossSitesTreeLikelihood::getLikelihoodForEachSiteForEachRateClass() const
{
	unsigned int nbSites   = getNumberOfSites();
	unsigned int nbClasses = getNumberOfClasses();
	VVdouble l(nbSites);
	for(unsigned int i = 0; i < nbSites; i++) {
		l[i].resize(nbClasses);
		for(unsigned int j = 0; j < nbClasses; j++)
			l[i][j] = getLikelihoodForASiteForARateClass(i, j);
	}
	return l;
}

/******************************************************************************/

double AbstractDiscreteRatesAcrossSitesTreeLikelihood::getLikelihoodForASiteForAState(unsigned int site, int state) const
{
	unsigned int nbClasses = getNumberOfClasses();
	double l = 0;
	for(unsigned int i = 0; i < nbClasses; i++) {
		l += getLikelihoodForASiteForARateClassForAState(site, i, state) * _rateDistribution -> getProbability(i);
	}
	return l;
}

/******************************************************************************/

double AbstractDiscreteRatesAcrossSitesTreeLikelihood::getLogLikelihoodForASiteForAState(unsigned int site, int state) const
{
	unsigned int nbClasses = getNumberOfClasses();
	double l = 0;
	for(unsigned int i = 0; i < nbClasses; i++) {
		l += getLikelihoodForASiteForARateClassForAState(site, i, state) * _rateDistribution -> getProbability(i);
	}
	//if(l <= 0.) cerr << "WARNING!!! Negative likelihood." << endl;
	return log(l);
}

/******************************************************************************/

VVdouble AbstractDiscreteRatesAcrossSitesTreeLikelihood::getLogLikelihoodForEachSiteForEachRateClass() const
{
	unsigned int nbSites   = getNumberOfSites();
	unsigned int nbClasses = getNumberOfClasses();
	VVdouble l(nbSites);
	for(unsigned int i = 0; i < nbSites; i++) {
		l[i] = Vdouble(nbClasses);
		for(unsigned int j = 0; j < nbClasses; j++)
			l[i][j] = getLogLikelihoodForASiteForARateClass(i, j);
	}
	return l;
}

/******************************************************************************/

VVVdouble AbstractDiscreteRatesAcrossSitesTreeLikelihood::getLikelihoodForEachSiteForEachRateClassForEachState() const
{
	unsigned int nbSites   = getNumberOfSites();
	unsigned int nbClasses = getNumberOfClasses();
	unsigned int nbStates  = getNumberOfStates();
	VVVdouble l(nbSites);
	for(unsigned int i = 0; i < nbSites; i++) {
		l[i].resize(nbClasses);
		for(unsigned int j = 0; j < nbClasses; j++) {
			l[i][j].resize(nbStates);
			for(int x = 0; x < nbStates; x++) {
				l[i][j][x] = getLikelihoodForASiteForARateClassForAState(i, j, x);
			}
		}
	}
	return l;
}

/******************************************************************************/

VVVdouble AbstractDiscreteRatesAcrossSitesTreeLikelihood::getLogLikelihoodForEachSiteForEachRateClassForEachState() const
{
	unsigned int nbSites   = getNumberOfSites();
	unsigned int nbClasses = getNumberOfClasses();
	unsigned int nbStates  = getNumberOfStates();
	VVVdouble l(nbSites);
	for(unsigned int i = 0; i < nbSites; i++) {
		l[i].resize(nbClasses);
		for(unsigned int j = 0; j < nbClasses; j++) {
			l[i][j].resize(nbStates);
			for(int x = 0; x < nbStates; x++) {
				l[i][j][x] = getLogLikelihoodForASiteForARateClassForAState(i, j, x);
			}
		}
	}
	return l;
}

/*******************************************************************************/

VVdouble AbstractDiscreteRatesAcrossSitesTreeLikelihood::getPosteriorProbabilitiesOfEachRate() const
{
	unsigned int nbSites   = getNumberOfSites();
	unsigned int nbClasses = getNumberOfClasses();
	VVdouble pb = getLikelihoodForEachSiteForEachRateClass();
	Vdouble  l  = getLikelihoodForEachSite();
	for(unsigned int i = 0; i < nbSites; i++) {
		for(unsigned int j = 0; j < nbClasses; j++) pb[i][j] = pb[i][j] * _rateDistribution -> getProbability(j) / l[i]; 
	}
	return pb;
}
	
/******************************************************************************/	

Vdouble AbstractDiscreteRatesAcrossSitesTreeLikelihood::getPosteriorRateOfEachSite() const
{
	unsigned int nbSites   = getNumberOfSites();
	unsigned int nbClasses = getNumberOfClasses();
	VVdouble lr = getLikelihoodForEachSiteForEachRateClass();
	Vdouble  l  = getLikelihoodForEachSite();
	Vdouble rates(nbSites, 0.);
	for(unsigned int i = 0; i < nbSites; i++) {
		for(unsigned int j = 0; j < nbClasses; j++) {
			rates[i] += (lr[i][j] / l[i]) * _rateDistribution -> getProbability(j) *  _rateDistribution -> getCategory(j);
		}
	}
	return rates;
}

/******************************************************************************/

Vint AbstractDiscreteRatesAcrossSitesTreeLikelihood::getRateClassWithMaxPostProbOfEachSite() const
{
	unsigned int nbSites   = getNumberOfSites();
	VVdouble l = getLikelihoodForEachSiteForEachRateClass();
	Vint classes(nbSites);
	for(unsigned int i = 0; i < nbSites; i++) classes[i] = posmax<double>(l[i]);
	return classes;
}

/******************************************************************************/

Vdouble AbstractDiscreteRatesAcrossSitesTreeLikelihood::getRateWithMaxPostProbOfEachSite() const
{
	unsigned int nbSites   = getNumberOfSites();
	VVdouble l = getLikelihoodForEachSiteForEachRateClass();
	Vdouble rates(nbSites);
	for(unsigned int i = 0; i < nbSites; i++) {
		rates[i] = _rateDistribution -> getCategory(posmax<double>(l[i]));
	}
	return rates;
}

/******************************************************************************/

