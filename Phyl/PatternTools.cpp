//
// File: PatternTools.h
// Created by: Julien Dutheil <julien.dutheil@ens-lyon.fr>
// Created on: Thu Mar 20 13:36:54 2003
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

#include "PatternTools.h"

#include "TreeTools.h"

// From SeqLib:
#include <Seq/VectorSiteContainer.h>
#include <Seq/SiteTools.h>

// From the STL:
#include <iostream>
#include <map>
#include <algorithm>

using namespace std;

/******************************************************************************/

SiteContainer * PatternTools::getSequenceSubset(const SiteContainer & sequenceSet, const Node & node) throw (Exception)
{
	VectorSiteContainer * sequenceSubset = new VectorSiteContainer(sequenceSet.getAlphabet());
	vector<const Node *> leaves = TreeTools::getLeaves(node);
	for(vector<const Node *>::iterator i = leaves.begin(); i < leaves.end(); i++) {
		const Sequence * newSeq = sequenceSet.getSequence((* i) -> getName());
		if(newSeq == NULL) throw SequenceNotFoundException("PatternToolsERROR: leaf name not found in sequence file: ", (* i) -> getName());
		sequenceSubset -> addSequence(* newSeq);
	}
	return sequenceSubset;
}

/******************************************************************************/

SiteContainer * PatternTools::getSequenceSubset(const SiteContainer & sequenceSet, const vector<string> & names) throw (Exception)
{
	VectorSiteContainer * sequenceSubset = new VectorSiteContainer(sequenceSet.getAlphabet());
	for(unsigned int i = 0; i < names.size(); i++) {
		const Sequence * newSeq = sequenceSet.getSequence(names[i]);
		if(newSeq == NULL) throw SequenceNotFoundException("PatternTools ERROR: name not found in sequence file: ", names[i]);
		sequenceSubset -> addSequence(* newSeq);
	}
	return sequenceSubset;
}

/******************************************************************************/

SiteContainer * PatternTools::shrinkSiteSet(const SiteContainer & siteSet) throw (Exception)
{
	if(siteSet.getNumberOfSites() == 0) throw Exception("PatternTools::shrinkSiteSet siteSet is void.");
	vector<const Site *> sites;
	for(unsigned int i = 0; i < siteSet.getNumberOfSites(); i++) {
		const Site * currentSite = siteSet.getSite(i);
		bool siteExists = false;
		for(unsigned int j = 0; !siteExists && j < sites.size(); j++) {
			if(SiteTools::areSitesIdentical(* currentSite, * sites[j])) siteExists = true;
		}
		if(!siteExists)	sites.push_back(currentSite);
	}
	SiteContainer * result = new VectorSiteContainer(sites, siteSet.getAlphabet());
	result -> setSequencesNames(siteSet.getSequencesNames(), false);
	return result;
}

/******************************************************************************/

 // O(n2)
const Pattern PatternTools::countSites(const SiteContainer & siteSet)
{
	Pattern pattern;
	unsigned int currentPos = 0;
	for(unsigned int i = 0; i < siteSet.getNumberOfSites(); i++) {
		pattern.indices.push_back(currentPos);
		const Site * currentSite = siteSet.getSite(i);
		bool siteExists = false;
		for(unsigned int j = 0; !siteExists && j < pattern.sites.size(); j++) {
			if(SiteTools::areSitesIdentical(* currentSite, * pattern.sites[j])) {
				siteExists = true;
				pattern.weights[j]++;
			}
		}
		if(!siteExists)	{
			pattern.sites.push_back(currentSite);
			pattern.weights.push_back(1);
			currentPos++;
		}
	}
	pattern.names = siteSet.getSequencesNames();
	return pattern;
}


/******************************************************************************/
 // o(n.log(n))
/*
const Pattern PatternTools::countSites(const SiteContainer & siteSet)
{
	unsigned int nbSites = siteSet.getNumberOfSites();
	vector<SortableSite *> ss(nbSites);
	for(unsigned int i = 0; i < nbSites; i++) {
		const Site * currentSite = siteSet.getSite(i);
		SortableSite * ssi = new SortableSite();
		ss[i] = ssi;
		ssi -> siteS = currentSite -> toString();
		ssi -> siteP = currentSite;
		ssi -> originalPosition = i;
	}

	// Quick sort according to site contents:
	sort(ss.begin(), ss.end(), SSComparator());
	
	// Now build patterns:
	Pattern pattern;
	SortableSite * ss0 = ss[0];
	const Site * previousSite = ss0 -> siteP;
	pattern.indices.resize(nbSites);
	pattern.indices[ss0 -> originalPosition] = 1;
	pattern.sites.push_back(previousSite);
	pattern.weights.push_back(1);
	
	unsigned int currentPos = 0;
	for(unsigned int i = 1; i < nbSites; i++) {
		SortableSite * ssi = ss[i];
		const Site * currentSite = ssi -> siteP;
		pattern.indices[ssi -> originalPosition] = currentPos;
		bool siteExists = SiteTools::areSitesIdentical(* currentSite, * previousSite);
		if(siteExists) {
			pattern.weights[currentPos]++;
		} else {
			pattern.sites.push_back(currentSite);
			pattern.weights.push_back(1);
			currentPos++;
		}
		delete ss[i - 1];
		previousSite = currentSite;
		pattern.indices.push_back(currentPos);
	}
	delete ss[nbSites - 1];
	pattern.names = siteSet.getSequencesNames();
	return pattern;
}
*/

/******************************************************************************/

/*
const Pattern PatternTools::countSites(const SiteContainer & siteSet)
{
	unsigned int nbSites = siteSet.getNumberOfSites();
	Pattern pattern;
	map<string, unsigned int> ss;
	unsigned int currentPos = 0;
	for(unsigned int i = 0; i < nbSites; i++) {
		pattern.indices.push_back(currentPos);
		const Site * currentSite = siteSet.getSite(i);
		unsigned int * c = & ss[currentSite -> toString()];
		if(*c > 0) {
			pattern.weights[*c]++;
		} else {
			pattern.sites.push_back(currentSite);
			pattern.weights.push_back(1);
			*c=currentPos;
			currentPos++;
		}
	}
	
	pattern.names = siteSet.getSequencesNames();
	return pattern;
}
*/

/******************************************************************************/

const vector<unsigned int> PatternTools::getWeights(const Pattern & pattern)
{
	return pattern.weights;
}

/******************************************************************************/

const vector<unsigned int> PatternTools::getIndices(const Pattern & pattern)
{
	return pattern.indices;
}

/******************************************************************************/


SiteContainer * PatternTools::getSites(const Pattern & pattern, const Alphabet * alpha)
{
	SiteContainer * sites = new VectorSiteContainer(pattern.sites, alpha);
	sites -> setSequencesNames(pattern.names, false);
	return sites;
}

/******************************************************************************/

Vint PatternTools::getIndexes(const SiteContainer & sequences1, const SiteContainer & sequences2) {
	int nbSites = sequences1.getNumberOfSites(); 
	Vint indexes(nbSites);
	for(int i = 0; i < nbSites; i++) {
		//For each site in sequences1,
		indexes[i] = -1;
		const Site * site = sequences1.getSite(i);
		for(unsigned int j = 0; j < sequences2.getNumberOfSites(); j++) {
			if(SiteTools::areSitesIdentical(* site, * sequences2.getSite(j))) {
				indexes[i] = j;
				break;
			}
		}
	}
	return indexes;
}

/******************************************************************************/
