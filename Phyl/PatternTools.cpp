//
// File: PatternTools.h
// Created by: Julien Dutheil <julien.dutheil@ens-lyon.fr>
// Created on: Thu Mar 20 13:36:54 2003
//

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
#include "PatternTools.h"

#include "TreeTools.h"

// From SeqLib:
#include <Seq/VectorSiteContainer.h>
#include <Seq/SiteTools.h>

// From the STL:
#include <iostream>

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

const Pattern PatternTools::countSites(const SiteContainer & siteSet)
{
	Pattern pattern;
	for(unsigned int i = 0; i < siteSet.getNumberOfSites(); i++) {
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
		}
	}
	return pattern;
}

/******************************************************************************/

const Vdouble PatternTools::getWeights(const Pattern & pattern)
{
	return pattern.weights;
}

/******************************************************************************/

SiteContainer * PatternTools::getSites(const Pattern & pattern, const Alphabet * alpha)
{
	return new VectorSiteContainer(pattern.sites, alpha);
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
