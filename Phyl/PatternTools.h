//
// File: PatternTools.h
// Created by: Julien Dutheil <julien.dutheil@ens-lyon.fr>
// Created on: Thu Mar 20 13:36:53 2003
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
 
#ifndef _PATTERNTOOLS_H_
#define _PATTERNTOOLS_H_

#include "Tree.h"

// From NumCalc:
#include <NumCalc/VectorTools.h>

//From SeqLib:
#include <Seq/SiteContainer.h>
#include <Seq/Site.h>

// From the STL:
#include <map>

struct Pattern
{
	vector<const Site *> sites;
	Vdouble weights;
};

/** 
 * const Site * points toward a unique site,
 * int is the number of sites identical to this sites.
 */


class PatternTools
{
	public:
		static SiteContainer * getSequenceSubset(const SiteContainer & sequenceSet, const Node & node) throw (Exception);
		static SiteContainer * shrinkSiteSet(const SiteContainer & sequenceSet) throw (Exception);
		static const Pattern countSites(const SiteContainer & sequences);
		static const Vdouble getWeights(const Pattern & pattern);
		static SiteContainer * getSites(const Pattern & pattern, const Alphabet * alpha);
		/**
	     * This method looks for the occurence of each site in sequences1 in sequences2 and send the
	     * position of the first occurence, or -1 if not found.
	     */
		static Vint getIndexes(const SiteContainer & sequences1, const SiteContainer & sequences2);
};


#endif	//_PATTERNTOOLS_H_
