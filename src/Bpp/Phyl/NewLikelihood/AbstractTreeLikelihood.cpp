//
// File: AbstractTreeLikelihood.cpp
// Created by: Julien Dutheil
// Created on: Fri Oct 17 17:57:21 2003
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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

#include "AbstractTreeLikelihood.h"

using namespace std;
using namespace bpp;
using namespace newlik;

/******************************************************************************/

AbstractTreeLikelihood::AbstractTreeLikelihood(
    SubstitutionProcess* process,
    bool verbose):
  AbstractParametrizable(""),
  data_(0),
  process_(process),
  computeFirstOrderDerivatives_(true),
  computeSecondOrderDerivatives_(true),
  initialized_(false),
  verbose_(verbose),
  nbSites_(0),
  nbDistinctSites_(0),
  nbStates_(process->getNumberOfStates()),
  nbClasses_(process->getNumberOfClasses())
{
  //initialize parameters:
  addParameters_(process_->getParameters());
}

/******************************************************************************/

Vdouble AbstractTreeLikelihood::getLikelihoodForEachSite() const
{
  Vdouble l(getNumberOfSites());
  for (unsigned int i = 0; i < l.size(); ++i)
    l[i] = getLikelihoodForASite(i);
  return l;
}

/******************************************************************************/

VVdouble AbstractTreeLikelihood::getLikelihoodForEachSiteForEachState() const
{
  VVdouble l(getNumberOfSites());
  for (unsigned int i = 0; i < l.size(); ++i)
    {
      Vdouble* l_i = &l[i];
      l_i->resize(getNumberOfStates());
		for(unsigned int x = 0; x < l_i->size(); ++x)
                  {
                    (* l_i)[x] = getLikelihoodForASiteForAState(i, x);
                  }
    }
  return l;
}

/******************************************************************************/

VVdouble AbstractTreeLikelihood::getLikelihoodForEachSiteForEachClass() const
{
  VVdouble l(getNumberOfSites());
  for (unsigned int i = 0; i < l.size(); ++i)
    {
      Vdouble* l_i = &l[i];
      l_i->resize(getNumberOfClasses());
      for (unsigned int c = 0; c < l_i->size(); ++c)
        {
          (* l_i)[c] = getLikelihoodForASiteForAClass(i, c);
        }
    }
  return l;
}

/******************************************************************************/

VVVdouble AbstractTreeLikelihood::getLikelihoodForEachSiteForEachClassForEachState() const
{
  VVVdouble l(getNumberOfSites());
  for (unsigned int i = 0; i < l.size(); ++i)
    {
      VVdouble* l_i = &l[i];
      l_i->resize(getNumberOfClasses());
      for (unsigned int c = 0; c < l_i->size(); ++c)
        {
          Vdouble* l_ic = &(*l_i)[c];
          l_ic->resize(getNumberOfStates());
          for (unsigned int x = 0; x < l_ic->size(); ++x)
            {
              (* l_ic)[x] = getLikelihoodForASiteForAClassForAState(i, c, x);
            }
        }
    }
  return l;
  
}

/******************************************************************************/

VVdouble AbstractTreeLikelihood::getPosteriorProbabilitiesOfEachClass() const
{
  size_t nbSites   = getNumberOfSites();
  size_t nbClasses = getNumberOfClasses();
  VVdouble pb = getLikelihoodForEachSiteForEachClass();
  Vdouble l = getLikelihoodForEachSite();
  for (size_t i = 0; i < nbSites; ++i)
  {
    for (size_t j = 0; j < nbClasses; ++j)
    {
      pb[i][j] = pb[i][j] * process_->getProbabilityForModel(j) / l[i];
    }
  }
  return pb;
}

/******************************************************************************/

vector<size_t> AbstractTreeLikelihood::getClassWithMaxPostProbOfEachSite() const
{
  size_t nbSites = getNumberOfSites();
  VVdouble l = getLikelihoodForEachSiteForEachClass();
  vector<size_t> classes(nbSites);
  for (size_t i = 0; i < nbSites; ++i)
  {
    classes[i] = VectorTools::whichMax<double>(l[i]);
  }
  return classes;
}

/******************************************************************************/

void AbstractTreeLikelihood::resetLikelihoodArray(VVVdouble& likelihoodArray)
{
  size_t nbSites   = likelihoodArray.size();
  size_t nbClasses = likelihoodArray[0].size();
  size_t nbStates  = likelihoodArray[0][0].size();
  for (size_t i = 0; i < nbSites; ++i)
  {
    for (size_t c = 0; c < nbClasses; ++c)
    {
      for (size_t s = 0; s < nbStates; ++s)
      {
        likelihoodArray[i][c][s] = 1.;
      }
    }
  }
}

/******************************************************************************/

void AbstractTreeLikelihood::displayLikelihoodArray(const VVVdouble& likelihoodArray)
{
  size_t nbSites   = likelihoodArray.size();
  size_t nbClasses = likelihoodArray[0].size();
  size_t nbStates  = likelihoodArray[0][0].size();
  for (size_t i = 0; i < nbSites; ++i)
  {
    cout << "Site " << i << ":" << endl;
    for (size_t c = 0; c < nbClasses; ++c)
    {
      cout << "Rate class " << c;
      for (size_t s = 0; s < nbStates; ++s)
      {
        cout << "\t" << likelihoodArray[i][c][s];
      }
      cout << endl;
    }
    cout << endl;
  }
}

/******************************************************************************/

