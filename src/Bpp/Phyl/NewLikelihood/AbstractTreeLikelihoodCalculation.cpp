//
// File: AbstractTreeLikelihoodCalculation.cpp
// Created by: Julien Dutheil
// Created on: Tue July 23 10:50 2013
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

#include "AbstractTreeLikelihoodCalculation.h"

using namespace bpp;
using namespace std;
using namespace newlik;

/******************************************************************************/

double AbstractTreeLikelihoodCalculation::getLogLikelihood() const
{
  vector<double> la(nbSites_);
  for (size_t i = 0; i < nbSites_; ++i)
  {
    la[i] = log(getLikelihoodForASite(i));
  }

  sort(la.begin(), la.end());
  double ll = 0;
  for (size_t i = nbSites_; i > 0; i--)
  {
    ll += la[i - 1];
  }
  return ll;
}

/******************************************************************************/

double AbstractTreeLikelihoodCalculation::getDLogLikelihood() const
{
  // Derivative of the sum is the sum of derivatives:
  vector<double> dla(nbSites_);
  for (size_t i = 0; i < nbSites_; ++i)
  {
    dla[i] = getDLogLikelihoodForASite(i);
  }
  
  sort(dla.begin(), dla.end());
  double dl = 0;
  for (size_t i = nbSites_; i > 0; --i)
  {
    dl += dla[i - 1];
  }

  return dl;
}

/******************************************************************************/

double AbstractTreeLikelihoodCalculation::getD2LogLikelihood() const
{
  // Derivative of the sum is the sum of derivatives:
  double dl = 0;
  for (size_t i = 0; i < nbSites_; ++i)
  {
    dl += getD2LogLikelihoodForASite(i);
  }

  return dl;
}


/******************************************************************************/

void AbstractTreeLikelihoodCalculation::displayLikelihoodArray(const VVVdouble& likelihoodArray)
{
  size_t nbClasses = likelihoodArray.size();
  size_t nbSites    = likelihoodArray[0].size();
  size_t nbStates  = likelihoodArray[0][0].size();
  for (size_t c = 0; c < nbClasses; ++c)
  {
    cout << "Model class " << c;
    for (size_t i = 0; i < nbSites; ++i)
    {
      cout << "Site " << i << ":" << endl;
      for (size_t s = 0; s < nbStates; ++s)
      {
        cout << "\t" << likelihoodArray[c][i][s];
      }
      cout << endl;
    }
    cout << endl;
  }
}

/******************************************************************************/

