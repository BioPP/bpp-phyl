// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "AbstractTreeLikelihood.h"

using namespace bpp;

/******************************************************************************/

Vdouble AbstractTreeLikelihood::getLikelihoodPerSite() const
{
  Vdouble l(getNumberOfSites());
  for (size_t i = 0; i < l.size(); i++)
  {
    l[i] = getLikelihoodForASite(i);
  }
  return l;
}

/******************************************************************************/

Vdouble AbstractTreeLikelihood::getLogLikelihoodPerSite() const
{
  Vdouble l(getNumberOfSites());
  for (size_t i = 0; i < l.size(); i++)
  {
    l[i] = getLogLikelihoodForASite(i);
  }
  return l;
}

/******************************************************************************/

VVdouble AbstractTreeLikelihood::getLikelihoodPerSitePerState() const
{
  VVdouble l(getNumberOfSites());
  for (size_t i = 0; i < l.size(); i++)
  {
    Vdouble* l_i = &l[i];
    l_i->resize(getNumberOfStates());
    for (size_t x = 0; x < l_i->size(); x++)
    {
      (*l_i)[x] = getLikelihoodForASiteForAState(i, static_cast<int>(x));
    }
  }
  return l;
}

/******************************************************************************/

VVdouble AbstractTreeLikelihood::getLogLikelihoodPerSitePerState() const
{
  VVdouble l(getNumberOfSites());
  for (size_t i = 0; i < l.size(); i++)
  {
    Vdouble* l_i = &l[i];
    l_i->resize(getNumberOfStates());
    for (size_t x = 0; x < l_i->size(); x++)
    {
      (*l_i)[x] = getLogLikelihoodForASiteForAState(i, static_cast<int>(x));
    }
  }
  return l;
}

/******************************************************************************/
