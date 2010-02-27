//
// File: DRNonHomogeneousMixedTreeLikelihood.cpp
// Created by: Laurent Gueguen
//

/*
   Copyright or © or Copr. CNRS, (November 16, 2004)

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

#include "DRNonHomogeneousMixedTreeLikelihood.h"


// From the STL:
#include <iostream>

#include <math.h>
#include "PatternTools.h"

#include <NumCalc/VectorTools.h>

#include <Utils/ApplicationTools.h>

using namespace bpp;
using namespace std;

DRNonHomogeneousMixedTreeLikelihood::DRNonHomogeneousMixedTreeLikelihood(
  const Tree& tree,
  SubstitutionModelSet* modelSet,
  DiscreteDistribution* rDist,
  bool verbose) throw (Exception) :
  DRNonHomogeneousTreeLikelihood(tree, modelSet, rDist, verbose),
  treeLikelihoodsContainer_(), internParam_()
{
  if (!modelSet->isFullySetUpFor(tree))
    throw Exception("RNonHomogeneousMixedTreeLikelihood(constructor). Model set is not fully specified.");

  map<int, int> mapmodels;
  unsigned int ttmodels = 1;
  unsigned int nbmodels = modelSet->getNumberOfModels();
  MixedSubstitutionModel* pmsm;
  SubstitutionModel* psm;
  for (unsigned int i = 0; i < nbmodels; i++)
  {
    if ((pmsm = dynamic_cast< MixedSubstitutionModel*>(modelSet->getModel(i))) != 0)
      mapmodels[i] = pmsm->getNumberOfModels();
    else
      mapmodels[i] = 1;
    ttmodels *= mapmodels[i];
  }

  SubstitutionModelSet* psms;
  int s;
  vector<string> vn;

  for (unsigned int i = 0; i < ttmodels; i++)
  {
    s = i;

    psms = new SubstitutionModelSet(modelSet->getAlphabet(),modelSet->getRootFrequenciesSet()->clone());

    for (unsigned int j = 0; j < nbmodels; j++)
    {
      if (mapmodels[j] == 1)
        psms->addModel(modelSet->getModel(j),modelSet->getNodesWithModel(j),
                       modelSet->getModel(j)->getParameters().getParameterNames());
      else
      {
        pmsm = dynamic_cast< MixedSubstitutionModel*>(modelSet->getModel(j));
        psm = pmsm->getNModel(s % mapmodels[j]);
        psms->addModel(psm,modelSet->getNodesWithModel(j),
                       psm->getParameters().getParameterNames());

        vn = psms->getModelParameters(j).getParameterNames();
        for (unsigned int i2 = 0; i2 < vn.size(); i2++)
        {
        Parameter p = psm->getParameter(psm->getParameterNameWithoutNamespace(psms->getParameterModelName(vn[i2])));
          if (!internParam_.hasParameter(vn[i2]))
            internParam_.addParameter(
                Parameter(vn[i2], p.getValue(), p.getConstraint()->clone(), true));
        }

        s /= mapmodels[j];
      }
    }
    treeLikelihoodsContainer_.push_back(
        new DRNonHomogeneousTreeLikelihood(tree, psms, rDist, false));
  }
}

DRNonHomogeneousMixedTreeLikelihood::DRNonHomogeneousMixedTreeLikelihood(
  const Tree& tree,
  const SiteContainer& data,
  SubstitutionModelSet* modelSet,
  DiscreteDistribution* rDist,
  bool verbose) throw (Exception) :
  DRNonHomogeneousTreeLikelihood(tree, modelSet,rDist, verbose),
  treeLikelihoodsContainer_(), internParam_()
{
  if (!modelSet->isFullySetUpFor(tree))
    throw Exception("RNonHomogeneousMixedTreeLikelihood(constructor). Model set is not fully specified.");

  map<int, int> mapmodels;
  unsigned int ttmodels = 1;
  unsigned int nbmodels = modelSet->getNumberOfModels();
  MixedSubstitutionModel* pmsm;
  SubstitutionModel* psm;
  for (unsigned int i = 0; i < nbmodels; i++)
  {
    if ((pmsm = dynamic_cast< MixedSubstitutionModel*>(modelSet->getModel(i))) != 0)
      mapmodels[i] = pmsm->getNumberOfModels();
    else
      mapmodels[i] = 1;
    ttmodels *= mapmodels[i];
  }

  SubstitutionModelSet* psms;
  int s;
  vector<string> vn;

  for (unsigned int i = 0; i < ttmodels; i++)
  {
    s = i;

    psms = new SubstitutionModelSet(modelSet->getAlphabet(),modelSet->getRootFrequenciesSet()->clone());

    for (unsigned int j = 0; j < nbmodels; j++)
    {
      if (mapmodels[j] == 1)
        psms->addModel(modelSet->getModel(j),modelSet->getNodesWithModel(j),
                       modelSet->getModel(j)->getParameters().getParameterNames());
      else
      {
        pmsm = dynamic_cast< MixedSubstitutionModel*>(modelSet->getModel(j));
        psm = pmsm->getNModel(s % mapmodels[j]);
        psms->addModel(psm,modelSet->getNodesWithModel(j),
                       psm->getParameters().getParameterNames());

        vn = psms->getModelParameters(j).getParameterNames();
        for (unsigned int i2 = 0; i2 < vn.size(); i2++)
        {
          Parameter p = psm->getParameter(psm->getParameterNameWithoutNamespace(psms->getParameterModelName(vn[i2])));
          if (!internParam_.hasParameter(vn[i2]))
            internParam_.addParameter(
                Parameter(vn[i2], p.getValue(), p.getConstraint()->clone(), true));
        }

        s /= mapmodels[j];
      }
    }
    treeLikelihoodsContainer_.push_back(
        new DRNonHomogeneousTreeLikelihood(tree, psms, rDist, false));
  }

  setData(data);
}


DRNonHomogeneousMixedTreeLikelihood & DRNonHomogeneousMixedTreeLikelihood::operator=(const DRNonHomogeneousMixedTreeLikelihood& lik)
{
  DRNonHomogeneousTreeLikelihood::operator=(lik);

  internParam_ = lik.internParam_;
  
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
   treeLikelihoodsContainer_.push_back(lik.treeLikelihoodsContainer_[i]->clone());
  }

  return *this;
}

DRNonHomogeneousMixedTreeLikelihood::DRNonHomogeneousMixedTreeLikelihood(const DRNonHomogeneousMixedTreeLikelihood& lik) :
  DRNonHomogeneousTreeLikelihood(lik),
  treeLikelihoodsContainer_(lik.treeLikelihoodsContainer_.size()), internParam_(lik.internParam_)

{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i] = lik.treeLikelihoodsContainer_[i]->clone();
  }
}

DRNonHomogeneousMixedTreeLikelihood::~DRNonHomogeneousMixedTreeLikelihood()
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    delete treeLikelihoodsContainer_[i];
  }
}


void DRNonHomogeneousMixedTreeLikelihood::initialize() throw (Exception)
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->initialize();
  }
  DRNonHomogeneousTreeLikelihood::initialize();
}

void DRNonHomogeneousMixedTreeLikelihood::setData(const SiteContainer& sites) throw (Exception)
{
   DRNonHomogeneousTreeLikelihood::setData(sites);

  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->setData(sites);
  }
}


double DRNonHomogeneousMixedTreeLikelihood::getLikelihood() const
{
  return exp(getLogLikelihood());
}

double DRNonHomogeneousMixedTreeLikelihood::getLogLikelihood() const
{
   vector<double> reslog;

  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
   reslog.push_back(treeLikelihoodsContainer_[i]->getLogLikelihood());
  }

  return VectorTools::logmeanexp(reslog);
}

double DRNonHomogeneousMixedTreeLikelihood::getLikelihoodForASite(unsigned int site) const
{
  return exp(getLogLikelihoodForASite(site));
}

double DRNonHomogeneousMixedTreeLikelihood::getLogLikelihoodForASite(unsigned int site) const
{
   vector<double> reslog;
   unsigned int i;

  for (i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
   reslog.push_back(treeLikelihoodsContainer_[i]->getLogLikelihoodForASite(site));
  }

  return VectorTools::logmeanexp(reslog);
}

double DRNonHomogeneousMixedTreeLikelihood::getLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const
{
  return exp(getLogLikelihoodForASiteForARateClass(site, rateClass));
}

double DRNonHomogeneousMixedTreeLikelihood::getLogLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const
{
   vector<double> reslog;
   unsigned int i;

  for (i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
   reslog.push_back(treeLikelihoodsContainer_[i]->getLogLikelihoodForASiteForARateClass(site, rateClass));
  }

  return VectorTools::logmeanexp(reslog);
}

double DRNonHomogeneousMixedTreeLikelihood::getLikelihoodForASiteForARateClassForAState(unsigned int site, unsigned int rateClass, int state) const
{
  return exp(getLogLikelihoodForASiteForARateClassForAState(site, rateClass, state));
}

double DRNonHomogeneousMixedTreeLikelihood::getLogLikelihoodForASiteForARateClassForAState(unsigned int site, unsigned int rateClass, int state) const
{
   vector<double> reslog;
   unsigned int i;

  for (i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
   reslog.push_back(treeLikelihoodsContainer_[i]->getLogLikelihoodForASiteForARateClassForAState(site, rateClass, state));
  }

  return VectorTools::logmeanexp(reslog);
}

void DRNonHomogeneousMixedTreeLikelihood::fireParameterChanged(const ParameterList& params)
{
  applyParameters();

  unsigned int s;
  const SubstitutionModel* psm;
  MixedSubstitutionModel* pm;
  SubstitutionModelSet* modelSet = getSubstitutionModelSet();
  SubstitutionModelSet* psms;

  ParameterList par;
  vector<string> vp;

  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    s = i;
    for (unsigned int j = 0; j < modelSet->getNumberOfModels(); j++)
    {
      pm = dynamic_cast<MixedSubstitutionModel*>(modelSet->getModel(j));

      if (pm != NULL)
      {
        psm = pm->getNModel(s % pm->getNumberOfModels());

        psms = treeLikelihoodsContainer_[i]->getSubstitutionModelSet();
        par = psms->getNodeParameters();
        vp = par.getParameterNames();

        for (unsigned int i2 = 0; i2 < vp.size(); i2++)
        {
          if (internParam_.hasParameter(vp[i2]))
            internParam_.setParameterValue(vp[i2],
                                           psm->getParameterValue(psm->getParameterNameWithoutNamespace(psms->getParameterModelName(vp[i2]))));
        }


        treeLikelihoodsContainer_[i]->matchParametersValues(internParam_);

        s /= pm->getNumberOfModels();
      }
    }
    treeLikelihoodsContainer_[i]->matchParametersValues(getParameters());
  }
  minusLogLik_ = -getLogLikelihood();
}

void DRNonHomogeneousMixedTreeLikelihood::computeTreeDLikelihoods()
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeTreeDLikelihoods();
  }
}

double DRNonHomogeneousMixedTreeLikelihood::getFirstOrderDerivative(const string& variable) const throw (Exception)
{
   vector<double> rescontainer;
   unsigned int i;

  for (i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
   rescontainer.push_back(treeLikelihoodsContainer_[i]->getFirstOrderDerivative(variable));
  }

  return VectorTools::mean<double, double>(rescontainer);
}

void DRNonHomogeneousMixedTreeLikelihood::computeTreeD2LikelihoodAtNode(const Node* node)
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeTreeD2LikelihoodAtNode(node);
  }
}

void DRNonHomogeneousMixedTreeLikelihood::computeTreeD2Likelihoods()
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeTreeD2Likelihoods();
  }
}

double DRNonHomogeneousMixedTreeLikelihood::getSecondOrderDerivative(const string& variable) const throw (Exception)
{
   vector<double> rescontainer;
   unsigned int i;

  for (i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
   rescontainer.push_back(treeLikelihoodsContainer_[i]->getSecondOrderDerivative(variable));
  }

  return VectorTools::mean<double, double>(rescontainer);
}

void DRNonHomogeneousMixedTreeLikelihood::resetLikelihoodArrays(const Node* node)
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->resetLikelihoodArrays(node);
  }
}

void DRNonHomogeneousMixedTreeLikelihood::computeTreeLikelihood()
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeTreeLikelihood();
  }
}

void DRNonHomogeneousMixedTreeLikelihood::computeSubtreeLikelihoodPostfix(const Node* node)
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeSubtreeLikelihoodPostfix(node);
  }
}

void DRNonHomogeneousMixedTreeLikelihood::computeSubtreeLikelihoodPrefix(const Node* node)
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeSubtreeLikelihoodPostfix(node);
  }
}

void DRNonHomogeneousMixedTreeLikelihood::computeRootLikelihood()
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeRootLikelihood();
  }
}

void DRNonHomogeneousMixedTreeLikelihood::computeLikelihoodAtNode(const Node * node, VVVdouble& likelihoodArray) const
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeLikelihoodAtNode(node, likelihoodArray);
  }
}

void DRNonHomogeneousMixedTreeLikelihood::computeLikelihoodFromArrays(
  const vector<const VVVdouble*>& iLik,
  const vector<const VVVdouble*>& tProb,
  VVVdouble& oLik,
  unsigned int nbNodes,
  unsigned int nbDistinctSites,
  unsigned int nbClasses,
  unsigned int nbStates,
  bool reset)
{
  cerr << "Invalid use of static method DRNonHomogeneousMixedTreeLikelihood::computeLikelihoodFromArrays" << endl;
  exit(0);
}

void DRNonHomogeneousMixedTreeLikelihood::computeLikelihoodFromArrays(
  const std::vector<const VVVdouble*>& iLik,
  const std::vector<const VVVdouble*>& tProb,
  const VVVdouble* iLikR,
  const VVVdouble* tProbR,
  VVVdouble& oLik,
  unsigned int nbNodes,
  unsigned int nbDistinctSites,
  unsigned int nbClasses,
  unsigned int nbStates,
  bool reset)
{
  cerr << "Invalid use of static method DRNonHomogeneousMixedTreeLikelihood::computeLikelihoodFromArrays" << endl;
  exit(0);
}

void DRNonHomogeneousMixedTreeLikelihood::displayLikelihood(const Node* node)
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->displayLikelihood(node);
  }
}
