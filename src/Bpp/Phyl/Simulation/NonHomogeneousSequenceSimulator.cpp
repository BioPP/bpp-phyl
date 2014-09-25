//
// File: NonHomogeneousSequenceSimulator.cpp
//       (previously SequenceSimulator.cpp, then HomogeneousSequenceSimulator.cpp)
// Created by: Julien Dutheil
//             Bastien Boussau
// Created on: Wed Feb  4 16:30:51 2004
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

#include "NonHomogeneousSequenceSimulator.h"
#include "../Model/SubstitutionModelSetTools.h"

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>

// From SeqLib:
#include <Bpp/Seq/Container/VectorSiteContainer.h>

using namespace bpp;
using namespace std;

/******************************************************************************/

NonHomogeneousSequenceSimulator::NonHomogeneousSequenceSimulator(
  const SubstitutionModelSet* modelSet,
  const DiscreteDistribution* rate,
  const Tree* tree) throw (Exception) :
  modelSet_(modelSet),
  alphabet_(modelSet_->getAlphabet()),
  supportedStates_(modelSet_->getAlphabetStates()),
  rate_(rate),
  templateTree_(tree),
  tree_(*tree),
  ownModelSet_(false),
  leaves_(tree_.getLeaves()),
  seqNames_(),
  nbNodes_(),
  nbClasses_(rate_->getNumberOfCategories()),
  nbStates_(modelSet_->getNumberOfStates()),
  continuousRates_(false)
{
  if (!modelSet->isFullySetUpFor(*tree))
    throw Exception("NonHomogeneousSequenceSimulator(constructor). Model set is not fully specified.");
  init();
}

/******************************************************************************/

NonHomogeneousSequenceSimulator::NonHomogeneousSequenceSimulator(
  const SubstitutionModel* model,
  const DiscreteDistribution* rate,
  const Tree* tree) :
  modelSet_(0),
  alphabet_(model->getAlphabet()),
  supportedStates_(model->getAlphabetStates()),
  rate_(rate),
  templateTree_(tree),
  tree_(*tree),
  ownModelSet_(true),
  leaves_(tree_.getLeaves()),
  seqNames_(),
  nbNodes_(),
  nbClasses_(rate_->getNumberOfCategories()),
  nbStates_(model->getNumberOfStates()),
  continuousRates_(false)
{
  FixedFrequenciesSet* fSet = new FixedFrequenciesSet(model->getStateMap().clone(), model->getFrequencies());
  fSet->setNamespace("anc.");
  modelSet_ = SubstitutionModelSetTools::createHomogeneousModelSet(dynamic_cast<SubstitutionModel*>(model->clone()), fSet, templateTree_);
  init();
}

/******************************************************************************/

void NonHomogeneousSequenceSimulator::init()
{
  seqNames_.resize(leaves_.size());
  for (size_t i = 0; i < seqNames_.size(); i++)
  {
    seqNames_[i] = leaves_[i]->getName();
  }
  // Initialize cumulative pxy:
  vector<SNode*> nodes = tree_.getNodes();
  nodes.pop_back(); // remove root
  nbNodes_ = nodes.size();

  for (size_t i = 0; i < nodes.size(); i++)
  {
    SNode* node = nodes[i];
    node->getInfos().model = modelSet_->getModelForNode(node->getId());
    double d = node->getDistanceToFather();
    VVVdouble* cumpxy_node_ = &node->getInfos().cumpxy;
    cumpxy_node_->resize(nbClasses_);
    for (size_t c = 0; c < nbClasses_; c++)
    {
      VVdouble* cumpxy_node_c_ = &(*cumpxy_node_)[c];
      cumpxy_node_c_->resize(nbStates_);
      RowMatrix<double> P = node->getInfos().model->getPij_t(d * rate_->getCategory(c));
      for (size_t x = 0; x < nbStates_; x++)
      {
        Vdouble* cumpxy_node_c_x_ = &(*cumpxy_node_c_)[x];
        cumpxy_node_c_x_->resize(nbStates_);
        (*cumpxy_node_c_x_)[0] = P(x, 0);
        for (size_t y = 1; y < nbStates_; y++)
        {
          (*cumpxy_node_c_x_)[y] = (*cumpxy_node_c_x_)[y - 1] + P(x, y);
        }
      }
    }
  }
}

/******************************************************************************/

Site* NonHomogeneousSequenceSimulator::simulateSite() const
{
  // Draw an initial state randomly according to equilibrum frequencies:
  size_t initialStateIndex = 0;
  double r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
  double cumprob = 0;
  vector<double> freqs = modelSet_->getRootFrequencies();
  for (size_t i = 0; i < nbStates_; i++)
  {
    cumprob += freqs[i];
    if (r <= cumprob)
    {
      initialStateIndex = i;
      break;
    }
  }
  return simulateSite(initialStateIndex);
}

/******************************************************************************/

Site* NonHomogeneousSequenceSimulator::simulateSite(size_t ancestralStateIndex) const
{
  if (continuousRates_)
  {
    // Draw a random rate:
    double rate = rate_->randC();
    // Make this state evolve:
    return simulateSite(ancestralStateIndex, rate);
  }
  else
  {
    // Draw a random rate:
    size_t rateClass = RandomTools::giveIntRandomNumberBetweenZeroAndEntry<size_t>(rate_->getNumberOfCategories());
    // Make this state evolve:
    return simulateSite(ancestralStateIndex, rateClass);
  }
}

/******************************************************************************/

Site* NonHomogeneousSequenceSimulator::simulateSite(size_t ancestralStateIndex, size_t rateClass) const
{
  // Launch recursion:
  SNode* root = tree_.getRootNode();
  root->getInfos().state = ancestralStateIndex;
  for (size_t i = 0; i < root->getNumberOfSons(); ++i)
  {
    evolveInternal(root->getSon(i), rateClass);
  }
  // Now create a Site object:
  Vint site(leaves_.size());
  for (size_t i = 0; i < leaves_.size(); ++i)
  {
    site[i] = leaves_[i]->getInfos().model->getAlphabetStateAsInt(leaves_[i]->getInfos().state);
  }
  return new Site(site, alphabet_);
}

/******************************************************************************/

Site* NonHomogeneousSequenceSimulator::simulateSite(size_t ancestralStateIndex, double rate) const
{
  // Launch recursion:
  SNode* root = tree_.getRootNode();
  root->getInfos().state = ancestralStateIndex;
  for (size_t i = 0; i < root->getNumberOfSons(); i++)
  {
    evolveInternal(root->getSon(i), rate);
  }
  // Now create a Site object:
  Vint site(leaves_.size());
  for (size_t i = 0; i < leaves_.size(); i++)
  {
    site[i] = leaves_[i]->getInfos().model->getAlphabetStateAsInt(leaves_[i]->getInfos().state);
  }
  return new Site(site, alphabet_);
}

/******************************************************************************/

Site* NonHomogeneousSequenceSimulator::simulateSite(double rate) const
{
  // Draw an initial state randomly according to equilibrum frequencies:
  size_t ancestralStateIndex = 0;
  double r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
  double cumprob = 0;
  vector<double> freqs = modelSet_->getRootFrequencies();
  for (size_t i = 0; i < nbStates_; i++)
  {
    cumprob += freqs[i];
    if (r <= cumprob)
    {
      ancestralStateIndex = i;
      break;
    }
  }
  // Make this state evolve:
  return simulateSite(ancestralStateIndex, rate);
}

/******************************************************************************/

SiteContainer* NonHomogeneousSequenceSimulator::simulate(size_t numberOfSites) const
{
  vector<size_t> ancestralStateIndices(numberOfSites, 0);
  for (size_t j = 0; j < numberOfSites; j++)
  {
    double r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
    double cumprob = 0;
    vector<double> freqs = modelSet_->getRootFrequencies();
    for (size_t i = 0; i < nbStates_; i++)
    {
      cumprob += freqs[i];
      if (r <= cumprob)
      {
        ancestralStateIndices[j] = i;
        break;
      }
    }
  }
  if (continuousRates_)
  {
    VectorSiteContainer* sites = new VectorSiteContainer(seqNames_.size(), alphabet_);
    sites->setSequencesNames(seqNames_);
    for (size_t j = 0; j < numberOfSites; j++)
    {
      Site* site = simulateSite();
      site->setPosition(static_cast<int>(j));
      sites->addSite(*site);
      delete site;
    }
    return sites;
  }
  else
  {
    // More efficient to do site this way:
    // Draw random rates:
    vector<size_t> rateClasses(numberOfSites);
    size_t nCat = rate_->getNumberOfCategories();
    for (size_t j = 0; j < numberOfSites; j++)
    {
      rateClasses[j] = RandomTools::giveIntRandomNumberBetweenZeroAndEntry<size_t>(nCat);
    }
    // Make these states evolve:
    SiteContainer* sites = multipleEvolve(ancestralStateIndices, rateClasses);
    return sites;
  }
}

/******************************************************************************/

RASiteSimulationResult* NonHomogeneousSequenceSimulator::dSimulateSite() const
{
  // Draw an initial state randomly according to equilibrum frequencies:
  size_t ancestralStateIndex = 0;
  double r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
  double cumprob = 0;
  vector<double> freqs = modelSet_->getRootFrequencies();
  for (size_t i = 0; i < nbStates_; i++)
  {
    cumprob += freqs[i];
    if (r <= cumprob)
    {
      ancestralStateIndex = i;
      break;
    }
  }

  return dSimulateSite(ancestralStateIndex);
}

/******************************************************************************/

RASiteSimulationResult* NonHomogeneousSequenceSimulator::dSimulateSite(size_t ancestralStateIndex) const
{
  // Draw a random rate:
  if (continuousRates_)
  {
    double rate = rate_->randC();
    return dSimulateSite(ancestralStateIndex, rate);
  }
  else
  {
    size_t rateClass = RandomTools::giveIntRandomNumberBetweenZeroAndEntry<size_t>(rate_->getNumberOfCategories());
    return dSimulateSite(ancestralStateIndex, rateClass);
    // NB: this is more efficient than dSimulate(initialState, rDist_->rand())
  }
}

/******************************************************************************/

RASiteSimulationResult* NonHomogeneousSequenceSimulator::dSimulateSite(size_t ancestralStateIndex, double rate) const
{
  // Make this state evolve:
  RASiteSimulationResult* hssr = new RASiteSimulationResult(templateTree_, modelSet_->getAlphabet(), ancestralStateIndex, rate);
  dEvolve(ancestralStateIndex, rate, *hssr);
  return hssr;
}

/******************************************************************************/

RASiteSimulationResult* NonHomogeneousSequenceSimulator::dSimulateSite(size_t ancestralStateIndex, size_t rateClass) const
{
  return dSimulateSite(ancestralStateIndex, rate_->getCategory(rateClass));
}

/******************************************************************************/

RASiteSimulationResult* NonHomogeneousSequenceSimulator::dSimulateSite(double rate) const
{
  // Draw an initial state randomly according to equilibrum frequencies:
  size_t ancestralStateIndex = 0;
  double r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
  double cumprob = 0;
  vector<double> freqs = modelSet_->getRootFrequencies();
  for (size_t i = 0; i < nbStates_; i++)
  {
    cumprob += freqs[i];
    if (r <= cumprob)
    {
      ancestralStateIndex = i;
      break;
    }
  }
  return dSimulateSite(ancestralStateIndex, rate);
}

/******************************************************************************/

size_t NonHomogeneousSequenceSimulator::evolve(const SNode* node, size_t initialStateIndex, size_t rateClass) const
{
  const Vdouble* cumpxy_node_c_x_ = &node->getInfos().cumpxy[rateClass][initialStateIndex];
  double rand = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
  for (size_t y = 0; y < nbStates_; y++)
  {
    if (rand < (*cumpxy_node_c_x_)[y]) return y;
  }
  throw Exception("HomogeneousSequenceSimulator::evolve. The impossible happened! rand = " + TextTools::toString(rand) + ".");
}

/******************************************************************************/

size_t NonHomogeneousSequenceSimulator::evolve(const SNode* node, size_t initialStateIndex, double rate) const
{
  double cumpxy = 0;
  double rand = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
  double l = rate * node->getDistanceToFather();
  const SubstitutionModel* model = node->getInfos().model;
  for (size_t y = 0; y < nbStates_; y++)
  {
    cumpxy += model->Pij_t(initialStateIndex, y, l);
    if (rand < cumpxy) return y;
  }
  MatrixTools::print(model->getPij_t(l));
  throw Exception("HomogeneousSequenceSimulator::evolve. The impossible happened! rand = " + TextTools::toString(rand) + ".");
}

/******************************************************************************/

void NonHomogeneousSequenceSimulator::multipleEvolve(
    const SNode* node,
    const std::vector<size_t>& initialStateIndices,
    const vector<size_t>& rateClasses,
    std::vector<size_t>& finalStateIndices) const
{
  const VVVdouble* cumpxy_node_ = &node->getInfos().cumpxy;
  for (size_t i = 0; i < initialStateIndices.size(); i++)
  {
    const Vdouble* cumpxy_node_c_x_ = &(*cumpxy_node_)[rateClasses[i]][initialStateIndices[i]];
    double rand = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
    for (size_t y = 0; y < nbStates_; y++)
    {
      if (rand < (*cumpxy_node_c_x_)[y])
      {
        finalStateIndices[i] = y;
        break;
      }
    }
  }
}

/******************************************************************************/

void NonHomogeneousSequenceSimulator::evolveInternal(SNode* node, size_t rateClass) const
{
  if (!node->hasFather())
  {
    cerr << "DEBUG: NonHomogeneousSequenceSimulator::evolveInternal. Forbidden call of method on root node." << endl;
    return;
  }
  node->getInfos().state = evolve(node, node->getFather()->getInfos().state, rateClass);
  for (size_t i = 0; i < node->getNumberOfSons(); i++)
  {
    evolveInternal(node->getSon(i), rateClass);
  }
}

/******************************************************************************/

void NonHomogeneousSequenceSimulator::evolveInternal(SNode* node, double rate) const
{
  if (!node->hasFather())
  {
    cerr << "DEBUG: NonHomogeneousSequenceSimulator::evolveInternal. Forbidden call of method on root node." << endl;
    return;
  }
  node->getInfos().state = evolve(node, node->getFather()->getInfos().state, rate);
  for (size_t i = 0; i < node->getNumberOfSons(); i++)
  {
    evolveInternal(node->getSon(i), rate);
  }
}

/******************************************************************************/

void NonHomogeneousSequenceSimulator::multipleEvolveInternal(SNode* node, const vector<size_t>& rateClasses) const
{
  if (!node->hasFather())
  {
    cerr << "DEBUG: NonHomogeneousSequenceSimulator::multipleEvolveInternal. Forbidden call of method on root node." << endl;
    return;
  }
  const vector<size_t>* initialStates = &node->getFather()->getInfos().states;
  size_t n = initialStates->size();
  node->getInfos().states.resize(n); // allocation.
  multipleEvolve(node, node->getFather()->getInfos().states, rateClasses, node->getInfos().states);
  for (size_t i = 0; i < node->getNumberOfSons(); i++)
  {
    multipleEvolveInternal(node->getSon(i), rateClasses);
  }
}

/******************************************************************************/

SiteContainer* NonHomogeneousSequenceSimulator::multipleEvolve(
    const std::vector<size_t>& initialStateIndices,
    const vector<size_t>& rateClasses) const
{
  // Launch recursion:
  SNode* root = tree_.getRootNode();
  root->getInfos().states = initialStateIndices;
  for (size_t i = 0; i < root->getNumberOfSons(); i++)
  {
    multipleEvolveInternal(root->getSon(i), rateClasses);
  }
  // Now create a SiteContainer object:
  AlignedSequenceContainer* sites = new AlignedSequenceContainer(alphabet_);
  size_t n = leaves_.size();
  size_t nbSites = initialStateIndices.size();
  const SubstitutionModel* model = 0;
  for (size_t i = 0; i < n; i++)
  {
    vector<int> content(nbSites);
    vector<size_t>& states = leaves_[i]->getInfos().states;
    model = leaves_[i]->getInfos().model;
    for (size_t j = 0; j < nbSites; j++)
    {
      content[j] = model->getAlphabetStateAsInt(states[j]);
    }
    sites->addSequence(BasicSequence(leaves_[i]->getName(), content, alphabet_), false);
  }
  return sites;
}

/******************************************************************************/

void NonHomogeneousSequenceSimulator::dEvolve(size_t initialState, double rate, RASiteSimulationResult& rassr) const
{
  // Launch recursion:
  SNode* root = tree_.getRootNode();
  root->getInfos().state = initialState;
  for (size_t i = 0; i < root->getNumberOfSons(); i++)
  {
    dEvolveInternal(root->getSon(i), rate, rassr);
  }
}

/******************************************************************************/

void NonHomogeneousSequenceSimulator::dEvolveInternal(SNode* node, double rate, RASiteSimulationResult& rassr) const
{
  if (!node->hasFather())
  {
    cerr << "DEBUG: NonHomogeneousSequenceSimulator::evolveInternal. Forbidden call of method on root node." << endl;
    return;
  }
  SimpleMutationProcess process(node->getInfos().model);
  MutationPath mp = process.detailedEvolve(node->getFather()->getInfos().state, node->getDistanceToFather() * rate);
  node->getInfos().state = mp.getFinalState();

  // Now append infos in rassr:
  rassr.addNode(node->getId(), mp);

  // Now jump to son nodes:
  for (size_t i = 0; i < node->getNumberOfSons(); i++)
  {
    dEvolveInternal(node->getSon(i), rate, rassr);
  }
}

/******************************************************************************/

