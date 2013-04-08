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
  FixedFrequenciesSet* fSet = new FixedFrequenciesSet(model->getAlphabet(), model->getFrequencies());
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
Site* NonHomogeneousSequenceSimulator::simulate() const
{
  // Draw an initial state randomly according to equilibrum frequencies:
  int initialState = 0;
  double r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
  double cumprob = 0;
  vector<double> freqs = modelSet_->getRootFrequencies();
  for (size_t i = 0; i < nbStates_; i++)
  {
    cumprob += freqs[i];
    if (r <= cumprob)
    {
      initialState = static_cast<int>(i);
      break;
    }
  }
  return simulate(initialState);
}

/******************************************************************************/
Site* NonHomogeneousSequenceSimulator::simulate(int initialState) const
{
  if (continuousRates_)
  {
    // Draw a random rate:
    double rate = rate_->randC();
    // Make this state evolve:
    return simulate(initialState, rate);
  }
  else
  {
    // Draw a random rate:
    size_t rateClass = static_cast<size_t>(RandomTools::giveIntRandomNumberBetweenZeroAndEntry(static_cast<int>(rate_->getNumberOfCategories())));
    // Make this state evolve:
    return simulate(initialState, rateClass);
  }
}

/******************************************************************************/
Site* NonHomogeneousSequenceSimulator::simulate(int initialState, size_t rateClass) const
{
  // Launch recursion:
  SNode* root = tree_.getRootNode();
  root->getInfos().state = initialState;
  for (size_t i = 0; i < root->getNumberOfSons(); i++)
  {
    evolveInternal(root->getSon(i), rateClass);
  }
  // Now create a Site object:
  Vint site(leaves_.size());
  for (size_t i = 0; i < leaves_.size(); i++)
  {
    site[i] = leaves_[i]->getInfos().model->getAlphabetChar(leaves_[i]->getInfos().state);
  }
  return new Site(site, alphabet_);
}

/******************************************************************************/
Site* NonHomogeneousSequenceSimulator::simulate(int initialState, double rate) const
{
  // Launch recursion:
  SNode* root = tree_.getRootNode();
  root->getInfos().state = initialState;
  for (size_t i = 0; i < root->getNumberOfSons(); i++)
  {
    evolveInternal(root->getSon(i), rate);
  }
  // Now create a Site object:
  Vint site(leaves_.size());
  for (size_t i = 0; i < leaves_.size(); i++)
  {
    site[i] = leaves_[i]->getInfos().model->getAlphabetChar(leaves_[i]->getInfos().state);
  }
  return new Site(site, alphabet_);
}

/******************************************************************************/
Site* NonHomogeneousSequenceSimulator::simulate(double rate) const
{
  // Draw an initial state randomly according to equilibrum frequencies:
  int initialState = 0;
  double r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
  double cumprob = 0;
  vector<double> freqs = modelSet_->getRootFrequencies();
  for (size_t i = 0; i < nbStates_; i++)
  {
    cumprob += freqs[i];
    if (r <= cumprob)
    {
      initialState = (int)i;
      break;
    }
  }
  // Make this state evolve:
  return simulate(initialState, rate);
}

/******************************************************************************/
SiteContainer* NonHomogeneousSequenceSimulator::simulate(size_t numberOfSites) const
{
  Vint initialStates(numberOfSites, 0);
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
        initialStates[j] = (int)i;
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
      Site* site = simulate();
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
      rateClasses[j] = static_cast<size_t>(RandomTools::giveIntRandomNumberBetweenZeroAndEntry(static_cast<int>(nCat)));
    }
    // Make these states evolve:
    SiteContainer* sites = multipleEvolve(initialStates, rateClasses);
    return sites;
  }
}

/******************************************************************************/
RASiteSimulationResult* NonHomogeneousSequenceSimulator::dSimulate() const
{
  // Draw an initial state randomly according to equilibrum frequencies:
  int initialState = 0;
  double r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
  double cumprob = 0;
  vector<double> freqs = modelSet_->getRootFrequencies();
  for (size_t i = 0; i < nbStates_; i++)
  {
    cumprob += freqs[i];
    if (r <= cumprob)
    {
      initialState = static_cast<int>(i);
      break;
    }
  }

  return dSimulate(initialState);
}

/******************************************************************************/
RASiteSimulationResult* NonHomogeneousSequenceSimulator::dSimulate(int initialState) const
{
  // Draw a random rate:
  if (continuousRates_)
  {
    double rate = rate_->randC();
    return dSimulate(initialState, rate);
  }
  else
  {
    size_t rateClass = static_cast<size_t>(RandomTools::giveIntRandomNumberBetweenZeroAndEntry(static_cast<int>(rate_->getNumberOfCategories())));
    return dSimulate(initialState, rateClass);
    // NB: this is more efficient than dSimulate(initialState, rDist_->rand())
  }
}

/******************************************************************************/
RASiteSimulationResult* NonHomogeneousSequenceSimulator::dSimulate(int initialState, double rate) const
{
  // Make this state evolve:
  RASiteSimulationResult* hssr = new RASiteSimulationResult(templateTree_, modelSet_->getAlphabet(), initialState, rate);
  dEvolve(initialState, rate, *hssr);
  return hssr;
}

/******************************************************************************/
RASiteSimulationResult* NonHomogeneousSequenceSimulator::dSimulate(int initialState, size_t rateClass) const
{
  return dSimulate(initialState, rate_->getCategory(rateClass));
}

/******************************************************************************/
RASiteSimulationResult* NonHomogeneousSequenceSimulator::dSimulate(double rate) const
{
  // Draw an initial state randomly according to equilibrum frequencies:
  int initialState = 0;
  double r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
  double cumprob = 0;
  vector<double> freqs = modelSet_->getRootFrequencies();
  for (size_t i = 0; i < nbStates_; i++)
  {
    cumprob += freqs[i];
    if (r <= cumprob)
    {
      initialState = static_cast<int>(i);
      break;
    }
  }
  return dSimulate(initialState, rate);
}

/******************************************************************************/
int NonHomogeneousSequenceSimulator::evolve(const SNode* node, int initialState, size_t rateClass) const
{
  const Vdouble* cumpxy_node_c_x_ = &node->getInfos().cumpxy[rateClass][initialState];
  double rand = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
  for (int y = 0; y < static_cast<int>(nbStates_); y++)
  {
    if (rand < (*cumpxy_node_c_x_)[y]) return y;
  }
  cerr << "DEBUG: This message should never happen! (HomogeneousSequenceSimulator::evolve)" << endl;
  cout << "   rand = " << rand << endl;
  return -1;
}

/******************************************************************************/
int NonHomogeneousSequenceSimulator::evolve(const SNode* node, int initialState, double rate) const
{
  double cumpxy = 0;
  double rand = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
  double l = rate * node->getDistanceToFather();
  const SubstitutionModel* model = node->getInfos().model;
  for (int y = 0; y < static_cast<int>(nbStates_); y++)
  {
    cumpxy += model->Pij_t(initialState, y, l);
    if (rand < cumpxy) return y;
  }
  cerr << "DEBUG: This message should never happen! (NonHomogeneousSequenceSimulator::evolve)" << endl;
  cout << "  rand = " << rand << endl;
  cout << "cumpxy = " << cumpxy << endl;
  MatrixTools::print(model->getPij_t(l));
  return -1;
}

/******************************************************************************/
void NonHomogeneousSequenceSimulator::multipleEvolve(const SNode* node, const Vint& initialStates, const vector<size_t>& rateClasses, Vint& finalStates) const
{
  const VVVdouble* cumpxy_node_ = &node->getInfos().cumpxy;
  for (size_t i = 0; i < initialStates.size(); i++)
  {
    const Vdouble* cumpxy_node_c_x_ = &(*cumpxy_node_)[rateClasses[i]][initialStates[i]];
    double rand = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
    for (size_t y = 0; y < nbStates_; y++)
    {
      if (rand < (*cumpxy_node_c_x_)[y])
      {
        finalStates[i] = (int)y;
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
  const vector<int>* initialStates = &node->getFather()->getInfos().states;
  size_t n = initialStates->size();
  node->getInfos().states.resize(n); // allocation.
  multipleEvolve(node, node->getFather()->getInfos().states, rateClasses, node->getInfos().states);
  for (size_t i = 0; i < node->getNumberOfSons(); i++)
  {
    multipleEvolveInternal(node->getSon(i), rateClasses);
  }
}

/******************************************************************************/
SiteContainer* NonHomogeneousSequenceSimulator::multipleEvolve(const Vint& initialStates, const vector<size_t>& rateClasses) const
{
  // Launch recursion:
  SNode* root = tree_.getRootNode();
  root->getInfos().states = initialStates;
  for (size_t i = 0; i < root->getNumberOfSons(); i++)
  {
    multipleEvolveInternal(root->getSon(i), rateClasses);
  }
  // Now create a SiteContainer object:
  AlignedSequenceContainer* sites = new AlignedSequenceContainer(alphabet_);
  size_t n = leaves_.size();
  size_t nbSites = initialStates.size();
  const SubstitutionModel* model = 0;
  for (size_t i = 0; i < n; i++)
  {
    vector<int> content(nbSites);
    vector<int>* states = &leaves_[i]->getInfos().states;
    model = leaves_[i]->getInfos().model;
    for (size_t j = 0; j < nbSites; j++)
    {
      content[j] = model->getAlphabetChar((*states)[j]);
    }
    sites->addSequence(BasicSequence(leaves_[i]->getName(), content, alphabet_), false);
  }
  return sites;
}

/******************************************************************************/
void NonHomogeneousSequenceSimulator::dEvolve(int initialState, double rate, RASiteSimulationResult& rassr) const
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

