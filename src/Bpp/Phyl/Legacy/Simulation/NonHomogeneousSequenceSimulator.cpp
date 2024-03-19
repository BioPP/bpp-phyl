// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "NonHomogeneousSequenceSimulator.h"
#include "../Model/SubstitutionModelSetTools.h"
#include "../../Tree/PhyloTreeTools.h"

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>

// From bpp-seq:
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/AlignedSequenceContainer.h>

using namespace bpp;
using namespace std;

/******************************************************************************/

NonHomogeneousSequenceSimulator::NonHomogeneousSequenceSimulator(
    std::shared_ptr<const SubstitutionModelSet> modelSet,
    std::shared_ptr<const DiscreteDistributionInterface> rate,
    std::shared_ptr<const Tree> tree) :
  modelSet_(modelSet),
  alphabet_(modelSet_->getAlphabet()),
  supportedStates_(modelSet_->getAlphabetStates()),
  rate_(rate),
  templateTree_(tree),
  tree_(*tree),
  phyloTree_(make_shared<const ParametrizablePhyloTree>(*PhyloTreeTools::buildFromTreeTemplate(*tree))),
  leaves_(tree_.getLeaves()),
  seqNames_(),
  nbNodes_(),
  nbClasses_(rate_->getNumberOfCategories()),
  nbStates_(modelSet_->getNumberOfStates()),
  continuousRates_(false),
  outputInternalSequences_(false)
{
  if (!modelSet->isFullySetUpFor(*tree))
    throw Exception("NonHomogeneousSequenceSimulator(constructor). Model set is not fully specified.");
  init();
}

/******************************************************************************/

NonHomogeneousSequenceSimulator::NonHomogeneousSequenceSimulator(
    std::shared_ptr<const TransitionModelInterface> model,
    std::shared_ptr<const DiscreteDistributionInterface> rate,
    std::shared_ptr<const Tree> tree) :
  modelSet_(0),
  alphabet_(model->getAlphabet()),
  supportedStates_(model->getAlphabetStates()),
  rate_(rate),
  templateTree_(tree),
  tree_(*tree),
  phyloTree_(make_shared<const ParametrizablePhyloTree>(*PhyloTreeTools::buildFromTreeTemplate(*tree))),
  leaves_(tree_.getLeaves()),
  seqNames_(),
  nbNodes_(),
  nbClasses_(rate_->getNumberOfCategories()),
  nbStates_(model->getNumberOfStates()),
  continuousRates_(false),
  outputInternalSequences_(false)
{
  auto fSet = make_shared<FixedFrequencySet>(model->getStateMap(), model->getFrequencies());
  fSet->setNamespace("anc.");
  modelSet_ = SubstitutionModelSetTools::createHomogeneousModelSet(shared_ptr<TransitionModelInterface>(model->clone()), fSet, *templateTree_);
  init();
}

/******************************************************************************/

void NonHomogeneousSequenceSimulator::init()
{
  vector<SNode*> nodes = tree_.getNodes();

  if (outputInternalSequences_)
  {
    seqNames_.resize(nodes.size());
    for (size_t i = 0; i < seqNames_.size(); i++)
    {
      if (nodes[i]->isLeaf())
      {
        seqNames_[i] = nodes[i]->getName();
      }
      else
      {
        seqNames_[i] = TextTools::toString( nodes[i]->getId() );
      }
    }
  }
  else
  {
    seqNames_.resize(leaves_.size());
    for (size_t i = 0; i < seqNames_.size(); i++)
    {
      seqNames_[i] = leaves_[i]->getName();
    }
  }
  // Initialize cumulative pxy:
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

unique_ptr<Site> NonHomogeneousSequenceSimulator::simulateSite() const
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

unique_ptr<Site> NonHomogeneousSequenceSimulator::simulateSite(size_t ancestralStateIndex) const
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

std::unique_ptr<Site> NonHomogeneousSequenceSimulator::simulateSite(size_t ancestralStateIndex, size_t rateClass) const
{
  // Launch recursion:
  SNode* root = tree_.getRootNode();
  root->getInfos().state = ancestralStateIndex;
  for (size_t i = 0; i < root->getNumberOfSons(); ++i)
  {
    evolveInternal(root->getSon(i), rateClass);
  }
  // Now create a Site object:
  size_t n = nbNodes_ + 1;
  if (!outputInternalSequences_)
  {
    n = leaves_.size();
  }
  Vint site(n);

  if (outputInternalSequences_)
  {
    vector<SNode*> nodes = tree_.getNodes();
    for (size_t i = 0; i < n; i++)
    {
      if (i == n - 1)   // We take the model of node i-1 because the root has no model
      {
        site[i] = nodes[i - 1]->getInfos().model->getAlphabetStateAsInt(nodes[i]->getInfos().state);
      }
      else
      {
        site[i] = nodes[i]->getInfos().model->getAlphabetStateAsInt(nodes[i]->getInfos().state);
      }
    }
  }
  else
  {
    for (size_t i = 0; i < leaves_.size(); ++i)
    {
      site[i] = leaves_[i]->getInfos().model->getAlphabetStateAsInt(leaves_[i]->getInfos().state);
    }
  }
  return make_unique<Site>(site, alphabet_);
}

/******************************************************************************/

std::unique_ptr<Site> NonHomogeneousSequenceSimulator::simulateSite(size_t ancestralStateIndex, double rate) const
{
  // Launch recursion:
  SNode* root = tree_.getRootNode();
  root->getInfos().state = ancestralStateIndex;
  for (size_t i = 0; i < root->getNumberOfSons(); i++)
  {
    evolveInternal(root->getSon(i), rate);
  }
  // Now create a Site object:
  size_t n = nbNodes_ + 1;
  if (!outputInternalSequences_)
  {
    n = leaves_.size();
  }
  Vint site(n);

  if (outputInternalSequences_)
  {
    vector<SNode*> nodes = tree_.getNodes();
    for (size_t i = 0; i < n; i++)
    {
      if (i == n - 1)   // We take the model of node i-1 because the root has no model
      {
        site[i] = nodes[i - 1]->getInfos().model->getAlphabetStateAsInt(nodes[i]->getInfos().state);
      }
      else
      {
        site[i] = nodes[i]->getInfos().model->getAlphabetStateAsInt(nodes[i]->getInfos().state);
      }
    }
  }
  else
  {
    for (size_t i = 0; i < leaves_.size(); i++)
    {
      site[i] = leaves_[i]->getInfos().model->getAlphabetStateAsInt(leaves_[i]->getInfos().state);
    }
  }
  return make_unique<Site>(site, alphabet_);
}

/******************************************************************************/

unique_ptr<Site> NonHomogeneousSequenceSimulator::simulateSite(double rate) const
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

unique_ptr<SiteContainerInterface> NonHomogeneousSequenceSimulator::simulate(size_t numberOfSites) const
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
    auto sites = make_unique<VectorSiteContainer>(seqNames_.size(), alphabet_);
    sites->setSequenceNames(seqNames_, true);
    for (size_t j = 0; j < numberOfSites; j++)
    {
      auto site = simulateSite();
      site->setCoordinate(static_cast<int>(j));
      sites->addSite(site);
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
    return multipleEvolve(ancestralStateIndices, rateClasses);
  }
}

/******************************************************************************/

std::unique_ptr<SiteSimulationResult> NonHomogeneousSequenceSimulator::dSimulateSite() const
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

std::unique_ptr<SiteSimulationResult> NonHomogeneousSequenceSimulator::dSimulateSite(size_t ancestralStateIndex) const
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

unique_ptr<SiteSimulationResult> NonHomogeneousSequenceSimulator::dSimulateSite(size_t ancestralStateIndex, double rate) const
{
  // Make this state evolve:
  auto hssr = make_unique<RASiteSimulationResult>(phyloTree_, modelSet_->getStateMap(), ancestralStateIndex, rate);
  dEvolve(ancestralStateIndex, rate, *hssr);
  return std::move(hssr);
}

/******************************************************************************/

unique_ptr<SiteSimulationResult> NonHomogeneousSequenceSimulator::dSimulateSite(size_t ancestralStateIndex, size_t rateClass) const
{
  return dSimulateSite(ancestralStateIndex, rate_->getCategory(rateClass));
}

/******************************************************************************/

unique_ptr<SiteSimulationResult> NonHomogeneousSequenceSimulator::dSimulateSite(double rate) const
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
    if (rand < (*cumpxy_node_c_x_)[y])
      return y;
  }
  throw Exception("HomogeneousSequenceSimulator::evolve. The impossible happened! rand = " + TextTools::toString(rand) + ".");
}

/******************************************************************************/

size_t NonHomogeneousSequenceSimulator::evolve(const SNode* node, size_t initialStateIndex, double rate) const
{
  double cumpxy = 0;
  double rand = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
  double l = rate * node->getDistanceToFather();
  auto model = node->getInfos().model;
  for (size_t y = 0; y < nbStates_; y++)
  {
    cumpxy += model->Pij_t(initialStateIndex, y, l);
    if (rand < cumpxy)
      return y;
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

unique_ptr<SiteContainerInterface> NonHomogeneousSequenceSimulator::multipleEvolve(
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
  auto sites = make_unique<AlignedSequenceContainer>(alphabet_);
  size_t nbSites = initialStateIndices.size();
  shared_ptr<const TransitionModelInterface> model = nullptr;
  if (outputInternalSequences_)
  {
    vector<SNode*> nodes = tree_.getNodes();
    size_t n = nbNodes_ + 1;
    for (size_t i = 0; i < n; i++)
    {
      vector<int> content(nbSites);
      vector<size_t>& states = nodes[i]->getInfos().states;
      if (i == n - 1)   // If at the root, there is no model, so we take the model of node n-1.
      {
        model = nodes[i - 1]->getInfos().model;
      }
      else
      {
        model = nodes[i]->getInfos().model;
      }
      for (size_t j = 0; j < nbSites; j++)
      {
        content[j] = model->getAlphabetStateAsInt(states[j]);
      }
      if (nodes[i]->isLeaf())
      {
        auto seq = make_unique<Sequence>(nodes[i]->getName(), content, alphabet_);
        sites->addSequence(nodes[i]->getName(), seq);
      }
      else
      {
        auto seq = make_unique<Sequence>(TextTools::toString(nodes[i]->getId()), content, alphabet_);
        sites->addSequence(TextTools::toString(nodes[i]->getId()), seq);
      }
    }
  }
  else
  {
    size_t n = leaves_.size();
    for (size_t i = 0; i < n; i++)
    {
      vector<int> content(nbSites);
      vector<size_t>& states = leaves_[i]->getInfos().states;
      model = leaves_[i]->getInfos().model;
      for (size_t j = 0; j < nbSites; j++)
      {
        content[j] = model->getAlphabetStateAsInt(states[j]);
      }
      auto seq = make_unique<Sequence>(leaves_[i]->getName(), content, alphabet_);
      sites->addSequence(leaves_[i]->getName(), seq);
    }
  }
  return std::move(sites);
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
  auto tm = node->getInfos().model;
  if (!dynamic_pointer_cast<const SubstitutionModelInterface>(tm))
    throw Exception("NonHomogeneousSequenceSimulator::dEvolveInternal : detailed simulation not possible for non-markovian model");

  SimpleMutationProcess process(dynamic_pointer_cast<const SubstitutionModelInterface>(tm));

  MutationPath mp = process.detailedEvolve(node->getFather()->getInfos().state, node->getDistanceToFather() * rate);
  node->getInfos().state = mp.getFinalState();

  // Now append infos in rassr:
  rassr.addNode(static_cast<unsigned int>(node->getId()), mp);

  // Now jump to son nodes:
  for (size_t i = 0; i < node->getNumberOfSons(); i++)
  {
    dEvolveInternal(node->getSon(i), rate, rassr);
  }
}

/******************************************************************************/

void NonHomogeneousSequenceSimulator::outputInternalSequences(bool yn)
{
  outputInternalSequences_ = yn;
  if (outputInternalSequences_)
  {
    vector<SNode*> nodes = tree_.getNodes();
    seqNames_.resize(nodes.size());
    for (size_t i = 0; i < seqNames_.size(); i++)
    {
      if (nodes[i]->isLeaf())
      {
        seqNames_[i] = nodes[i]->getName();
      }
      else
      {
        seqNames_[i] = TextTools::toString( nodes[i]->getId() );
      }
    }
  }
  else
  {
    seqNames_.resize(leaves_.size());
    for (size_t i = 0; i < seqNames_.size(); i++)
    {
      seqNames_[i] = leaves_[i]->getName();
    }
  }
}

/******************************************************************************/
