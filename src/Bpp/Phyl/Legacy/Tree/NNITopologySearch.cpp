// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Text/TextTools.h>

#include "../Likelihood/NNIHomogeneousTreeLikelihood.h"
#include "NNITopologySearch.h"

using namespace bpp;

// From the STL:
#include <cmath>

using namespace std;

const string NNITopologySearch::FAST   = "Fast";
const string NNITopologySearch::BETTER = "Better";
const string NNITopologySearch::PHYML  = "PhyML";

void NNITopologySearch::notifyAllPerformed(const TopologyChangeEvent& event)
{
  searchableTree_->topologyChangePerformed(event);
  for (size_t i = 0; i < topoListeners_.size(); i++)
  {
    topoListeners_[i]->topologyChangePerformed(event);
  }
}

void NNITopologySearch::notifyAllTested(const TopologyChangeEvent& event)
{
  searchableTree_->topologyChangeTested(event);
  for (size_t i = 0; i < topoListeners_.size(); i++)
  {
    topoListeners_[i]->topologyChangeTested(event);
  }
}

void NNITopologySearch::notifyAllSuccessful(const TopologyChangeEvent& event)
{
  searchableTree_->topologyChangeSuccessful(event);
  for (size_t i = 0; i < topoListeners_.size(); i++)
  {
    topoListeners_[i]->topologyChangeSuccessful(event);
  }
}

void NNITopologySearch::search()
{
  if (algorithm_ == FAST)
    searchFast();
  else if (algorithm_ == BETTER)
    searchBetter();
  else if (algorithm_ == PHYML)
    searchPhyML();
  else
    throw Exception("Unknown NNI algorithm: " + algorithm_ + ".\n");
}

void NNITopologySearch::searchFast()
{
  bool test = true;
  do
  {
    TreeTemplate<Node> tree(searchableTree_->topology());
    vector<Node*> nodes = tree.getNodes();

    vector<Node*> nodesSub = nodes;
    for (size_t i = nodesSub.size(); i > 0; i--)
    {
      // !!! must not reach i==0 because of size_t
      if (!(nodesSub[i - 1]->hasFather()))
        nodesSub.erase(nodesSub.begin() + static_cast<ptrdiff_t>(i - 1)); // Remove root node.
      else if (!(nodesSub[i - 1]->getFather()->hasFather()))
        nodesSub.erase(nodesSub.begin() + static_cast<ptrdiff_t>(i - 1)); // Remove son of root node.
    }

    // Test all NNIs:
    test = false;
    for (size_t i = 0; !test && i < nodesSub.size(); i++)
    {
      Node* node = nodesSub[i];
      double diff = searchableTree_->testNNI(node->getId());
      if (verbose_ >= 3)
      {
        ApplicationTools::displayResult("   Testing node " + TextTools::toString(node->getId())
            + " at " + TextTools::toString(node->getFather()->getId()),
            TextTools::toString(diff));
      }

      if (diff < 0.)
      { // Good NNI found...
        if (verbose_ >= 2)
        {
          ApplicationTools::displayResult("   Swapping node " + TextTools::toString(node->getId())
              + " at " + TextTools::toString(node->getFather()->getId()),
              TextTools::toString(diff));
        }
        searchableTree_->doNNI(node->getId());
        // Notify:
        notifyAllPerformed(TopologyChangeEvent());
        test = true;

        if (verbose_ >= 1)
          ApplicationTools::displayResult("   Current value", TextTools::toString(searchableTree_->getTopologyValue(), 10));
      }
    }
  }
  while (test);
}

void NNITopologySearch::searchBetter()
{
  bool test = true;
  do
  {
    TreeTemplate<Node> tree(searchableTree_->topology());
    vector<Node*> nodes = tree.getNodes();

    if (verbose_ >= 3)
      ApplicationTools::displayTask("Test all possible NNIs...");

    vector<Node*> nodesSub = nodes;
    for (size_t i = nodesSub.size(); i > 0; i--)
    { // !!! must not reach i==0 because of size_t
      if (!(nodesSub[i - 1]->hasFather()))
        nodesSub.erase(nodesSub.begin() + static_cast<ptrdiff_t>(i - 1)); // Remove root node.
      else if (!(nodesSub[i - 1]->getFather()->hasFather()))
        nodesSub.erase(nodesSub.begin() + static_cast<ptrdiff_t>(i - 1)); // Remove son of root node.
    }

    // Test all NNIs:
    vector<Node*> improving;
    vector<double> improvement;
    if (verbose_ >= 2 && ApplicationTools::message)
      ApplicationTools::message->endLine();
    for (size_t i = 0; i < nodesSub.size(); i++)
    {
      Node* node = nodesSub[i];
      double diff = searchableTree_->testNNI(node->getId());
      if (verbose_ >= 3)
      {
        ApplicationTools::displayResult("   Testing node " + TextTools::toString(node->getId())
            + " at " + TextTools::toString(node->getFather()->getId()),
            TextTools::toString(diff));
      }

      if (diff < 0.)
      {
        improving.push_back(node);
        improvement.push_back(diff);
      }
    }
    if (verbose_ >= 3)
      ApplicationTools::displayTaskDone();
    test = improving.size() > 0;
    if (test)
    {
      size_t nodeMin = VectorTools::whichMin(improvement);
      Node* node = improving[nodeMin];
      if (verbose_ >= 2)
        ApplicationTools::displayResult("   Swapping node " + TextTools::toString(node->getId())
            + " at " + TextTools::toString(node->getFather()->getId()),
            TextTools::toString(improvement[nodeMin]));
      searchableTree_->doNNI(node->getId());

      // Notify:
      notifyAllPerformed(TopologyChangeEvent());

      if (verbose_ >= 1)
        ApplicationTools::displayResult("   Current value", TextTools::toString(searchableTree_->getTopologyValue(), 10));
    }
  }
  while (test);
}

void NNITopologySearch::searchPhyML()
{
  bool test = true;
  do
  {
    if (verbose_ >= 3)
      ApplicationTools::displayTask("Test all possible NNIs...");
    TreeTemplate<Node> tree(searchableTree_->topology());
    vector<Node*> nodes = tree.getNodes();
    vector<Node*> nodesSub = nodes;
    for (size_t i = nodesSub.size(); i > 0; i--)
    {
      // !!! must not reach i==0 because of size_t
      if (!(nodesSub[i - 1]->hasFather()))
        nodesSub.erase(nodesSub.begin() + static_cast<ptrdiff_t>(i - 1)); // Remove root node.
      else if (!(nodesSub[i - 1]->getFather()->hasFather()))
        nodesSub.erase(nodesSub.begin() + static_cast<ptrdiff_t>(i - 1)); // Remove son of root node.
    }

    // Test all NNIs:
    vector<int> improving;
    vector<Node*> improvingNodes;
    vector<double> improvement;
    if (verbose_ >= 2 && ApplicationTools::message)
      ApplicationTools::message->endLine();
    for (size_t i = 0; i < nodesSub.size(); i++)
    {
      Node* node = nodesSub[i];
      double diff = searchableTree_->testNNI(node->getId());
      if (verbose_ >= 3)
      {
        ApplicationTools::displayResult("   Testing node " + TextTools::toString(node->getId())
            + " at " + TextTools::toString(node->getFather()->getId()),
            TextTools::toString(diff));
      }

      if (diff < 0.)
      {
        bool ok = true;
        // Must test for incompatible NNIs...
        for (size_t j = improving.size(); j > 0; j--)
        {
          if (improvingNodes[j - 1]->getFather()->getId() == node->getFather()->getId()
              || improvingNodes[j - 1]->getFather()->getFather()->getId() == node->getFather()->getFather()->getId()
              || improvingNodes[j - 1]->getId() == node->getFather()->getId()                           || improvingNodes[j - 1]->getFather()->getId() == node->getId()
              || improvingNodes[j - 1]->getId() == node->getFather()->getFather()->getId()              || improvingNodes[j - 1]->getFather()->getFather()->getId() == node->getId()
              || improvingNodes[j - 1]->getFather()->getId() == node->getFather()->getFather()->getId() || improvingNodes[j - 1]->getFather()->getFather()->getId() == node->getFather()->getId())
          {
            // These are incompatible NNIs. We only keep the best:
            if (diff < improvement[j - 1])
            { // Erase previous node
              improvingNodes.erase(improvingNodes.begin() + static_cast<ptrdiff_t>(j - 1));
              improving.erase(improving.begin() + static_cast<ptrdiff_t>(j - 1));
              improvement.erase(improvement.begin() + static_cast<ptrdiff_t>(j - 1));
            } // Otherwise forget about this NNI.
            else
            {
              ok = false;
            }
          }
        }
        if (ok)
        { // We add this NNI to the list,
          // by decreasing improvement:
          size_t pos = improvement.size();
          for (size_t j = 0; j < improvement.size(); j++)
          {
            if (diff < improvement[j])
            {
              pos = j; break;
            }
          }
          if (pos < improvement.size())
          {
            improvingNodes.insert(improvingNodes.begin() + static_cast<ptrdiff_t>(pos), node);
            improving.insert(improving.begin() + static_cast<ptrdiff_t>(pos), node->getId());
            improvement.insert(improvement.begin() + static_cast<ptrdiff_t>(pos), diff);
          }
          else
          {
            improvingNodes.insert(improvingNodes.end(), node);
            improving.insert(improving.end(), node->getId());
            improvement.insert(improvement.end(), diff);
          }
        }
      }
    }
    // This array is no more useful.
    // Moreover, if a backward movement is performed,
    // the underlying node will not exist anymore...
    improvingNodes.clear();
    if (verbose_ >= 3)
      ApplicationTools::displayTaskDone();
    test = improving.size() > 0;
    if (test)
    {
      double currentValue = searchableTree_->getTopologyValue();
      bool test2 = true;
      // Make a backup copy:
      NNISearchable* backup = dynamic_cast<NNISearchable*>(searchableTree_->clone());
      do
      {
        if (verbose_ >= 1)
          ApplicationTools::displayMessage("Trying to perform " + TextTools::toString(improving.size()) + " NNI(s).");
        for (size_t i = 0; i < improving.size(); i++)
        {
          int nodeId = improving[i];
          if (verbose_ >= 2)
          {
            ApplicationTools::displayResult(string("   Swapping node ") + TextTools::toString(nodeId)
                + string(" at ") + TextTools::toString(searchableTree_->topology().getFatherId(nodeId)),
                TextTools::toString(improvement[i]));
          }
          searchableTree_->doNNI(improving[i]);
        }

        // Notify:
        notifyAllTested(TopologyChangeEvent());
        if (verbose_ >= 1)
          ApplicationTools::displayResult("   Current value", TextTools::toString(searchableTree_->getTopologyValue(), 10));
        if (searchableTree_->getTopologyValue() >= currentValue)
        {
          // No improvement!
          // Restore backup:
          searchableTree_.reset(backup->clone());
          if (verbose_ >= 1)
          {
            ApplicationTools::displayResult("Score >= current score! Moving backward", TextTools::toString(searchableTree_->getTopologyValue()));
          }
          // And try doing half of the movements:
          if (improving.size() == 1)
          {
            // Problem! This should have worked!!!
            throw Exception("NNITopologySearch::searchPhyML. Error, no improving NNI!\n This may be due to a change in parameters between testNNI and doNNI. Check your code!");
          }
          size_t n = (size_t)ceil((double)improving.size() / 2.);
          improving.erase(improving.begin() + static_cast<ptrdiff_t>(n), improving.end());
          improvement.erase(improvement.begin() + static_cast<ptrdiff_t>(n), improvement.end());
        }
        else
        {
          test2 = false;
        }
      }
      while (test2);
      delete backup;
      // Notify:
      notifyAllSuccessful(TopologyChangeEvent());
    }
  }
  while (test);
}
