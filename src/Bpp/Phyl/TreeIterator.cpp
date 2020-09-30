#include "TreeIterator.h"

// imports from bpp-core
#include <Bpp/BppBoolean.h>
using namespace bpp;
using namespace std;

<<<<<<< HEAD
=======
#define LAST_VISITED_SON "last_visited_son"
# define IS_VISITED "is_visited"

/******************************************************************************/

TreeIterator::~TreeIterator()
{
  vector<Node*> nodes = tree_.getNodes();
  for (size_t i = 0; i < nodes.size(); ++i)
  {
    try
    {
      nodes[i]->deleteNodeProperty(LAST_VISITED_SON);
    }
    catch (const std::exception& e)
    {}
    try
    {
      nodes[i]->deleteNodeProperty(IS_VISITED);
    }
    catch (const std::exception& e)
    {}
  }
}

/******************************************************************************/

void TreeIterator::setNodeStatus(Node* node, bool visited)
{
  BppBoolean* visitedProperty = new BppBoolean(visited);
  node->setNodeProperty(IS_VISITED, *visitedProperty);
  delete visitedProperty;
}

>>>>>>> upstream/master

/******************************************************************************/

void TreeIterator::init()
{
<<<<<<< HEAD
    vector <const Node*> nodes = tree_.getNodes();
    for (size_t i=0; i<nodes.size(); ++i)
    {
        nodeToVisited_[nodes[i]->getId()] = false;
        nodeToSonVisited_[nodes[i]->getId()] = false;
    }
=======
  vector<Node*> nodes = tree_.getNodes();
  for (size_t i = 0; i < nodes.size(); ++i)
  {
    setLastVisitedSon(nodes[i], -1);
    setNodeStatus(nodes[i], false);
  }
}

/******************************************************************************/

Node* TreeIterator::begin()
{
  setNodeStatus(curNode_, true);
  if (curNode_->hasFather())
  {
    setLastVisitedSon(curNode_->getFather(), getLastVisitedSon(curNode_->getFather()) + 1);
  }
  return curNode_;
>>>>>>> upstream/master
}

/******************************************************************************/

const Node* TreeIterator::begin()
{
<<<<<<< HEAD
    nodeToVisited_[currNode_->getId()] = true;
    if (currNode_->hasFather())
    {
        nodeToSonVisited_[currNode_->getFatherId()] = true;
        for (size_t i=0; i<currNode_->getFather()->getNumberOfSons(); ++i) {
            if (currNode_ == currNode_->getFather()->getSon(i))
                nodeToLastVisitedSonIndex_[currNode_->getFatherId()] = i;
        }
    }
    return currNode_;
=======
  vector<Node*> nodes = tree_.getNodes();
  for (size_t i = 0; i < nodes.size(); ++i)
  {
    try
    {
      nodes[i]->deleteNodeProperty(LAST_VISITED_SON);
    }
    catch (const std::exception& e)
    {}
    try
    {
      nodes[i]->deleteNodeProperty(IS_VISITED);
    }
    catch (const std::exception& e)
    {}
  }
>>>>>>> upstream/master
}

/******************************************************************************/

TreeIterator& TreeIterator::operator++()    // prefix ++
{
  this->next();
  return *this;
}

/******************************************************************************/

<<<<<<< HEAD
const Node* PostOrderTreeIterator::getLeftMostPredecessor(const Node* startNode)
{
    size_t nextPos = 0;
    if (nodeToSonVisited_[startNode->getId()])
        nextPos = nodeToLastVisitedSonIndex_[startNode->getId()] + 1;
    // if startNode has no predecessors ->  move to him
    if (nextPos >= startNode->getNumberOfSons())
    {
        return startNode;
    }
    const Node* node = startNode->getSon(nextPos);
	while (node->getNumberOfSons() > 0)
    {
		node = node->getSon(0);
	}
	return node;
=======
int TreeIterator::getLastVisitedSon(Node* node)
{
  BppInteger* posProperty = dynamic_cast<BppInteger*>(node->getNodeProperty(LAST_VISITED_SON));
  int posValue = posProperty->getValue();
  return posValue;
}

/******************************************************************************/

void TreeIterator::setLastVisitedSon(Node* node, int visitedSon)
{
  BppInteger* visitedProperty = new BppInteger(visitedSon);
  node->setNodeProperty(LAST_VISITED_SON, *visitedProperty);
  delete visitedProperty;
}

/******************************************************************************/

Node* PostOrderTreeIterator::getLeftMostPredessesor(Node* startNode)
{
  int nextPos = getLastVisitedSon(startNode) + 1;
  // if startNode has no predessesors ->  move to him
  if (nextPos >= static_cast<int>(startNode->getNumberOfSons()))
  {
    return startNode;
  }
  Node* next = startNode->getSon(nextPos);
  Node* node = next;
  while (node->getNumberOfSons() > 0)
  {
    node = node->getSon(0);
  }
  return node;
>>>>>>> upstream/master
}

/******************************************************************************/

<<<<<<< HEAD
const Node* PostOrderTreeIterator::next()
{   // order: (left, right, parent)

    // stop condition: currNode_ is root
    if (!currNode_->hasFather())
    {
        return NULL;
    }

    // by the time you visit currNode_, all the nodes in its subtree were already visited
    size_t numOfBrothers =  currNode_->getFather()->getNumberOfSons(); // get the number of brothers of currNode_
    size_t lastVisitedBrotherPos = nodeToLastVisitedSonIndex_[currNode_->getFatherId()]; // since currNode is already visited, its father must have at least one visited son (which is currNode_)
    // if all the brothers were already visited -> now visit the father
    if (lastVisitedBrotherPos == numOfBrothers-1)
    {
        currNode_ = currNode_->getFather();
    }
    // else -> visit the leftmost predecessor next brother on the list
    else
    {
        size_t nextSiblingPos = lastVisitedBrotherPos + 1;
        currNode_ = currNode_->getFather()->getSon(nextSiblingPos); // the next brother is not visited yet if its has a non-visited leftmost predecessor
        currNode_ = getLeftMostPredecessor(currNode_);
    }

    // update the returned node to be visited and its father's last visited son, if needed
    nodeToVisited_[currNode_->getId()] = true;
    if (currNode_->hasFather())
    {
        if (!nodeToSonVisited_[currNode_->getFatherId()])
        {
            nodeToSonVisited_[currNode_->getFatherId()] = true;
            nodeToLastVisitedSonIndex_[currNode_->getFatherId()] = 0;
        }
        else
        {
            nodeToLastVisitedSonIndex_[currNode_->getFatherId()] = nodeToLastVisitedSonIndex_[currNode_->getFatherId()] + 1;
        }
    }
    return currNode_;
=======
Node* PostOrderTreeIterator::next()
{   // order: (left, right, parent)
    // stop condition: curNode_ is root
  if (!curNode_->hasFather())
  {
    clearProperties();
    return NULL;
  }
  // by the time you visit curNode_, all the nodes in its subtree were already visited
  int numOfBrothers =  static_cast<int>(curNode_->getFather()->getNumberOfSons()); // get the number of brothers of curNode_
  int lastVisitedBrotherPos = getLastVisitedSon(curNode_->getFather());
  // if all the brothers were already visited -> now visit the father
  if (lastVisitedBrotherPos == numOfBrothers - 1)
  {
    curNode_ = curNode_->getFather();
  }
  // else -> visit the leaftmost presdessesor next brother on the list
  else
  {
    int nextSiblingPos = getLastVisitedSon(curNode_->getFather()) + 1;
    curNode_ = curNode_->getFather()->getSon(nextSiblingPos);
    curNode_ = getLeftMostPredessesor(curNode_);
  }
  // update the returned node to be visited and its father's last visited son, if needed
  if (curNode_->hasFather())
  {
    setLastVisitedSon(curNode_->getFather(), getLastVisitedSon(curNode_->getFather()) + 1);
  }
  setNodeStatus(curNode_, true);
  return curNode_;
>>>>>>> upstream/master
}

/******************************************************************************/

<<<<<<< HEAD
const Node* PreOrderTreeIterator::next()
{   // order: (parent, left, right)

    
    vector<int> leafIds = tree_.getLeavesId();
    bool hasVisitedSons = nodeToSonVisited_[currNode_->getId()];
    size_t numOfSons = currNode_->getNumberOfSons();

    // stop condition: the node is the rightmost leaf of the tree
    if (currNode_->getId() == leafIds[leafIds.size()-1])
    {
        return NULL;
    }

    // by the time you visit currNode_, all its ancestors and left brothers (if exist) were already visited
    if (!hasVisitedSons && numOfSons > 0)
    {
        currNode_ = currNode_->getSon(0); // this somehow modifies the value of nodeToSonVisited_[currNode_->getFatherId()] to true... how?
    }
    // as long as there are still sons to visit -> visit the leftmost unvisited son
    else if (hasVisitedSons && nodeToLastVisitedSonIndex_[currNode_->getId()] < numOfSons-1)
    {
        // if this point is reached, the current node already has visited sons, or it has no sons at all and we are done
        currNode_ = currNode_->getSon(nodeToLastVisitedSonIndex_[currNode_->getId()]+1);
    // else -> traverse to the leftmost brother which is right to currNode_ (also occurs when currNode has no sons
    }
    else
    {
       currNode_ = currNode_->getFather();
       size_t lastVisitedSonPos = nodeToLastVisitedSonIndex_[currNode_->getId()]; // the father of the original currNode_ must have at least one visited child which is the original currNode_
       numOfSons = currNode_->getNumberOfSons();
       while (lastVisitedSonPos == numOfSons-1)
       {
             currNode_ = currNode_->getFather(); // the father node must have visited sons as the currNode is a visited son of his
            lastVisitedSonPos = nodeToLastVisitedSonIndex_[currNode_->getId()];
            numOfSons = currNode_->getNumberOfSons();
       }
        currNode_ = currNode_->getSon(lastVisitedSonPos+1);    // need to remove +1??
    }
    // update the returned node to be visited and its father's last visited son, if needed
    if (currNode_->hasFather())
    {
        if (nodeToSonVisited_[currNode_->getFatherId()])
            nodeToLastVisitedSonIndex_[currNode_->getFatherId()] = nodeToLastVisitedSonIndex_[currNode_->getFatherId()] + 1;
        else
        {
            nodeToSonVisited_[currNode_->getFatherId()] = true;
            nodeToLastVisitedSonIndex_[currNode_->getFatherId()] = 0;
        }
    }
    nodeToVisited_[currNode_->getId()] = true;
    return currNode_;
=======
Node* PreOrderTreeIterator::next()
{   // order: (parent, left, right)
    // stop condition: the node is the righmost leaf of the tree
  vector<int> leafIds = tree_.getLeavesId();
  if (curNode_->getId() == leafIds[leafIds.size() - 1])
  {
    clearProperties();
    return NULL;
  }
  // by the time you visit curNode_, all its ancestors and left brothers (if exist) were already visited
  int numOfSons =  static_cast<int>(curNode_->getNumberOfSons());
  int lastVisitedSonPos = getLastVisitedSon(curNode_);
  // as long as there are still sons to visit -> visit the lesftmost unvisited son
  if (lastVisitedSonPos < numOfSons - 1)
  {
    curNode_ = curNode_->getSon(lastVisitedSonPos + 1);
    // else -> traverse to the leftmost brother which is right to curNode_
  }
  else
  {
    curNode_ = curNode_->getFather();
    lastVisitedSonPos = getLastVisitedSon(curNode_);
    numOfSons = static_cast<int>(curNode_->getNumberOfSons());
    while (lastVisitedSonPos == numOfSons - 1)
    {
      curNode_ = curNode_->getFather();
      lastVisitedSonPos = getLastVisitedSon(curNode_);
      numOfSons = static_cast<int>(curNode_->getNumberOfSons());
    }
    curNode_ = curNode_->getSon(lastVisitedSonPos + 1);
  }
  // update the returned node to be visited and its father's last visited son, if needed
  if (curNode_->hasFather())
  {
    setLastVisitedSon(curNode_->getFather(), getLastVisitedSon(curNode_->getFather()) + 1);
  }
  setNodeStatus(curNode_, true);
  return curNode_;
>>>>>>> upstream/master
}

/******************************************************************************/

const Node* InOrderTreeIterator::doStep(const Node* node)
{
<<<<<<< HEAD
    // if the node has unvisited left sons -> visit the leftmost unvisited son
    size_t lastVisitedSon;
    if (nodeToSonVisited_[node->getId()])
        lastVisitedSon = nodeToLastVisitedSonIndex_[node->getId()];
    size_t numOfSons = node->getNumberOfSons();
    bool is_visited =  nodeToVisited_[node->getId()];
   
    // if the node has unvisited left sons - go to the leftmost unvisited son
    if (!nodeToSonVisited_[node->getId()] && numOfSons > 0)
    {
        return node->getSon(0);
 
    }
    else if (nodeToSonVisited_[node->getId()] && lastVisitedSon < static_cast<size_t>(floor(numOfSons/2)-1))
    {
        return node->getSon(lastVisitedSon+1);

    }
    
    // if the node last visited its last left son or is a not yet visited leaf / parent of a single node -> go to it
    else if (numOfSons > 1 && lastVisitedSon == static_cast<size_t>(floor(numOfSons/2)-1) && !is_visited)
    {
        return node;
    }

    // if the node has unvisited right sons -> go to the leftmost unvisited right son
    else if (nodeToSonVisited_[node->getId()] && lastVisitedSon < numOfSons-1)
    {
        return node->getSon(lastVisitedSon+1);
    }

    // else - the entire subtree of the node has been scanned - move to its father
    else 
    {
        return node->getFather();
    }
=======
  // if the node has unvisited left sons -> visit the leftmost unvisited son
  int lastVisitedSon = getLastVisitedSon(node);
  int numOfSons = static_cast<int>(node->getNumberOfSons());
  int nextSon;
  bool is_visited =  dynamic_cast<BppBoolean*>(node->getNodeProperty(IS_VISITED))->getValue();
  // if the node has unvisited left sons - go to the leftmost unvisited son
  if (lastVisitedSon < floor(numOfSons / 2) - 1)
  {
    nextSon = lastVisitedSon + 1;
    return node->getSon(nextSon);
  }
  // if the node last visited its last left son or is a leaf / parent of a single node -> go to it
  else if (node->getNumberOfSons() > 1 && lastVisitedSon == floor(numOfSons / 2) - 1 && !is_visited)
  {
    return node;
  }
  // if the node has unvisited right sons -> go to the leftmost unviited right son
  else if (lastVisitedSon < numOfSons - 1)
  {
    nextSon = lastVisitedSon + 1;
    return node->getSon(nextSon);
  }
  // else - the entire subtree of the node has been scanned - move to its father
  else
  {
    return node->getFather();
  }
>>>>>>> upstream/master
}


/******************************************************************************/

const Node* InOrderTreeIterator::next()
{   // order: (left (0,..,n/2-1), parent, right (n/2,....,n))
<<<<<<< HEAD

    vector<int> leafIds = tree_.getLeavesId();

    // stop condition: the node is the rightmost leaf of the tree
    if (currNode_->getId() == leafIds[leafIds.size()-1])
    {
        return NULL;
    }

    // while curNode still has unvisited left sons -> do another step
    while (nodeToVisited_[currNode_->getId()]
     || (!nodeToSonVisited_[currNode_->getId()] && !currNode_->isLeaf()) || (nodeToSonVisited_[currNode_->getId()] && nodeToLastVisitedSonIndex_[currNode_->getId()] < static_cast<size_t>(floor(currNode_->getNumberOfSons()/2)-1) && !currNode_->isLeaf()))
    {
        currNode_ = doStep(currNode_);
    }

    // update the returned node to be visited and its father's last visited son, if needed
    if (currNode_->hasFather())
    {
        if (nodeToSonVisited_[currNode_->getFatherId()])
            nodeToLastVisitedSonIndex_[currNode_->getFatherId()] = nodeToLastVisitedSonIndex_[currNode_->getFatherId()] + 1;
        else
        {
            nodeToSonVisited_[currNode_->getFatherId()] = true;
            nodeToLastVisitedSonIndex_[currNode_->getFatherId()] = 0;
        }
    }
    nodeToVisited_[currNode_->getId()] = true;
    return currNode_;
=======
    // stop condition: the node is the righmost leaf of the tree
  vector<int> leafIds = tree_.getLeavesId();
  if (curNode_->getId() == leafIds[leafIds.size() - 1])
  {
    clearProperties();
    return NULL;
  }
  bool is_visited = dynamic_cast<BppBoolean*>(curNode_->getNodeProperty(IS_VISITED))->getValue();
  while (is_visited || (getLastVisitedSon(curNode_) < floor(curNode_->getNumberOfSons() / 2) - 1 && !curNode_->isLeaf()))  // while curNode still has unvisited left sons -> do another step
  {
    curNode_ = doStep(curNode_);
    is_visited = dynamic_cast<BppBoolean*>(curNode_->getNodeProperty(IS_VISITED))->getValue();
  }
  // update the returned node to be visited and its father's last visited son, if needed
  if (curNode_->hasFather())
  {
    setLastVisitedSon(curNode_->getFather(), getLastVisitedSon(curNode_->getFather()) + 1);
  }
  setNodeStatus(curNode_, true);
  return curNode_;
>>>>>>> upstream/master
}
