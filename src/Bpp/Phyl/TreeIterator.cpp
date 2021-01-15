#include "TreeIterator.h"

// imports from bpp-core
using namespace bpp;
using namespace std;


/******************************************************************************/

void TreeIterator::init()
{
    vector <const Node*> nodes = tree_.getNodes();
    for (size_t i=0; i<nodes.size(); ++i)
    {
        nodeToVisited_[nodes[i]->getId()] = false;
        nodeToSonVisited_[nodes[i]->getId()] = false;
    }
}

/******************************************************************************/

const Node* TreeIterator::begin()
{
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
}

/******************************************************************************/

TreeIterator& TreeIterator::operator++()    // prefix ++
{
    this->next();
    return *this;
}

/******************************************************************************/

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
}

/******************************************************************************/

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
}

/******************************************************************************/

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
}

/******************************************************************************/

const Node* InOrderTreeIterator::doStep(const Node* node)
{
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
}


/******************************************************************************/

const Node* InOrderTreeIterator::next()
{   // order: (left (0,..,n/2-1), parent, right (n/2,....,n))

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
}
