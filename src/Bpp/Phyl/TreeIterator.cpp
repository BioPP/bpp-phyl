#include "TreeIterator.h"

// imports from bpp-core
#include <Bpp/BppBoolean.h>

#include <iostream>
#include <fstream>

using namespace bpp;

using namespace std;

#define LAST_VISITED_SON "last_visited_son"
# define IS_VISITED "is_visited"

/******************************************************************************/

TreeIterator::~TreeIterator()
{
    vector <Node*> nodes = tree_.getNodes();
    for (size_t i=0; i<nodes.size(); ++i) 
    {
        try
        {
            nodes[i]->deleteNodeProperty(LAST_VISITED_SON);
        }
        catch(const std::exception& e) {}
        try
        {
            nodes[i]->deleteNodeProperty(IS_VISITED);
        }
        catch(const std::exception& e) {}
        
    }
}

/******************************************************************************/

void TreeIterator::setNodeStatus(Node* node, bool visited)
{
    BppBoolean* visitedProperty = new BppBoolean(visited);
    node->setNodeProperty(IS_VISITED, *visitedProperty);
    delete visitedProperty; 
}


/******************************************************************************/

void TreeIterator::init()
{
    vector <Node*> nodes = tree_.getNodes();
    for (size_t i=0; i<nodes.size(); ++i) 
    {
        setLastVisitedSon(nodes[i],-1);
        setNodeStatus(nodes[i],false);
    }
}

/******************************************************************************/

Node* TreeIterator::begin()
{
    setNodeStatus(curNode_, true);
    if (curNode_->hasFather())
    {
        setLastVisitedSon(curNode_->getFather(), getLastVisitedSon(curNode_->getFather())+1);
    }
    return curNode_;
}

/******************************************************************************/

void TreeIterator::clearProperties()
{
    vector <Node*> nodes = tree_.getNodes();
    for (size_t i=0; i<nodes.size(); ++i) 
    {
        try
        {
            nodes[i]->deleteNodeProperty(LAST_VISITED_SON);
        }
        catch(const std::exception& e) {}
        try
        {
            nodes[i]->deleteNodeProperty(IS_VISITED);
        }
        catch(const std::exception& e) {}
        
    }
}

/******************************************************************************/

TreeIterator& TreeIterator::operator++()    // prefix ++
{
    this->next();
    return *this;
}

/******************************************************************************/

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
    int nextPos = getLastVisitedSon(startNode)+1;
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
}

/******************************************************************************/

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
    if (lastVisitedBrotherPos == numOfBrothers-1)
    {
        curNode_ = curNode_->getFather();
    }
    // else -> visit the leaftmost presdessesor next brother on the list
    else 
    {
        int nextSiblingPos = getLastVisitedSon(curNode_->getFather())+1;
        curNode_ = curNode_->getFather()->getSon(nextSiblingPos);
        curNode_ = getLeftMostPredessesor(curNode_);
    }                                        
    // update the returned node to be visited and its father's last visited son, if needed
    if (curNode_->hasFather())
    {
        setLastVisitedSon(curNode_->getFather(), getLastVisitedSon(curNode_->getFather())+1);
    }
    setNodeStatus(curNode_, true);
    return curNode_;
}

/******************************************************************************/

Node* PreOrderTreeIterator::next() 
{   // order: (parent, left, right)
    // stop condition: the node is the righmost leaf of the tree 
    vector<int> leafIds = tree_.getLeavesId();
    if (curNode_->getId() == leafIds[leafIds.size()-1])
    {
        clearProperties();
        return NULL;
    }
    // by the time you visit curNode_, all its ancestors and left brothers (if exist) were already visited
    int numOfSons =  static_cast<int>(curNode_->getNumberOfSons());
    int lastVisitedSonPos = getLastVisitedSon(curNode_);
    // as long as there are still sons to visit -> visit the lesftmost unvisited son
    if (lastVisitedSonPos < numOfSons-1)
    {
        curNode_ = curNode_->getSon(lastVisitedSonPos+1);       
    // else -> traverse to the leftmost brother which is right to curNode_
    } 
    else 
    {
       curNode_ = curNode_->getFather();
       lastVisitedSonPos = getLastVisitedSon(curNode_);
       numOfSons = static_cast<int>(curNode_->getNumberOfSons());
       while (lastVisitedSonPos == numOfSons-1)
       {
            curNode_ = curNode_->getFather();
            lastVisitedSonPos = getLastVisitedSon(curNode_);
            numOfSons = static_cast<int>(curNode_->getNumberOfSons());
       }
        curNode_ = curNode_->getSon(lastVisitedSonPos+1);    
    }
    // update the returned node to be visited and its father's last visited son, if needed
    if (curNode_->hasFather())
    {
        setLastVisitedSon(curNode_->getFather(), getLastVisitedSon(curNode_->getFather())+1);
    }
    setNodeStatus(curNode_, true);
    return curNode_;
}

/******************************************************************************/

Node* InOrderTreeIterator::doStep(Node* node)
{
    // if the node has unvisited left sons -> visit the leftmost unvisited son
    int lastVisitedSon = getLastVisitedSon(node);
    int numOfSons = static_cast<int>(node->getNumberOfSons());
    int nextSon;
    bool is_visited =  dynamic_cast<BppBoolean*>(node->getNodeProperty(IS_VISITED))->getValue();
    // if the node has unvisited left sons - go to the leftmost unvisited son
    if (lastVisitedSon < floor(numOfSons/2)-1)
    {
        nextSon = lastVisitedSon+1;
        return node->getSon(nextSon);
    }
    // if the node last visited its last left son or is a leaf / parent of a single node -> go to it
    else if (node->getNumberOfSons() > 1 && lastVisitedSon == floor(numOfSons/2)-1 && !is_visited)
    {
        return node;
    }
    // if the node has unvisited right sons -> go to the leftmost unviited right son
    else if (lastVisitedSon < numOfSons-1)
    {
        nextSon = lastVisitedSon+1;
        return node->getSon(nextSon);
    }
    // else - the entire subtree of the node has been scanned - move to its father
    else 
    {
        return node->getFather();
    }
}

/******************************************************************************/

Node* InOrderTreeIterator::next()
{   // order: (left (0,..,n/2-1), parent, right (n/2,....,n))
    // stop condition: the node is the righmost leaf of the tree 
    vector<int> leafIds = tree_.getLeavesId();
    if (curNode_->getId() == leafIds[leafIds.size()-1])
    {
        clearProperties();
        return NULL;
    }
    bool is_visited = dynamic_cast<BppBoolean*>(curNode_->getNodeProperty(IS_VISITED))->getValue();
    while (is_visited || (getLastVisitedSon(curNode_) < floor(curNode_->getNumberOfSons()/2)-1 && !curNode_->isLeaf()))  // while curNode still has unvisited left sons -> do another step
    {
        curNode_ = doStep(curNode_);
        is_visited = dynamic_cast<BppBoolean*>(curNode_->getNodeProperty(IS_VISITED))->getValue();
    }
    // update the returned node to be visited and its father's last visited son, if needed
    if (curNode_->hasFather())
    {
        setLastVisitedSon(curNode_->getFather(), getLastVisitedSon(curNode_->getFather())+1);
    }
    setNodeStatus(curNode_, true);
    return curNode_;
}
