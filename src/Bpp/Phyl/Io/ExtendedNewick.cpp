// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/BppString.h>
#include <Bpp/Numeric/Number.h>
#include <Bpp/Text/NestedStringTokenizer.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Text/TextTools.h>

#include "../Tree/PhyloBranch.h"
#include "../Tree/PhyloNode.h"
#include "../Tree/PhyloDAG.h"
#include "ExtendedNewick.h"

using namespace bpp;

// From the STL:
#include <iostream>
#include <fstream>

using namespace std;

/******************************************************************************/

const string ExtendedNewick::getFormatName() const { return "ExtendedNewick"; }

/******************************************************************************/

const string ExtendedNewick::getFormatDescription() const
{
  return string("Extended Newick Format. ");
}

/**********************************************************/
/*  INPUT */
/**********************************************************/


unique_ptr<PhyloDAG> ExtendedNewick::readPhyloDAG(istream& in) const
{
  // Checking the existence of specified file
  if (!in)
  {
    throw IOException ("ExtendedNewick::readPhyloDAG: failed to read from stream");
  }

  // We concatenate all line in file till we reach the ending semi colon:
  string temp, description; // Initialization
  // Main loop : for all file lines
  while (getline(in, temp, '\n'))
  {
    string::size_type index = temp.find(";");
    if (index != string::npos)
    {
      description += temp.substr(0, index + 1);
      break;
    }
    else
      description += temp;
  }

  if (allowComments_)
    description = TextTools::removeSubstrings(description, '[', ']');
  if (TextTools::isEmpty(description))
    throw IOException("ExtendedNewick::read: no dag was found!");
  return parenthesisToPhyloDAG(description, verbose_);
}


/******************************************************************************/

void ExtendedNewick::readPhyloDAGs(istream& in, vector<unique_ptr<PhyloDAG>>& dags) const
{
  // Checking the existence of specified file
  if (!in)
  {
    throw IOException ("ExtendedNewick::readPhyloDAGs(vector): failed to read from stream");
  }

  // Main loop : for all file lines
  string temp, description; // Initialization
  string::size_type index;
  // We concatenate all line in file till we reach the ending semi colon:
  while (getline(in, temp, '\n'))
  {
    index = temp.find(";");
    if (index != string::npos)
    {
      description += temp.substr(0, index + 1);
      if (allowComments_)
        description = TextTools::removeSubstrings(description, '[', ']');
      dags.push_back(parenthesisToPhyloDAG(description, verbose_));
      description = temp.substr(index + 1);
    }
    else
      description += temp;
  }
  // In case the file is empty, the method will not add any neww dag to the vector.
}

/***************************************/

IODAG::Element ExtendedNewick::getElement(const string& elt) const
{
  IODAG::Element element;
  element.length    = ""; // default
  element.annotation = ""; // default
  element.isLeaf    = false; // default

  size_t colonIndex;
  bool hasColon = false;
  for (colonIndex = elt.size(); colonIndex > 0 && elt[colonIndex] != ')'; colonIndex--)
  {
    if (elt[colonIndex] == ':')
    {
      hasColon = true;
      break;
    }
  }
  try
  {
    string elt2;
    if (hasColon)
    {
      // this is an element with length:
      elt2 = elt.substr(0, colonIndex);
      element.length = TextTools::removeSurroundingWhiteSpaces(elt.substr(colonIndex + 1));
    }
    else
    {
      // this is an element without length;
      elt2 = elt;
    }

    string::size_type lastP = elt2.rfind(')');
    string::size_type firstP = elt2.find('(');
    if (firstP == string::npos)
    {
      // This is a leaf:
      element.content = elt2;
      element.annotation = elt2;
      element.isLeaf = true;
    }
    else
    {
      // This is a node:
      if (lastP < firstP)
        throw IOException("ExtendedNewick::getElement(). Invalid format: bad closing parenthesis in " + elt2);
      element.content = TextTools::removeSurroundingWhiteSpaces(elt2.substr(firstP + 1, lastP - firstP - 1));
      string annot = TextTools::removeSurroundingWhiteSpaces(elt2.substr(lastP + 1));
      if (!TextTools::isEmpty(annot))
      {
        element.annotation = annot;
      }
    }
  }
  catch (exception& e)
  {
    throw IOException("Bad dag description: " + elt);
  }
  return element;
}

/************************************************************/

shared_ptr<PhyloNode>  ExtendedNewick::parenthesisToNode(PhyloDAG& dag, std::shared_ptr<PhyloNode>  father, const std::string& description, unsigned int& nodeCounter, unsigned int& branchCounter, std::map<std::string, std::shared_ptr<PhyloNode> >& mapEvent, bool withId, bool verbose) const
{
//  cerr << "NODE: " << description << endl;
  IODAG::Element elt = getElement(description);

  // Is the node a connecting one?

  string annot = elt.annotation;
  size_t poshash = annot.find("#");

  shared_ptr<PhyloNode> node;

  
  // Check Event:
  if (poshash != string::npos)
  {
    string evId = annot.substr(poshash+1);
    string label = annot.substr(0, poshash);

    if (mapEvent.find(evId)!=mapEvent.end())
      node=mapEvent[evId];
    else
    {
      node = std::make_shared<PhyloNode>(label);
      if (evId[0]=='H')
      {
        auto event = std::make_shared<NodeEvent>(NodeEvent::hybridizationEvent);
        node->setProperty("event", *event);
      }
      mapEvent[evId]=node;
      dag.createNode(node);
    }
  }
  else
  {
    node = std::make_shared<PhyloNode>(annot);
    dag.createNode(node);
  }

  shared_ptr<PhyloBranch> branch(father ? new PhyloBranch() : 0);
  
  if (father)
  {
    dag.link(father, node, branch);
    
    if (!TextTools::isEmpty(elt.length))
      branch->setLength(TextTools::toDouble(elt.length));
  }

  
  if (annot.size()!=0)
  {
    if (withId)
    {
      auto id = static_cast<PhyloDAG::NodeIndex>(TextTools::toInt(elt.annotation));
      dag.setNodeIndex(node, id);
    }
  }

  NestedStringTokenizer nt(elt.content, "(", ")", ",");
  vector<string> elements;
  while (nt.hasMoreToken())
  {
    elements.push_back(nt.nextToken());
  }

  if (elt.isLeaf)
  {
    // This is a leaf:
    string name = TextTools::removeSurroundingWhiteSpaces(elements[0]);
    if (withId)
    {
      StringTokenizer st(name, "_", true, true);
      ostringstream realName;
      for (size_t i = 0; i < st.numberOfRemainingTokens() - 1; ++i)
      {
        if (i != 0)
          realName << "_";

        realName << st.getToken(i);
      }
      node->setName(realName.str());
      dag.setNodeIndex(node, static_cast<PhyloDAG::NodeIndex>(
            TextTools::toInt(st.getToken(st.numberOfRemainingTokens() - 1))));
    }
    else
      node->setName(name);
  }
  else
  {
    // This is a node:
    for (size_t i = 0; i < elements.size(); i++)
      parenthesisToNode(dag, node, elements[i], nodeCounter, branchCounter, mapEvent, withId, verbose);

  }

  if (!withId)
  {
    if (!dag.hasNodeIndex(node))
    {
      dag.setNodeIndex(node, nodeCounter);
      nodeCounter++;
    }
    
    if (branch){
      dag.setEdgeIndex(branch, branchCounter);
      branchCounter++;
    }
  }

  if (verbose)
    ApplicationTools::displayUnlimitedGauge(nodeCounter);
  return node;
}

/******************************************************************************/

unique_ptr<PhyloDAG> ExtendedNewick::parenthesisToPhyloDAG(const string& description, bool withId, bool verbose) const
{
  string::size_type semi = description.rfind(';');
  if (semi == string::npos)
    throw Exception("ExtendedNewick::parenthesisToPhyloDAG(). Bad format: no semi-colon found.");
  string content = description.substr(0, semi);
  unsigned int nodeCounter = 0;
  unsigned int branchCounter = 0;
  map<std::string, shared_ptr<PhyloNode> > mapEvent;
  
  auto dag = make_unique<PhyloDAG>();
  shared_ptr<PhyloNode> root = parenthesisToNode(*dag, 0, content, nodeCounter, branchCounter, mapEvent, withId, verbose);
  dag->rootAt(root);
  if (verbose)
  {
    (*ApplicationTools::message) << " nodes loaded.";
    ApplicationTools::message->endLine();
  }

  return dag;
}

/**********************************************************/
/*  OUTPUT */
/**********************************************************/

void ExtendedNewick::write_(const PhyloDAG& dag, ostream& out) const
{
  // Checking the existence of specified file, and possibility to open it in write mode
  if (!out)
  {
    throw IOException ("ExtendedNewick::writePhyloDAG: failed to write to stream");
  }
  out << dagToParenthesis(dag, writeId_);
}


/******************************************************************************/

void ExtendedNewick::write_(const vector<const PhyloDAG*>& dags, ostream& out) const
{
  // Checking the existence of specified file, and possibility to open it in write mode
  if (!out)
  {
    throw IOException ("ExtendedNewick::write: failed to write to stream");
  }
  for (unsigned int i = 0; i < dags.size(); i++)
    out << dagToParenthesis(*dags[i], writeId_);

}

/******************************************************************************/

string ExtendedNewick::edgeToParenthesis(const PhyloDAG& dag, const std::shared_ptr<PhyloBranch> edge,   std::vector<std::shared_ptr<PhyloNode>>& writtenNodes, bool writeId) const
{
  ostringstream s;
  shared_ptr<PhyloNode> node = dag.getSon(edge);

  if (std::find(writtenNodes.begin(), writtenNodes.end(), node)!=writtenNodes.end())
  {
    s << node->getName();
    if (edge->hasLength())
      s << ":" << edge->getLength();
    return s.str();
  }
  
  if (dag.getNumberOfSons(node) != 0)
  {
    s << "(";

    vector<shared_ptr<PhyloBranch>> vEdges = dag.getOutgoingEdges(node);

    for (vector<shared_ptr<PhyloBranch>>::const_iterator it = vEdges.begin(); it != vEdges.end(); it++)
    {
      if (it != vEdges.begin())
        s << ",";

      s << edgeToParenthesis(dag, *it, writtenNodes,  writeId);
    }
    s << ")";
  }
  s << node->getName();

  if (writeId)
  {
    if (dag.isLeaf(node))
      s << "_";
    s << dag.getNodeIndex(node);
  }

  if (edge->hasLength())
    s << ":" << edge->getLength();

  writtenNodes.push_back(node);
  
  return s.str();
}

/******************************************************************************/

string ExtendedNewick::dagToParenthesis(const PhyloDAG& dag, bool writeId) const
{
  ostringstream s;
  s << "(";

  shared_ptr<PhyloNode>  root = dag.getRoot();

  std::vector<shared_ptr<PhyloNode>> writtenNodes;
  
  std::vector<shared_ptr<PhyloBranch>> rEdges = dag.getOutgoingEdges(root);

  if (dag.isRooted())
  {
    for (size_t i = 0; i < rEdges.size(); ++i)
    {
      if (i != 0)
        s << ",";
      s << edgeToParenthesis(dag, rEdges[i], writtenNodes, writeId);
    }
  }
  else
  {
    s << root->getName();

    for (size_t i = 0; i < rEdges.size(); ++i)
    {
      if (i != 0)
        s << ",";
      s << edgeToParenthesis(dag, rEdges[i], writtenNodes, writeId);
    }
  }

  s << ")";

  s << ";" << endl;

  return s.str();
}
