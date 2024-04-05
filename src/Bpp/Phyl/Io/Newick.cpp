// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/BppString.h>
#include <Bpp/Numeric/Number.h>
#include <Bpp/Text/NestedStringTokenizer.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/TextTools.h>

#include "../Tree/PhyloBranch.h"
#include "../Tree/PhyloNode.h"
#include "../Tree/PhyloTree.h"
#include "../Tree/Tree.h"
#include "../Tree/TreeTemplate.h"
#include "../Tree/TreeTemplateTools.h"
#include "Newick.h"

using namespace bpp;

// From the STL:
#include <iostream>
#include <fstream>

using namespace std;

/******************************************************************************/

const string Newick::getFormatName() const { return "Newick"; }

/******************************************************************************/

const string Newick::getFormatDescription() const
{
  return string("New hampshire parenthesis format. ") +
         "See http://evolution.genetics.washington.edu/phylip/newicktree.html for more info.";
}

/**********************************************************/
/*  INPUT */
/**********************************************************/

unique_ptr<TreeTemplate<Node>> Newick::readTreeTemplate(istream& in) const
{
  // Checking the existence of specified file
  if (!in)
  {
    throw IOException ("Newick::read: failed to read from stream");
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
    throw IOException("Newick::read: no tree was found!");
  return TreeTemplateTools::parenthesisToTree(description, useBootstrap_, bootstrapPropertyName_, false, verbose_);
}

/*********************************************************************************/

unique_ptr<PhyloTree> Newick::readPhyloTree(istream& in) const
{
  // Checking the existence of specified file
  if (!in)
  {
    throw IOException ("Newick::readPhyloTree: failed to read from stream");
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
    throw IOException("Newick::read: no tree was found!");
  return parenthesisToPhyloTree(description, useBootstrap_, bootstrapPropertyName_, false, verbose_);
}

/******************************************************************************/

void Newick::readTrees(istream& in, vector<unique_ptr<Tree>>& trees) const
{
  // Checking the existence of specified file
  if (!in)
  {
    throw IOException ("Newick::readTrees(vector): failed to read from stream");
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
      trees.push_back(TreeTemplateTools::parenthesisToTree(description, useBootstrap_, bootstrapPropertyName_, false, verbose_));
      description = temp.substr(index + 1);
    }
    else
      description += temp;
  }
  // In case the file is empty, the method will not add any neww tree to the vector.
}

/******************************************************************************/

void Newick::readPhyloTrees(istream& in, vector<unique_ptr<PhyloTree>>& trees) const
{
  // Checking the existence of specified file
  if (!in)
  {
    throw IOException ("Newick::readTrees(vector): failed to read from stream");
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
      trees.push_back(parenthesisToPhyloTree(description, useBootstrap_, bootstrapPropertyName_, false, verbose_));
      description = temp.substr(index + 1);
    }
    else
      description += temp;
  }
  // In case the file is empty, the method will not add any neww tree to the vector.
}

/***************************************/

IOTree::Element Newick::getElement(const string& elt) const
{
  IOTree::Element element;
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
      element.isLeaf = true;
    }
    else
    {
      // This is a node:
      if (lastP < firstP)
        throw IOException("Newick::getElement(). Invalid format: bad closing parenthesis in " + elt2);
      element.content = TextTools::removeSurroundingWhiteSpaces(elt2.substr(firstP + 1, lastP - firstP - 1));
      string bootstrap = TextTools::removeSurroundingWhiteSpaces(elt2.substr(lastP + 1));
      // cout << "ELEMENT: BOOTSTRAP: " << bootstrap << endl;
      if (!TextTools::isEmpty(bootstrap))
      {
        element.annotation = bootstrap;
      }
    }
  }
  catch (exception& e)
  {
    throw IOException("Bad tree description: " + elt);
  }
  return element;
}

/************************************************************/

shared_ptr<PhyloNode>  Newick::parenthesisToNode(PhyloTree& tree, shared_ptr<PhyloNode>  father, const string& description, unsigned int& nodeCounter, bool bootstrap, const string& propertyName, bool withId, bool verbose) const
{
//  cout << "NODE: " << description << endl;
  IOTree::Element elt = getElement(description);

  // New node:
  std::shared_ptr<PhyloNode> node(new PhyloNode());

  shared_ptr<PhyloBranch> branch(father ? new PhyloBranch() : 0);

  if (father)
  {
    tree.createNode(father, node, branch);

    if (!TextTools::isEmpty(elt.length))
      branch->setLength(TextTools::toDouble(elt.length));
  }
  else
    tree.createNode(node);

  if (!TextTools::isEmpty(elt.annotation))
  {
    if (withId)
    {
      auto id = static_cast<PhyloTree::NodeIndex>(TextTools::toInt(elt.annotation));
      tree.setNodeIndex(node, id);
      if (branch)
        tree.setEdgeIndex(branch, id);
    }
    else
    {
      if (bootstrap)
      {
        if (branch)
          branch->setProperty("bootstrap", Number<double>(TextTools::toDouble(elt.annotation)));
        // cout << "NODE: BOOTSTRAP: " << * elt.bootstrap << endl;
      }
      else
      {
        if (branch)
          branch->setProperty(propertyName, BppString(elt.annotation));
      }
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
        {
          realName << "_";
        }
        realName << st.getToken(i);
      }
      node->setName(realName.str());
      tree.setNodeIndex(node, static_cast<PhyloTree::NodeIndex>(
            TextTools::toInt(st.getToken(st.numberOfRemainingTokens() - 1))));
      if (branch)
        tree.setEdgeIndex(branch, static_cast<PhyloTree::NodeIndex>(
              TextTools::toInt(st.getToken(st.numberOfRemainingTokens() - 1))));
    }
    else
      node->setName(name);
  }
  else
  {
    // This is a node:
    for (size_t i = 0; i < elements.size(); i++)
    {
      //    cout << "NODE: SUBNODE: " << i << ", " << elements[i] << endl;
      parenthesisToNode(tree, node, elements[i], nodeCounter, bootstrap, propertyName, withId, verbose);
    }
  }

  if (!withId)
  {
    tree.setNodeIndex(node, nodeCounter);
    if (branch)
      tree.setEdgeIndex(branch, nodeCounter);
  }

  nodeCounter++;
  if (verbose)
    ApplicationTools::displayUnlimitedGauge(nodeCounter);
  return node;
}

/******************************************************************************/

unique_ptr<PhyloTree> Newick::parenthesisToPhyloTree(const string& description, bool bootstrap, const string& propertyName, bool withId, bool verbose) const
{
  string::size_type semi = description.rfind(';');
  if (semi == string::npos)
    throw Exception("Newick::parenthesisToTree(). Bad format: no semi-colon found.");
  string content = description.substr(0, semi);
  unsigned int nodeCounter = 0;
  auto tree = make_unique<PhyloTree>();
  shared_ptr<PhyloNode> root = parenthesisToNode(*tree, 0, content, nodeCounter, bootstrap, propertyName, withId, verbose);
  tree->rootAt(root);
  if (verbose)
  {
    (*ApplicationTools::message) << " nodes loaded.";
    ApplicationTools::message->endLine();
  }

  return tree;
}

/**********************************************************/
/*  OUTPUT */
/**********************************************************/

void Newick::write_(const Tree& tree, ostream& out) const
{
  // Checking the existence of specified file, and possibility to open it in write mode
  if (!out)
  {
    throw IOException ("Newick::writeTree: failed to write to stream");
  }
  if (useBootstrap_)
  {
    out << TreeTools::treeToParenthesis(tree, writeId_);
  }
  else
  {
    out << TreeTools::treeToParenthesis(tree, false, bootstrapPropertyName_);
  }
}

void Newick::write_(const PhyloTree& tree, ostream& out) const
{
  // Checking the existence of specified file, and possibility to open it in write mode
  if (!out)
  {
    throw IOException ("Newick::writeTree: failed to write to stream");
  }
  if (useBootstrap_)
  {
    out << treeToParenthesis(tree, writeId_);
  }
  else
  {
    out << treeToParenthesis(tree, false, bootstrapPropertyName_);
  }
}

/******************************************************************************/

template<class N>
void Newick::write_(const TreeTemplate<N>& tree, ostream& out) const
{
  // Checking the existence of specified file, and possibility to open it in write mode
  if (!out)
  {
    throw IOException ("Newick::writeTree: failed to write to stream");
  }
  if (useBootstrap_)
  {
    out << TreeTemplateTools::treeToParenthesis(tree, writeId_);
  }
  else
  {
    out << TreeTemplateTools::treeToParenthesis(tree, false, bootstrapPropertyName_);
  }
}


/******************************************************************************/

void Newick::write_(const vector<const Tree*>& trees, ostream& out) const
{
  // Checking the existence of specified file, and possibility to open it in write mode
  if (!out)
  {
    throw IOException ("Newick::write: failed to write to stream");
  }
  for (unsigned int i = 0; i < trees.size(); i++)
  {
    if (useBootstrap_)
    {
      out << TreeTools::treeToParenthesis(*trees[i], writeId_);
    }
    else
    {
      out << TreeTools::treeToParenthesis(*trees[i], false, bootstrapPropertyName_);
    }
  }
}

/******************************************************************************/

template<class N>
void Newick::write_(const vector<TreeTemplate<N>*>& trees, ostream& out) const
{
  // Checking the existence of specified file, and possibility to open it in write mode
  if (!out)
  {
    throw IOException ("Newick::write: failed to write to stream");
  }
  for (unsigned int i = 0; i < trees.size(); i++)
  {
    if (useBootstrap_)
    {
      out << TreeTemplateTools::treeToParenthesis(*trees[i], writeId_);
    }
    else
    {
      out << TreeTemplateTools::treeToParenthesis(*trees[i], false, bootstrapPropertyName_);
    }
  }
}

/******************************************************************************/

void Newick::write_(const vector<const PhyloTree*>& trees, ostream& out) const
{
  // Checking the existence of specified file, and possibility to open it in write mode
  if (!out)
  {
    throw IOException ("Newick::write: failed to write to stream");
  }
  for (unsigned int i = 0; i < trees.size(); i++)
  {
    if (useBootstrap_)
    {
      out << treeToParenthesis(*trees[i], writeId_);
    }
    else
    {
      out << treeToParenthesis(*trees[i], false, bootstrapPropertyName_);
    }
  }
}

/******************************************************************************/

string Newick::nodeToParenthesis(const PhyloTree& tree, const std::shared_ptr<PhyloNode> node, bool writeId) const
{
  ostringstream s;
  shared_ptr<PhyloBranch> branch = tree.hasFather(node) ? tree.getEdgeToFather(node) : 0;

  if (tree.getNumberOfSons(node) == 0)
  {
    s << node->getName();
  }
  else
  {
    s << "(";

    vector<shared_ptr<PhyloNode>> vSons = tree.getSons(node);

    for (vector<shared_ptr<PhyloNode>>::const_iterator it = vSons.begin(); it != vSons.end(); it++)
    {
      if (it != vSons.begin())
        s << ",";

      s << nodeToParenthesis(tree, *it);
    }

    s << ")";
  }

  if (writeId)
  {
    if (tree.isLeaf(node))
      s << "_";
    s << tree.getNodeIndex(node);
  }
  else
  {
    if (branch && branch->hasProperty("bootstrap"))
      s << (dynamic_cast<const Number<double>*>(branch->getProperty("bootstrap"))->getValue());
  }

  if (branch && branch->hasLength())
    s << ":" << branch->getLength();
  return s.str();
}

/******************************************************************************/

string Newick::nodeToParenthesis(const PhyloTree& tree, const std::shared_ptr<PhyloNode> node, bool bootstrap, const string& propertyName) const
{
  ostringstream s;
  shared_ptr<PhyloBranch> branch = tree.hasFather(node) ? tree.getEdgeToFather(node) : 0;

  if (tree.getNumberOfSons(node) == 0)
  {
    s << node->getName();
  }
  else
  {
    s << "(";

    vector<shared_ptr<PhyloNode>> vSons = tree.getSons(node);

    for (vector<shared_ptr<PhyloNode>>::const_iterator it = vSons.begin(); it != vSons.end(); it++)
    {
      if (it != vSons.begin())
        s << ",";

      s << nodeToParenthesis(tree, *it, bootstrap, propertyName);
    }

    s << ")";

    if (branch)
    {
      if (bootstrap)
      {
        if (branch->hasProperty("bootstrap"))
          s << (dynamic_cast<const Number<double>*>(branch->getProperty("bootstrap"))->getValue());
      }
      else
      {
        if (node->hasProperty(propertyName))
          s << *(dynamic_cast<const BppString*>(node->getProperty(propertyName)));
      }
    }
  }

  if (branch && branch->hasLength())
    s << ":" << branch->getLength();

  return s.str();
}

/******************************************************************************/

string Newick::treeToParenthesis(const PhyloTree& tree, bool writeId) const
{
  ostringstream s;
  s << "(";

  shared_ptr<PhyloNode>  root = tree.getRoot();

  std::vector<shared_ptr<PhyloNode>> rSons = tree.getSons(root);

  if (tree.isRooted())
  {
    for (size_t i = 0; i < rSons.size(); ++i)
    {
      if (i != 0)
        s << ",";
      s << nodeToParenthesis(tree, rSons[i], writeId);
    }
  }
  else
  {
    s << root->getName();

    for (size_t i = 0; i < rSons.size(); ++i)
    {
      if (i != 0)
        s << ",";
      s << nodeToParenthesis(tree, rSons[i], writeId);
    }
  }

  s << ")";

  const shared_ptr<PhyloBranch> branch = tree.hasFather(root) ? tree.getEdgeToFather(root) : 0;

  if (branch && branch->hasLength())
    s << ":" << branch->getLength();
  s << ";" << endl;

  return s.str();
}

/******************************************************************************/

string Newick::treeToParenthesis(const PhyloTree& tree, bool bootstrap, const string& propertyName) const
{
  ostringstream s;
  s << "(";

  shared_ptr<PhyloNode>  root = tree.getRoot();

  std::vector<shared_ptr<PhyloNode>> rSons = tree.getSons(root);

  if (tree.isRooted())
  {
    for (size_t i = 0; i < rSons.size(); ++i)
    {
      if (i != 0)
        s << ",";
      s << nodeToParenthesis(tree, rSons[i], bootstrap, propertyName);
    }
  }
  else
  {
    s << root->getName();

    for (size_t i = 0; i < rSons.size(); ++i)
    {
      if (i != 0)
        s << ",";
      s << nodeToParenthesis(tree, rSons[i], bootstrap, propertyName);
    }
  }

  s << ")";

  shared_ptr<PhyloBranch> branch = tree.hasFather(root) ? tree.getEdgeToFather(root) : 0;

  if (branch)
  {
    if (bootstrap)
    {
      if (branch->hasProperty("bootstrap"))
        s << (dynamic_cast<const Number<double>*>(branch->getProperty("bootstrap"))->getValue());
    }
    else
    {
      if (branch->hasProperty(propertyName))
        s << *(dynamic_cast<const BppString*>(branch->getProperty(propertyName)));
    }
  }

  s << ";" << endl;
  return s.str();
}
