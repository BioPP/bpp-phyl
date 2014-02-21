//
// File: Nhx.cpp
// Created by: Bastien Boussau
// Created on: Thu Oct 19 11:06:03 2010
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

#include "Nhx.h"
#include "../Tree/Tree.h"
#include "../Tree/TreeTemplate.h"

//From bpp-core:
#include <Bpp/Text/TextTools.h>
#include <Bpp/BppString.h>
#include <Bpp/BppBoolean.h>
#include <Bpp/Numeric/Number.h>

using namespace bpp;

// From the STL:
#include <iostream>
#include <fstream>

using namespace std;

/******************************************************************************/

Nhx::Nhx(bool useTagsAsPptNames):
  supportedProperties_(),
  useTagsAsPropertyNames_(useTagsAsPptNames),
  hasIds_(false)
{
  registerProperty(Property("Gene name", "GN", false, 0));
  registerProperty(Property("Sequence accession", "AC", false, 0));
  registerProperty(Property("Node ID", "ND", false, 0));
  registerProperty(Property(TreeTools::BOOTSTRAP, "B", true, 2));
  registerProperty(Property("Event", "Ev", true, 0));
  registerProperty(Property("EC number", "E", false, 0));
  registerProperty(Property("Function", "Fu", false, 0));
  registerProperty(Property("Domain structure", "DS", false, 0));
  registerProperty(Property("Species name", "S", false, 0));
  registerProperty(Property("Taxonomy ID", "T", false, 1));
  registerProperty(Property("Width of parent branch", "W", true, 1));
  registerProperty(Property("Color of parent branch", "C", true, 0));
  registerProperty(Property("Collapse", "C", false, 3));
  registerProperty(Property("Custom", "XB", true, 0));
  registerProperty(Property("Custom", "XN", false, 0));
  registerProperty(Property("Orthologous", "O", false, 1));
  registerProperty(Property("Subtree neighbors", "SN", false, 1));
  registerProperty(Property("Super orthologous", "SO", false, 1));
}

/******************************************************************************/

const string Nhx::getFormatName() const { return "Nhx"; }

/******************************************************************************/

const string Nhx::getFormatDescription() const
{
  return string("New Hampshire eXtended parenthesis format. ") +
    "See http://www.phylosoft.org/NHX/ for more info.";
}

/******************************************************************************/

#if defined(NO_VIRTUAL_COV)
    Tree*
#else
    TreeTemplate<Node>* 
#endif
Nhx::read(istream& in) const throw (Exception)
{
  // Checking the existence of specified file
  if (! in) { throw IOException ("Nhx ::read: failed to read from stream"); }
  
  //We concatenate all line in file till we reach the ending semi colon:
  string temp, description;// Initialization
  // Main loop : for all file lines
  while (! in.eof())
  {
    getline(in, temp, '\n');  // Copy current line in temporary string
    string::size_type index = temp.find(";");
    if (index != string::npos)
    {
      description += temp.substr(0, index + 1);
      break;
    }
    else description += temp;
  }
  vector<string> beginnings, endings;
  beginnings.push_back("[&&NHX:");
  description = TextTools::removeSubstrings(description, '[', ']', beginnings, endings);
  return parenthesisToTree(description);
}

/******************************************************************************/

void Nhx::write_(const Tree& tree, ostream& out) const throw (Exception)
{
  // Checking the existence of specified file, and possibility to open it in write mode
  if (! out) { throw IOException ("Nhx::writeTree: failed to write to stream"); }
    out << treeToParenthesis(tree);
}

/******************************************************************************/

template<class N>
void Nhx::write_(const TreeTemplate<N>& tree, ostream& out) const throw (Exception)
{
  // Checking the existence of specified file, and possibility to open it in write mode
  if (! out) { throw IOException ("Nhx::writeTree: failed to write to stream"); }
    out << treeToParenthesis(tree);
}

/******************************************************************************/

void Nhx::read(istream& in, vector<Tree*>& trees) const throw (Exception)
{
  // Checking the existence of specified file
  if (! in) { throw IOException ("Nhx::read: failed to read from stream"); }
  
  // Main loop : for all file lines
  string temp, description;// Initialization
  string::size_type index;
  vector <string > beginnings, endings;
  beginnings.push_back("[&&NHX:");
  while (!in.eof())
  {
    //We concatenate all line in file till we reach the ending semi colon:
    while (!in.eof())
    {
      getline(in, temp, '\n');  // Copy current line in temporary string
      index = temp.find(";");
      if (index != string::npos)
      {
        description += temp.substr(0, index + 1);
        description = TextTools::removeSubstrings(description, '[', ']', beginnings, endings);
        trees.push_back(parenthesisToTree(description));
        description = temp.substr(index + 1);
      }
      else description += temp;
    }
  }
}

/******************************************************************************/

void Nhx::write_(const vector<const Tree*>& trees, ostream& out) const throw (Exception)
{
  // Checking the existence of specified file, and possibility to open it in write mode
  if (! out) { throw IOException ("Nhx::write: failed to write to stream"); }
  for(unsigned int i = 0; i < trees.size(); i++)
  {
    if (dynamic_cast<const TreeTemplate<Node>* >(trees[i]))
      out << treeToParenthesis(*dynamic_cast<const TreeTemplate<Node>* >(trees[i]));
  }
}

/******************************************************************************/

template<class N>
void Nhx::write_(const vector<TreeTemplate<N>*>& trees, ostream& out) const throw (Exception)
{
  // Checking the existence of specified file, and possibility to open
  // it in write mode
  if (! out) { throw IOException ("Nhx::write: failed to write to stream"); }
  for(unsigned int i = 0; i < trees.size(); i++)
  {
    out << treeToParenthesis(*trees[i]);
  }
}

/******************************************************************************/

Nhx::Element Nhx::getElement(const string& elt) const throw (IOException)
{
  Element element;
  element.length     = ""; //default
  element.annotation = ""; //default
  element.isLeaf     = false; // default
 
  //cout << "ELT=" << elt << endl;
  size_t lastP = elt.rfind(")"), firstP;
  size_t beginAnno = string::npos;
  if (lastP == string::npos)
    beginAnno = elt.rfind("[&&NHX:");
  else
    beginAnno = elt.find("[&&NHX:", lastP + 1);
  string elementWithoutAnnotation;
  if (beginAnno != string::npos) {
    size_t endAnno = elt.find("]", beginAnno + 7);
    element.annotation = elt.substr(beginAnno + 7, endAnno - beginAnno - 7);
    elementWithoutAnnotation = elt.substr(0, beginAnno);
  } else {
    element.annotation = "";
    elementWithoutAnnotation = elt;
  } 
  //cout << "ANNO=" << element.annotation << endl;
  //cout << "ELT =" << elementWithoutAnnotation << endl;

  size_t colonIndex;
  bool hasColon = false;
  for (colonIndex = elementWithoutAnnotation.size() - 1; colonIndex > 0 && elementWithoutAnnotation[colonIndex] != ')' && !hasColon; --colonIndex)
  {
    if (elementWithoutAnnotation[colonIndex] == ':')
    {
      hasColon = true;
    }
  }
  try
  {
    string elt2;
    if (hasColon)
    {
      //this is an element with length:
      elt2 = elementWithoutAnnotation.substr(0, colonIndex + 1);
      element.length = elementWithoutAnnotation.substr(colonIndex + 2);
    }
    else
    {
      //this is an element without length;
      elt2 = elementWithoutAnnotation;
    }
  
    lastP = elt2.rfind(')');
    firstP = elt2.find('(');
    if (firstP == string::npos)
    {
      //This is a leaf:
      element.content = elt2;
      element.isLeaf  = true;
    }
    else
    {
      //This is a node:
      if (lastP < firstP) throw IOException("Nhx::getElement(). Invalid format: bad closing parenthesis in " + elt2);
      element.content = elt2.substr(firstP + 1, lastP - firstP - 1);
    }
  }
  catch (exception& e)
  {
    throw IOException("Bad tree description: " + elt);
  }
  //cout << endl;
  //cout << "CONTENT:" << endl << element.content << endl;
  //cout << endl;
  //cout << "ANNOTATION:" << endl << element.annotation << endl;
  //cout << endl;
 
  return element;
}  

/******************************************************************************/


Node* Nhx::parenthesisToNode(const string& description) const
{
  //cout << "NODE: " << description << endl;
  Element elt = getElement(description);
  
  //New node:
  Node* node = new Node();
  if (!TextTools::isEmpty(elt.length))
  {
    node->setDistanceToFather(TextTools::toDouble(elt.length));
  }
  if (!TextTools::isEmpty(elt.annotation))
  {
    bool hasId = setNodeProperties(*node, elt.annotation);
    hasIds_ |= hasId;
    if (hasIds_ && !hasId)
      throw Exception("Nhx::parenthesisToNode. At least one node is missing an id (ND tag).");
  }
 
  NestedStringTokenizer nt(elt.content, "(", ")", ",");
  vector<string> elements;
  while (nt.hasMoreToken())
  {
    elements.push_back(nt.nextToken());
  }
  if (elt.isLeaf)
  {
    //This is a leaf:
    string name = TextTools::removeSurroundingWhiteSpaces(elements[0]);
    node->setName(name);
  }
  else
  {
    //This is a node:
    for (size_t i = 0; i < elements.size(); ++i)
    {
      //cout << "NODE: SUBNODE: " << i << ", " << elements[i] << endl;
      Node* son = parenthesisToNode(elements[i]);
      node->addSon(son);
    }
  }
  return node;
}

/******************************************************************************/

TreeTemplate<Node>* Nhx::parenthesisToTree(const string& description) const throw (Exception) 
{
  hasIds_ = false;
  string::size_type semi = description.rfind(';');
  if (semi == string::npos)
    throw Exception("Nhx::parenthesisToTree(). Bad format: no semi-colon found.");
  string content = description.substr(0, semi);
  Node* node = parenthesisToNode(content);
  TreeTemplate<Node>* tree = new TreeTemplate<Node>();
  tree->setRootNode(node);
  if (!hasIds_)
  {
    tree->resetNodesId();
  }
  return tree;
}

/******************************************************************************/

string Nhx::propertyToString_(const Clonable* pptObject, short type) throw (Exception)
{
  if (type == 0) {
    const BppString* castedPptObject = dynamic_cast<const BppString*>(pptObject);
    if (castedPptObject)
      return castedPptObject->toSTL();
    else
      throw Exception("Nhx::propertyToString_. Unvalid property type, should be of class BppString.");
  } else if (type == 1) {
    const Number<int>* castedPptObject = dynamic_cast<const Number<int>*>(pptObject);
    if (castedPptObject)
      return TextTools::toString(castedPptObject->getValue());
    else
      throw Exception("Nhx::propertyToString_. Unvalid property type, should be of class Number<int>.");
  } else if (type == 2) {
    const Number<double>* castedPptObject = dynamic_cast<const Number<double>*>(pptObject);
    if (castedPptObject)
      return TextTools::toString(castedPptObject->getValue());
    else
      throw Exception("Nhx::propertyToString_. Unvalid property type, should be of class Number<double>.");
  } else if (type == 3) {
    const BppBoolean* castedPptObject = dynamic_cast<const BppBoolean*>(pptObject);
    if (castedPptObject)
      return TextTools::toString(castedPptObject->getValue());
    else
      throw Exception("Nhx::propertyToString_. Unvalid property type, should be of class BppBoolean.");
  } else {
    throw Exception("Nhx::propertyToString_. Unsupported type: " + TextTools::toString(type));
  }
}

/******************************************************************************/

Clonable* Nhx::stringToProperty_(const string& pptDesc, short type) throw (Exception)
{
  if (type == 0) {
    return new BppString(pptDesc);
  } else if (type == 1) {
    return new Number<int>(TextTools::toInt(pptDesc));
  } else if (type == 2) {
    return new Number<double>(TextTools::toDouble(pptDesc));
  } else if (type == 3) {
    return new BppBoolean(TextTools::to<bool>(pptDesc));
  } else {
    throw Exception("Nhx::stringToProperty_. Unsupported type: " + TextTools::toString(type));
  }
}


/******************************************************************************/

string Nhx::propertiesToParenthesis(const Node& node) const
{
  ostringstream s;
  s << "[&&NHX";
  for (set<Property>::iterator it = supportedProperties_.begin(); it != supportedProperties_.end(); ++it) {
    string ppt = (useTagsAsPropertyNames_ ? it->tag : it->name);
    if (it->onBranch) {
      if (node.hasBranchProperty(ppt)) {
        const Clonable* pptObject = node.getBranchProperty(ppt);
        s << ":" << it->tag << "=" << propertyToString_(pptObject, it->type);
      }
    } else {
      if (node.hasNodeProperty(ppt)) {
        const Clonable* pptObject = node.getNodeProperty(ppt);
        s << ":" << it->tag << "=" << propertyToString_(pptObject, it->type);
      }
    }
  }
  //If no special node id is provided, we output the one from the tree:
  if (!node.hasNodeProperty(useTagsAsPropertyNames_ ? "ND" : "Node ID"))
  {
    s << ":ND="<<TextTools::toString(node.getId());
  }
  s << "]";
  return s.str();  
}

/******************************************************************************/

string Nhx::nodeToParenthesis(const Node& node) const
{
  ostringstream s;
  if (node.isLeaf())
  {
    s << node.getName();
  }
  else
  {
    s << "(";
    s << nodeToParenthesis(* node[0]);
    for (int i = 1; i < static_cast<int>(node.getNumberOfSons()); i++)
    {
      s << "," << nodeToParenthesis(*node[i]);
    }
    s << ")";
  }
  if (node.hasDistanceToFather()) s << ":" << node.getDistanceToFather();
  s << propertiesToParenthesis(node);
  return s.str();  
}

/******************************************************************************/

string Nhx::treeToParenthesis(const TreeTemplate<Node>& tree) const
{
  ostringstream s;
  s << "(";

  const Node* node = tree.getRootNode();

  if (node->isLeaf())
  {
    s << node->getName();
    for (size_t i = 0; i < node->getNumberOfSons(); ++i)
    {
      s << "," << nodeToParenthesis(*node->getSon(i));
    }
  }
  else
  {
    s << nodeToParenthesis(* node->getSon(0));
    for (size_t i = 1; i < node->getNumberOfSons(); ++i)
    {
      s << "," << nodeToParenthesis(*node->getSon(i));
    }
  }

  s << ")" ;
  if (node->hasDistanceToFather()) s << ":" << node->getDistanceToFather();
  s << propertiesToParenthesis(*node);
  s << ";" << endl;
  return s.str();  
}

/******************************************************************************/

bool Nhx::setNodeProperties(Node& node, const string properties) const
{
  string propsDesc = TextTools::removeChar(properties, ']');
  StringTokenizer st(propsDesc, ":", true, true);
  map<string, string> props; 
  while (st.hasMoreToken())
  {
    string token = st.nextToken();
    if (TextTools::hasSubstring(token, "=")) {
      StringTokenizer pt(token, "=", true, true);
      string tag = pt.nextToken();
      string value = pt.nextToken();
      props[tag] = value;
    }
  }

  for (set<Property>::iterator it = supportedProperties_.begin(); it != supportedProperties_.end(); ++it) {
    if (props.find(it->tag) != props.end()) {
      //Property found
      string ppt = (useTagsAsPropertyNames_ ? it->tag : it->name);
      if (it->onBranch) {
        node.setBranchProperty(ppt, *auto_ptr<Clonable>(stringToProperty_(props[it->tag], it->type)));
      } else {
        node.setNodeProperty(ppt, *auto_ptr<Clonable>(stringToProperty_(props[it->tag], it->type)));
      }
    }
  }
     
  //If the ND tag is present and is decimal, we use it has the node id:
  bool hasId = false;
  
  if (props.find("ND") != props.end()) {
    string prop = props["ND"];
    if (TextTools::isDecimalNumber(prop))
    {
      node.setId(TextTools::toInt(prop));
      hasId = true;   
    }
  }
  return hasId;
}

/******************************************************************************/

void Nhx::changeTagsToNames(Node& node) const {
  for (set<Property>::iterator it = supportedProperties_.begin(); it != supportedProperties_.end(); ++it) {
    if (it->onBranch) {
      if (node.hasBranchProperty(it->tag)) {
        node.setBranchProperty(it->name, *node.getBranchProperty(it->tag));
        node.deleteBranchProperty(it->tag);
      }
    } else {
      if (node.hasNodeProperty(it->tag)) {
        node.setNodeProperty(it->name, *node.getNodeProperty(it->tag));
        node.deleteNodeProperty(it->tag);
      }
    }
  }
  for (unsigned int i = 0; i < node.getNumberOfSons(); ++i)
    changeTagsToNames(*node.getSon(i));
}

/******************************************************************************/

void Nhx::changeNamesToTags(Node& node) const {
  for (set<Property>::iterator it = supportedProperties_.begin(); it != supportedProperties_.end(); ++it) {
    if (it->onBranch) {
      if (node.hasBranchProperty(it->name)) {
        node.setBranchProperty(it->tag, *node.getBranchProperty(it->name));
        node.deleteBranchProperty(it->name);
      }
    } else {
      if (node.hasNodeProperty(it->name)) {
        node.setNodeProperty(it->tag, *node.getNodeProperty(it->name));
        node.deleteNodeProperty(it->name);
      }
    }
  }
  for (unsigned int i = 0; i < node.getNumberOfSons(); ++i)
    changeNamesToTags(*node.getSon(i));
}

/******************************************************************************/

