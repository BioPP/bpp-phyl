//
// File: NHX.cpp
// Created by: Bastien Boussau
// Created on: Thu Oct 19 11:06:03 2010
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004)

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

#include "NHX.h"
#include "../Tree.h"
#include "../TreeTemplate.h"

#include <Bpp/Text/TextTools.h>

using namespace bpp;

// From the STL:
#include <iostream>
#include <fstream>

using namespace std;

/******************************************************************************/

const string NHX::getFormatName() const { return "NHX"; }

/******************************************************************************/

const string NHX::getFormatDescription() const
{
	return string("New Hampshire eXtended parenthesis format. ") +
		"See http://www.phylosoft.org/NHX/ for more info.";
}

/******************************************************************************/

#if defined(NO_VIRTUAL_COV)
		Tree*
#else
		TreeTemplate<Node> * 
#endif
NHX::read(istream& in) const throw (Exception)
{
	// Checking the existence of specified file
	if (! in) { throw IOException ("NHX ::read: failed to read from stream"); }
	
	//We concatenate all line in file till we reach the ending semi colon:
	string temp, description;// Initialization
	// Main loop : for all file lines
	while (! in.eof())
  {
		getline(in, temp, '\n');  // Copy current line in temporary string
    string::size_type index = temp.find(";");
		if(index != string::npos)
    {
			description += temp.substr(0, index + 1);
			break;
		}
    else description += temp;
	}
  vector <string > beginnings, endings;
  beginnings.push_back("[&&NHX:");
	description = TextTools::removeSubstrings(description, '[', ']', beginnings, endings);
	return parenthesisToTree(description);
}

/******************************************************************************/

void NHX::write_(const Tree& tree, ostream& out) const throw (Exception)
{
	// Checking the existence of specified file, and possibility to open it in write mode
	if (! out) { throw IOException ("NHX::writeTree: failed to write to stream"); }
    out << treeToParenthesis(tree);
}

/******************************************************************************/

template<class N>
void NHX::write_(const TreeTemplate<N>& tree, ostream& out) const throw (Exception)
{
	// Checking the existence of specified file, and possibility to open it in write mode
	if (! out) { throw IOException ("NHX::writeTree: failed to write to stream"); }
    out << treeToParenthesis(tree);
}

/******************************************************************************/

void NHX::read(istream& in, vector<Tree*>& trees) const throw (Exception)
{
	// Checking the existence of specified file
	if (! in) { throw IOException ("NHX::read: failed to read from stream"); }
	
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

void NHX::write_(const vector<Tree*>& trees, ostream& out) const throw (Exception)
{
	// Checking the existence of specified file, and possibility to open it in write mode
	if (! out) { throw IOException ("NHX::write: failed to write to stream"); }
  for(unsigned int i = 0; i < trees.size(); i++)
  {
      out << treeToParenthesis(*trees[i]);
  }
}

/******************************************************************************/

template<class N>
void NHX::write_(const vector<TreeTemplate<N>*>& trees, ostream& out) const throw (Exception)
{
	// Checking the existence of specified file, and possibility to open it in write mode
	if (! out) { throw IOException ("NHX::write: failed to write to stream"); }
  for(unsigned int i = 0; i < trees.size(); i++)
  {
      out << treeToParenthesis(*trees[i]);
  }
}

/******************************************************************************/

NHX::Element NHX::getElement(const string& elt) const throw (IOException)
{
  Element element;
  element.length    = ""; //default
  element.annotation = ""; //default
    
  StringTokenizer st(elt, "[&&NHX:", true, true);
  
  StringTokenizer st2(st.getToken(st.numberOfRemainingTokens()-1), "]", true, true);
  element.annotation = st2.getToken(0);
  
  
  unsigned int colonIndex;
  bool hasColon = false;
  string::size_type lastAnnot = elt.rfind("[&&NHX");
  string elementWithoutAnnotation=elt.substr(0,lastAnnot-1);
  for (colonIndex = elementWithoutAnnotation.size(); colonIndex > 0 && elementWithoutAnnotation[colonIndex] != ')'; colonIndex--)
    {
    if (elementWithoutAnnotation[colonIndex] == ':')
      {
      hasColon = true;
      break;
      }
    }
  try
  {
  string elt2;
  if(hasColon)
    {
    //this is an element with length:
    elt2 = elementWithoutAnnotation.substr(0, colonIndex);
  //  if (st2.numberOfRemainingTokens()>1)
    element.length = elt.substr(colonIndex + 1);
    }
  else
    {
    //this is an element without length;
   // elt2 = st.getToken(0);
    //if (st2.numberOfRemainingTokens()>1)
      elt2 = elementWithoutAnnotation;
    }
  
  string::size_type lastP = elt2.rfind(')');
  string::size_type firstP = elt2.find('(');
  if(firstP == string::npos)
    {
    //This is a leaf:
    element.content = elt2;
    }
  else
    {
    //This is a node:
    if(lastP < firstP) throw IOException("NHX::getElement(). Invalid format: bad closing parenthesis in " + elt2);
    element.content = elt2.substr(firstP + 1, lastP - firstP - 1);
    }
  }
  catch(exception e)
  {
  throw IOException("Bad tree description: " + elt);
  }
  return element;
}  

/******************************************************************************/


Node* NHX::parenthesisToNode(const string& description) const
{
  //cout << "NODE: " << description << endl;
  Element elt = getElement(description);
  //New node:
  Node* node = new Node();
  if(!TextTools::isEmpty(elt.length))
    {
    node->setDistanceToFather(TextTools::toDouble(elt.length));
    }
  if(!TextTools::isEmpty(elt.annotation))
    {
    setNodeProperties(*node, elt.annotation);
    }
  
  NestedStringTokenizer nt(elt.content, "(", ")", ",");
  vector<string> elements;
  while (nt.hasMoreToken())
    {
    elements.push_back(nt.nextToken());
    }
  if(elements.size() == 1)
    {
    //This is a leaf:
    
    string name = TextTools::removeSurroundingWhiteSpaces(elements[0]);
    node->setName(name);
    }
  else
    {
    //This is a node:
    for(unsigned int i = 0; i < elements.size(); i++)
      {
      //cout << "NODE: SUBNODE: " << i << ", " << elements[i] << endl;
      Node* son = parenthesisToNode(elements[i]);
      node->addSon(son);
      }
    }
  return node;
}

/******************************************************************************/

TreeTemplate<Node>* NHX::parenthesisToTree(const string& description) const throw (Exception) 
{
  string::size_type lastP  = description.rfind(')');
  bool hasId = false;
  if (lastP == string::npos)
    throw Exception("NHX::parenthesisToTree(). Bad format: no closing parenthesis found.");
  string::size_type firstP = description.find('(');
  if (firstP == string::npos)
    throw Exception("NHX::parenthesisToTree(). Bad format: no opening parenthesis found.");
  string::size_type semi = description.rfind(';');
  if (semi == string::npos)
    throw Exception("NHX::parenthesisToTree(). Bad format: no semi-colon found.");
  if (lastP <= firstP)
    throw Exception("NHX::parenthesisToTree(). Bad format: closing parenthesis before opening parenthesis.");
  string content = description.substr(firstP + 1, lastP - firstP - 1);
  string element = (semi == string::npos) ? description.substr(lastP + 1) : description.substr(lastP + 1, semi - lastP - 1);
  //cout << "TREE: " << content << endl;
  //New root node:
  Node* node = new Node();
  NestedStringTokenizer nt(content,"(", ")", ",");
  vector<string> elements;
  while (nt.hasMoreToken())
    {
    elements.push_back(nt.nextToken());
    }
  if (elements.size() == 1)
    {
    //This is a leaf:
    StringTokenizer st(elements[0], "[&&NHX", true, true);
    StringTokenizer st2(st.getToken(0), ":", true, true);
    if (st2.getToken(0) != "") {
      node->setName(st2.getToken(0));
    }
    node->setDistanceToFather(TextTools::toDouble(st2.getToken(1)));
    hasId = setNodeProperties(*node, st.getToken(1));
    }
  else
    {
    //This is a node:
    for(unsigned int i = 0; i < elements.size(); i++)
      {
      Node* son = parenthesisToNode(elements[i]);
      node->addSon(son);
      }
    if(! TextTools::isEmpty(element))
      {
      StringTokenizer st(element, "[&&NHX", true, true);
      StringTokenizer st2(st.getToken(0), ":", true, true);
      if (st2.getToken(0) != "") 
        {
        node->setName(st2.getToken(0));
        }
      if (st2.numberOfRemainingTokens () >= 1) 
        {
        if (st2.numberOfRemainingTokens ()>1 )
          node->setDistanceToFather(TextTools::toDouble(st2.getToken(1)));       
      }
      hasId = setNodeProperties(*node, st.getToken(1));
      }
    }
  TreeTemplate<Node>* tree = new TreeTemplate<Node>();
  tree->setRootNode(node);
  if (!hasId)
    {
    tree->resetNodesId();
    }
  return tree;
}

/******************************************************************************/

string NHX::propertiesToParenthesis(const Node& node) const
{
  ostringstream s;
  s << "[&&NHX";
  vector <string> bProps = node.getBranchPropertyNames();
  for (unsigned int i= 0 ; i <bProps.size(); i++)
    {
    if (bProps[i]==TreeTools::BOOTSTRAP)
      {
      s << ":"<<"B"<<"="<< (dynamic_cast<const Number<double> *>(node.getBranchProperty(bProps[i]))->getValue());
      }
    else
      {
      s << ":"<<bProps[i]<<"="<< *(dynamic_cast<const BppString*>(node.getBranchProperty(bProps[i])));
      }
    }
  vector <string> nProps = node.getNodePropertyNames();
  for (unsigned int i= 0 ; i <nProps.size(); i++)
    {
    s << ":"<<nProps[i]<<"="<<*(dynamic_cast<const BppString*>(node.getNodeProperty(nProps[i])));
    }
  if (!node.hasNodeProperty ("ND"))
    {
    s << ":ND="<<TextTools::toString(node.getId());
    }
  s << "]";
  return s.str();  
}

/******************************************************************************/

string NHX::nodeToParenthesis(const Node& node) const
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
    for (unsigned int i = 1; i < node.getNumberOfSons(); i++)
      {
      s << "," << nodeToParenthesis(* node[i]);
      }
    s << ")";
        }
  if (node.hasDistanceToFather()) s << ":" << node.getDistanceToFather();
  s << propertiesToParenthesis(node);
  return s.str();  
}

/******************************************************************************/

string NHX::treeToParenthesis(const TreeTemplate<Node>& tree) const
{
  ostringstream s;
  s << "(";

  const Node* node = tree.getRootNode();

  if (node->isLeaf())
    {

    s << node->getName();

    for (unsigned int i = 0; i < node->getNumberOfSons(); ++i)
      {
      s << "," << nodeToParenthesis(*node->getSon(i));
      }
    }
  else
    {
    s << nodeToParenthesis(* node->getSon(0));
    for(unsigned int i = 1; i < node->getNumberOfSons(); ++i)
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

bool NHX::setNodeProperties(Node& node, const string properties) const
{
  string prop = TextTools::removeChar(properties, ']');
  StringTokenizer st(prop, ":", true, true);
  vector <string> props ; 
  bool hasId = false;
  while (st.hasMoreToken())
    {
    props.push_back(st.nextToken());
    }
  //This is a node:

  for(unsigned int i = 0; i < props.size(); i++)
    {
    if (TextTools::hasSubstring(st.getToken(i), "="))
      {
      
      StringTokenizer pt(st.getToken(i), "=", true, true);
      
      BppString property (pt.getToken(1));
      
      //Set branch properties, like support values ("B"), presence or absence of a duplication ("D"), of an event ("Ev"), 
      //width of parent branch ("W"), color of parent branch ("C"), custom data ("XB")
      if (pt.getToken(0)=="D" || pt.getToken(0)=="Ev" || pt.getToken(0)=="W" || pt.getToken(0)=="C" || pt.getToken(0)=="XB")
        {
        node.setBranchProperty (pt.getToken(0), property);        
        }
      else if (pt.getToken(0)=="B" )
        {        
          node.setBranchProperty (TreeTools::BOOTSTRAP, Number<double> (TextTools::toDouble(pt.getToken(1))));        
        }
      else if ((pt.getToken(0)=="ND") )//&& ((double) ((int) property.getValue())) == property.getValue())
        {        
          if (TextTools::isDecimalNumber(property.toSTL()))
            {          
              node.setId (TextTools::toInt(property.toSTL()));          
              hasId = true;   
            }
          else 
            {          
              node.setNodeProperty ("ND", property);           
            }
        }
      else {
        node.setNodeProperty (pt.getToken(0), property);   
      }
      }
    }
  return hasId;
}

/******************************************************************************/


