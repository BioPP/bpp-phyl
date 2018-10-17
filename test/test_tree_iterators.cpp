// Test for TreeIterator implementation

// From the STL:
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

// From bpp-core:
#include <Bpp/Text/TextTools.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/AutoParameter.h>


// From bpp-phyl
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/TreeIterator.h>

using namespace bpp;

void testPostOderIterator()
{
        // binary tree
        cout << "Test PostOrderIterator" << endl;
        cout << "Test on simple binary tree" << endl;
        TreeTemplate<Node>* binPostOrderTree = TreeTemplateTools::parenthesisToTree("((A:1,B:1):1,D:1);");
        vector <Node*> nodes = binPostOrderTree->getNodes();
        for (size_t i=0; i < nodes.size(); ++i)
        {
            Node* node = nodes[i];
            if (!node->isLeaf() && !node->hasFather())      // node is root
            {
                node->setName("E"); 
            }
            else if (!node->isLeaf() && node->hasFather())  // node is internal but not the root
            {
                node->setName("C");
            }
        }
        cout << "Visiting nodes of ((A,B)C,D)E; in post-order -> printed nodes names should match lexicographic order" << endl;
        PostOrderTreeIterator* binPostOrderTreeIt = new PostOrderTreeIterator(*binPostOrderTree);
        for (Node* node = binPostOrderTreeIt->begin(); node != binPostOrderTreeIt->end(); node = binPostOrderTreeIt->next())
        {
             cout << node->getName() << endl;
        }
        cout << "\n" << endl;
        delete(binPostOrderTreeIt);
        delete(binPostOrderTree);

        // non binary tree
        cout << "Test on non-binary tree" << endl;
        TreeTemplate<Node>* nonBinPostOrderTree = TreeTemplateTools::parenthesisToTree("(((A:1):1,(C:1,D:1,E:1):1):1,H:1);");
        nodes = nonBinPostOrderTree->getNodes();
        for (size_t i=0; i < nodes.size(); ++i) // who is node 7?
        {
            Node* node = nodes[i];
            if (!node->isLeaf() && !node->hasFather())      // node is root
            {
                node->setName("I"); 
            }
            else if (!node->isLeaf() && node->hasFather())  // node is internal but not the root
            {
                if (node->getSon(0)->getName() == "A") {
                    node->setName("B");
                } else if (node->getSon(0)->getName() == "C") {
                    node->setName("F");
                } else if (node->getSon(0)->getName() == "B") {
                    node->setName("G"); // this is the 3rd node in nodes -which indicates it doesn't match the post order
                }
                
            }
        }

        cout << "Visiting nodes of ((A)B,(C,D,E)F)G,H)I; in post-order -> printed nodes names should match lexicographic order" << endl;
        PostOrderTreeIterator* nonBinPostOrderTreeIt = new PostOrderTreeIterator(*nonBinPostOrderTree);
        for (Node* node = nonBinPostOrderTreeIt->begin(); node != nonBinPostOrderTreeIt->end(); node = nonBinPostOrderTreeIt->next())
        {
             cout << node->getName() << endl;
        }
        cout << "\n" << endl;
        delete(nonBinPostOrderTreeIt);
        delete(nonBinPostOrderTree);      
}

void testPreOrderIterator()
{
        // binary tree
        cout << "Test PreOrderIterator" << endl;
        cout << "Test on simple binary tree" << endl;
        TreeTemplate<Node>* binPreOrderTree  = TreeTemplateTools::parenthesisToTree("((C:1,D:1):1,E:1);");
        vector <Node*> nodes = binPreOrderTree->getNodes();
        for (size_t i=0; i < nodes.size(); ++i)
        {
            Node* node = nodes[i];
            if (!node->isLeaf() && !node->hasFather()) // node is root
            {
                node->setName("A"); 
            }
            else if (!node->isLeaf() && node->hasFather()) // node is internal but not the root
            {
                node->setName("B");
            }
        }
        cout << "Visiting nodes of ((C,D)B,E)A; in pre-order -> printed nodes names should match lexicographic order" << endl;
        PreOrderTreeIterator* binPreOrderTreeIt = new PreOrderTreeIterator(*binPreOrderTree);
        for (Node* node = binPreOrderTreeIt->begin(); node != binPreOrderTreeIt->end(); node = binPreOrderTreeIt->next())
        {
             cout << node->getName() << endl;
        }
        cout << "\n" << endl;
        delete(binPreOrderTreeIt);
        delete(binPreOrderTree);

        // non binary tree
        cout << "Test on non-binary tree" << endl;
        TreeTemplate<Node>* nonBinPreOrderTree  = TreeTemplateTools::parenthesisToTree("(((D:1):1,(F:1,G:1,H:1):1):1,I:1);");
        nodes = nonBinPreOrderTree->getNodes();
        for (size_t i=0; i < nodes.size(); ++i)
        {
            Node* node = nodes[i];
            if (!node->isLeaf() && !node->hasFather()) {
                node->setName("A"); 
            } else if (!node->isLeaf() && node->getSon(0)->getName() == "D") {
                node->setName("C");
            } else if (!node->isLeaf() && node->getSon(0)->getName() == "F") {
                node->setName("E");
            } else if (!node->isLeaf() && node->getSon(0)->getName() == "C") {
                node->setName("B");
            } 
        }

        cout << "Visiting nodes of (((D)C,(F,G,H)E)B,I)A; in pre-order -> printed nodes names should match lexicographic order" << endl;
        PreOrderTreeIterator* nonBinPreOrderTreeIt = new PreOrderTreeIterator(*nonBinPreOrderTree);
        for (Node* node = nonBinPreOrderTreeIt->begin(); node != nonBinPreOrderTreeIt->end(); node = nonBinPreOrderTreeIt->next())
        {
             cout << node->getName() << endl;
        }
        cout << "\n" << endl;
        delete(nonBinPreOrderTreeIt);
        delete(nonBinPreOrderTree);
}

void testInOderIterator()
{
        // binary tree
        cout << "Test InOrderIterator" << endl;
        cout << "Test on simple binary tree" << endl;
        TreeTemplate<Node>* binInOrderTree   = TreeTemplateTools::parenthesisToTree("((A:1,C:1):1,E:1);");
        vector <Node*> nodes = binInOrderTree->getNodes();
        for (size_t i=0; i < nodes.size(); ++i)
        {
            Node* node = nodes[i];
            if (!node->isLeaf() && !node->hasFather()) // node is root
            {
                node->setName("D"); 
            }
            else if (!node->isLeaf() && node->hasFather()) // node is internal but not the root
            {
                node->setName("B");
            }
        }
        cout << "Visiting nodes of ((A,C)B,E)D; in in-order -> printed nodes names should match lexicographic order" << endl;
        InOrderTreeIterator* binInOrderTreeIt = new InOrderTreeIterator(*binInOrderTree);
        for (Node* node = binInOrderTreeIt->begin(); node != binInOrderTreeIt->end(); node = binInOrderTreeIt->next())
        {
             cout << node->getName() << endl;
        }
        cout << "\n" << endl;
        delete(binInOrderTreeIt);
        delete(binInOrderTree);

        // non binary tree
        cout << "Test on non-binary tree" << endl;
        TreeTemplate<Node>* nonBinInOrderTree   = TreeTemplateTools::parenthesisToTree("(((A:1):1,(D:1,F:1,G:1):1):1,(I:1,J:1,L:1,M:1):1);");
        nodes = nonBinInOrderTree->getNodes();
        for (size_t i=0; i < nodes.size(); ++i)
        {
            Node* node = nodes[i];
            if (!node->isLeaf() && !node->hasFather()) {
                node->setName("H"); 
            } else if (!node->isLeaf() && node->getSon(0)->getName() == "A") {
                node->setName("B");
            } else if (!node->isLeaf() && node->getSon(0)->getName() == "D") {
                node->setName("E");
            } else if (!node->isLeaf() && node->getSon(0)->getName() == "B") {
                node->setName("C");
            } else if (!node->isLeaf() && node->getSon(0)->getName() == "I") {
                node->setName("K");
            } 
        }
        cout << "Visiting nodes of (((A)B,(D,F,G)E)C,(I,J,L,M)K)H; in in-order -> printed nodes names should match lexicographic order" << endl;
        InOrderTreeIterator* nonBinInOrderTreeIt = new InOrderTreeIterator(*nonBinInOrderTree);
        for (Node* node = nonBinInOrderTreeIt->begin(); node != nonBinInOrderTreeIt->end(); node = nonBinInOrderTreeIt->next())
        {
             cout << node->getName() << endl;
        }
        cout << "\n" << endl;
        delete(nonBinInOrderTreeIt);
        delete(nonBinInOrderTree);
}

int main() 
{
    try
    {
        // Bio++ doesn't handle internal node names when reading the tree, so first I need to set the internal nodes names in each case

        testPostOderIterator();

        testPreOrderIterator();

        testInOderIterator();
        
    }
    catch (exception & e)
    {
        cout << e.what() << endl;
        return 1;
    }
    return 0;
}
