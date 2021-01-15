#include <Bpp/Seq/Alphabet/ChromosomeAlphabet.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Io/AbstractISequence.h>
#include <Bpp/Seq/Io/ISequence.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/TreeTemplateTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/ChromosomeSubstitutionModel.h>
#include <Bpp/Seq/Container/VectorSequenceContainer.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <string>
#include <vector>
#include <iostream>

using namespace bpp;
using namespace std;
//unsigned int _MAX_CHR_NUM= 15;


VectorSiteContainer* resizeAlphabetForSequenceContainer(VectorSequenceContainer* vs, unsigned int max_chr_number);
VectorSiteContainer* getSequenceData(const std :: string &path, ChromosomeAlphabet* alpha, unsigned int global_max, unsigned int* maxObservedChrNum, unsigned int* minObservedChrNum);
void rescale_tree(TreeTemplate<Node>* tree, double scale_tree_factor, unsigned int ChrRange);
bool checkPijtMatrix(RowMatrix<double>& mat);

/******************************************************************************/
VectorSiteContainer* getSequenceData(const std :: string &path, ChromosomeAlphabet* alpha, unsigned int global_max, unsigned int* maxObservedChrNum, unsigned int* minObservedChrNum){
    Fasta fasta;
    VectorSequenceContainer* initial_set_of_sequences = fasta.readSequences(path, alpha);
    size_t number_of_sequences = initial_set_of_sequences->getNumberOfSequences();
    vector <string> sequence_names = initial_set_of_sequences->getSequencesNames();
    //unsigned int max_number_of_chr = 1; //the minimal number of chromosomes cannot be zero
    for (size_t i = 0; i < number_of_sequences; i++){
    //Sequence* seq_pointer = new BasicSequence(alpha);
        BasicSequence seq = initial_set_of_sequences->getSequence(sequence_names[i]);
        int sequence_content = seq.getValue(0);
        if ((sequence_content == -1) | (sequence_content == static_cast<int>(global_max)+1)){
            continue;
        }
        if ((unsigned int) sequence_content > *maxObservedChrNum){
            *maxObservedChrNum = sequence_content;
        }
        if ((unsigned int) sequence_content < *minObservedChrNum){
            *minObservedChrNum = sequence_content;
        }
    }
    //*maxObservedChrNum += 10;
    VectorSiteContainer* vsc = resizeAlphabetForSequenceContainer(initial_set_of_sequences, *maxObservedChrNum);
    delete initial_set_of_sequences;
    return vsc;

}



/******************************************************************************/
VectorSiteContainer* resizeAlphabetForSequenceContainer(VectorSequenceContainer* vs, unsigned int max_chr_number){

  size_t number_of_sequences = vs->getNumberOfSequences();
  vector <string> sequence_names = vs->getSequencesNames();
  ChromosomeAlphabet* new_alphabet = new ChromosomeAlphabet(1,max_chr_number);
  VectorSiteContainer* resized_alphabet_site_container = new VectorSiteContainer(new_alphabet);
  for (size_t i = 0; i < number_of_sequences; i++){
    BasicSequence seq = vs->getSequence(sequence_names[i]);
    BasicSequence new_seq = BasicSequence(seq.getName(), seq.getChar(0), new_alphabet);
    resized_alphabet_site_container->addSequence(new_seq);

  }
  return resized_alphabet_site_container;

}
/******************************************************************************/

void rescale_tree(TreeTemplate<Node>* tree, double scale_tree_factor, unsigned int chrRange){
    string tree_str = TreeTemplateTools::treeToParenthesis(*tree);
    std :: cout << tree_str << endl;
    bool rooted = tree->isRooted();
    if (!rooted){
        throw UnrootedTreeException("The given input tree is unrooted. Tree must be rooted!", tree);
    }
    if (scale_tree_factor == 1.0){
        return;
    }else{
        //tree must be rescaled
        double treeLength = tree->getTotalLength();
        if (scale_tree_factor == 999){
            //size_t number_of_chr_states = alpha->getNumberOfStates()-1; //not including 'X'
            //the number of chromosome changes observed in the data is the minimum bound
            scale_tree_factor = (double)chrRange/treeLength;

        }
        tree->scaleTree(scale_tree_factor);


    }
}

/******************************************************************************/

bool checkPijtMatrix(RowMatrix <double>& mat){
    for (size_t i = 0; i < mat.getNumberOfRows(); i++){
        double sum_of_cols = 0.0;
        for (size_t j = 0; j < mat.getNumberOfColumns(); j++){
            sum_of_cols += mat(i,j);

        }
        cout <<"Matrix row sum is: "<< sum_of_cols <<endl;
        if (fabs(sum_of_cols - 1.0) > 0.0001){
            return false;
        }
    }
    return true;

}

int main(){
    string path_tree = "/home/anats/phd_project/tree.newick";
    //string path_tree = "/home/anats/phd_project/unrooted_tree.newick";
    string path_seq = "/home/anats/phd_project/example.fasta";
    unsigned int global_max = 25;
    unsigned int maxObservedChrNum = 1;
    unsigned int minObservedChrNum = 500;

    ChromosomeAlphabet* alpha_initial = new ChromosomeAlphabet(1,global_max);
    VectorSiteContainer* vsc = getSequenceData(path_seq, alpha_initial, global_max, &maxObservedChrNum, &minObservedChrNum);
    const Alphabet* alpha = vsc->getAlphabet();
    Alphabet* new_alpha = alpha->clone();
    ChromosomeAlphabet* chr_alpha = dynamic_cast <ChromosomeAlphabet*>(new_alpha);
    //ChromosomeAlphabet* alphc = new ChromosomeAlphabet(20);
    double scale_tree_factor = 999;
    Newick newick;
    TreeTemplate<Node>* tree = newick.readTree(path_tree);
    rescale_tree(tree, scale_tree_factor, maxObservedChrNum-minObservedChrNum);
    ChromosomeSubstitutionModel* chr_model = new ChromosomeSubstitutionModel(chr_alpha, 2, 1, 3, 1.3, IgnoreParam, IgnoreParam, IgnoreParam, IgnoreParam, IgnoreParam, maxObservedChrNum-minObservedChrNum, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, ChromosomeSubstitutionModel::LINEAR);
    RowMatrix<double> Qij;
    MatrixTools::copy(chr_model->getGenerator(), Qij);
    MatrixTools::print(Qij);

    RowMatrix <double> pijt;
    RowMatrix <double> pijt_2;
    RowMatrix <double> pijt_3;


    std::vector <double> branches = tree->getBranchLengths();
    for (int i = 0; i < (int)branches.size(); i++){
        double branchLength = branches[i];
        std::cout <<"branch length: "<< branchLength << endl;
        MatrixTools::copy(chr_model->getPij_t(branchLength), pijt);
        MatrixTools::print(pijt);
        MatrixTools::copy(chr_model->getPij_t_func2(branchLength), pijt_2);
        MatrixTools::print(pijt_2);   
        MatrixTools::copy(chr_model->getPij_t_func4(branchLength), pijt_3);  
        MatrixTools::print(pijt_3);
        bool converged = chr_model->checkIfReachedConvergence(pijt, pijt_2);
        cout << converged << endl;
        converged = chr_model->checkIfReachedConvergence(pijt, pijt_3);
        cout << converged << endl;
        converged = chr_model->checkIfReachedConvergence(pijt_2, pijt_3);
        cout << converged << endl;
        bool prob = checkPijtMatrix(pijt);
        cout << "check if probability matrix" <<endl;
        cout << prob<<"good!"<<endl;
    }
    

    



    //const RowMatrix <double> dPijt = chr_model->getdPij_dt(3.6);
    //MatrixTools::print(dPijt);
    //const RowMatrix <double> d2Pijt = chr_model->getd2Pij_dt2(3.6);
    //MatrixTools::print(d2Pijt);
    
    //delete tree;
    //delete chr_model;


}