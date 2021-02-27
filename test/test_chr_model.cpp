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


//VectorSiteContainer* resizeAlphabetForSequenceContainer(VectorSequenceContainer* vs, unsigned int max_chr_number);
//VectorSiteContainer* getSequenceData(const std :: string &path, ChromosomeAlphabet* alpha, unsigned int global_max, unsigned int* maxObservedChrNum, unsigned int* minObservedChrNum);
//void rescale_tree(TreeTemplate<Node>* tree, double scale_tree_factor, unsigned int ChrRange);
bool checkPijtMatrix(RowMatrix<double>& mat);

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


    //const NucleicAlphabet* alphabet = &AlphabetTools::DNA_ALPHABET;

   // VectorSiteContainer sites(alphabet);
    //sites.addSequence(BasicSequence("A", "ATCCAGACATGCCGGGACTTTGCAGAGAAGGAGTTGTTTCCCATTGCAGCCCAGGTGGATAAGGAACAGC", alphabet));
    //sites.addSequence(BasicSequence("B", "CGTCAGACATGCCGTGACTTTGCCGAGAAGGAGTTGGTCCCCATTGCGGCCCAGCTGGACAGGGAGCATC", alphabet));
    //sites.addSequence(BasicSequence("C", "GGTCAGACATGCCGGGAATTTGCTGAAAAGGAGCTGGTTCCCATTGCAGCCCAGGTAGACAAGGAGCATC", alphabet));
    //sites.addSequence(BasicSequence("D", "TTCCAGACATGCCGGGACTTTACCGAGAAGGAGTTGTTTTCCATTGCAGCCCAGGTGGATAAGGAACATC", alphabet));

    const ChromosomeAlphabet* alphabet = new ChromosomeAlphabet(1,5);
    ChromosomeSubstitutionModel* model1 = new ChromosomeSubstitutionModel(alphabet, 2, 1, 3, 1.3, 0.0123, 0.35, IgnoreParam, IgnoreParam, IgnoreParam, 4, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, ChromosomeSubstitutionModel::LINEAR);
    ChromosomeSubstitutionModel* model2 = new ChromosomeSubstitutionModel(alphabet, 2, 1, 3, 1.3, 0.0123, 0.35, IgnoreParam, IgnoreParam, IgnoreParam, 4, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, ChromosomeSubstitutionModel::LINEAR, true);
    RowMatrix<double> Qij;
    MatrixTools::copy(model1->getGenerator(), Qij);
    MatrixTools::print(Qij);
    std::vector <double> t;
    t.push_back(0.05829885627);
    t.push_back(0.6344500368);
    t.push_back(0.6344500368);
    t.push_back(0.02345);
    t.push_back(0.001);

    RowMatrix <double> pijt;
    RowMatrix <double> pijdt;
    RowMatrix <double> pijtdt2;

    MatrixXef pijt_ef;
    MatrixXef pijdt_ef;
    MatrixXef pijtdt2_ef;

    for (size_t i = 0; i < t.size(); i++){
        std::cout << "********************************" << endl;
        std::cout << "********************************" << endl;
        std::cout << "Starting with t = " << t[i] << endl;
        pijt = model1->getPijt_test(t[i]);
        pijt_ef = model2->getPijtEf_test(t[i]);
        pijdt = model1->getdPij_dt(t[i]);
        pijdt_ef = model2->getdPijEf_dt(t[i]);
        pijtdt2 = model1->getd2Pij_dt2(t[i]);
        pijtdt2_ef = model2->getd2PijEf_dt2(t[i]);

        Eigen::MatrixXd pijt_db = ExtendedFloatTools::convertBppMatrixToEigenDB(pijt);
        Eigen::MatrixXd pijt_ef_db = ExtendedFloatTools::convertMatToDouble(pijt_ef);

        Eigen::MatrixXd pijdt_db = ExtendedFloatTools::convertBppMatrixToEigenDB(pijdt);
        Eigen::MatrixXd pijdt_ef_db = ExtendedFloatTools::convertMatToDouble(pijdt_ef);

        Eigen::MatrixXd pijtdt2_db = ExtendedFloatTools::convertBppMatrixToEigenDB(pijtdt2);
        Eigen::MatrixXd pijtdt2_ef_db = ExtendedFloatTools::convertMatToDouble(pijtdt2_ef);


        std::cout << "***" << endl;
        std::cout << "Pijt double:" << endl;
        std::cout << pijt_db << endl;
        std::cout << "***" << endl;
        std::cout << "Pijt Extended float:" << endl;
        std::cout << pijt_ef_db << endl;
        std::cout << "******************************" << endl;
        std::cout << "Pijdt double:" << endl;
        std::cout << pijdt_db << endl;
        std::cout << "***" << endl;
        std::cout << "Pijdt Extended float:" << endl;
        std::cout << pijdt_ef_db << endl;
        std::cout << "******************************" << endl;
        std::cout << "pijtd2t2 double:" << endl;
        std::cout << pijtdt2_db << endl;
        std::cout << "***" << endl;
        std::cout << "Pijd2t2 Extended float:" << endl;
        std::cout << pijtdt2_ef_db << endl;
        std::cout << "******************************" << endl;

    }
 
    delete model1;
    delete model2;
    delete alphabet;


    // std::vector <double> branches = tree->getBranchLengths();
    // for (int i = 0; i < (int)branches.size(); i++){
    //     double branchLength = branches[i];
    //     std::cout <<"branch length: "<< branchLength << endl;
    //     MatrixTools::copy(chr_model->getPij_t(branchLength), pijt);
    //     MatrixTools::print(pijt);
    //     MatrixTools::copy(chr_model->getPij_t_func2(branchLength), pijt_2);
    //     MatrixTools::print(pijt_2);   
    //     MatrixTools::copy(chr_model->getPij_t_func4(branchLength), pijt_3);  
    //     MatrixTools::print(pijt_3);
    //     bool converged = chr_model->checkIfReachedConvergence(pijt, pijt_2);
    //     cout << converged << endl;
    //     converged = chr_model->checkIfReachedConvergence(pijt, pijt_3);
    //     cout << converged << endl;
    //     converged = chr_model->checkIfReachedConvergence(pijt_2, pijt_3);
    //     cout << converged << endl;
    //     bool prob = checkPijtMatrix(pijt);
    //     cout << "check if probability matrix" <<endl;
    //     cout << prob<<"good!"<<endl;
    //}
    

    



    //const RowMatrix <double> dPijt = chr_model->getdPij_dt(3.6);
    //MatrixTools::print(dPijt);
    //const RowMatrix <double> d2Pijt = chr_model->getd2Pij_dt2(3.6);
    //MatrixTools::print(d2Pijt);
    
    //delete tree;
    //delete chr_model;


}