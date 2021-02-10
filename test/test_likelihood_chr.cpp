//from bpp-core
#include <Bpp/Version.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/App/BppApplication.h>
//#include <Bpp/Numeric/AutoParameter.h>
//#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
//#include <Bpp/Numeric/Matrix/MatrixTools.h>
//#include <Bpp/Numeric/Random/RandomTools.h>
//#include <Bpp/Text/StringTokenizer.h>



//from bpp-seq
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Alphabet/ChromosomeAlphabet.h>
//#include <Bpp/Seq/Container/VectorSequenceContainer.h>
//#include <Bpp/Seq/Container/VectorSiteContainer.h>
//#include <Bpp/Seq/Container/SiteContainerTools.h>
//#include <Bpp/Seq/Io/AbstractISequence.h>
//#include <Bpp/Seq/Io/ISequence.h>
//#include <Bpp/Seq/Io/chrFasta.h>
//#include <Bpp/Seq/SiteTools.h>
//#include <Bpp/Seq/App/SequenceApplicationTools.h>

//from bpp-phyl
//#include <Bpp/Phyl/TreeTemplate.h>
//#include <Bpp/Phyl/TreeTemplateTools.h>
//#include <Bpp/Phyl/Io/Newick.h>
//#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
//#include <Bpp/Phyl/Likelihood/DRNonHomogeneousTreeLikelihood.h>
//#include <Bpp/Phyl/Likelihood/MLAncestralStateReconstruction.h>
//#include <Bpp/Phyl/Likelihood/MarginalNonRevAncestralStateReconstruction.h>
//#include <Bpp/Phyl/Likelihood/ChromosomeNumberOptimizer.h>
//#include <Bpp/Phyl/Mapping/ComputeChangesExpectations.h>
//#include <Bpp/Phyl/Mapping/ComputeChromosomeTransitionsExp.h>
//#include <Bpp/Phyl/Model/ChromosomeSubstitutionModel.h>
//#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
//#include <Bpp/Phyl/Model/SubstitutionModelSet.h>
//#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/App/ChromEvolOptions.h>
#include <Bpp/Phyl/App/ChromosomeNumberMng.h>


//standard libraries
#include <string>
#include <vector>
#include <iostream>
#include <time.h>

using namespace bpp;
using namespace std;




int main(int args, char **argv) {

    if (args == 1){
        std::cout << "No arguments provided"<<endl;
        return 0;
    }
    try{
        time_t t1;
        time(&t1);
        time_t t2;
        BppApplication ChromEvol(args, argv, "ChromEvol");
        ChromEvolOptions::initAllParameters(ChromEvol);
        ChromosomeNumberMng* mng = new ChromosomeNumberMng();
        if (!ChromEvolOptions::simulateData_){
            mng->getCharacterData(ChromEvolOptions::characterFilePath_);
        }
        mng->getTree(ChromEvolOptions::treeFilePath_, ChromEvolOptions::treeLength_);       
        //mng->runChromEvol();
        mng->runChromEvol();
        time(&t2);
        std::cout << "****** Max allowed chromosome number: "<< ChromEvolOptions::maxChrNum_ <<endl;
        std::cout <<"Total running time is: "<< static_cast<int>(t2-t1) <<endl;
        delete mng;

    }
    catch (exception& e)
    {
        cout << e.what() << endl;
        return 1;
    }

    return 0;
}
