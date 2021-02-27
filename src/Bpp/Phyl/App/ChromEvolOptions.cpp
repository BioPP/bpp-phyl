
#include "ChromEvolOptions.h"
using namespace bpp;
using namespace std;

string ChromEvolOptions::treeFilePath_;
string ChromEvolOptions::characterFilePath_;
int ChromEvolOptions::maxChrNum_;
int ChromEvolOptions::minChrNum_;
double ChromEvolOptions::branchMul_;
std::vector <unsigned int> ChromEvolOptions::OptPointsNum_;
std::vector <unsigned int> ChromEvolOptions::OptIterNum_;
double ChromEvolOptions::constGain_;
double ChromEvolOptions::constLoss_;
double ChromEvolOptions::constDupl_;
double ChromEvolOptions::constDemiDupl_;
double ChromEvolOptions::gainR_;
double ChromEvolOptions::lossR_;
int ChromEvolOptions::baseNum_;
double ChromEvolOptions::baseNumR_;
double ChromEvolOptions::duplR_;
double ChromEvolOptions::tolerance_;
unsigned int ChromEvolOptions::maxIterations_;
bool ChromEvolOptions::maxParsimonyBound_;
unsigned int ChromEvolOptions::maxAlpha_;
unsigned int ChromEvolOptions::minAlpha_;
int ChromEvolOptions::BrentBracketing_;
bool ChromEvolOptions::standardOptimization_;
string ChromEvolOptions::optimizationMethod_;
int ChromEvolOptions::seed_;
std::vector <double> ChromEvolOptions::probsForMixedOptimization_;
string ChromEvolOptions::rootFreqs_;
string ChromEvolOptions::fixedFrequenciesFilePath_;
ChromosomeSubstitutionModel::rateChangeFunc ChromEvolOptions::rateChangeType_;
//bool ChromEvolOptions::optimizeBaseNumber_;
string ChromEvolOptions::baseNumOptimizationMethod_;
std::vector<unsigned int> ChromEvolOptions::fixedParams_;
int ChromEvolOptions::NumOfSimulations_;
int ChromEvolOptions::jumpTypeMethod_;
bool ChromEvolOptions::simulateData_;
int ChromEvolOptions::numOfDataToSimulate_;
string ChromEvolOptions::resultsPathDir_;
int ChromEvolOptions::maxBaseNumTransition_;
double ChromEvolOptions::treeLength_;
int ChromEvolOptions::maxNumOfTrials_;
uint ChromEvolOptions::transitionMatFactor_;

/*************************************************************************/
void ChromEvolOptions::initAllParameters(BppApplication& ChromEvol){
    initDefaultParameters();
    initParametersFromFile(ChromEvol);

}
/*************************************************************************/
void ChromEvolOptions::initDefaultParameters(){
    maxAlpha_ = 500;
    minAlpha_ = 1;
    maxChrNum_ = -10;
    minChrNum_ = 1;
    maxIterations_ = 5;
    tolerance_ = 0.01;
    branchMul_ = 999;
    constGain_ = IgnoreParam;
    constLoss_ = IgnoreParam;
    constDupl_ = IgnoreParam;
    constDemiDupl_= IgnoreParam;
    gainR_ = IgnoreParam;
    lossR_ = IgnoreParam;
    baseNum_ = IgnoreParam;
    baseNumR_ = IgnoreParam;
    duplR_ = IgnoreParam;
    maxParsimonyBound_ = false;
    standardOptimization_ = false;
    BrentBracketing_ = 2;
    optimizationMethod_ = "Brent";
    seed_ = 0;
    rootFreqs_ = "weighted";
    rateChangeType_ = ChromosomeSubstitutionModel::LINEAR;
    //optimizeBaseNumber_ = false;
    baseNumOptimizationMethod_ = "Brent";
    NumOfSimulations_ = 10000;
    jumpTypeMethod_ = 0;
    simulateData_ = false;
    numOfDataToSimulate_ = 1;
    maxBaseNumTransition_ = 18;
    treeLength_ = 0;
    maxNumOfTrials_ = 100;
    transitionMatFactor_ = 1;



}
/*************************************************************************/
void ChromEvolOptions::setFixedParams(std::vector<unsigned int> fixedParams){
    
    //if (optimizeBaseNumber_){
    for (size_t i = 0; i < ChromosomeSubstitutionModel::NUM_OF_CHR_PARAMS; i++){
        switch (i)
        {
        case ChromosomeSubstitutionModel::BASENUM:
            if (baseNum_ != IgnoreParam){
                fixedParams_.push_back(fixedParams[ChromosomeSubstitutionModel::BASENUM]);      
            }
            break;
        case ChromosomeSubstitutionModel::BASENUMR:
            if (baseNumR_ != IgnoreParam){
                fixedParams_.push_back(fixedParams[ChromosomeSubstitutionModel::BASENUMR]);
            }
            break;
        case ChromosomeSubstitutionModel::DUPL:
            if (constDupl_ != IgnoreParam){
                fixedParams_.push_back(fixedParams[ChromosomeSubstitutionModel::DUPL]);
            }
            break;
        case ChromosomeSubstitutionModel::LOSS:
            if (constLoss_ != IgnoreParam){
                fixedParams_.push_back(fixedParams[ChromosomeSubstitutionModel::LOSS]);
            }
            break;
        case ChromosomeSubstitutionModel::GAIN:
            if (constGain_ != IgnoreParam){
                fixedParams_.push_back(fixedParams[ChromosomeSubstitutionModel::GAIN]);
            }
            break;
        case ChromosomeSubstitutionModel::LOSSR:
            if (lossR_ != IgnoreParam){
                fixedParams_.push_back(fixedParams[ChromosomeSubstitutionModel::LOSSR]);
            }
            break;
        case ChromosomeSubstitutionModel::GAINR:
            if (gainR_ != IgnoreParam){
                fixedParams_.push_back(fixedParams[ChromosomeSubstitutionModel::GAINR]);
            }
            break;
        case ChromosomeSubstitutionModel::DUPLR:
            if (duplR_ != IgnoreParam){
                fixedParams_.push_back(fixedParams[ChromosomeSubstitutionModel::DUPLR]);
            }
            break;
        case ChromosomeSubstitutionModel::DEMIDUPL:
            if ((constDemiDupl_ != IgnoreParam) && (constDemiDupl_ != DemiEqualDupl)){
                fixedParams_.push_back(fixedParams[ChromosomeSubstitutionModel::DEMIDUPL]);
            }
            break;
       
        default:
            throw Exception("ChromEvolOptions::setFixedParams(): Invalid rate type!");
            break;
        }

    }
}
/*************************************************************************/
void ChromEvolOptions::initParametersFromFile(BppApplication& ChromEvol){
    maxChrNum_ = ApplicationTools::getIntParameter("_maxChrNum", ChromEvol.getParams(), maxChrNum_, "", true, 0);
    minChrNum_ = ApplicationTools::getIntParameter("_minChrNum", ChromEvol.getParams(), minChrNum_, "", true, 0);
    seed_ = ApplicationTools::getIntParameter("_seed", ChromEvol.getParams(), seed_, "", true, 0);
    simulateData_ = ApplicationTools::getBooleanParameter("_simulateData", ChromEvol.getParams(), simulateData_, "", true, 0);
    if (simulateData_){
        characterFilePath_ = ApplicationTools::getAFilePath("_dataFile", ChromEvol.getParams(), false, true, "", true, "none", 1);
    }else{
        characterFilePath_ = ApplicationTools::getAFilePath("_dataFile", ChromEvol.getParams(), true, true, "", true, "none", 1);
    }
    treeFilePath_ = ApplicationTools::getAFilePath("_treeFile", ChromEvol.getParams(), true, true, "", true, "none", 1);
    branchMul_ = ApplicationTools::getDoubleParameter("_branchMul", ChromEvol.getParams(), branchMul_, "", true, 0);
    maxIterations_ = (unsigned int)ApplicationTools::getIntParameter("_maxOptimizationItarations", ChromEvol.getParams(), maxIterations_, "", true, 0);
    tolerance_ = ApplicationTools::getDoubleParameter("_tolParamOptimization", ChromEvol.getParams(), tolerance_, "", true, 0);
    constGain_ = ApplicationTools::getDoubleParameter("_gainConstR", ChromEvol.getParams(), constGain_, "", true, 0);
    constLoss_ = ApplicationTools::getDoubleParameter("_lossConstR", ChromEvol.getParams(), constLoss_, "", true, 0);
    constDupl_ = ApplicationTools::getDoubleParameter("_duplConstR", ChromEvol.getParams(), constDupl_, "", true, 0);
    constDemiDupl_ = ApplicationTools::getDoubleParameter("_demiPloidyR", ChromEvol.getParams(), constDemiDupl_, "", true, 0);
    gainR_ = ApplicationTools::getDoubleParameter("_gainR", ChromEvol.getParams(), gainR_, "", true, 0);
    lossR_ = ApplicationTools::getDoubleParameter("_lossR", ChromEvol.getParams(), lossR_, "", true, 0);
    baseNum_ = ApplicationTools::getIntParameter("_baseNum", ChromEvol.getParams(), baseNum_, "", true, 0);
    baseNumR_ = ApplicationTools::getDoubleParameter("_baseNumR", ChromEvol.getParams(), baseNumR_, "", true, 0);
    duplR_ = ApplicationTools::getDoubleParameter("_duplR", ChromEvol.getParams(), duplR_, "", true, 0);
    maxParsimonyBound_ = ApplicationTools::getBooleanParameter("_maxParsimonyBound", ChromEvol.getParams(), maxParsimonyBound_, "", true, 0);
    standardOptimization_ = ApplicationTools::getBooleanParameter("_standardOptimization", ChromEvol.getParams(), standardOptimization_, "", true, 0);
    BrentBracketing_ = ApplicationTools::getIntParameter("_BrentBracketing", ChromEvol.getParams(), BrentBracketing_, "", true, 0);
    optimizationMethod_ = ApplicationTools::getStringParameter("_optimizationMethod", ChromEvol.getParams(), optimizationMethod_, "", true, 0);
    string defaultValForOptPointsNum = "10,3,1";
    string defaultValForOptIterNum = "0,2,5";
    string defaultValForProbsForMixedOpt = "1,0";
    OptPointsNum_ = ApplicationTools::getVectorParameter<unsigned int>("_optimizePointsNum", ChromEvol.getParams(), ',', defaultValForOptPointsNum, "", true, 0);
    OptIterNum_ = ApplicationTools::getVectorParameter<unsigned int>("_optimizeIterNum", ChromEvol.getParams(), ',', defaultValForOptIterNum, "", true, 0);
    probsForMixedOptimization_ = ApplicationTools::getVectorParameter<double>("_probsForMixedOptimization", ChromEvol.getParams(), ',', defaultValForProbsForMixedOpt, "", true, 0);
    fixedFrequenciesFilePath_ = ApplicationTools::getAFilePath("_fixedFrequenciesFilePath", ChromEvol.getParams(), false, true, "", true, "none", 0);
    rootFreqs_ = ApplicationTools::getStringParameter("_rootFreqs", ChromEvol.getParams(), rootFreqs_, "", true, 0);
    int rateChangeType = ApplicationTools::getIntParameter("_rateChangeType", ChromEvol.getParams(), rateChangeType_, "", true, 0);
    if (rateChangeType){
        rateChangeType_ = ChromosomeSubstitutionModel::EXP;
    }else{
        rateChangeType_ = ChromosomeSubstitutionModel::LINEAR;
    }
    //optimizeBaseNumber_ = ApplicationTools::getBooleanParameter("_optimizeBaseNumber", ChromEvol.getParams(), optimizeBaseNumber_, "", true, 0);
    baseNumOptimizationMethod_ = ApplicationTools::getStringParameter("_baseNumOptimizationMethod", ChromEvol.getParams(), baseNumOptimizationMethod_, "", true, 0);
    string defaultValForFixedParams = "0,0,0,0,0,0,0,0,0";
    std::vector<unsigned int> fixedParams = ApplicationTools::getVectorParameter<unsigned int>("_fixedParams", ChromEvol.getParams(), ',', defaultValForFixedParams, "", true, 0);
    setFixedParams(fixedParams);
    NumOfSimulations_ = ApplicationTools::getIntParameter("_NumOfSimulations", ChromEvol.getParams(), NumOfSimulations_, "", true, 0);
    jumpTypeMethod_ = ApplicationTools::getIntParameter("_jumpTypeMethod", ChromEvol.getParams(), jumpTypeMethod_, "", true, 0);
    numOfDataToSimulate_ = ApplicationTools::getIntParameter("_numOfDataToSimulate", ChromEvol.getParams(), numOfDataToSimulate_, "", true, 0);
    resultsPathDir_ = ApplicationTools::getAFilePath("_resultsPathDir", ChromEvol.getParams(), false, true, "", true, "none", 0);
    maxBaseNumTransition_ = ApplicationTools::getIntParameter("_maxBaseNumTransition", ChromEvol.getParams(), maxBaseNumTransition_, "", true, 0);
    treeLength_ = ApplicationTools::getDoubleParameter("_treeLength", ChromEvol.getParams(), treeLength_, "", true, 0);
    maxNumOfTrials_ = ApplicationTools::getIntParameter("_maxNumOfTrials", ChromEvol.getParams(), maxNumOfTrials_, "", true, 0);
    transitionMatFactor_ = (uint)(ApplicationTools::getIntParameter("_transitionMatFactor", ChromEvol.getParams(), transitionMatFactor_, "", true, 0));

}
/*************************************************************************/
void ChromEvolOptions::initVectorOfChrNumParameters(vector<double>& paramVector){
    for (size_t i = 0; i < ChromosomeSubstitutionModel::NUM_OF_CHR_PARAMS; i++){
        switch(i){
            case ChromosomeSubstitutionModel::BASENUM:
                paramVector.push_back(baseNum_);
                break;
            case ChromosomeSubstitutionModel::BASENUMR:
                paramVector.push_back(baseNumR_);
                break;
            case ChromosomeSubstitutionModel::DUPL:
                paramVector.push_back(constDupl_);
                break;
            case ChromosomeSubstitutionModel::LOSS:
                paramVector.push_back(constLoss_);
                break;
            case ChromosomeSubstitutionModel::GAIN:
                paramVector.push_back(constGain_);
                break;
            case ChromosomeSubstitutionModel::DEMIDUPL:
                paramVector.push_back(constDemiDupl_);
                break;
            case ChromosomeSubstitutionModel::LOSSR:
                paramVector.push_back(lossR_);
                break;
            case ChromosomeSubstitutionModel::GAINR:
                paramVector.push_back(gainR_);
                break;
            case ChromosomeSubstitutionModel::DUPLR:
                paramVector.push_back(duplR_);
                break;
            default:
                throw Exception("ChromEvolOptions::initVectorOfChrNumParameters(): Invalid rate type!");
                break;
        }

    }   

}
