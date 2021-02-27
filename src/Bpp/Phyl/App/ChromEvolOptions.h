#ifndef _CHROMEVOLOPTIONS_H_
#define _CHROMEVOLOPTIONS_H_

// From bpp-core:

#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>


// From bpp-seq:
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From bpp-phyl:
#include <Bpp/Phyl/Model/ChromosomeSubstitutionModel.h>
#include "PhylogeneticsApplicationTools.h"



//standard libraries
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
namespace bpp{

class ChromEvolOptions
{

public:
    static void initAllParameters(BppApplication& ChromEvol);
    static void initVectorOfChrNumParameters(vector<double>& paramVector);
    virtual ~ChromEvolOptions(){};
    
public:
    static string treeFilePath_;
    static string characterFilePath_;
    static int maxChrNum_;
    static int minChrNum_;
    static double branchMul_;
    static std::vector <unsigned int> OptPointsNum_;
    static std::vector <unsigned int> OptIterNum_;
    static double constGain_;
    static double constLoss_;
    static double constDupl_;
    static double constDemiDupl_;
    static double gainR_;
    static double lossR_;
    static int baseNum_;
    static double baseNumR_;
    static double duplR_;
    static double tolerance_;
    static unsigned int maxIterations_;
    static bool maxParsimonyBound_;
    static bool standardOptimization_;
    static int BrentBracketing_;
    static string optimizationMethod_;
    static unsigned int maxAlpha_;
    static unsigned int minAlpha_;
    static int seed_;
    static std::vector <double> probsForMixedOptimization_;
    static string rootFreqs_;
    static string fixedFrequenciesFilePath_;
    static ChromosomeSubstitutionModel::rateChangeFunc rateChangeType_;
    //static bool optimizeBaseNumber_;
    static string baseNumOptimizationMethod_;
    static std::vector<unsigned int> fixedParams_; //1 if parameter should be fixed. The order corresponds to the one in the model definition.
    static int NumOfSimulations_;
    static int jumpTypeMethod_;
    static bool simulateData_;
    static int numOfDataToSimulate_;
    static string resultsPathDir_;
    static int maxBaseNumTransition_; // needed for the simulator, since there is no data to infer it!
    static double treeLength_;
    static int maxNumOfTrials_; // to test the severity of the underflow problems
    static uint transitionMatFactor_;

private:
    static void initDefaultParameters();
    static void initParametersFromFile(BppApplication& ChromEvol);
    static void setFixedParams(std::vector<unsigned int> fixedParams);

};

}


#endif  // _CHROMEVOLOPTIONS_H_