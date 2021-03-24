#include "ChromosomeNumberOptimizer.h"
using namespace bpp;

void ChromosomeNumberOptimizer::initModels(vector<double> modelParams, double parsimonyBound, ChromosomeSubstitutionModel::rateChangeFunc rateChange, int seed, unsigned int numberOfModels, const string& fixedRootFreqPath, vector<unsigned int>& fixedParams){
    //optimizeBaseNumber_ = optimizeBaseNumber;
    fixedParams_ = fixedParams;
    map <int, double> setOfFixedParams;
    ChromosomeSubstitutionModel::getSetOfFixedParameters(modelParams, fixedParams_, setOfFixedParams);
    optimizeBaseNumber_ = setOfFixedParams.count(ChromosomeSubstitutionModel::BASENUM) == 0;
    vectorOfLikelohoods_.reserve(numberOfModels);
    vectorOfContexts_.reserve(numberOfModels);
    //vector <int> nodeIds = tree_->getNodesId();
    //nodeIds.pop_back();
    if (seed != 0){
        RandomTools::setSeed(static_cast<long>(seed));
    }
    for (size_t n = 0; n < numberOfModels; n++){
        DiscreteDistribution* rdist = new GammaDiscreteRateDistribution(1, 1.0);
        std::shared_ptr<ChromosomeSubstitutionModel> chrModel;
        if (n == 0){
            chrModel = std::make_shared<ChromosomeSubstitutionModel>(alphabet_, modelParams, baseNumberUpperBound_, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, rateChange);
            //chrModel = new ChromosomeSubstitutionModel(alphabet_, modelParams, baseNumberUpperBound_, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, rateChange);//initModel(alpha, chrRange);
        }else{
            //chrModel = initRandomModel(alpha, chrRange, parsimonyBound * (double)n);
            chrModel = std::shared_ptr<ChromosomeSubstitutionModel>(ChromosomeSubstitutionModel::initRandomModel(alphabet_, modelParams, baseNumberUpperBound_, ChromosomeSubstitutionModel::ROOT_LL, rateChange, fixedParams_, parsimonyBound * (double)n));
    
        }

        SingleProcessPhyloLikelihood* lik = getLikelihoodFunction(tree_, vsc_, chrModel, rdist, fixedRootFreqPath);
        
        if (std::isnan(lik->getValue())){
            std::cout << "value is nan"<<endl;
        }
        int countNumOfTrials = 0;
        
        while (((std::isinf(lik->getValue())) || (std::isnan(lik->getValue())))||(lik->getValue() < 0))
        {
            if (countNumOfTrials >= ChromEvolOptions::maxNumOfTrials_){
                break;
            }
            //deleteTreeLikAssociatedAttributes(lik);
            delete lik;
            vectorOfContexts_.pop_back();
            rdist = new GammaDiscreteRateDistribution(1, 1.0);
            chrModel = std::shared_ptr<ChromosomeSubstitutionModel>(ChromosomeSubstitutionModel::initRandomModel(alphabet_, modelParams, baseNumberUpperBound_, ChromosomeSubstitutionModel::ROOT_LL, rateChange, fixedParams_, parsimonyBound * (double)n));
            lik = getLikelihoodFunction(tree_, vsc_, chrModel, rdist, fixedRootFreqPath);          
            countNumOfTrials ++;

        }
            
        

        //initializing the likelihood instance
        //lik.initialize();
        vectorOfLikelohoods_.push_back(lik);//add to vector of likelihoods
        
    }

}
/****************************************************************************/
SingleProcessPhyloLikelihood* ChromosomeNumberOptimizer::getLikelihoodFunction(const PhyloTree* tree, const VectorSiteContainer* vsc, std::shared_ptr<ChromosomeSubstitutionModel> &chrModel, DiscreteDistribution* rdist, const string& fixedRootFreqPath){
    // bool calculateDerivatives = true;
    // if (ChromEvolOptions::optimizationMethod_ == "Brent"){
    //     calculateDerivatives  = false;
    // }
    
    bool weightedRootFreqs;
    std::shared_ptr<SubstitutionModel> model(static_pointer_cast<SubstitutionModel>(chrModel)->clone());
    //unsigned int nbStates = (unsigned int)(model->getNumberOfStates());
    //unsigned int nbEdges = (unsigned int)(tree_->getBranchLengths().size());
    unsigned int factor = ChromEvolOptions::transitionMatFactor_;
    NonHomogeneousSubstitutionProcess* subProSim;
    ParametrizablePhyloTree parTree(*tree_);
    if (fixedRootFreqPath != "none"){
        vector <double> rootFreqs = setFixedRootFrequencies(ChromEvolOptions::fixedFrequenciesFilePath_, chrModel);
        std::shared_ptr<FixedFrequencySet> rootFreqsFixed = std::make_shared<FixedFrequencySet>(std::shared_ptr<const StateMap>(new CanonicalStateMap(chrModel->getStateMap(), false)), rootFreqs);
        std::shared_ptr<FrequencySet> rootFrequencies = static_pointer_cast<FrequencySet>(rootFreqsFixed);
        subProSim= NonHomogeneousSubstitutionProcess::createHomogeneousSubstitutionProcess(model, rdist, parTree.clone(), shared_ptr<FrequencySet>(rootFrequencies->clone()));
        weightedRootFreqs = false;
        
    }else{
        subProSim= NonHomogeneousSubstitutionProcess::createHomogeneousSubstitutionProcess(model, rdist, parTree.clone());
        weightedRootFreqs = true;
        
    }
    
    SubstitutionProcess* nsubPro=subProSim->clone();
    Context context;
    vectorOfContexts_.push_back(context);

    auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(vectorOfContexts_[vectorOfContexts_.size()-1], *vsc_->clone(), *nsubPro, factor, weightedRootFreqs);
    lik->setFactor(factor);
    
    SingleProcessPhyloLikelihood* ntl = new SingleProcessPhyloLikelihood(vectorOfContexts_[vectorOfContexts_.size()-1], lik, lik->getParameters());
    delete subProSim;
    return ntl;
}
/****************************************************************************/
vector <double> ChromosomeNumberOptimizer::setFixedRootFrequencies(const std::string &path, std::shared_ptr<ChromosomeSubstitutionModel> chrModel){
    ifstream stream;
    stream.open(path.c_str());
    vector <double> freqs;
    vector <string> lines = FileTools::putStreamIntoVectorOfStrings(stream);
    stream.close();
    for (size_t i = 0; i < lines.size(); i++){
        string freq_i_str = TextTools::removeSurroundingWhiteSpaces(lines[i]);
        if (freq_i_str == ""){
            continue;
        }
        double freq_i = TextTools::toDouble(freq_i_str);
        if (static_cast<unsigned int>(freqs.size()) >= chrModel->getNumberOfStates()){
            if (freq_i > 0){
                throw Exception("Invalid fixed frequencies file!");
            }

        }else{
            freqs.push_back(freq_i);
        }
        
    }
    size_t nbStates = chrModel->getNumberOfStates();
    if (freqs.size() < nbStates){
        for (size_t s = freqs.size(); s < nbStates; s++){
            freqs.push_back(0);
        }
        
    }
    if (nbStates != freqs.size()){
        throw Exception("Invalid fixed frequencies file!");
    }
    return freqs;
}

/****************************************************************************/

void ChromosomeNumberOptimizer::optimize()
{

    unsigned int totalNumOfEvaluations = 0;
    unsigned int numOfEvaluations;
    unsigned int numOfEvaluationsPerCycle;
    vector <unsigned int> baseNumCandidates;

    // If base number is one of the parameters
    if ((baseNumOptimizationMethod_ != "Brent") && (optimizeBaseNumber_)){
        fillVectorOfBaseNumCandidates(baseNumCandidates, lowerBoundBaseNumber, baseNumberUpperBound_);

    }

    //Go over each cycle
    for (size_t i = 0; i < numOfIterations_.size(); i++){
        numOfEvaluationsPerCycle = 0;
        clearVectorOfLikelihoods(numOfPoints_[i]);
        cout <<"##################################" << endl;
        cout << "*********  cycle "<< i <<"  **************"<<endl;     
        //Go over each point at cycle i 
        for (size_t j = 0; j < numOfPoints_[i]; j++){
            numOfEvaluations = 0;
            std::cout << "Starting cycle with Point #" << j <<"...."<<endl;;
            printLikParameters(vectorOfLikelohoods_[j], 0);
            //If the number of optimization iterations is larger than zero, optimize the number of times as specified
            if (numOfIterations_[i] > 0){
                numOfEvaluations = optimizeModelParameters(vectorOfLikelohoods_[j], tolerance_, numOfIterations_[i], baseNumCandidates);
                          
            }
            std:: cout << "Number of evaluations per point is : " << numOfEvaluations << endl;
            numOfEvaluationsPerCycle += numOfEvaluations;
            std:: cout <<"*****************************" << endl;            
        }
        totalNumOfEvaluations += numOfEvaluationsPerCycle;
        //sort the vector of likelihoods, such that the worst likelihood is at the end
        sort(vectorOfLikelohoods_.begin(), vectorOfLikelohoods_.end(), compareLikValues);
        printLikelihoodVectorValues(vectorOfLikelohoods_);
        
    }
    const string outPath = (ChromEvolOptions::resultsPathDir_ == "none") ? (ChromEvolOptions::resultsPathDir_) : (ChromEvolOptions::resultsPathDir_ + "//" + "likelihood.txt");
    const string outPathFreq = (ChromEvolOptions::resultsPathDir_ == "none") ? (ChromEvolOptions::resultsPathDir_) : (ChromEvolOptions::resultsPathDir_ + "//" + "inferred_rootFreq.txt");
   
    printRootFrequencies(vectorOfLikelohoods_[0], outPathFreq);
    cout <<"*****  Final Optimized -logL *********"  <<endl;
    printLikParameters(vectorOfLikelohoods_[0], 1, outPath);
    std:: cout << "final number of evaluations is : " << totalNumOfEvaluations << endl;
}

/********************************************************************************/
void ChromosomeNumberOptimizer::clearVectorOfLikelihoods(size_t new_size){
    while(vectorOfLikelohoods_.size() > new_size){
        //deleteTreeLikAssociatedAttributes(vectorOfLikelohoods_[vectorOfLikelohoods_.size()-1]);
        SingleProcessPhyloLikelihood* lik_to_del = vectorOfLikelohoods_.back(); 
        vectorOfLikelohoods_.pop_back();
        delete lik_to_del;
    }
}
/*********************************************************************************/
// void ChromosomeNumberOptimizer::deleteTreeLikAssociatedAttributes(SingleProcessPhyloLikelihood &lik){
//     const SubstitutionModelSet* modelSet = lik.getSubstitutionModelSet();
//     const DiscreteDistribution* rateDist = lik.getRateDistribution();
//     delete modelSet;
//     delete rateDist;
// }
/***********************************************************************************/
bool ChromosomeNumberOptimizer::compareLikValues(SingleProcessPhyloLikelihood* lik1, SingleProcessPhyloLikelihood* lik2){
    return (lik1->getValue() < lik2->getValue());
}
/***********************************************************************************/
void ChromosomeNumberOptimizer::printLikParameters(SingleProcessPhyloLikelihood* lik, unsigned int optimized, const string filePath) const{
    //double res = lik.getLikelihood();
    ofstream outFile;
    if (filePath != "none"){
        outFile.open(filePath);
    }
    if (optimized == 0){
        std:: cout << "Initial likelihood is : "<< lik->getValue() << endl;
    }else{
        std:: cout << "Optimized likelihood is : "<< lik->getValue() << endl;
        if (filePath != "none"){
            outFile << "Final optimized likelihood is: "<< lik->getValue() << endl;
        }
    }
    
    std:: cout << "Parameters are:" << endl;
    if (filePath != "none"){
        outFile << "Optimized parameters are:"<<endl;
    }
    ParameterList substitutionModelParams = lik->getSubstitutionModelParameters();
    std::vector<std::string> paramsNames = substitutionModelParams.getParameterNames();
    for (int i = 0; i < (int)(paramsNames.size()); i++){
        if (paramsNames[i] == "Chromosome.baseNum_1"){
            std::cout << paramsNames[i] << "value is "<< (int)(lik->getLikelihoodCalculation()->getParameter(paramsNames[i]).getValue()) <<endl;
            if (filePath != "none"){
                outFile <<  paramsNames[i] << "value is "<< (int)(lik->getLikelihoodCalculation()->getParameter(paramsNames[i]).getValue()) <<endl;
            }
        }else{
            std::cout << paramsNames[i] << "value is "<< lik->getLikelihoodCalculation()->getParameter(paramsNames[i]).getValue() <<endl;
            if (filePath != "none"){
                outFile << paramsNames[i] << "value is "<< lik->getLikelihoodCalculation()->getParameter(paramsNames[i]).getValue() <<endl;
            }
        }
        
    }
    if (filePath != "none"){
        outFile.close();
    }
    std::cout <<"***"<<endl;

}
/*************************************************************************************/
void ChromosomeNumberOptimizer::printLikelihoodVectorValues(std::vector <SingleProcessPhyloLikelihood*> lik_vec) const{
    std :: cout <<"The likelihoods at the end of cycle are :"<<endl;
    for (size_t i = 0; i < lik_vec.size(); i++){
        std :: cout << lik_vec[i]->getValue() << endl;
    }
}

/******************************************************************************/
void ChromosomeNumberOptimizer::printRootFrequencies(SingleProcessPhyloLikelihood* lik, const string filePath) const{
    ofstream outFile;
    if (filePath != "none"){
        outFile.open(filePath);
    }
    
    ValueRef <RowLik> rootFreqVector = lik->getLikelihoodCalculationSingleProcess()->getRootFreqs();
    for (size_t s = 0; s < (size_t)rootFreqVector->getTargetValue().size(); s ++){
        cout << "F[" << s + alphabet_->getMin() << "] = " << rootFreqVector.get()->getTargetValue()[s] << endl;
        if (filePath != "none"){
            outFile << "F[" << s + alphabet_->getMin() << "] = " << rootFreqVector.get()->getTargetValue()[s] << endl;
        }
    }
    if (filePath != "none"){
        outFile.close();
    }

}

/***********************************************************************************/
void ChromosomeNumberOptimizer::fillVectorOfBaseNumCandidates(std::vector <unsigned int> &baseNumCandidates, unsigned int lowerBound, unsigned int upperBound) const{
    if (baseNumOptimizationMethod_ == "Ranges"){
        getAllPossibleChrRanges(baseNumCandidates);

    }
    else if ((baseNumOptimizationMethod_ == "Sequential") || (baseNumCandidates.size() == 0)){

        for (unsigned int chr = (unsigned int)lowerBound; chr <= upperBound; chr++){
            baseNumCandidates.push_back(chr);
        }

    }

}
/***************************************************************************************/
void ChromosomeNumberOptimizer::getAllPossibleChrRanges(std::vector <unsigned int> &baseNumCandidates) const{
    size_t numOfSequences = vsc_->getNumberOfSequences();
    unsigned int minRange = 0;
    vector <string> sequenceNames = vsc_->getSequencesNames();
    for (size_t i = 0; i < numOfSequences; i++){
        if (i == numOfSequences-1){
            continue;
        }
        BasicSequence seq1 = vsc_->getSequence(sequenceNames[i]);
        int chrNum1 = seq1.getValue(0);
        if (chrNum1 == -1){
            continue;
        }
        for (size_t j = i + 1; j < numOfSequences; j++){
            BasicSequence seq2 = vsc_->getSequence(sequenceNames[j]);
            int chrNum2 = seq2.getValue(0);
            if (chrNum2 == -1){
                continue;
            }
            unsigned int chrRange = (unsigned int)(abs(chrNum1 - chrNum2));
            if (chrRange == 0 || chrRange == 1){
                continue;
            }
            else if (chrRange == 2){
                continue;
            }
            if (!std::count(baseNumCandidates.begin(), baseNumCandidates.end(), chrRange)){
                if (minRange == 0){
                    minRange = chrRange;
                }else{
                    if (chrRange < minRange){
                        minRange = chrRange;
                    }
                }
                baseNumCandidates.push_back(chrRange);

            }

        }
    }
    if (minRange > 3){
        for (unsigned int i = 3; i < minRange; i++){
            baseNumCandidates.push_back(i);
        }

    }

}
/**********************************************************************************/
unsigned int ChromosomeNumberOptimizer::optimizeModelParameters(SingleProcessPhyloLikelihood* tl, double tol, unsigned int maxNumOfIterations, std::vector <unsigned int> &baseNumCandidates){
    unsigned int numOfEvaluations = 0;
    // if (standardOptimization_){
    //     double prevLogLik;
    //     ParameterList params;
    //     for (size_t i = 0; i < maxNumOfIterations; i++){
    //         std::cout << "Iteration #"<<i <<endl;
    //         prevLogLik = tl->getValue();
    //         params = tl->getSubstitutionModelParameters();
    //         numOfEvaluations += OptimizationTools::optimizeNumericalParameters(tl, params, 0, 1, tol, 2, ApplicationTools::message.get(), ApplicationTools::message.get(), false, 0, OptimizationTools::OPTIMIZATION_NEWTON, OptimizationTools::OPTIMIZATION_BRENT, (unsigned int)(BrentBracketing_));
    //         printLikParameters(*tl, 1);
    //         if (abs(tl->getValue() - prevLogLik) < tol){
    //             break;
    //         }
    //     }
    //     std::cout <<"..."<<endl;
    //}else{
    if (typeOfOptimizer_ == "Brent"){
        numOfEvaluations += optimizeModelParametersOneDimension(tl, tol, maxNumOfIterations, baseNumCandidates);
    }else if (typeOfOptimizer_ == "gradient"){

        numOfEvaluations += optimizeMultiDimensions(tl, tol, maxNumOfIterations);

    }else{
        numOfEvaluations += useMixedOptimizers(tl, tol, maxNumOfIterations, baseNumCandidates);
    }
        
    //}
    
    return numOfEvaluations;
    
}
/****************************************************************************************/
unsigned int ChromosomeNumberOptimizer::optimizeMultiDimensions(SingleProcessPhyloLikelihood* tl, double tol, unsigned int maxNumOfIterations, bool mixed, unsigned int currentIterNum){
    DerivableSecondOrder* f = tl;
    unique_ptr<AbstractNumericalDerivative> fnum;
    fnum.reset(new TwoPointsNumericalDerivative(f));
    fnum->setInterval(0.0000001);
    ConjugateGradientMultiDimensions* optimizer = new ConjugateGradientMultiDimensions(fnum.get());
    ParameterList tmp = tl->getSubstitutionModelParameters();
    fnum->setParametersToDerivate(tmp.getParameterNames());
    optimizer->setVerbose(1);
    optimizer->setProfiler(0);
    optimizer->setMessageHandler(0);
    optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    optimizer->getStopCondition()->setTolerance(tol* 0.1);
    optimizer->setMaximumNumberOfEvaluations(1000);

    unsigned int numOfEvaluations = 0;
    double currentLikelihood = tl->getValue();
    double prevLikelihood;
    for (size_t i = 0; i < maxNumOfIterations; i++){
        if(mixed){
            std::cout << "Iteration #"<< currentIterNum <<endl;

        }else{
            std::cout << "Iteration #"<< i <<endl;
        }
        
        ParameterList paramsFull = tl->getSubstitutionModelParameters();
        std::vector <string> paramsNames = getNonFixedParams(fixedParams_, paramsFull);
        ParameterList params = tl->getParameters().createSubList(paramsNames);
        std::shared_ptr<IntervalConstraint> interval = make_shared<IntervalConstraint>(lowerBoundOfRateParam + 0.0000000001, upperBoundOfRateParam, true, true);
        for (size_t j = 0; j < params.size(); j++){
            params[j].setConstraint(interval);
        }
        prevLikelihood = currentLikelihood;
        optimizer->init(params);
        currentLikelihood = optimizer->optimize();
        printLikParameters(tl, 1);
        if (abs(prevLikelihood-currentLikelihood) < tol){
            break;
        }
        
        
    }
    numOfEvaluations += optimizer->getNumberOfEvaluations();
    if (!mixed){
        std::cout <<"..."<<endl;
    }
    //std::cout << "The final number of evaluations is: "<< numOfEvaluations << endl;
    delete optimizer;
    return numOfEvaluations;

}
/*******************************************************************************/

unsigned int ChromosomeNumberOptimizer::useMixedOptimizers(SingleProcessPhyloLikelihood* tl, double tol, unsigned int maxNumOfIterations, vector <unsigned int> &baseNumCandidates){
    std::vector<size_t> optimization = RandomTools::randMultinomial(maxNumOfIterations, probsForMixedOptimization_);
    unsigned int numOfEvaluations = 0;
    for (size_t i = 0; i < maxNumOfIterations; i++){
        double prevLikelihood = tl->getValue();
        if (optimization[i] == 0){
            std::cout << "Optimizing with Brent" <<endl;
            numOfEvaluations += optimizeModelParametersOneDimension(tl, tol, 1, baseNumCandidates, true, (unsigned int)i);
        }else{
            std::cout << "Optimizing with Gradient Descent" <<endl;
            numOfEvaluations += optimizeMultiDimensions(tl, tol, 1, true, (unsigned int)i);
        }
        double currentLikValue = tl->getValue();
        if (abs(prevLikelihood-currentLikValue) < tol){
            break;
        }


    }
    return numOfEvaluations;

}
/*******************************************************************************/
unsigned int ChromosomeNumberOptimizer::optimizeModelParametersOneDimension(SingleProcessPhyloLikelihood* tl, double tol, unsigned int maxNumOfIterations, std::vector <unsigned int> &baseNumCandidates, bool mixed, unsigned curentIterNum){

    // Initialize optimizer
    map <string,pair<string, bool>> paramPairs;
    constructParamPairsMap(paramPairs);
    DerivableSecondOrder* f = tl;
    BrentOneDimension* optimizer = new BrentOneDimension(f);
    optimizer->setVerbose(1);
    //optimizer->setProfiler(ApplicationTools::message.get());
    //optimizer->setMessageHandler(ApplicationTools::message.get());
    optimizer->setProfiler(0);
    optimizer->setMessageHandler(0);
    optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    optimizer->setMaximumNumberOfEvaluations(100);
    std::cout <<"max chromosome number: " << alphabet_->getMax() << endl;
    if (BrentBracketing_ == 1){
        optimizer->setBracketing(BrentOneDimension::BRACKET_INWARD);

    }else if (BrentBracketing_ == 2){
        optimizer->setBracketing(BrentOneDimension::BRACKET_SIMPLE);
    }else{
        optimizer->setBracketing(BrentOneDimension::BRACKET_OUTWARD);
    }
    
    double currentLikelihood = tl->getValue();
    double prevLikelihood;
    unsigned int numOfEvaluations = 0;
    //unsigned int numOfBaseNumEval = 0;

    for (size_t i = 0; i < maxNumOfIterations; i++){
        if (mixed){
            std::cout << "Iteration #"<<curentIterNum <<endl;

        }else{
            std::cout << "Iteration #"<<i <<endl;
        }
        ParameterList params = tl->getParameters();
        //ParameterList params = optimizer->getFunction()->getParameters();
        ParameterList substitutionModelParams = tl->getSubstitutionModelParameters();
        vector <string> paramsNames = substitutionModelParams.getParameterNames();        
        size_t nbParams = substitutionModelParams.size();
        prevLikelihood = currentLikelihood;
        
        for (size_t j = 0; j < nbParams; j ++){
            //numOfLikEvaluations = tl->getNumberOfLikelihoodEvaluations();
            //ParameterList params = tl->getSubstitutionModelParameters();
            if (fixedParams_[j]){
                continue;
            }
            const string nameOfParam = substitutionModelParams[j].getName();
            Parameter param = params.getParameter(nameOfParam);
            cout <<"Parameter name is: "<< nameOfParam << endl;
            string paramNameInModel = findParameterNameInModel(nameOfParam);
            const ChromosomeSubstitutionModel* model = dynamic_cast<const ChromosomeSubstitutionModel*>(dynamic_cast<const SubstitutionModel*>(tl->getLikelihoodCalculationSingleProcess()->getSubstitutionProcess().getModel(1)));
            //const ChromosomeSubstitutionModel* model = dynamic_cast<const ChromosomeSubstitutionModel*>(substitutionModel);

            //model->setBoundsForEquivalentParameter(param, paramNameInModel);
            //model->checkParametersBounds();
            double lowerBound = dynamic_pointer_cast<IntervalConstraint>(tl->getLikelihoodCalculation()->getParameter(param.getName()).getConstraint())->getLowerBound();
            double upperBound = dynamic_pointer_cast<IntervalConstraint>(tl->getLikelihoodCalculation()->getParameter(param.getName()).getConstraint())->getUpperBound();
            setNewBounds(tl->getLikelihoodCalculation()->getParameters(), param, paramPairs, &lowerBound, model);
 
            //model->checkParametersBounds();
            if ((baseNumOptimizationMethod_ != "Brent") && (param.getName() == "Chromosome.baseNum_1")){
                if (optimizeBaseNumber_){
                    optimizeBaseNum(tl, j, baseNumCandidates, &currentLikelihood, lowerBound, upperBound);
                    std::cout << "parameter value after optimization "<< tl->getLikelihoodCalculation()->getParameter(param.getName()).getValue() << endl;
                    continue;

                }
            }        
            if ((i == 1) & (maxNumOfIterations > 2)){
                optimizer->getStopCondition()->setTolerance(tol* 2);
            }else{
                optimizer->getStopCondition()->setTolerance(tol);
            }
            if (param.getName() != "Chromosome.baseNum_1"){
                optimizer->setInitialInterval(lowerBound + 1e-10, upperBound);
            }else{
                optimizer->setInitialInterval(lowerBound, upperBound);
            }            
            optimizer->init(params.createSubList(param.getName()));
            currentLikelihood = optimizer->optimize();
            std::cout <<"Parameter value after optimization: "<< tl->getLikelihoodCalculation()->getParameter(param.getName()).getValue() <<endl;
            std::cout << "***"<<endl;                        
        }
        printLikParameters(tl, 1);
        
        if (abs(prevLikelihood-currentLikelihood) < tol){
            break;
        }
        numOfEvaluations += optimizer->getNumberOfEvaluations();
       
    }
    //std::cout << "Number of likelihood evaluations per parameter is "<< numOfBaseNumEval << endl;
    if (!mixed){
        std::cout <<"..."<<endl;
    }
    delete optimizer;
    return numOfEvaluations;
}
/*******************************************************************************/
void ChromosomeNumberOptimizer::setNewBounds(const ParameterList params, Parameter &param, map<string, pair<string, bool>> &paramPairsMap, double* lowerBound, const ChromosomeSubstitutionModel* model){
    if (paramPairsMap.find(param.getName()) == paramPairsMap.end()){
        return;
    }
    vector<string> parametersNames = params.getParameterNames();
    pair <string, bool> matchedParamAndType = paramPairsMap[param.getName()];
    if (!(std::count(parametersNames.begin(), parametersNames.end(), matchedParamAndType.first))){
        return;
    }
    Parameter matchedParam = params.getParameter(matchedParamAndType.first);
    double valueOfMatched = matchedParam.getValue();
    std::shared_ptr<IntervalConstraint> interval = dynamic_pointer_cast<IntervalConstraint>(param.getConstraint());
    if (matchedParamAndType.second){
        if (!(ChromEvolOptions::rateChangeType_)){
            *lowerBound = std::max(lowerBoundOfRateParam, -valueOfMatched*(model->getMax()-1));
            interval->setLowerBound(*lowerBound, interval->strictLowerBound());
        }
        

    }else{
        if (!(ChromEvolOptions::rateChangeType_)){
            *lowerBound = -valueOfMatched/(model->getMax()-1);
            interval->setLowerBound(*lowerBound, interval->strictLowerBound());

        }
        
    }
    
}
/*******************************************************************************/
void ChromosomeNumberOptimizer::constructParamPairsMap(map<string, pair<string, bool>> &paramPairsMap){
    paramPairsMap["Chromosome.gain_1"] = pair<string, bool>("Chromosome.gainR_1", true);
    paramPairsMap["Chromosome.gainR_1"] = pair<string, bool>("Chromosome.gain_1", false);
    paramPairsMap["Chromosome.loss_1"] = pair<string, bool>("Chromosome.lossR_1", true);
    paramPairsMap["Chromosome.lossR_1"] = pair<string, bool>("Chromosome.loss_1", false);
    paramPairsMap["Chromosome.dupl_1"] = pair<string, bool>("Chromosome.duplR_1", true);
    paramPairsMap["Chromosome.duplR_1"] = pair<string, bool>("Chromosome.dupl_1", false);
}
/*******************************************************************************/
std::vector <string> ChromosomeNumberOptimizer::getNonFixedParams(std::vector <unsigned int> fixedParams, ParameterList &allParams) const{
    std:: vector <string> paramsNames;
    for (size_t i = 0; i < allParams.size(); i++){
        if (!(fixedParams[i])){
            paramsNames.push_back(allParams[i].getName());
        }

    }
    return paramsNames;
}
/***********************************************************************************/
string ChromosomeNumberOptimizer::findParameterNameInModel(string fullParameterName) const{
    StringTokenizer st(fullParameterName, "._", false, false);
    string paramName;
    while(st.hasMoreToken()){
        string token = st.nextToken();
        if (token == "Chromosome"){
            continue;
        }else if (token == "1"){
            continue;
        }else{
            paramName = token;
            break;
        }
    }
    return paramName;
}
/***************************************************************************************/
void ChromosomeNumberOptimizer::optimizeBaseNum(SingleProcessPhyloLikelihood* tl, size_t index, std::vector <unsigned int> baseNumCandidates, double* currentLikelihood, double lowerBound, double upperBound){

    Function* func = tl;
    ParameterList params = tl->getSubstitutionModelParameters();
    size_t best_i = (size_t)(params[index].getValue());
    //double f_value = func->f(params);
    double f_value = *currentLikelihood;
    
    for (size_t i = 0; i < baseNumCandidates.size(); i++){
        unsigned int baseNum = baseNumCandidates[i];
        params[index].setValue((double)baseNum);
        //ParameterList params = tl->getSubstitutionModelParameters();
        double f_i = func->f(params);
        if (f_i < f_value){
            best_i = baseNum;
            f_value = f_i;
        }
    }
    params[index].setValue((double)best_i);
    func->f(params);
    //param.setValue((double)best_i);
    *currentLikelihood = f_value;

}