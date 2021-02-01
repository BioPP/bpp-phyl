#include "ChromosomeNumberMng.h"
#include "ChromEvolOptions.h"

using namespace bpp;

void ChromosomeNumberMng::getCharacterData (const string& path){
    ChromosomeAlphabet* alphaInitial = new ChromosomeAlphabet(ChromEvolOptions::minAlpha_, ChromEvolOptions::maxAlpha_);
    VectorSequenceContainer* initialSetOfSequences = chrFasta::readSequencesFromFile(path, alphaInitial);
    size_t numOfSequences = initialSetOfSequences->getNumberOfSequences();
    vector <string> sequenceNames = initialSetOfSequences->getSequencesNames();

    unsigned int maxNumberOfChr = 1; //the minimal number of chromosomes cannot be zero
    unsigned int minNumOfChr = ChromEvolOptions::maxAlpha_;

    std::vector <int> UniqueCharacterStates;
    cout<<"vector size is "<< UniqueCharacterStates.size()<<endl;
    for (size_t i = 0; i < numOfSequences; i++){
        BasicSequence seq = initialSetOfSequences->getSequence(sequenceNames[i]);
        int character = seq.getValue(0);
        if (character == -1){
            continue;
        }
        if (character == static_cast<int>(ChromEvolOptions::maxAlpha_)+1){
            continue;
        }
        // if it is a composite state
        if (character > static_cast<int>(ChromEvolOptions::maxAlpha_) +1){
            const std::vector<int> compositeCharacters = alphaInitial->getSetOfStatesForAComposite(character);
            for (size_t j = 0; j < compositeCharacters.size(); j++){
                if ((unsigned int) compositeCharacters[j] > maxNumberOfChr){
                    maxNumberOfChr = compositeCharacters[j];
                }
                if ((unsigned int) compositeCharacters[j] < minNumOfChr){
                    minNumOfChr = compositeCharacters[j];
                }
                
            }
            continue;
        }

        if (!std::count(UniqueCharacterStates.begin(), UniqueCharacterStates.end(), character)){
            UniqueCharacterStates.push_back(character);
        }
        if ((unsigned int) character > maxNumberOfChr){
            maxNumberOfChr = character;
        }
        if ((unsigned int) character < minNumOfChr){
            minNumOfChr = character;
        }

    }
    numberOfUniqueStates_ = (unsigned int)UniqueCharacterStates.size() + alphaInitial->getNumberOfCompositeStates();
    chrRange_ = maxNumberOfChr - minNumOfChr;
    if (ChromEvolOptions::baseNum_ != IgnoreParam){
        if (ChromEvolOptions::baseNum_ > (int)chrRange_){
            chrRange_ = ChromEvolOptions::baseNum_ + 1;
        }
    }
    cout <<"Number of unique states is " << numberOfUniqueStates_ <<endl;

    setMaxChrNum(maxNumberOfChr);
    setMinChrNum(minNumOfChr);

    vsc_ = resizeAlphabetForSequenceContainer(initialSetOfSequences, alphaInitial);
    delete initialSetOfSequences;
    delete alphaInitial;
    return;
}
// /*******************************************************************************************************************/
VectorSiteContainer* ChromosomeNumberMng::resizeAlphabetForSequenceContainer(VectorSequenceContainer* vsc, ChromosomeAlphabet* alphaInitial){
    size_t numOfSequences = vsc->getNumberOfSequences();
    vector <string> sequenceNames = vsc->getSequencesNames();
    alphabet_ = new ChromosomeAlphabet(ChromEvolOptions::minChrNum_,ChromEvolOptions::maxChrNum_);
        // fill with composite values
    if (alphaInitial->getNumberOfCompositeStates() > 0){
        const std::map <int, std::map<int, double>> compositeStates = alphaInitial->getCompositeStatesMap();
        std::map <int, std::map<int, double>>::const_iterator it = compositeStates.begin();
        while (it != compositeStates.end()){
            int compositeState = it->first;
            std::string charComposite = alphaInitial->intToChar(compositeState);
            alphabet_->setCompositeState(charComposite);
            it++;
        }
    }
    VectorSiteContainer* resized_alphabet_site_container = new VectorSiteContainer(alphabet_);
    for (size_t i = 0; i < numOfSequences; i++){
        BasicSequence seq = vsc->getSequence(sequenceNames[i]);
        BasicSequence new_seq = BasicSequence(seq.getName(), seq.getChar(0), alphabet_);
        resized_alphabet_site_container->addSequence(new_seq);

    }
    return resized_alphabet_site_container;
}
/*******************************************************************************************/
void ChromosomeNumberMng::setMaxChrNum(unsigned int maxNumberOfChr){
    if (ChromEvolOptions::maxChrNum_ < 0){
        ChromEvolOptions::maxChrNum_ = maxNumberOfChr + abs(ChromEvolOptions::maxChrNum_);
    }else{
        if ((int)maxNumberOfChr > ChromEvolOptions::maxChrNum_){
            ChromEvolOptions::maxChrNum_ = maxNumberOfChr;
        }
    }

}
/****************************************************************************/
void ChromosomeNumberMng::setMinChrNum(unsigned int minNumberOfChr){
    if (ChromEvolOptions::minChrNum_ < 0){
        ChromEvolOptions::minChrNum_ = minNumberOfChr - abs(ChromEvolOptions::minChrNum_);
    }else{
        if ((int)minNumberOfChr < ChromEvolOptions::minChrNum_){
            ChromEvolOptions::minChrNum_ = minNumberOfChr;
        }

    }
}
/********************************************************************************************/

void ChromosomeNumberMng::getTree(const string& path, double treeLength){
    Newick reader;
    tree_ = reader.readPTree(path);
    double treeLengthToScale = (treeLength > 0) ? treeLength : (double) numberOfUniqueStates_;
    rescale_tree(tree_, treeLengthToScale);
    return;

}
/****************************************************************************/
void ChromosomeNumberMng::rescale_tree(PhyloTree* tree, double chrRange){
    double scale_tree_factor = 1.0;
    //string tree_str = TreeTemplateTools::treeToParenthesis(*tree);
    //std :: cout << tree_str << endl;
    bool rooted = tree->isRooted();
    if (!rooted){
        //throw UnrootedTreeException("The given input tree is unrooted. Tree must be rooted!", tree);
        throw Exception("The given input tree is unrooted. Tree must be rooted!\n");
    }
    if (ChromEvolOptions::branchMul_ == 1.0){
        return;
    }else{
        //tree must be rescaled
        double treeLength = tree->getTotalLength();
        if (ChromEvolOptions::branchMul_ == 999){
            scale_tree_factor = chrRange/treeLength;
        }
        tree->scaleTree(scale_tree_factor);

    }

}
/*****************************************************************************************/
/* void ChromosomeNumberMng::getMaxParsimonyUpperBound(double* parsimonyBound) const{

    TreeTemplate<Node>* treeForParsimonyBound = tree_->clone();
    DRTreeParsimonyScore maxParsimonyObject = DRTreeParsimonyScore(*treeForParsimonyBound, *vsc_);
    *parsimonyBound = (maxParsimonyObject.getScore())/(tree_->getTotalLength());
    delete treeForParsimonyBound;
    return;   

} */
/*****************************************************************************************/
/* ChromosomeNumberOptimizer* ChromosomeNumberMng::optimizeLikelihoodMultiStartPoints() const{
    
    vector<double> modelParams;
    modelParams.reserve(ChromosomeSubstitutionModel::NUM_OF_CHR_PARAMS);
    ChromEvolOptions::initVectorOfChrNumParameters(modelParams);
    double parsimonyBound = 0;
    if (ChromEvolOptions::maxParsimonyBound_){
        getMaxParsimonyUpperBound(&parsimonyBound);
    }
    bool calculateDerivatives = true;
    if (ChromEvolOptions::optimizationMethod_ == "Brent"){
        calculateDerivatives  = false;
    }
    unsigned int maxBaseNumTransition = (ChromEvolOptions::simulateData_) ? ChromEvolOptions::maxBaseNumTransition_ : chrRange_;
    ChromosomeNumberOptimizer* opt = new ChromosomeNumberOptimizer(tree_, alphabet_, vsc_, maxBaseNumTransition);
    opt->initModels(modelParams, parsimonyBound, ChromEvolOptions::rateChangeType_, ChromEvolOptions::seed_, ChromEvolOptions::OptPointsNum_[0], calculateDerivatives, ChromEvolOptions::fixedFrequenciesFilePath_, ChromEvolOptions::fixedParams_);

    //initialize all the optimization specific parameters
    opt->initOptimizer(ChromEvolOptions::OptPointsNum_, ChromEvolOptions::OptIterNum_, ChromEvolOptions::optimizationMethod_, ChromEvolOptions::baseNumOptimizationMethod_,
        ChromEvolOptions::tolerance_, ChromEvolOptions::standardOptimization_, ChromEvolOptions::BrentBracketing_, 
        ChromEvolOptions::probsForMixedOptimization_);
    //optimize models
    opt->optimize();
    // it is safe to delete the chrOptimizer, because the destructor doesn't delete nothing associated with the vector of likelihoods
    return opt;
       
} */
/******************************************************************************************************/
/* void ChromosomeNumberMng::getJointMLAncestralReconstruction(DRNonHomogeneousTreeLikelihood* lik) const{
    std::cout << "ML Joint Ancestral Reconstruction"<<endl;
    TransitionModel* model = lik->getSubstitutionModelSet()->getModel(0);
    std::vector<double> rootFreqs = lik->getRootFrequencies(0);
    std::map<int, VVVdouble> Pijt;
    TreeTemplate<Node> tree = lik->getTree();
    vector <int> nodesIds = tree.getNodesId();
    for (size_t n = 0; n < nodesIds.size(); n++){
        Pijt[nodesIds[n]] = lik->getTransitionProbabilitiesPerRateClass(nodesIds[n], 0);
    }
    MLAncestralStateReconstruction* ancr = new MLAncestralStateReconstruction(lik, model, rootFreqs, &Pijt);
    ancr->computeJointLikelihood();
    std::map<int, std::vector<size_t> > ancestors = ancr->getAllAncestralStates();
    if (ChromEvolOptions::resultsPathDir_ == "none"){
        printTreeWithStates(lik->getTree(), ancestors, ChromEvolOptions::resultsPathDir_);
    }else{
        const string outFilePath = ChromEvolOptions::resultsPathDir_ + "//" + "MLAncestralReconstruction.tree";
        printTreeWithStates(lik->getTree(), ancestors, outFilePath);

    }
    

    delete ancr;

} */
/**************************************************************************************************************/
/* std::map<int, std::map<size_t, VVdouble>> ChromosomeNumberMng::getMarginalAncestralReconstruction(DRNonHomogeneousTreeLikelihood* lik) const{
    std::cout << "Marginal Ancestral Reconstruction"<<endl;
    MarginalNonRevAncestralStateReconstruction* ancr = new MarginalNonRevAncestralStateReconstruction(lik);
    ancr->computePosteriorProbabilitiesOfNodesForEachStatePerSite();
    std::map<int, std::vector<size_t> > ancestors = ancr->getAllAncestralStates();
    std::map<int, map<size_t, std::vector<double>>>* probs = ancr->getPosteriorProbForAllNodesAndStatesPerSite();
    if (ChromEvolOptions::resultsPathDir_ == "none"){
        printTreeWithStates(lik->getTree(), ancestors, ChromEvolOptions::resultsPathDir_, probs);
    }else{
        const string outFilePath = ChromEvolOptions::resultsPathDir_ +"//"+"MarginalAncestralReconstruction.tree";
        printTreeWithStates(lik->getTree(), ancestors, outFilePath, probs);

    }
    
    std::map<int, std::map<size_t, VVdouble>> jointProbabilitiesFatherSon = ancr->getAllJointFatherNodeProbabilities();
    vector<double> rootPosterior = ancr->getRootPosteriorProb();
    printPosteriorProbNodes(jointProbabilitiesFatherSon, rootPosterior);
    delete ancr;
    return jointProbabilitiesFatherSon;

} */
/**********************************************************************************************************************/
/* void ChromosomeNumberMng::computeExpectations(DRNonHomogeneousTreeLikelihood* lik, std::map<int, std::map<size_t, VVdouble>>& jointProbabilitiesFatherSon, int numOfSimulations) const{

    ComputeChromosomeTransitionsExp* expCalculator = new ComputeChromosomeTransitionsExp(lik, jointProbabilitiesFatherSon, ChromEvolOptions::jumpTypeMethod_);
    expCalculator->runSimulations(numOfSimulations);
    expCalculator->computeExpectationPerType();
    if (ChromEvolOptions::resultsPathDir_ == "none"){
        expCalculator->printResults();
    }else{
        const string outFilePath = ChromEvolOptions::resultsPathDir_+"//"+ "expectations.txt";
        //const string outFilePathForNonAccounted = ChromEvolOptions::resultsPathDir_+"//"+ "exp_nonAccounted_branches.txt";
        const string outFilePathHeuristics = ChromEvolOptions::resultsPathDir_+"//"+ "expectations_second_round.txt";
        const string outTreePath = ChromEvolOptions::resultsPathDir_+"//"+ "exp.tree";
        expCalculator->printResults(outFilePath);
        expCalculator->runHeuristics();
        expCalculator->printResults(outFilePathHeuristics);
        TreeTemplate<Node>* expTree = expCalculator->getResultTree();
        string tree_str = printTree(*expTree);
        delete expTree;
        ofstream treeFile;
        treeFile.open(outTreePath);
        treeFile << tree_str;
        treeFile.close();
    }
    

    //delete
    delete expCalculator;
} */
/***********************************************************************************/
void ChromosomeNumberMng::runTest(){

    vector<double> modelParams;
    modelParams.reserve(ChromosomeSubstitutionModel::NUM_OF_CHR_PARAMS);
    ChromEvolOptions::initVectorOfChrNumParameters(modelParams);
    ParametrizablePhyloTree parTree(*tree_);
    unsigned int maxBaseNumTransition = chrRange_;
    std::shared_ptr<SubstitutionModel> chrModel = std::make_shared<ChromosomeSubstitutionModel>(alphabet_, modelParams, maxBaseNumTransition, ChromosomeSubstitutionModel::rootFreqType::FIXED, ChromosomeSubstitutionModel::LINEAR);
    DiscreteDistribution* rdist =  new GammaDiscreteRateDistribution(1, 1.0);
    vector <double> rootFreqs = getRootFrequencies(ChromEvolOptions::fixedFrequenciesFilePath_);
    std::shared_ptr<FixedFrequencySet> rootFreqsFixed = std::make_shared<FixedFrequencySet>(std::shared_ptr<const StateMap>(new CanonicalStateMap(chrModel.get()->getStateMap(), false)), rootFreqs);
    std::shared_ptr<FrequencySet> rootFrequencies = static_pointer_cast<FrequencySet>(rootFreqsFixed);
    std::vector<std::string> globalParameterNames;
    globalParameterNames.push_back("Chromosome.baseNum");
    globalParameterNames.push_back("Chromosome.baseNumR");
    globalParameterNames.push_back("Chromosome.dupl");
    globalParameterNames.push_back("Chromosome.loss");
    globalParameterNames.push_back("Chromosome.gain");
    //globalParameterNames.push_back("Chromosome.lossR");
    //globalParameterNames.push_back("Chromosome.gainR");
    //globalParameterNames.push_back("Chromosome.duplR");
    globalParameterNames.push_back("Chromosome.demi");
    NonHomogeneousSubstitutionProcess* subProSim= NonHomogeneousSubstitutionProcess::createNonHomogeneousSubstitutionProcess(chrModel, rdist, parTree.clone(), rootFrequencies, globalParameterNames);
    SubstitutionProcess* nsubPro=subProSim->clone();
    Context context;
    auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(context, *vsc_->clone(), *nsubPro, true);
    
    SingleProcessPhyloLikelihood ntl(context, lik, lik->getParameters());
    cout << setprecision(10) << "NewTL init: "  << ntl.getValue()  << endl;



}
/***********************************************************************************/
vector <double> ChromosomeNumberMng::getRootFrequencies(const string &path){
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
        if (static_cast<unsigned int> (freqs.size()) >= alphabet_->getSize()){
            if (freq_i > 0){
                throw Exception ("Invalid root frequencies file!\n");
            }

        }else{
            freqs.push_back(freq_i);
        }
    }
    return freqs;

}
/***********************************************************************************/
/* void ChromosomeNumberMng::runChromEvol(){
    if (ChromEvolOptions::simulateData_){
        //simulate data using a tree and a set of model parameters
        RandomTools::setSeed(static_cast<long>(ChromEvolOptions::seed_));
        simulateData();
        if (ChromEvolOptions::numOfDataToSimulate_ > 1){
            return;
        }

    }
    // optimize likelihood
    ChromosomeNumberOptimizer* chrOptimizer = optimizeLikelihoodMultiStartPoints();
    std::vector <DRNonHomogeneousTreeLikelihood> lik_vec = chrOptimizer->getVectorOfLikelihoods();
    // get joint ML ancestral reconstruction
    getJointMLAncestralReconstruction(&lik_vec[0]);
    //get Marginal ML ancestral reconstruction, and with the help of them- calculate expectations of transitions
    std::map<int, std::map<size_t, VVdouble>>  jointProbabilitiesFatherSon = getMarginalAncestralReconstruction(&lik_vec[0]);
    
    //compute expectations
    computeExpectations(&lik_vec[0], jointProbabilitiesFatherSon, ChromEvolOptions::NumOfSimulations_);
    //Clear the vector of likelihoods entirely
    delete chrOptimizer;


} */
/**************************************************************************************/
/* void ChromosomeNumberMng::printPosteriorProbNodes(std::map<int, std::map<size_t, VVdouble>>& jointProbabilitiesFatherSon, vector<double>& rootPosterior) const{
    if (ChromEvolOptions::resultsPathDir_ == "none"){
        throw Exception("Error in ChromosomeNumberMng::printPosteriorProbNodes(): No results file path!\n");
    }
    const string outPath = ChromEvolOptions::resultsPathDir_+"//"+ "ancestorsProbs.txt";
    ofstream outFile;
    outFile.open(outPath);
    vector<int> nodesIds = tree_->getNodesId();
    outFile <<"NODE\t";
    for (size_t i = 0; i < alphabet_->getSize(); i++){
        (i < alphabet_->getSize()-1) ? (outFile << (int)i + alphabet_->getMin() << "\t") : (outFile << (int)i + alphabet_->getMin() <<"\n");
    }
    for (size_t n = 0; n < nodesIds.size(); n++){
        int nodeId = nodesIds[n];
        string nodeName;
        if (tree_->isLeaf(nodeId)){
            nodeName = tree_->getNodeName(nodeId);

        }else{
            nodeName = "N" + std::to_string(nodeId);
        }
        outFile << nodeName <<"\t";
        if (tree_->getRootId() == nodeId){
            for (size_t state = 0; state < alphabet_->getSize(); state++){
                (state < alphabet_->getSize()-1) ? (outFile << rootPosterior[state] << "\t") : (outFile << rootPosterior[state] <<"\n");
            }
            continue;
        }
        
        for (size_t son = 0; son < alphabet_->getSize(); son++){
            double posteriorProb = 0;
            for (size_t father = 0; father < alphabet_->getSize(); father++){
                posteriorProb += jointProbabilitiesFatherSon[nodeId][0][son][father];

            }
            (son == alphabet_->getSize()-1) ? (outFile << posteriorProb <<"\n") : (outFile << posteriorProb <<"\t");          
        }
        
    }
    outFile.close();

} */
/**************************************************************************************/
/* void ChromosomeNumberMng::printSimulatedEvoPath(TreeTemplate<Node> tree, const string outPath, RASiteSimulationResult* simResult) const{
    ofstream outFile;
    outFile.open(outPath);
    size_t totalNumTransitions = 0;
    vector<int> nodesIds = tree.getNodesId();
    for (size_t n = 0; n < nodesIds.size(); n++){
        if (tree.isRoot(nodesIds[n])){
            outFile << "N-" + to_string(nodesIds[n]) << endl;
            outFile <<"\tThe root state is: "<< ((int)(simResult->getRootAncestralState())+ alphabet_->getMin()) <<endl;
        }else{
            if (tree.isLeaf(nodesIds[n])){
                outFile << tree.getNodeName(nodesIds[n]) << endl;
            }else{
                outFile << "N-" + to_string(nodesIds[n]) <<endl;

            }
            MutationPath mutPath = simResult->getMutationPath(nodesIds[n]);
            vector<size_t> states = mutPath.getStates();
            vector<double> times = mutPath.getTimes();
            totalNumTransitions += static_cast<int>(times.size());
            for (size_t i = 0; i < states.size(); i++){
                outFile <<"\tt = "<<times[i] << " to state = "<< ((int)(states[i]) + alphabet_->getMin()) << endl;

            }
            outFile <<"# Number of transitions per branch: "<< times.size() <<endl;   
            
        }
        
        outFile <<"*************************************"<<endl;
        
    }
    outFile <<"Total number of transitions is: "<< totalNumTransitions << endl;
    outFile.close();

} */
/**************************************************************************************/
/* void ChromosomeNumberMng::printTreeWithStates(TreeTemplate<Node> tree, std::map<int, std::vector<size_t> > ancestors, const string &filePath, std::map<int, map<size_t, std::vector<double>>>* probs) const{
    map <string, double> mapOfNodeNameProb;
    vector<int> nodesIds = tree.getNodesId();
    for (size_t n= 0; n < nodesIds.size(); n++){
        for (size_t i = 0; i < ancestors[nodesIds[0]].size(); i++){
            size_t state = ancestors[nodesIds[n]][i] + alphabet_->getMin();
            string prevName;
            
            if (tree.isLeaf(nodesIds[n])){
                prevName = tree.getNodeName(nodesIds[n]);
                const string newName = (prevName + "-"+ std::to_string(state));
                tree.setNodeName(nodesIds[n], newName);
                if (probs){
                    double postProb = (*probs)[nodesIds[n]][i][ancestors[nodesIds[n]][i]];
                    cout <<newName<< " state index: " << state << " : " << postProb <<endl;
                    mapOfNodeNameProb[newName] = postProb;

                }else{
                    cout <<newName<< " state index: " << state <<endl;
                }
            }else{
                prevName = "N" + std::to_string(nodesIds[n]);
                const string newName = (prevName + "-"+ std::to_string(state));
                tree.setNodeName(nodesIds[n], newName);
                if (probs){
                    double postProb = (*probs)[nodesIds[n]][i][ancestors[nodesIds[n]][i]];
                    cout <<newName<< " state index: " << state << " : " << postProb <<endl;
                    mapOfNodeNameProb[newName] = postProb;

                }else{
                    cout <<newName<< " state index: " << state <<endl;
                }
            }
            
            
        }
        
    }
    string tree_str;
    if (probs){
        tree_str = printTree(tree, &mapOfNodeNameProb);
    }
    else{
        tree_str = printTree(tree);
    }
    cout << tree_str << endl;
    if (filePath != "none"){
       ofstream outFile;
       outFile.open(filePath);
       outFile << tree_str << endl;
       outFile.close();
    }

} */
/****************************************************************************************/
/* string ChromosomeNumberMng::printTree(const TreeTemplate<Node>& tree, map<string, double>* mapNameProb)
{
  ostringstream s;
  s << "(";
  const Node* node = tree.getRootNode();
  if (node->isLeaf() && node->hasName()) // In case we have a tree like ((A:1.0)); where the root node is an unamed leaf!
  {
    s << node->getName();
    for (size_t i = 0; i < node->getNumberOfSons(); ++i)
    {
      s << "," << nodeToParenthesis(*node->getSon(i), mapNameProb);
    }
  }
  else
  {
    s << nodeToParenthesis(*node->getSon(0), mapNameProb);
    for (size_t i = 1; i < node->getNumberOfSons(); ++i)
    {
      s << "," << nodeToParenthesis(*node->getSon(i), mapNameProb);
    }
  }
  s << ")";
  s << tree.getRootNode()->getName();
  if (! mapNameProb){
    if (node->hasDistanceToFather()){
        s << ":" << node->getDistanceToFather();

    }
        
  }else{
      s << ":" << (*mapNameProb)[node->getName()];
  }

 
  
  s << ";" << endl;
  return s.str();
} */

/******************************************************************************/
/* string ChromosomeNumberMng::nodeToParenthesis(const Node& node, map<string, double>* mapNameProb)
{
  ostringstream s;
  if (node.isLeaf())
  {
    s << node.getName();
  }
  else
  {
    s << "(";
    s << nodeToParenthesis(*node[0], mapNameProb);
    for (int i = 1; i < static_cast<int>(node.getNumberOfSons()); i++)
    {
      s << "," << nodeToParenthesis(*node[i], mapNameProb);
    }
    s << ")";
  }
  if (! node.isLeaf()){
      s << node.getName();
  }
  if (node.hasDistanceToFather()){
    if (!mapNameProb){
        s << ":" << node.getDistanceToFather();
    }else{
        s << ":" << (*mapNameProb)[node.getName()];

    }
  }

  return s.str();
} */
/*********************************************************************************/
/* void ChromosomeNumberMng::simulateData(){
    if ((ChromEvolOptions::minChrNum_ <= 0) || (ChromEvolOptions::maxChrNum_ < 0)){
        throw Exception("ERROR!!! ChromosomeNumberMng::simulateData(): minimum and maximum chromsome number should be positive!");
    }
    if (ChromEvolOptions::maxChrNum_ <= ChromEvolOptions::minChrNum_){
        throw Exception("ERROR!!! ChromosomeNumberMng::simulateData(): maximum chromsome number should be larger than minimum chromosome number!");
    }
    alphabet_ = new ChromosomeAlphabet(ChromEvolOptions::minChrNum_,ChromEvolOptions::maxChrNum_);
    vector<double> modelParams;
    modelParams.reserve(ChromosomeSubstitutionModel::NUM_OF_CHR_PARAMS);
    ChromEvolOptions::initVectorOfChrNumParameters(modelParams);
    DiscreteDistribution* rdist = new GammaDiscreteRateDistribution(1, 1.0);
    SubstitutionModelSet* modelSet = new SubstitutionModelSet(alphabet_);
    ChromosomeSubstitutionModel* chrModel = new ChromosomeSubstitutionModel(alphabet_, modelParams, (unsigned int)ChromEvolOptions::maxBaseNumTransition_, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL,  ChromEvolOptions::rateChangeType_);
    vector <int> nodeIds = tree_->getNodesId();
    nodeIds.pop_back();
    modelSet->addModel(chrModel, nodeIds);
    if (ChromEvolOptions::fixedFrequenciesFilePath_ == "none"){
        throw Exception("ERROR!!! ChromosomeNumberMng::simulateData(): You need to path the file of fixed root frequencies!!");
    }
    ChromosomeNumberOptimizer::setFixedRootFrequencies(ChromEvolOptions::fixedFrequenciesFilePath_, modelSet);

    for (size_t i = 0; i < (size_t)ChromEvolOptions::numOfDataToSimulate_; i++){
        NonHomogeneousSequenceSimulator* sim = new NonHomogeneousSequenceSimulator(modelSet, rdist, tree_);
        RASiteSimulationResult* simResult = sim->dSimulateSite();
        vector <size_t> leavesStates = simResult->getFinalStates();
        vector<string> leavesNames = simResult->getLeaveNames();
        printSimulatedData(leavesStates, leavesNames, i);
        printSimulatedDataAndAncestors(simResult);
        if (ChromEvolOptions::resultsPathDir_ != "none"){
            printSimulatedEvoPath(*tree_, ChromEvolOptions::resultsPathDir_ +"//"+ "simulatedEvolutionPaths.txt", simResult);
        }
        delete simResult;
        delete sim;

    }
    delete modelSet;
    delete rdist;

} */
/*******************************************************************************/
/* void ChromosomeNumberMng::printSimulatedData(vector<size_t> leavesStates, vector<string> leavesNames, size_t iter){
    cout << "Simulated data #" << iter << endl;
    for (size_t i = 0; i < leavesNames.size(); i++){
        cout << leavesNames[i] << " "<< leavesStates[i] + alphabet_->getMin() <<endl;
    }
    cout << "******************************"<<endl;
    
    if (ChromEvolOptions::resultsPathDir_ != "none"){
        //create vector site container object and save fasta file.
        VectorSiteContainer* simulatedData = new VectorSiteContainer(alphabet_);
        for (size_t i = 0; i < leavesNames.size(); i++){
            int state = (int)leavesStates[i] + alphabet_->getMin();
            BasicSequence seq = BasicSequence(leavesNames[i], alphabet_->intToChar(state), static_cast <const Alphabet*>(alphabet_));
            simulatedData->addSequence(seq);
        }
        vsc_ = simulatedData;
        string pathForSimulatedData = ChromEvolOptions::resultsPathDir_ + "//"+ "chr_counts"+ to_string(static_cast<int>(iter)) +".fasta";
        Fasta fasta;
        fasta.writeSequences(pathForSimulatedData, *simulatedData);

    }


    
} */
/****************************************************************************/
/* void ChromosomeNumberMng::printSimulatedDataAndAncestors(RASiteSimulationResult* simResult) const{
    std::map<int, std::vector<size_t> > ancestors;
    vector<int> nodesIds = tree_->getNodesId();
    for (size_t i = 0; i < nodesIds.size(); i++){
        vector<size_t> nodesStates;
        if (tree_->isRoot(nodesIds[i])){
            nodesStates.push_back(simResult->getRootAncestralState());
        }else{
            nodesStates.push_back(simResult->getAncestralState(nodesIds[i]));
        }
        ancestors[nodesIds[i]] = nodesStates;
    }
    if (ChromEvolOptions::resultsPathDir_ == "none"){
        printTreeWithStates(*tree_, ancestors, ChromEvolOptions::resultsPathDir_);
    }else{
        const string outFilePath = ChromEvolOptions::resultsPathDir_ +"//"+ "simulatedDataAncestors.tree";
        printTreeWithStates(*tree_, ancestors, outFilePath);
    }
  
} */