// From bpp-core:
#include <Bpp/Version.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Numeric/Function/BrentOneDimension.h>
#include <Bpp/Numeric/Function/PowellMultiDimensions.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Io/BppOSequenceReaderFormat.h>
#include <Bpp/Seq/Io/BppOAlignmentReaderFormat.h>
#include <Bpp/Seq/Io/BppOSequenceWriterFormat.h>
#include <Bpp/Seq/Io/BppOAlignmentWriterFormat.h>


// From bpp-phyl:
#include <Bpp/Phyl/Tree/PhyloTree.h>
#include <Bpp/Phyl/Tree/PhyloNode.h>
#include <Bpp/Phyl/Legacy/Likelihood/RHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Legacy/Likelihood/DRHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Legacy/Likelihood/RNonHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Legacy/Likelihood/DRNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Legacy/Likelihood/RASTools.h>
#include <Bpp/Phyl/Legacy/PatternTools.h>
#include <Bpp/Phyl/Legacy/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Legacy/OptimizationTools.h>
#include <Bpp/Phyl/Model/TwoParameterBinarySubstitutionModel.h>
#include <Bpp/Phyl/Model/Protein/CoalaCore.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/FrequencySet/MvaFrequencySet.h>
#include <Bpp/Phyl/Model/FrequencySet/FrequencySet.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Io/BppOFrequencySetFormat.h>
#include <Bpp/Phyl/Mapping/StochasticMapping.h>
#include <Bpp/Phyl/Legacy/Likelihood/JointLikelihoodFunction.h>

// From the STL:
#include <iostream>
#include <iomanip>
#include <limits>
#include <map>
#include <vector>
//#include <boost/math/distributions/chi_squared.hpp>


using namespace bpp;
using namespace std;

typedef vector<vector<double>> VVDouble;
typedef vector<double> VDouble;
typedef unsigned int uint;

/******************************************************************************/
/**************************** Auxiliary functions *****************************/
/******************************************************************************/

shared_ptr<const CodonAlphabet> getCodonAlphabet()
{
  map<string, string> alphabetSettings;
  alphabetSettings["alphabet"] = "Codon(letter=DNA)";
  shared_ptr<const Alphabet> alphabet = SequenceApplicationTools::getAlphabet(alphabetSettings, "", false); 
  auto codonAlphabet = dynamic_pointer_cast<const CodonAlphabet>(alphabet);
  return codonAlphabet;
}

/******************************************************************************/

shared_ptr<TransitionModelInterface> setCharacterModel(BppApplication* bppml,
    const VectorSiteContainer& charData,
    shared_ptr<const BinaryAlphabet> alphabet,
    shared_ptr<Tree> tree,
    shared_ptr<DRTreeParsimonyScore> mpData)
{
  // create the model
  // TO DO: ADD HERE PROCESSING OF INITIAL CHARACTER MODEL PARAMETERS AND PADD TO CONSTRUCTOR
  // extract the user initial value of k for potential later use
  double init_mu = ApplicationTools::getDoubleParameter("character_model.mu", bppml->getParams(), 1);
  double init_pi0 = ApplicationTools::getDoubleParameter("character_model.pi0", bppml->getParams(), 0.5);
  auto model = make_shared<TwoParameterBinarySubstitutionModel>(alphabet, init_mu, init_pi0);

  // compute the maximum parsimony score and set the lower and upper bounds on mu (=rate) as mp/tree_size, 2*mp/tree_size 
  VDouble treeBranches = tree->getBranchLengths();
  double treeSize = 0;
  for (size_t i = 0; i < treeBranches.size(); ++i)
  {
    treeSize += treeBranches[i];
  }
  double characterMuLb = mpData->getScore()/treeSize;
  double characterMuUb = 4*characterMuLb; // used 4*lb instead of 2*lb because factor of 2 is not enough (optimization converges to upper bound)
 
  // set the initial values of the model
  if (!ApplicationTools::getBooleanParameter("character_model.set_initial_parameters", bppml->getParams(), true, "", true, false))
  {
  // set the value of mu to be the middle of the interval
  model->setParameterValue(string("mu"),(characterMuLb + characterMuUb) / 2);
  // estimate the initial Frequencies as observedPseudoCount with pseudocount as 1 to avoid possible case of frequency = 0
  model->setFreqFromData(charData, 1); // the second argument stands for pesudocount 1
  }
  else
  {
    double mu = ApplicationTools::getDoubleParameter("character_model.mu", bppml->getParams(), 10);
    model->setParameterValue(string("mu"), mu);      
    double pi0 = ApplicationTools::getDoubleParameter("character_model.pi0", bppml->getParams(), 0.5);
    map<int,double>Frequencies;
    Frequencies[0] = pi0;
    Frequencies[1] = 1-pi0;
    model->setFreq(Frequencies);
  }

  return dynamic_pointer_cast<TransitionModelInterface>(model);
  
}

/******************************************************************************/

void giveNamesToInternalNodes(Tree* tree)
{
  TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(tree);
  vector<Node*> nodes = ttree->getNodes();
  for (size_t i=0; i<nodes.size(); ++i) {
    if (!nodes[i]->hasName())
      nodes[i]->setName("_baseInternal_" + TextTools::toString(nodes[i]->getId()));
  }  
}

/******************************************************************************/

void setMpPartition(BppApplication* bppml,
    shared_ptr<DRTreeParsimonyScore> mpData,
    const VectorSiteContainer& characterData,
    shared_ptr<TransitionModelInterface> characterModel,
    shared_ptr<Tree> tree)
{
//  mpData->computeSolution();
  const Tree& solution = mpData->topology();
  vector <const Node*> nodes = (dynamic_cast<const TreeTemplate<Node>&>(solution)).getNodes();
  
  // set assignment to branches
  string character0NodesIds, character1NodesIds = "";
  for (size_t i=0; i<nodes.size(); ++i)
  {
  if (!tree->isRoot(nodes[i]->getId()))
  {
    int nodeState = dynamic_cast<const BppInteger*>(nodes[i]->getNodeProperty("state"))->getValue();
    if (nodeState == 0)
    {
    character0NodesIds = character0NodesIds + TextTools::toString(nodes[i]->getId()) + ","; 
    } 
    else 
    {
    character1NodesIds = character1NodesIds + TextTools::toString(nodes[i]->getId()) + ","; 
    }
  }
  } 
  bppml->getParam("model1.nodes_id") = character0NodesIds.substr(0,character0NodesIds.length()-1);
  bppml->getParam("model2.nodes_id") = character1NodesIds.substr(0,character1NodesIds.length()-1);  
}

/******************************************************************************/

shared_ptr<MixedSubstitutionModelSet> setSequenceModel(BppApplication* bppml,
    shared_ptr<const VectorSiteContainer> codon_data,
    shared_ptr<const CodonAlphabet> codonAlphabet,
    shared_ptr<DRTreeParsimonyScore> mpData,
    shared_ptr<const VectorSiteContainer> characterData,
    shared_ptr<TransitionModelInterface> characterModel,
    shared_ptr<Tree> tree,
    shared_ptr<const GeneticCode> gCode)
{
  bppml->getParam("model1") = "RELAX(kappa=1,p=0.1,omega1=1.0,omega2=2.0,theta1=0.5,theta2=0.8,k=1,Frequencies=F3X4,initFreqs=observed,initFreqs.observedPseudoCount=1)";
  bppml->getParam("model2") = "RELAX(kappa=RELAX.kappa_1,p=RELAX.p_1,omega1=RELAX.omega1_1,omega2=RELAX.omega2_1,theta1=RELAX.theta1_1,theta2=RELAX.theta2_1,k=1,1_Full.theta=RELAX.1_Full.theta_1,1_Full.theta1=RELAX.1_Full.theta1_1,1_Full.theta2=RELAX.1_Full.theta2_1,2_Full.theta=RELAX.2_Full.theta_1,2_Full.theta1=RELAX.2_Full.theta1_1,2_Full.theta2=RELAX.2_Full.theta2_1,3_Full.theta=RELAX.3_Full.theta_1,3_Full.theta1=RELAX.3_Full.theta1_1,3_Full.theta2=RELAX.3_Full.theta2_1,Frequencies=F3X4,initFreqs=observed,initFreqs.observedPseudoCount=1)";
  bppml->getParam("nonhomogeneous")="general";
  bppml->getParam("nonhomogeneous.number_of_models") = "2";
  bppml->getParam("nonhomogeneous.stationarity") = "yes"; // constrain root Frequencies to be the same as stationary (since RELAX is a time reversible model, this should not cause issues)
  
  // set likelihood computation to restrict the same selective regime for each site along the phylogeny
  bppml->getParam("site.number_of_paths") = "2";                 // the 3rd path mapping omega3 in the branches under chatacter states 0 and 1 is imlies the the other two paths
  bppml->getParam("site.path1") = "model1[YN98.omega_1]&model2[YN98.omega_1]"; // map omega1 in the branches under character state 0 (=model1) to omega1 in the branches under character state 1 (=model2) 
  bppml->getParam("site.path2") = "model1[YN98.omega_2]&model2[YN98.omega_2]"; // do the same for omega2
  
//   string codeDesc = ApplicationTools::getStringParameter("genetic_code", bppml->getParams(), "Standard", "", true, true);
//   GeneticCode* gCode = SequenceApplicationTools::getGeneticCode(codonAlphabet->shareNucleicAlphabet(), codeDesc);
  
//   // set initial partition, based on maximum parsimony
//   setMpPartition(bppml, mpData, characterData, characterModel, tree); // the partition is set on tree
//   // create the set of models
  shared_ptr<SubstitutionModelSet> initialModelSet = PhylogeneticsApplicationToolsOld::getSubstitutionModelSet(codonAlphabet, gCode, codon_data, bppml->getParams());
  auto modelSet = dynamic_pointer_cast<MixedSubstitutionModelSet>(initialModelSet);
  return modelSet;
}

/******************************************************************************/
/*********************************** Main *************************************/
/******************************************************************************/

int main(int args, char** argv)
{
  try
  {

  /* process input from params file */
  int argNum = 1;
  char* argVals[1];
  argVals[0] = (char *)"";
  BppApplication bpp(argNum, argVals, "bpp");

  // process tree
  shared_ptr<TreeTemplate<Node>> ttree = TreeTemplateTools::parenthesisToTree("(((A:1,B:1):1,C:1):1,D:3);");
  auto tree = dynamic_pointer_cast<Tree>(ttree); 
  
  // process character data
  auto balpha = make_shared<const BinaryAlphabet>();
  auto balpha2 = dynamic_pointer_cast<const Alphabet>(balpha);
  auto charData = make_shared<VectorSiteContainer>(balpha);
  auto seq1A = make_unique<Sequence>("A", "0", balpha2);
  charData->addSequence("A", seq1A);
  auto seq1B = make_unique<Sequence>("B", "0", balpha2);
  charData->addSequence("B", seq1B);
  auto seq1C = make_unique<Sequence>("C", "0", balpha2);
  charData->addSequence("C", seq1C);
  auto seq1D = make_unique<Sequence>("D", "1", balpha2);
  charData->addSequence("D", seq1D);

  //compute the maximum parsimony for the purpose of setting bounds on the rate parameter of the character model and an intial tree partition for the starting point
  auto mpData = make_shared<DRTreeParsimonyScore>(ttree, charData); 

  // set the character model
  auto charModel = setCharacterModel(&bpp, *charData, balpha, tree, mpData);

  // process sequence data
  auto calpha = getCodonAlphabet();
  auto calpha2 = dynamic_pointer_cast<const Alphabet>(calpha);
  auto seqData = make_shared<VectorSiteContainer>(calpha);
  auto seq2A = make_unique<Sequence>("A", "AAATGGCTGTGCACGTCT", calpha2);
  seqData->addSequence("A", seq2A);
  auto seq2B = make_unique<Sequence>("B", "AACTGGATCTGCATGTCT", calpha2);
  seqData->addSequence("B", seq2B);
  auto seq2C = make_unique<Sequence>("C", "ATCTGGACGTGCACGTGT", calpha2);
  seqData->addSequence("C", seq2C);
  auto seq2D = make_unique<Sequence>("D", "CAACGGGAGTGCGCCTAT", calpha2);
  seqData->addSequence("D", seq2D);

  string codeDesc = ApplicationTools::getStringParameter("genetic_code", bpp.getParams(), "Standard", "", true, true);
  shared_ptr<GeneticCode> gCode = SequenceApplicationTools::getGeneticCode(calpha->getNucleicAlphabet(), codeDesc);
  // unique_ptr<GeneticCode> gCode;
  // gCode.reset(SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc));


  // set the sequence model
  auto seqModel = setSequenceModel(&bpp, seqData, calpha, mpData, charData, charModel, tree, gCode);

  // create joint likelihood function instance
  auto rDist = make_shared<ConstantRateDistribution>();
  auto jlf = make_shared<JointLikelihoodFunction>(&bpp, tree, charData, charModel, seqData, seqModel, rDist);
  jlf->setOptimizationScope(JointLikelihoodFunction::OptimizationScope(0)); // no optimization of sequece parameters is required at these stage
  
  // execute computation under the null hypothesis
  cout << "\ncomputing null model likelihood...\n" << endl;
  jlf->setHypothesis(JointLikelihoodFunction::Hypothesis(0));
  ParameterList params;
  jlf->fireParameterChanged(params);
  double nullLoglL_1 = jlf->getValue();

  // execute computation under the alternative hypothesis
  cout << "\ncomputing alternative model likelihood...\n" << endl;
  jlf->setHypothesis(JointLikelihoodFunction::Hypothesis(1));
  jlf->fireParameterChanged(params);
  double alternativeLoglL_1 = jlf->getValue();
  double kVal = jlf->getParameterValue("RELAX.k_2");
  if (kVal != 1 && nullLoglL_1 == alternativeLoglL_1)
  {
    cerr << "Error! the joint likelihood of the null and alternative hypotheses are the same despite k not being 1" << endl;
    return 1;
  }

  // alter the value of mu
  cout << "\n\n**** altering the value of mu: should affect both the alternative and the null, in different ways ****\n\n" << endl;
  jlf->setParameterValue("TwoParameterBinary.mu", 2);
  
  // execute computation under the null hypothesis
  cout << "\ncomputing null model likelihood...\n" << endl;
  jlf->setHypothesis(JointLikelihoodFunction::Hypothesis(0));
  jlf->fireParameterChanged(params);
  double nullLoglL_2 = jlf->getValue();
  if (nullLoglL_1 == nullLoglL_2)
  {
    cerr << "Error! the joint likelihood of the null model remains the same upon alteration of the value of mu\n" << endl;
    return 1;
  }

  // execute computation under the alternative hypothesis
  cout << "\ncomputing alternative model likelihood...\n" << endl;
  jlf->setHypothesis(JointLikelihoodFunction::Hypothesis(1)); // caused change in the value of mu - need to check how and fix the bug
  jlf->fireParameterChanged(params);
  double alternativeLoglL_2 = jlf->getValue();
  if (alternativeLoglL_1 == alternativeLoglL_2)
  {
    cerr << "Error! the joint likelihood of the alternative model remains the same upon alteration of the value of mu" << endl;
    return 1;
  }

  // alter the value of k
  cout << "\n\n**** altering the value of k: should only affect the alternative**** \n\n" << endl;
  jlf->setParameterValue("RELAX.k_2", 2);

  // execute computation under the null hypothesis
  cout << "\ncomputing null model likelihood...\n" << endl;
  jlf->setHypothesis(JointLikelihoodFunction::Hypothesis(0));
  jlf->fireParameterChanged(params);
  double nullLoglL_3 = jlf->getValue();
  if ((nullLoglL_2 - nullLoglL_3) > 0.0001)
  {
    cerr << "Error! the joint likelihood of the null model is altered upon modification of the value of k" << endl;
    return 1;
  }

  // execute computation under the alternative hypothesis
  cout << "\ncomputing alternative model likelihood...\n" << endl;
  jlf->setHypothesis(JointLikelihoodFunction::Hypothesis(1));
  jlf->fireParameterChanged(params);
  double alternativeLoglL_3 = jlf->getValue();
  if (alternativeLoglL_2 == alternativeLoglL_3)
  {
    cerr << "Error! the joint likelihood of the alternative model is not affected by change of the value of k" << endl;
    return 1;
  }

  // likelihood of a single site is NOT the multiplicity of the likelihood of the sequence site and the one of the character model, but we can regard it as it were for the purpose of sitewise comparison to the null TraitRELAX model
  // test likelihood computation by site
  double totalLikelihood = jlf->getLikelihood();
  vector<double> likelihoodBySite = jlf->getLikelihoodForEachSite();
  double multOverSites = 1.;
  for (size_t s=0; s<likelihoodBySite.size(); ++s)
  {
    multOverSites = multOverSites + likelihoodBySite[s]; // + because we are dealing with log liklihood values
  }
  if ((totalLikelihood-multOverSites) > 0.0001)
  {
    cerr << "Error! the joint likelihood computation by site returns a different sum from the overall joint likelihood" << endl;
		cerr << "joint likelihood = " << totalLikelihood << endl;
		cerr << " sum over all = " << multOverSites << endl;
    return 1;
  }

  }

  catch (exception& e)
  {
  cout << e.what() << endl;
  return 1;
  }
  return 0;
}  
