// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Phyl/Model/Codon/YNGP_M2.h>
#include <Bpp/Phyl/Model/Codon/RELAX.h>
#include <Bpp/Phyl/Model/FrequencySet/CodonFrequencySet.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Alphabet/CodonAlphabet.h>
#include <Bpp/Seq/GeneticCode/StandardGeneticCode.h>
#include <Bpp/Numeric/Function/ReparametrizationFunctionWrapper.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Phyl/Tree/TreeTemplate.h>
#include <Bpp/Phyl/Legacy/Likelihood/RHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Legacy/Likelihood/DRHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Legacy/Likelihood/RNonHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Legacy/Likelihood/DRNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Legacy/Likelihood/RASTools.h>
#include <Bpp/Phyl/PatternTools.h>
#include <Bpp/Phyl/Legacy/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Legacy/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <iostream>


using namespace bpp;
using namespace std;


void printModelParameters(const TreeLikelihoodInterface& tl)
{
  ParameterList parameters = tl.getParameters();
  for (size_t i = 0; i < parameters.size(); ++i)
  {
    ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
  }
  cout << "\n" << endl;
}


int main() 
{
/*    try
    {
        // process tree
        TreeTemplate<Node>* tree = TreeTemplateTools::parenthesisToTree("(((A:0.01, B:0.01):0.02,C:0.03):0.01,D:0.04);");
        Tree* ttree = dynamic_cast<Tree*>(tree);

        // process sequence data
        map<string, string> alphabetParams;
        alphabetParams["alphabet"] = "Codon(letter=DNA)";
        alphabetParams["genetic_code"] = "Standard";
        const Alphabet* nucAlphabet = SequenceApplicationTools::getAlphabet(alphabetParams, "", false);
        const CodonAlphabet* alphabet = dynamic_cast<const CodonAlphabet*>(nucAlphabet);
        unique_ptr<GeneticCode> gCode;
        gCode.reset(SequenceApplicationTools::getGeneticCode(alphabet->getNucleicAlphabet(), "Standard"));
        VectorSiteContainer sites(alphabet);
        sites.addSequence(BasicSequence("A", "AAATGGCTGTGCACGTCT", alphabet));
        sites.addSequence(BasicSequence("B", "AACTGGATCTGCATGTCT", alphabet));
        sites.addSequence(BasicSequence("C", "ATCTGGACGTGCACGTGT", alphabet));
        sites.addSequence(BasicSequence("D", "CAACGGGAGTGCGCCTAT", alphabet));

        // set partition A and feed it to the RELAX model with k=1
        map<string,string> params;
        params["model1"] = "RELAX(kappa=2.0,p=0.1,omega1=1.0,omega2=2.0,k=1.0,theta1=0.5,theta2=0.8,Frequency=F0)";
        params["model2"] = "RELAX(kappa=RELAX.kappa_1,p=RELAX.p_1,omega1=RELAX.omega1_1,omega2=RELAX.omega2_1,theta1=RELAX.theta1_1,theta2=RELAX.theta2_1,Frequency=F0,k=1.0)";
        params["nonhomogeneous"]="general";
        params["nonhomogeneous.number_of_models"] = "2";
        params["nonhomogeneous.stationarity"] = "yes";
        params["site.number_of_paths"] = "2";                               // the 3rd path mapping omega3 in the branches under chatacter states 0 and 1 is imlies the the other two paths
        params["site.path1"] = "model1[YN98.omega_1]&model2[YN98.omega_1]"; // map omega1 in the branches under character state 0 (=model1) to omega1 in the branches under character state 1 (=model2) 
        params["site.path2"] = "model1[YN98.omega_2]&model2[YN98.omega_2]"; // these to complement the path of omega2
        params["model1.nodes_id"] = "0";
        params["model2.nodes_id"] = "1,2,3,4,5";
        SubstitutionModelSet* RELAXModel = PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet, gCode.get(), dynamic_cast<const SiteContainer*>(&sites), params);
		MixedSubstitutionModelSet* RELAXModel_1 = dynamic_cast<MixedSubstitutionModelSet*>(RELAXModel);

        // create likelihood function
        ConstantRateDistribution* rdist = new ConstantRateDistribution();
        RNonHomogeneousMixedTreeLikelihood* RELAXTreeLikelihood_1 = new RNonHomogeneousMixedTreeLikelihood(*ttree, dynamic_cast<const SiteContainer&>(sites), RELAXModel_1, rdist, true, false);
        RELAXTreeLikelihood_1->initialize();
        double RELAXLogLikelihood_1 = -1*RELAXTreeLikelihood_1->getValue();

        // set partition 2 -> make sure likelihood has not changed
        params["model1.nodes_id"] = "1,2,3,4,5";
        params["model2.nodes_id"] = "0";
        MixedSubstitutionModelSet* RELAXModel_2 = dynamic_cast<MixedSubstitutionModelSet*>(PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet, gCode.get(), dynamic_cast<const SiteContainer*>(&sites), params));
        RNonHomogeneousMixedTreeLikelihood* RELAXTreeLikelihood_2 = new RNonHomogeneousMixedTreeLikelihood(*ttree, dynamic_cast<const SiteContainer&>(sites), RELAXModel_2, rdist, true, false);
        RELAXTreeLikelihood_2->initialize();
        double RELAXLogLikelihood_2 = -1*RELAXTreeLikelihood_2->getValue();
        if (abs(RELAXLogLikelihood_1 - RELAXLogLikelihood_2) > 0.001)
        {
            cout << "Error! different likelihood is computed in RELAX for different partitions when k=1" << endl;
            return 1;
        }

        // make sure you receive the same likelihood as YNGP_M2 (simple site model)
        map<string,string> m2params;
        m2params["model"] = "YNGP_M2(kappa=2.0,omega0=0.1,omega2=2.0,theta1=0.5,theta2=0.8,Frequency=F0)";
        m2params["nonhomogeneous"] = "no";
        TransitionModel* M2Model = PhylogeneticsApplicationTools::getTransitionModel(alphabet, gCode.get(), dynamic_cast<const SiteContainer*>(&sites), m2params);
        RHomogeneousMixedTreeLikelihood* M2TreeLikelihood = new RHomogeneousMixedTreeLikelihood(*ttree, dynamic_cast<const SiteContainer&>(sites), M2Model, rdist, true, false);
        M2TreeLikelihood->initialize();
        double M2LogLikelihood = -1*M2TreeLikelihood->getValue();
        if (abs(RELAXLogLikelihood_1 - M2LogLikelihood) > 0.001)
        {
            cout << "Error! RELAX when k=1 yields different likelihood than M2 model" << endl;
			cout << "RELAX Log Likelihood: " << RELAXLogLikelihood_1 << endl;
			printModelParameters(RELAXTreeLikelihood_1);
			cout << "M2 Log Likelihood: " << M2LogLikelihood << endl;
			printModelParameters(M2TreeLikelihood);
            //return 1;
        }
        
        // set k to 2 -> fit two YNGP_M2 copies with the induced omega values and make sure the smae likelihood is obtained
        RELAXTreeLikelihood_2->setParameterValue("RELAX.k_2", 2);
		RELAXTreeLikelihood_2->computeTreeLikelihood();
        // make sure that updating other parameters except for k is done sucessfully 
        double RELAXLogLikelihood_3 = -1*RELAXTreeLikelihood_2->getValue();

        params["model1"] = "YNGP_M2(kappa=2.0,omega0=0.1,omega2=2.0,theta1=0.5,theta2=0.8,Frequency=F0)";
        params["model2"] = "YNGP_M2(kappa=YNGP_M2.kappa_1,omega0=0.01,omega2=4,theta1=YNGP_M2.theta1_1,theta2=YNGP_M2.theta2_1,Frequency=F0)";
        MixedSubstitutionModelSet* DoubleM2Model = dynamic_cast<MixedSubstitutionModelSet*>(PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet, gCode.get(), dynamic_cast<const SiteContainer*>(&sites), params));
        RNonHomogeneousMixedTreeLikelihood* DoubleM2TreeLikelihood = new RNonHomogeneousMixedTreeLikelihood(*ttree, dynamic_cast<const SiteContainer&>(sites), DoubleM2Model, rdist, true, false);
        DoubleM2TreeLikelihood->initialize();
        double DoubleM2LogLikelihood = -1*DoubleM2TreeLikelihood->getValue(); 
        if (abs(RELAXLogLikelihood_3 - DoubleM2LogLikelihood) > 0.001)
        {
            cout << "Error! RELAX yields different likelihood from two copies of YNGP_M2 that produce the same BG and FG as RELAX" << endl;
			cout << "RELAX Log Likelihood: " << RELAXLogLikelihood_3 << endl;
			printModelParameters(RELAXTreeLikelihood_2);
			cout << "M2 Log Likelihood: " << DoubleM2LogLikelihood << endl;
			printModelParameters(DoubleM2TreeLikelihood);
            return 1;
        }  

        // free resources
        delete rdist;
        delete RELAXTreeLikelihood_1;
        delete RELAXModel_1;
        delete RELAXTreeLikelihood_2;
        delete RELAXModel_2;
        delete M2Model;
        delete M2TreeLikelihood;
        delete DoubleM2Model;
        delete DoubleM2TreeLikelihood;
    }
    catch (exception & e)
    {
        cout << e.what() << endl;
        return 1;
    }
*/
  return 0;
}
