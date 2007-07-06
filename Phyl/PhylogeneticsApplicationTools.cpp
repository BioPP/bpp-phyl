//
// File: PhylogeneticsApplicationTools.cpp
// Created by: Julien Dutheil
// Created on: Fri Oct 21 16:49 2005
// from old file ApplicationTools.cpp created on Sun Dec 14 09:36:26 2003
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004)

This software is a computer program whose purpose is to provide classes
for phylogenetic data analysis.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#include "PhylogeneticsApplicationTools.h"
#include "models"
#include "OptimizationTools.h"
#include "Tree.h"
#include "Newick.h"
#include "HomogeneousTreeLikelihood.h"

// From Utils:
#include <Utils/FileTools.h>
#include <Utils/TextTools.h>
#include <Utils/ApplicationTools.h>
#include <Utils/StringTokenizer.h>

// From NumCalc:
#include <NumCalc/ConstantDistribution.h>
#include <NumCalc/GammaDiscreteDistribution.h>
#include <NumCalc/optimizers>

// From the STL:
#include <fstream>
#include <iomanip>
using namespace std;

/******************************************************************************/

TreeTemplate<Node> * PhylogeneticsApplicationTools::getTree(
	map<string, string> & params,
	const string & suffix,
	bool suffixIsOptional,
	bool verbose) throw (Exception)
{
	string treeFilePath = ApplicationTools::getAFilePath("tree.file", params, true, true, suffix, suffixIsOptional);
	
	//Read the tree file:
	Newick newick(true);
	TreeTemplate<Node> * tree = dynamic_cast<TreeTemplate<Node> *>(newick.read(treeFilePath));
	if(verbose) ApplicationTools::displayResult("Tree file", treeFilePath);
	return tree;
}

/******************************************************************************/

void PhylogeneticsApplicationTools::printInputTreeHelp()
{
  if(!ApplicationTools::message) return;
	*ApplicationTools::message << "Input tree parameters:" << endl;
	*ApplicationTools::message << "tree.file                     | file where to write the tree" << endl;
  *ApplicationTools::message << "______________________________|___________________________________________" << endl;
}

/******************************************************************************/

SubstitutionModel * PhylogeneticsApplicationTools::getSubstitutionModel(
	const Alphabet * alphabet,
	const SiteContainer * data,
	map<string, string> & params,
	const string & suffix,
	bool suffixIsOptional,
	bool verbose) throw (Exception)
{
	SubstitutionModel * model = NULL;

	if(alphabet->getAlphabetType() == "DNA alphabet"
	|| alphabet->getAlphabetType() == "RNA alphabet")
  {
	  string modelName = ApplicationTools::getStringParameter("model", params, "JCnuc", suffix, suffixIsOptional);
		const NucleicAlphabet * alpha = dynamic_cast<const NucleicAlphabet *>(alphabet);
    if(verbose) ApplicationTools::displayResult("Substitution model", modelName);

		if(modelName == "GTR")
    {
			double piA = 0.25, piC = 0.25, piG = 0.25, piT = 0.25;
			double a = ApplicationTools::getDoubleParameter("a", params, 1, suffix, suffixIsOptional);
			double b = ApplicationTools::getDoubleParameter("b", params, 1, suffix, suffixIsOptional);
			double c = ApplicationTools::getDoubleParameter("c", params, 1, suffix, suffixIsOptional);
			double d = ApplicationTools::getDoubleParameter("d", params, 1, suffix, suffixIsOptional);
			double e = ApplicationTools::getDoubleParameter("e", params, 1, suffix, suffixIsOptional);
			bool useObsFreq = ApplicationTools::getBooleanParameter("model.use_observed_freq", params, false, suffix, suffixIsOptional);
			if(useObsFreq && data != NULL)
      {
				model = new GTR(alpha, a, b, c, d, e);
				dynamic_cast<GTR *>(model)->setFreqFromData(*data);
				piA = model->getParameterValue("piA");
				piC = model->getParameterValue("piC");
				piG = model->getParameterValue("piG");
				piT = model->getParameterValue("piT");
			}
      else
      {
				piA = ApplicationTools::getDoubleParameter("piA", params, 0.25, suffix, suffixIsOptional);
				piC = ApplicationTools::getDoubleParameter("piC", params, 0.25, suffix, suffixIsOptional);
				piG = ApplicationTools::getDoubleParameter("piG", params, 0.25, suffix, suffixIsOptional);
				piT = ApplicationTools::getDoubleParameter("piT", params, 0.25, suffix, suffixIsOptional);
				if( fabs(1-(piA + piT + piG + piC)) > 0.00000000000001 )
        {
					throw Exception("Equilibrium base frequencies must equal 1. Aborting...");
				}
				model = new GTR(alpha, a, b, c, d, e, piA, piC, piG, piT);
			}
			if(verbose)
      {
				ApplicationTools::displayResult("a"  , TextTools::toString(a));
				ApplicationTools::displayResult("b"  , TextTools::toString(b));
				ApplicationTools::displayResult("c"  , TextTools::toString(c));
				ApplicationTools::displayResult("d"  , TextTools::toString(d));
				ApplicationTools::displayResult("e"  , TextTools::toString(e));
				ApplicationTools::displayResult("piA", TextTools::toString(piA));
				ApplicationTools::displayResult("piC", TextTools::toString(piC));
				ApplicationTools::displayResult("piG", TextTools::toString(piG));
				ApplicationTools::displayResult("piT", TextTools::toString(piT));
			}
		}
    else if(modelName == "TN93")
    {
			double piA = 0.25, piC = 0.25, piG = 0.25, piT = 0.25;
			double kappa1 = ApplicationTools::getDoubleParameter("kappa1", params, 2, suffix, suffixIsOptional);
			double kappa2 = ApplicationTools::getDoubleParameter("kappa2", params, 2, suffix, suffixIsOptional);
			bool useObsFreq = ApplicationTools::getBooleanParameter("model.use_observed_freq", params, false, suffix, suffixIsOptional);
			if(useObsFreq && data != NULL)
      {
				model = new TN93(alpha, kappa1, kappa2);
				dynamic_cast<TN93 *>(model)->setFreqFromData(*data);
				piA = model->getParameterValue("piA");
				piC = model->getParameterValue("piC");
				piG = model->getParameterValue("piG");
				piT = model->getParameterValue("piT");
			}
      else
      {
				piA = ApplicationTools::getDoubleParameter("piA", params, 0.25, suffix, suffixIsOptional);
				piC = ApplicationTools::getDoubleParameter("piC", params, 0.25, suffix, suffixIsOptional);
				piG = ApplicationTools::getDoubleParameter("piG", params, 0.25, suffix, suffixIsOptional);
				piT = ApplicationTools::getDoubleParameter("piT", params, 0.25, suffix, suffixIsOptional);
				if( fabs(1-(piA + piT + piG + piC)) > 0.00000000000001 )
        {
					throw Exception("Equilibrium base frequencies must equal 1. Aborting...");
					exit(-1);
				}
				model = new TN93(alpha, kappa1, kappa2, piA, piC, piG, piT);
			}
			if(verbose)
      {
				ApplicationTools::displayResult("kappa1", TextTools::toString(kappa1));
				ApplicationTools::displayResult("kappa2", TextTools::toString(kappa2));
				ApplicationTools::displayResult("piA"   , TextTools::toString(piA));
				ApplicationTools::displayResult("piC"   , TextTools::toString(piC));
				ApplicationTools::displayResult("piG"   , TextTools::toString(piG));
				ApplicationTools::displayResult("piT"   , TextTools::toString(piT));
			}
		}
    else if(modelName == "HKY85")
    {
			double piA = 0.25, piC = 0.25, piG = 0.25, piT = 0.25;
			double kappa = ApplicationTools::getDoubleParameter("kappa", params, 2, suffix, suffixIsOptional);
			bool useObsFreq = ApplicationTools::getBooleanParameter("model.use_observed_freq", params, false, suffix, suffixIsOptional);
			if(useObsFreq && data != NULL)
      {
				model = new HKY85(alpha, kappa);
				dynamic_cast<HKY85 *>(model)->setFreqFromData(*data);
				piA = model->getParameterValue("piA");
				piC = model->getParameterValue("piC");
				piG = model->getParameterValue("piG");
				piT = model->getParameterValue("piT");
			}
      else
      {
				piA = ApplicationTools::getDoubleParameter("piA", params, 0.25, suffix, suffixIsOptional);
				piC = ApplicationTools::getDoubleParameter("piC", params, 0.25, suffix, suffixIsOptional);
				piG = ApplicationTools::getDoubleParameter("piG", params, 0.25, suffix, suffixIsOptional);
				piT = ApplicationTools::getDoubleParameter("piT", params, 0.25, suffix, suffixIsOptional);
				if( fabs(1-(piA + piT + piG + piC)) > 0.00000000000001 )
        {
					throw Exception("Equilibrium base frequencies must equal 1. Aborting...");
					exit(-1);
				}
				model = new HKY85(alpha, kappa, piA, piC, piG, piT);
			}
			if(verbose)
      {
				ApplicationTools::displayResult("kappa", TextTools::toString(kappa));
				ApplicationTools::displayResult("piA"  , TextTools::toString(piA));
				ApplicationTools::displayResult("piC"  , TextTools::toString(piC));
				ApplicationTools::displayResult("piG"  , TextTools::toString(piG));
				ApplicationTools::displayResult("piT"  , TextTools::toString(piT));
			}
		}
    else if(modelName == "F84")
    {
			double piA = 0.25, piC = 0.25, piG = 0.25, piT = 0.25;
			double kappa = ApplicationTools::getDoubleParameter("kappa", params, 2, suffix, suffixIsOptional);
			bool useObsFreq = ApplicationTools::getBooleanParameter("model.use_observed_freq", params, false, suffix, suffixIsOptional);
			if(useObsFreq && data != NULL)
      {
				model = new F84(alpha, kappa);
				dynamic_cast<F84 *>(model)->setFreqFromData(*data);
				piA = model->getParameterValue("piA");
				piC = model->getParameterValue("piC");
				piG = model->getParameterValue("piG");
				piT = model->getParameterValue("piT");
			}
      else
      {
				piA = ApplicationTools::getDoubleParameter("piA", params, 0.25, suffix, suffixIsOptional);
				piC = ApplicationTools::getDoubleParameter("piC", params, 0.25, suffix, suffixIsOptional);
				piG = ApplicationTools::getDoubleParameter("piG", params, 0.25, suffix, suffixIsOptional);
				piT = ApplicationTools::getDoubleParameter("piT", params, 0.25, suffix, suffixIsOptional);
				if( fabs(1-(piA + piT + piG + piC)) > 0.00000000000001 )
        {
					throw Exception("Equilibrium base frequencies must equal 1. Aborting...");
					exit(-1);
				}
				model = new F84(alpha, kappa, piA, piC, piG, piT);
			}
			if(verbose)
      {
				ApplicationTools::displayResult("kappa", TextTools::toString(kappa));
				ApplicationTools::displayResult("piA"  , TextTools::toString(piA));
				ApplicationTools::displayResult("piC"  , TextTools::toString(piC));
				ApplicationTools::displayResult("piG"  , TextTools::toString(piG));
				ApplicationTools::displayResult("piT"  , TextTools::toString(piT));
			}
		}
    else if(modelName == "T92")
    {
			double kappa = ApplicationTools::getDoubleParameter("kappa", params, 2, suffix, suffixIsOptional);
			double theta = 0.5;
			bool useObsFreq = ApplicationTools::getBooleanParameter("model.use_observed_freq", params, false, suffix, suffixIsOptional);
			if(useObsFreq && data != NULL)
      {
				model = new T92(alpha, kappa);
				dynamic_cast<T92 *>(model)->setFreqFromData(*data);
				theta = model->getParameterValue("theta");
			}
      else
      {
				theta = ApplicationTools::getDoubleParameter("theta", params, 0.5, suffix, suffixIsOptional);
				model = new T92(alpha, kappa, theta);
			}
			if(verbose)
      {
				ApplicationTools::displayResult("kappa", TextTools::toString(kappa));
				ApplicationTools::displayResult("theta", TextTools::toString(theta));
			}
		}
    else if(modelName == "K80")
    {
    	double kappa = ApplicationTools::getDoubleParameter("kappa", params, 2, suffix, suffixIsOptional);
			model = new K80(alpha, kappa);
  		if(verbose)
      {
				ApplicationTools::displayResult("kappa", TextTools::toString(kappa));
			}
  	}
    else if(modelName == "JCnuc")
    {
			model = new JCnuc(alpha);
		}
    else
    {
			throw Exception("Model '" + modelName + "' unknown.");
		}
	}
  else
  { 
    // Alphabet supposed to be proteic!
	  string modelName = ApplicationTools::getStringParameter("model", params, "JCprot", suffix, suffixIsOptional);
		const ProteicAlphabet * alpha = dynamic_cast<const ProteicAlphabet *>(alphabet);
		bool useObsFreq = ApplicationTools::getBooleanParameter("model.use_observed_freq", params, false, suffix, suffixIsOptional);
		if(modelName == "JCprot")
    {
			model = new JCprot(alpha);
		}
    else if(modelName == "DSO78")
    {
			model = new DSO78(alpha);
			if(useObsFreq && data != NULL) dynamic_cast<DSO78 *>(model)->setFreqFromData(*data);
		}
    else if(modelName == "JTT92")
    {
			model = new JTT92(alpha);
			if(useObsFreq && data != NULL) dynamic_cast<JTT92 *>(model)->setFreqFromData(*data);	
		}
    else if(modelName == "empirical")
    {
			string file = ApplicationTools::getAFilePath("model_empirical.file", params, true, true, suffix, true);
			model = new UserProteinSubstitutionModel(alpha, file);
			if(useObsFreq && data != NULL)
        dynamic_cast<UserProteinSubstitutionModel *>(model)->setFreqFromData(*data);	
	  }
    else
    {
			throw Exception("Model '" + modelName + "' unknown. Aborting...");
		}
 		if(verbose)
    {
			ApplicationTools::displayResult("Substitution model", modelName + (useObsFreq && (model != NULL) ? "-F" : ""));
		}
	}
  
	return model;
}

/******************************************************************************/

void PhylogeneticsApplicationTools::printSubstitutionModelHelp()
{
  if(!ApplicationTools::message) return;
  *ApplicationTools::message << "Substitution Model:" << endl;
	*ApplicationTools::message << "model               | Nucleotides (N): [JCnuc, K80, T92, F84, HKY85, TN93," << endl;
  *ApplicationTools::message << "                    | GTR]" << endl;
  *ApplicationTools::message << "                    | Proteins (P): [JCprot, DSO78, JTT92, empirical]" << endl;
  *ApplicationTools::message << "kappa               | kappa(N)  parameter in Q matrix" << endl;
  *ApplicationTools::message << "kappa1              | kappa1(N) parameter in Q matrix" << endl;
  *ApplicationTools::message << "kappa2              | kappa2(N) parameter in Q matrix" << endl;
  *ApplicationTools::message << "a,b,c,d,e,f         | GTR rates parameter in Q matrix" << endl;
  *ApplicationTools::message << "theta               | theta(N)  parameter in Q matrix" << endl;
  *ApplicationTools::message << "piA                 | piA(N)    parameter in Q matrix" << endl;
  *ApplicationTools::message << "piT                 | piT(N)    parameter in Q matrix" << endl;
  *ApplicationTools::message << "piC                 | piC(N)    parameter in Q matrix" << endl;
  *ApplicationTools::message << "piG                 | piG(N)    parameter in Q matrix" << endl;
	*ApplicationTools::message << "use_observed_freq   | (N,P) Tell if the observed frequencies must be used." << endl; 
	*ApplicationTools::message << "model_empirical.file| (P) The path toward data file to use (PAML format)." << endl; 
  *ApplicationTools::message << "____________________|_____________________________________________________" << endl;
}

/******************************************************************************/

MarkovModulatedSubstitutionModel * PhylogeneticsApplicationTools::getCovarionProcess(
  SubstitutionModel * model,
  DiscreteDistribution * rDist,
  map<string, string> & params,
	string suffix,
	bool suffixIsOptional,
	bool verbose) throw (Exception)
{
	string covarionName = ApplicationTools::getStringParameter("covarion", params, "none", suffix, suffixIsOptional);
  if(covarionName == "G2001")
  {
	  double nu = ApplicationTools::getDoubleParameter("covarion_G2001.nu", params, 1., suffix, suffixIsOptional);
		if(verbose)
    {
			ApplicationTools::displayResult("Covarion model" , covarionName);
			ApplicationTools::displayResult("nu"   , TextTools::toString(nu));
		}
    return new G2001(model, rDist, nu);
  }
  else if(covarionName == "TS98")
  {
	  double s1 = ApplicationTools::getDoubleParameter("covarion_TS98.s1", params, 1., suffix, suffixIsOptional);
	  double s2 = ApplicationTools::getDoubleParameter("covarion_TS98.s2", params, 1., suffix, suffixIsOptional);
		if(verbose)
    {
			ApplicationTools::displayResult("Covarion model" , covarionName);
			ApplicationTools::displayResult("s1"   , TextTools::toString(s1));
			ApplicationTools::displayResult("s2"   , TextTools::toString(s2));
		}
    return new TS98(model, s1, s2);
  }
  else
  {
	  throw Exception("Process unknown: " + covarionName + ".");
	}
}

/******************************************************************************/

void PhylogeneticsApplicationTools::printCovarionModelHelp()
{
  if(!ApplicationTools::message) return;
  *ApplicationTools::message << "Covarion Model:" << endl;
	*ApplicationTools::message << "covarion            | [none|G2001|TS98]" << endl;
  *ApplicationTools::message << "covarion_G2001.nu   | covarion rate parameter in G2001 model." << endl;
  *ApplicationTools::message << "covarion_TS98.s1    | covarion rate parameter in TS98 model." << endl;
  *ApplicationTools::message << "covarion_TS98.s2    | covarion rate parameter in TS98 model." << endl;
  *ApplicationTools::message << "____________________|_____________________________________________________" << endl;
}

/******************************************************************************/

DiscreteDistribution * PhylogeneticsApplicationTools::getRateDistribution(
	map<string, string> & params,
	string suffix,
	bool suffixIsOptional,
	bool verbose) throw (Exception)
{
	string distributionType = ApplicationTools::getStringParameter("rate_distribution", params, "constant", suffix, suffixIsOptional);
	DiscreteDistribution * rDist;
	if(distributionType == "constant")
  {
		rDist = new ConstantDistribution(1.);
		if(verbose)
    {
			ApplicationTools::displayResult("Rate distribution", distributionType);
		}
	}
  else if(distributionType == "gamma")
  {
		double alpha = ApplicationTools::getDoubleParameter("rate_distribution_gamma.alpha", params, 1., suffix, suffixIsOptional);
		int nbClasses = ApplicationTools::getIntParameter("rate_distribution.classes_number", params, 4, suffix, suffixIsOptional);

		if(alpha < 0)
    {
			throw Exception("Alpha parameter in gamma distribution of rates must be > 0, found " + TextTools::toString(alpha) + ".");
		}
    else
    {
			rDist = new GammaDiscreteDistribution(nbClasses, alpha);
			if(verbose)
      {
				ApplicationTools::displayResult("Rate distribution", distributionType);
				ApplicationTools::displayResult("Shape (alpha)", TextTools::toString((dynamic_cast<GammaDiscreteDistribution *>(rDist)) -> getParameterValue("alpha")));
				ApplicationTools::displayResult("Number of classes", TextTools::toString(rDist -> getNumberOfCategories()));
				for(unsigned int c = 0; c < rDist -> getNumberOfCategories(); c++)
        {
					ApplicationTools::displayResult("- Category " + TextTools::toString(c)
					+ "(rate = " + TextTools::toString(rDist -> getCategory(c))
					+ "), prob = ", TextTools::toString(rDist -> getProbability(c)));
				}
			}
		}
	}
  else
  {
		throw Exception("Distribution unknown: " + distributionType + ".");
	}

	return rDist;
}

/******************************************************************************/

void PhylogeneticsApplicationTools::printRateDistributionHelp()
{
  if(!ApplicationTools::message) return;
  *ApplicationTools::message << "Rate distribution parameters:" << endl;
 	*ApplicationTools::message << "rate_distribution               | uniform or gamma." << endl;
 	*ApplicationTools::message << "rate_distribution_gamma.alpha   | the gamma law's alpha parameter." << endl;
 	*ApplicationTools::message << "rate_distribution.classes_number| discrete approximation: number of" << endl;
  *ApplicationTools::message << "                                | categories (default to 4)." << endl;
  *ApplicationTools::message << "________________________________|_________________________________________" << endl;
}

/******************************************************************************/

TreeLikelihood * PhylogeneticsApplicationTools::optimizeParameters(
	TreeLikelihood * tl,
	map<string, string> & params,
	const string & suffix,
	bool suffixIsOptional,
	bool verbose)
	throw (Exception)
{
	bool optimize = ApplicationTools::getBooleanParameter("optimization", params, true, suffix, suffixIsOptional, false);
	if(!optimize) return tl;
	
	unsigned int optVerbose = ApplicationTools::getParameter<unsigned int>("optimization.verbose", params, 2, suffix, suffixIsOptional);
	
	string mhPath = ApplicationTools::getAFilePath("optimization.message_handler", params, false, false, suffix, suffixIsOptional);
	ostream * messageHandler = 
		(mhPath == "none") ? NULL :
			(mhPath == "std") ? &cout :
				new ofstream(mhPath.c_str(), ios::out);
	if(verbose) ApplicationTools::displayResult("Message handler", mhPath);

	string prPath = ApplicationTools::getAFilePath("optimization.profiler", params, false, false, suffix, suffixIsOptional);
	ostream * profiler = 
		(prPath == "none") ? NULL :
			(prPath == "std") ? &cout :
				new ofstream(prPath.c_str(), ios::out);
	if(profiler != NULL) (*profiler) << setprecision(20);
	if(verbose) ApplicationTools::displayResult("Profiler", prPath);

	bool scaleFirst = ApplicationTools::getBooleanParameter("optimization.scale_first", params, false, suffix, suffixIsOptional, false);
	if(scaleFirst)
  {
		// We scale the tree before optimizing each branch length separately:
		if(verbose) ApplicationTools::displayMessage("Scaling the tree before optimizing each branch length separately.");
		double tolerance = ApplicationTools::getDoubleParameter("optimization.scale_first.tolerance", params, .0001, suffix, suffixIsOptional, true);
		if(verbose) ApplicationTools::displayResult("Scaling tolerance", TextTools::toString(tolerance));
		int nbEvalMax = ApplicationTools::getIntParameter("optimization.scale_first.max_number_f_eval", params, 1000000, suffix, suffixIsOptional, true);
		if(verbose) ApplicationTools::displayResult("Scaling max # f eval", TextTools::toString(nbEvalMax));
		int n = OptimizationTools::optimizeTreeScale(
			tl,
			tolerance,
			nbEvalMax,
			messageHandler,
			profiler);
		if(verbose) ApplicationTools::displayMessage("Performed " + TextTools::toString(n) + " function evaluations.");
	}
	
	// Should I ignore some parameters?
	string paramListDesc = ApplicationTools::getStringParameter("optimization.ignore_parameter", params, "", suffix, suffixIsOptional, false);
	StringTokenizer st(paramListDesc, ",");
	while(st.hasMoreToken())
  {
		try
    {
      string param = st.nextToken();
      if(param == "BrLen")
      {
        vector<string> vs = tl->getBranchLengthsParameters().getParameterNames();
        for(unsigned int i = 0; i < vs.size(); i++)
        {
          dynamic_cast<AbstractHomogeneousTreeLikelihood *>(tl)->ignoreParameter(vs[i]);
        }
      }
      else dynamic_cast<AbstractHomogeneousTreeLikelihood *>(tl)->ignoreParameter(param);
		} 
    catch(ParameterNotFoundException pnfe)
    {
			ApplicationTools::displayError("Parameter '" + pnfe.getParameter() + "' not found, and so can't be ignored!");
		}
    catch(exception e)
    {
			ApplicationTools::displayError("DEBUB: ERROR!!! This functionality can only be used with HomogeneousTreeLikelihood for now.");
		}
	}
	
	unsigned int nbEvalMax = ApplicationTools::getParameter<unsigned int>("optimization.max_number_f_eval", params, 1000000, suffix, suffixIsOptional);
	if(verbose) ApplicationTools::displayResult("Max # ML evaluations", TextTools::toString(nbEvalMax));
	
	double tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance", params, .000001, suffix, suffixIsOptional);
	if(verbose) ApplicationTools::displayResult("Tolerance", TextTools::toString(tolerance));
	
  bool optimizeTopo = ApplicationTools::getBooleanParameter("optimization.topology", params, false, suffix, suffixIsOptional, false);
  if(verbose) ApplicationTools::displayResult("Optimize topology", optimizeTopo ? "yes" : "no");
  
  string method = ApplicationTools::getStringParameter("optimization.method", params, "DB", suffix, suffixIsOptional, false);
  string order  = ApplicationTools::getStringParameter("optimization.method.derivatives", params, "newton", suffix, suffixIsOptional, false);
  string optMethod;
  if(order == "gradient")
  {
    optMethod = OptimizationTools::OPTIMIZATION_GRADIENT;
  }
  else if(order == "newton")
  {
    optMethod = OptimizationTools::OPTIMIZATION_NEWTON;
  }
  else throw Exception("Unknown derivatives algorithm: '" + order + "'.");
  if(verbose) ApplicationTools::displayResult("Optimization method", method);
  if(verbose) ApplicationTools::displayResult("Algorithm used for derivable parameters", order);
	
  unsigned int n = 0;
  if(method == "DB")
  {
    //Uses Newton-Brent method:
	  
    unsigned int nstep = ApplicationTools::getParameter<unsigned int>("optimization.method_DB.nstep", params, 1, suffix, suffixIsOptional, false);
    if(optimizeTopo)
    {
      bool optNumFirst = ApplicationTools::getBooleanParameter("optimization.topology.numfirst", params, true, suffix, suffixIsOptional, false);
      unsigned int n   = ApplicationTools::getParameter<unsigned int>("optimization.topology.nstep", params, 1, "", true, false);
	    double tolBefore = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.before", params, 100, suffix, suffixIsOptional);
	    double tolDuring = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.during", params, 100, suffix, suffixIsOptional);
      tl = OptimizationTools::optimizeTreeNNI(
			    dynamic_cast<NNIHomogeneousTreeLikelihood *>(tl),
          optNumFirst, tolBefore, tolDuring, nbEvalMax, n, messageHandler, profiler, optVerbose, optMethod, nstep);
    }

	  if(verbose && nstep > 1) ApplicationTools::displayResult("# of precision steps", TextTools::toString(nstep));
    n = OptimizationTools::optimizeNumericalParameters(
		  dynamic_cast<AbstractHomogeneousTreeLikelihood *>(tl),
      NULL,
      nstep,
			tolerance,
			nbEvalMax,
			messageHandler,
			profiler,
			optVerbose,
      optMethod);		
  }
  else if(method == "fullD")
  {
    //Uses Newton-raphson alogrithm with numerical derivatives when required.
    
    if(optimizeTopo)
    {
      bool optNumFirst = ApplicationTools::getBooleanParameter("optimization.topology.numfirst", params, true, suffix, suffixIsOptional, false);
      unsigned int n   = ApplicationTools::getParameter<unsigned int>("optimization.topology.nstep", params, 1, "", true, false);
	    double tolBefore = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.before", params, 100, suffix, suffixIsOptional);
	    double tolDuring = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.during", params, 100, suffix, suffixIsOptional);
      tl = OptimizationTools::optimizeTreeNNI2(
			    dynamic_cast<NNIHomogeneousTreeLikelihood *>(tl),
          optNumFirst, tolBefore, tolDuring, nbEvalMax, n, messageHandler, profiler, optVerbose, optMethod);
    }

    n = OptimizationTools::optimizeNumericalParameters2(
			dynamic_cast<AbstractHomogeneousTreeLikelihood *>(tl),
      NULL,
			tolerance,
			nbEvalMax,
			messageHandler,
			profiler,
			optVerbose,
      optMethod);		   
  }
  else throw Exception("Unknown optimization method: " + method);

  string finalMethod = ApplicationTools::getStringParameter("optimization.final", params, "none", suffix, suffixIsOptional, false);
  Optimizer * finalOptimizer  = NULL;
  if(finalMethod == "none") {}
  else if(finalMethod == "simplex")
  {
    finalOptimizer = new DownhillSimplexMethod(tl);
  }
  else if(finalMethod == "powell")
  {
    finalOptimizer = new PowellMultiDimensions(tl);
  }
  else throw Exception("Unknown final optimization method: " + finalMethod);

  if(finalOptimizer)
  {
    if(verbose) ApplicationTools::displayResult("Final optimization step", finalMethod);
    finalOptimizer->setProfiler(profiler);
    finalOptimizer->setMessageHandler(messageHandler);
    finalOptimizer->setMaximumNumberOfEvaluations(nbEvalMax);
    finalOptimizer->getStopCondition()->setTolerance(tolerance);
    finalOptimizer->setVerbose(verbose);
    finalOptimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    finalOptimizer->init(tl->getParameters());
    finalOptimizer->optimize();
    n += finalOptimizer->getNumberOfEvaluations();
    delete finalOptimizer;
  }
  
	if(verbose) ApplicationTools::displayResult("Performed", TextTools::toString(n) + " function evaluations.");
  return tl;
}

/******************************************************************************/

void PhylogeneticsApplicationTools::optimizeParameters(
	ClockTreeLikelihood * tl,
	map<string, string> & params,
	const string & suffix,
	bool suffixIsOptional,
	bool verbose)
	throw (Exception)
{
	bool optimize = ApplicationTools::getBooleanParameter("optimization", params, true, suffix, suffixIsOptional, false);
	if(!optimize) return;
	
	unsigned int optVerbose = ApplicationTools::getParameter<unsigned int>("optimization.verbose", params, 2, suffix, suffixIsOptional);
	
	string mhPath = ApplicationTools::getAFilePath("optimization.message_handler", params, false, false, suffix, suffixIsOptional);
	ostream * messageHandler = 
		(mhPath == "none") ? NULL :
			(mhPath == "std") ? &cout :
				new ofstream(mhPath.c_str(), ios::out);
	if(verbose) ApplicationTools::displayResult("Message handler", mhPath);

	string prPath = ApplicationTools::getAFilePath("optimization.profiler", params, false, false, suffix, suffixIsOptional);
	ostream * profiler = 
		(prPath == "none") ? NULL :
			(prPath == "std") ? &cout :
				new ofstream(prPath.c_str(), ios::out);
	if(profiler != NULL) (*profiler) << setprecision(20);
	if(verbose) ApplicationTools::displayResult("Profiler", prPath);

	// Should I ignore some parameters?
	string paramListDesc = ApplicationTools::getStringParameter("optimization.ignore_parameter", params, "", suffix, suffixIsOptional, false);
	StringTokenizer st(paramListDesc, ",");
	while(st.hasMoreToken())
  {
		try
    {
      string param = st.nextToken();
      if(param == "BrLen")
      {
        vector<string> vs = tl->getBranchLengthsParameters().getParameterNames();
        for(unsigned int i = 0; i < vs.size(); i++)
        {
          dynamic_cast<AbstractHomogeneousTreeLikelihood *>(tl)->ignoreParameter(vs[i]);
        }
      }
      else dynamic_cast<AbstractHomogeneousTreeLikelihood *>(tl)->ignoreParameter(param);
		} 
    catch(ParameterNotFoundException pnfe)
    {
			ApplicationTools::displayError("Parameter '" + pnfe.getParameter() + "' not found, and so can't be ignored!");
		}
    catch(exception e)
    {
			ApplicationTools::displayError("DEBUB: ERROR!!! This functionality can only be used with HomogeneousTreeLikelihood for now.");
		}
	}
	
	unsigned int nbEvalMax = ApplicationTools::getParameter<unsigned int>("optimization.max_number_f_eval", params, 1000000, suffix, suffixIsOptional);
	if(verbose) ApplicationTools::displayResult("Max # ML evaluations", TextTools::toString(nbEvalMax));
	
	double tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance", params, .000001, suffix, suffixIsOptional);
	if(verbose) ApplicationTools::displayResult("Tolerance", TextTools::toString(tolerance));
	
  string method = ApplicationTools::getStringParameter("optimization.method", params, "DB", suffix, suffixIsOptional, false);
  string order  = ApplicationTools::getStringParameter("optimization.method.derivatives", params, "gradient", suffix, suffixIsOptional, false);
	string optMethod, derMethod;
  if(order == "gradient")
  {
    optMethod = OptimizationTools::OPTIMIZATION_GRADIENT;
  }
  else if(order == "newton")
  {
    optMethod = OptimizationTools::OPTIMIZATION_NEWTON;
  }
  else throw Exception("Option '" + order + "' is not known for 'optimization.method.derivatives'.");
  if(verbose) ApplicationTools::displayResult("Optimization method", method);
  if(verbose) ApplicationTools::displayResult("Algorithm used for derivable parameters", order);
	
  unsigned int n = 0;
  if(method == "DB")
  {
    //Uses Newton-Brent method:
	  unsigned int nstep = ApplicationTools::getParameter<unsigned int>("optimization.method_DB.nstep", params, 1, suffix, suffixIsOptional, false);
	  if(verbose && nstep > 1) ApplicationTools::displayResult("# of precision steps", TextTools::toString(nstep));
    n = OptimizationTools::optimizeNumericalParametersWithGlobalClock(
		  tl,
      NULL,
      nstep,
			tolerance,
			nbEvalMax,
			messageHandler,
			profiler,
			optVerbose,
      optMethod);		
  }
  else if(method == "fullD")
  {
    //Uses Newton-raphson alogrithm with numerical derivatives when required.
    n = OptimizationTools::optimizeNumericalParametersWithGlobalClock2(
			tl,
      NULL,
			tolerance,
			nbEvalMax,
			messageHandler,
			profiler,
			optVerbose,
      optMethod);		   
  }
  else throw Exception("Unknown optimization method: " + method);

  string finalMethod = ApplicationTools::getStringParameter("optimization.final", params, "none", suffix, suffixIsOptional, false);
  Optimizer * finalOptimizer  = NULL;
  if(finalMethod == "none") {}
  else if(finalMethod == "simplex")
  {
    finalOptimizer = new DownhillSimplexMethod(tl);
  }
  else if(finalMethod == "powell")
  {
    finalOptimizer = new PowellMultiDimensions(tl);
  }
  else throw Exception("Unknown final optimization method: " + finalMethod);

  if(finalOptimizer)
  {
    ApplicationTools::displayResult("Final optimization step", finalMethod);
    finalOptimizer->setProfiler(profiler);
    finalOptimizer->setMessageHandler(messageHandler);
    finalOptimizer->setMaximumNumberOfEvaluations(nbEvalMax);
    finalOptimizer->getStopCondition()->setTolerance(tolerance);
    finalOptimizer->setVerbose(verbose);
    finalOptimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    finalOptimizer->init(tl->getParameters());
    finalOptimizer->optimize();
    n += finalOptimizer->getNumberOfEvaluations();
    delete finalOptimizer;
  }
  
  if(verbose) ApplicationTools::displayResult("Performed", TextTools::toString(n) + " function evaluations.");
}

/******************************************************************************/

void PhylogeneticsApplicationTools::printOptimizationHelp(bool topo, bool clock)
{
  if(!ApplicationTools::message) return;
	*ApplicationTools::message << "optimization                   | [yes/no] optimize parameters?" << endl;
	*ApplicationTools::message << "optimization.method            | [DB/fullD] method to use" << endl;
	*ApplicationTools::message << "optimization.method.derivatives| [gradient/newton] use Conjugate gradient" << endl;
  *ApplicationTools::message << "                               | or Newton-Raphson" << endl;
	*ApplicationTools::message << "optimization.method_DB.step    | Number of progressive step to perform." << endl;
	*ApplicationTools::message << "optimization.final             | [none|simplex|powell] final step." << endl;
	*ApplicationTools::message << "optimization.verbose           | [0,1,2] level of verbose" << endl;
	*ApplicationTools::message << "optimization.message_handler   | [none, std ot file path] where to dislay" << endl;
  *ApplicationTools::message << "                               | optimization messages" << endl;
	*ApplicationTools::message << "                               | (if std, uses 'cout' to display messages)." << endl;
	*ApplicationTools::message << "optimization.profiler          | [none, std ot file path] where to display" << endl;
  *ApplicationTools::message << "                               | optimization steps (if std, uses 'cout'" << endl;
	*ApplicationTools::message << "                               | to display optimization steps)." << endl;
	*ApplicationTools::message << "optimization.tolerance         | [double] tolerance parameter for stopping" << endl;
  *ApplicationTools::message << "                               | the estimation." << endl;
	*ApplicationTools::message << "optimization.max_number_f_eval | [int] max. # of likelihood computations." << endl;
	*ApplicationTools::message << "optimization.ignore_parameter  | [list] parameters to ignore during optimization." << endl;
  if(!clock)
  {
	*ApplicationTools::message << "optimization.scale_first       | [yes, no] tell if a global scale" << endl;
  *ApplicationTools::message << "                               | optimization must be done prior to" << endl;
	*ApplicationTools::message << "                               | separate estimation of branch lengths." << endl;
	*ApplicationTools::message << "optimization.scale_first       | " << endl;
	*ApplicationTools::message << "                     .tolerance| [double] tolerance parameter for global" << endl;
  *ApplicationTools::message << "                               | scale optimization." << endl;
	*ApplicationTools::message << "             .max_number_f_eval| [int] maximum number of computation for" << endl;
  *ApplicationTools::message << "                               | global scale optimization." << endl;
  *ApplicationTools::message << "_______________________________|__________________________________________" << endl;
  }
  if(topo && !clock)
  {
	*ApplicationTools::message << "optimization.topology          | [yes/no] Optimize tree topology?" << endl;
	*ApplicationTools::message << "optimization.topology.algorithm| [nni] Topology movements to use." << endl;
	*ApplicationTools::message << "optimization.topology.nstep    | estimate numerical parameters every 'n'" << endl;
  *ApplicationTools::message << "                               | topology movement rounds." << endl;
	*ApplicationTools::message << "optimization.topology.numfirst | [yes/no] Optimize num. parameters first?" << endl;
	*ApplicationTools::message << "optimization.topology.tolerance| " << endl;
	*ApplicationTools::message << "                        .before| Tolerance for prior estimation." << endl;
	*ApplicationTools::message << "                        .during| Tolerance during estimation." << endl;
  }
  *ApplicationTools::message << "_______________________________|__________________________________________" << endl;
}

/******************************************************************************/

void PhylogeneticsApplicationTools::writeTree(
	const TreeTemplate<Node> & tree,
	map<string, string> & params,
	const string & suffix,
	bool verbose) throw (Exception)
{
	string file = ApplicationTools::getAFilePath("output.tree.file", params, true, false, suffix, false);
	Newick newick;
	newick.write(tree, file, true);
	if(verbose) ApplicationTools::displayResult("Wrote tree to file ", file);
}

/******************************************************************************/

void PhylogeneticsApplicationTools::printOutputTreeHelp()
{
  if(!ApplicationTools::message) return;
	*ApplicationTools::message << "Output tree parameters:" << endl;
	*ApplicationTools::message << "output.tree.file              | file where to write the tree" << endl;
  *ApplicationTools::message << "______________________________|___________________________________________" << endl;
}

/******************************************************************************/

