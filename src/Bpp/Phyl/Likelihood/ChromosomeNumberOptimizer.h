//
// File: ChromosomeNumberOptimizer.h
// Created by: Anat Shafir
// Created on: Wednesday September 2 15:05 2020
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004, 2005, 2006)

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
#ifndef _CHROMOSOMENUMBEROPTIMIZER_H_
#define _CHROMOSOMENUMBEROPTIMIZER_H_

//from bpp-core
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Numeric/Function/BrentOneDimension.h>
#include <Bpp/Numeric/Function/ConjugateGradientMultiDimensions.h>
#include <Bpp/Numeric/Function/AbstractNumericalDerivative.h>
#include <Bpp/Numeric/Function/TwoPointsNumericalDerivative.h>

//from bpp-seq
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Alphabet/ChromosomeAlphabet.h>
#include <Bpp/Seq/Container/VectorSequenceContainer.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>


//from bpp-phyl
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/Tree/PhyloTree.h>
#include <Bpp/Phyl/TreeTemplateTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/Likelihood/DRNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Model/ChromosomeSubstitutionModel.h>
#include <Bpp/Phyl/NewLikelihood/NonHomogeneousSubstitutionProcess.h>
#include <Bpp/Phyl/NewLikelihood/RateAcrossSitesSubstitutionProcess.h>
#include <Bpp/Phyl/NewLikelihood/DataFlow/LikelihoodCalculationSingleProcess.h>
#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Model/SubstitutionModelSet.h>
#include <Bpp/Phyl/App/ChromEvolOptions.h>
#include <Bpp/Phyl/OptimizationTools.h>
// From Seqlib:
#include <vector>
#include <map>
#include <utility>
#include <string>

using namespace std;
namespace bpp
{
    class ChromosomeNumberOptimizer{
    // a class which is used for ChromEvol to run the likelihood optimization procedures
    // with different options available in ChromEvol

        private:
            vector <SingleProcessPhyloLikelihood*> vectorOfLikelohoods_;
            vector <Context> vectorOfContexts_;
            const PhyloTree* tree_;
            const ChromosomeAlphabet* alphabet_;
            const VectorSiteContainer* vsc_;
            bool optimizeBaseNumber_;
            vector<unsigned int> numOfPoints_;
            vector<unsigned int> numOfIterations_;
            string typeOfOptimizer_;
            string baseNumOptimizationMethod_;
            unsigned int baseNumberUpperBound_;
            double tolerance_;
            bool standardOptimization_;
            int BrentBracketing_;
            vector <double> probsForMixedOptimization_;
            vector<unsigned int> fixedParams_;

        public:
            ChromosomeNumberOptimizer(
                const PhyloTree* tree,
                const ChromosomeAlphabet* alpha,
                const VectorSiteContainer* vsc,
                unsigned int baseNumberUpperBound):
                    vectorOfLikelohoods_(),
                    vectorOfContexts_(),
                    tree_(tree),
                    alphabet_(alpha),
                    vsc_(vsc),
                    optimizeBaseNumber_(),
                    numOfPoints_(),
                    numOfIterations_(),
                    typeOfOptimizer_(),
                    baseNumOptimizationMethod_(),
                    baseNumberUpperBound_(baseNumberUpperBound),
                    tolerance_(),
                    standardOptimization_(),
                    BrentBracketing_(),
                    probsForMixedOptimization_(),
                    fixedParams_()


            {}

            ChromosomeNumberOptimizer(const ChromosomeNumberOptimizer& opt):
                vectorOfLikelohoods_(opt.vectorOfLikelohoods_),
                vectorOfContexts_(opt.vectorOfContexts_),
                tree_ (opt.tree_),
                alphabet_(opt.alphabet_),
                vsc_(opt.vsc_),
                optimizeBaseNumber_(opt.optimizeBaseNumber_),
                numOfPoints_(opt.numOfPoints_),
                numOfIterations_(opt.numOfIterations_),
                typeOfOptimizer_(opt.typeOfOptimizer_),
                baseNumOptimizationMethod_(opt.baseNumOptimizationMethod_),
                baseNumberUpperBound_(opt.baseNumberUpperBound_),
                tolerance_(opt.tolerance_),
                standardOptimization_(opt.standardOptimization_),
                BrentBracketing_(opt.BrentBracketing_),
                probsForMixedOptimization_(opt.probsForMixedOptimization_),
                fixedParams_(opt.fixedParams_)

            {}
            ChromosomeNumberOptimizer& operator=(const ChromosomeNumberOptimizer& opt){
                vectorOfLikelohoods_ = opt.vectorOfLikelohoods_;
                vectorOfContexts_ = opt.vectorOfContexts_;
                tree_ = opt.tree_;
                alphabet_ = opt.alphabet_;
                vsc_ = opt.vsc_;
                optimizeBaseNumber_ = opt.optimizeBaseNumber_;
                numOfPoints_ = opt.numOfPoints_;
                numOfIterations_ = opt.numOfIterations_;
                typeOfOptimizer_ = opt.typeOfOptimizer_;
                baseNumOptimizationMethod_ = opt.baseNumOptimizationMethod_;
                baseNumberUpperBound_ = opt.baseNumberUpperBound_;
                tolerance_ = opt.tolerance_;
                standardOptimization_ = opt.standardOptimization_;
                BrentBracketing_ = opt.BrentBracketing_;
                probsForMixedOptimization_ = opt.probsForMixedOptimization_;
                fixedParams_ = opt.fixedParams_;
                return *this;
            }
            ChromosomeNumberOptimizer* clone() const { return new ChromosomeNumberOptimizer(*this); }
            virtual ~ChromosomeNumberOptimizer(){clearVectorOfLikelihoods(0);};
            //init models
            void initModels(vector<double> modelParams, double parsimonyBound, ChromosomeSubstitutionModel::rateChangeFunc rateChange, int seed, unsigned int numberOfModels, const string& fixedRootFreqPath, vector<unsigned int>& fixedParams);
            //initialize all the optimization specific members
            void initOptimizer(
                vector<unsigned int> numOfPoints,
                vector<unsigned int> numOfIterations,
                string typeOfOptimizer,
                string baseNumOptimizationMethod,
                double tolerance,
                bool standardOptimization,
                int BrentBracketing,
                vector <double>& probsForMixedOptimization)
            {
                numOfPoints_ = numOfPoints;
                numOfIterations_ = numOfIterations;
                typeOfOptimizer_ = typeOfOptimizer;
                baseNumOptimizationMethod_ = baseNumOptimizationMethod;
                tolerance_ = tolerance;
                standardOptimization_ = standardOptimization;
                BrentBracketing_ =BrentBracketing;
                probsForMixedOptimization_ = probsForMixedOptimization;
                

            }
            void optimize();
            vector<SingleProcessPhyloLikelihood*> getVectorOfLikelihoods(){return vectorOfLikelohoods_;}
            static vector <double> setFixedRootFrequencies(const std::string &path, std::shared_ptr<ChromosomeSubstitutionModel> chrModel);


        protected:
            // for model initiation
            SingleProcessPhyloLikelihood* getLikelihoodFunction(const PhyloTree* tree, const VectorSiteContainer* vsc, std::shared_ptr<ChromosomeSubstitutionModel> &chrModel, DiscreteDistribution* rdist, const string& fixedRootFreqPath);
            
            
            // //functions of optimization
            unsigned int optimizeModelParameters(SingleProcessPhyloLikelihood* tl, double tol, unsigned int maxNumOfIterations, vector<unsigned int> &baseNumCandidates);//, unsigned int inwardBracketing, bool standardOptimization);
            unsigned int optimizeModelParametersOneDimension(SingleProcessPhyloLikelihood* tl, double tol, unsigned int maxNumOfIterations, std::vector<unsigned int> &baseNumCandidates, bool mixed = false, unsigned int currentIterNum = 0);
            unsigned int optimizeMultiDimensions(SingleProcessPhyloLikelihood* tl, double tol, unsigned int maxNumOfIterations, bool mixed = false, unsigned int currentIterNum = 0);
            unsigned int useMixedOptimizers(SingleProcessPhyloLikelihood* tl, double tol, unsigned int maxNumOfIterations, vector <unsigned int> &baseNumCandidates);
            void optimizeBaseNum(SingleProcessPhyloLikelihood* tl, size_t index, std::vector <unsigned int> baseNumCandidates, double* currentLikelihood, double lowerBound, double upperBound);

            // // function working on the likelihoods vector object
            void clearVectorOfLikelihoods(size_t new_size);
            // //void deleteTreeLikAssociatedAttributes(SingleProcessPhyloLikelihood &lik);
            static bool compareLikValues(SingleProcessPhyloLikelihood* lik1, SingleProcessPhyloLikelihood* lik2);

            // // helper functions for optimization
            vector <string> getNonFixedParams(vector <unsigned int> fixedParams, ParameterList &allParams) const;
            void fillVectorOfBaseNumCandidates(vector <unsigned int> &baseNumCandidates, unsigned int lowerBound, unsigned int upperBound) const;
            void getAllPossibleChrRanges(vector <unsigned int> &baseNumCandidates) const;
            string findParameterNameInModel(string fullParameterName) const;
            void constructParamPairsMap(map<string, pair<string, bool>> &paramPairsMap);
            void setNewBounds(const ParameterList params, Parameter &param, map<string, pair<string, bool>> &paramPairsMap, double* lowerBound, const ChromosomeSubstitutionModel* model);

            //print functions
            void printLikParameters(SingleProcessPhyloLikelihood* lik, unsigned int optimized, const string path = "none") const;
            void printRootFrequencies(SingleProcessPhyloLikelihood* lik, const string path = "none") const;
            void printLikelihoodVectorValues(vector <SingleProcessPhyloLikelihood*> lik_vec) const;

    };
}
#endif // _CHROMOSOMENUMBEROPTIMIZER_H_