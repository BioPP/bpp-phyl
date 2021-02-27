//
// File: CromosomeSubstitutionModel.h
// Created by: Anat Shafir
// Created on: 2020
//


#ifndef _CHROMOSOMESUBSTITUTIONMODEL_H_
#define _CHROMOSOMESUBSTITUTIONMODEL_H_

#include "AbstractSubstitutionModel.h"
#include <Bpp/Seq/Alphabet/ChromosomeAlphabet.h>
#include <Bpp/Exceptions.h>
#include <Bpp/Phyl/NewLikelihood/DataFlow/ExtendedFloatTools.h>

#define lowerBoundOfRateParam 0.0
#define lowerBoundOfExpParam -3.0
#define lowerBoundBaseNumber 3
#define upperBoundOfRateParam 100.0
#define upperBoundLinearRateParam 5.0
#define upperBoundExpParam 4.6
#define IgnoreParam -999
#define DemiEqualDupl -2
#define EPSILON 2.22045e-016
using namespace std;
namespace bpp
{


class ChromosomeSubstitutionModel :
  public AbstractSubstitutionModel
{
public:
  enum rootFreqType {UNIFORM, ROOT_LL, STATIONARY, FIXED};
  enum rateChangeFunc {LINEAR = 0, EXP = 1};
  enum typeOfTransition {GAIN_T = 0, LOSS_T = 1, DUPL_T = 2, DEMIDUPL_T = 3, BASENUM_T = 4, MAXCHR_T = 5, NUMTYPES = 6, ILLEGAL = 7};
  enum paramType {BASENUM = 0, BASENUMR = 1, DUPL = 2, LOSS = 3, GAIN = 4, LOSSR = 5, GAINR = 6, DUPLR = 7, DEMIDUPL = 8, NUM_OF_CHR_PARAMS = 9};

private:
  double gain_;
  double loss_;
  double dupl_;
  double demiploidy_;
  double gainR_;
  double lossR_;
  int baseNum_;
  double baseNumR_;
  double duplR_;
  unsigned int maxChrRange_;
  rootFreqType freqType_;
  rateChangeFunc rateChangeFuncType_;
  int ChrMinNum_;
  int ChrMaxNum_;
  double firstNormQ_;
  mutable bool pijtCalledFromDeriv_;
  mutable MatrixXef pijtEf_;
  mutable MatrixXef dpijtEf_;
  mutable MatrixXef d2pijtEf_;
  mutable bool useExtendedFloatTools_;

  

protected:
  mutable std::vector< RowMatrix<double> > vPowExp_;
  mutable std::vector <MatrixXef> vPowExpEf_;



public:
  ChromosomeSubstitutionModel(const ChromosomeAlphabet* alpha, 
    double gain, 
    double loss, 
    double dupl, 
    double demi,
    double gainR,
    double lossR,
    int baseNum,
    double baseNumR,
    double duplR,
    unsigned int maxChrRange, 
    rootFreqType freqType,
    rateChangeFunc rateChangeType,
    bool useExtendedFloat = false);

  //constructor with vector of parameters
  ChromosomeSubstitutionModel(const ChromosomeAlphabet* alpha, 
    vector<double> modelParams,
    unsigned int maxChrRange,
    rootFreqType freqType,
    rateChangeFunc rateChangeType,
    bool useExtendedFloat = false);

  virtual ~ChromosomeSubstitutionModel() {}

  ChromosomeSubstitutionModel* clone() const { return new ChromosomeSubstitutionModel(*this); }

  
public:
  static ChromosomeSubstitutionModel* initRandomModel(
    const ChromosomeAlphabet* alpha,
    vector<double> initParams,
    unsigned int chrRange,
    rootFreqType rootFrequenciesType,
    rateChangeFunc rateChangeType,
    vector<unsigned int>& fixedParams,
    double parsimonyBound = 0,
    bool useExtendedFloat = false);



  const Matrix<double>& getPij_t    (double d) const;
  const MatrixXef& getPijEf_t(double t) const;
  const Matrix<double>& getPij_t_func2(double d) const;
  const Matrix<double>& getPijt_test(double d) const;
  const MatrixXef& getPijtEf_test(double d) const;
  const Matrix<double>& getPij_t_func3(double d) const;
  const Matrix<double>& getPij_t_func4(double d) const;
  const Matrix<double>& getdPij_dt  (double d) const;
  const MatrixXef& getdPijEf_dt  (double d) const;
  const Matrix<double>& getd2Pij_dt2(double d) const;
  const MatrixXef& getd2PijEf_dt2(double d) const;

  //void calculateExp_Qt(int pow, double s, size_t m, double v);
  std::string getName() const { return "Chromosome"; }
  void setFreq(std::map<int, double>& freqs);
  size_t getNumberOfStates() const { return size_; }
  int getMin() const {return ChrMinNum_;}
  int getMax() const {return ChrMaxNum_;}
  unsigned int getMaxChrRange() const {return maxChrRange_;}
  bool checkIfReachedConvergence(const Matrix<double>& pijt, const Matrix<double>& mt_prev) const;
  bool checkIfReachedConvergence(const MatrixXef& pijt, const MatrixXef& mt_prev) const;
  double getInitValue(size_t i, int state) const;
  int getBaseNumber() const {return baseNum_;}
  double getDemiDupl() const {return demiploidy_;}
  double getConstDupl () const{return dupl_;}
  double getChangeRateDupl() const {return duplR_;}
  double getConstGain() const {return gain_;}
  double getChangeRateGain() const {return gainR_;}
  double getConstLoss() const {return loss_;}
  double getChangeRateLoss() const {return lossR_;}
  double getBaseNumR() const {return baseNumR_;}
  double getRate (size_t state, double constRate, double changeRate) const;
  static void getSetOfFixedParameters(vector<double>& initParams, vector<unsigned int>& fixedParams, map<int, double>& setOfFixedParams);

  
  //These functions should be used from chromsome number optimizer
  void setBoundsForEquivalentParameter(Parameter &param, string parameterName) const;
  void checkParametersBounds() const;


protected:
  void calculatePijtUsingEigenValues(double t) const;
  static void getRandomParameter(paramType type, double initParamValue, vector<double>& randomParams, double upperBound, double upperBoundLinear, double upperBoundExp, rateChangeFunc rateFunc, int maxChrNum, unsigned int chrRange, map<int, double>& setOfFixedParameters);
  void updateParameters();
  void updateLinearParameters();
  void updateExpParameters();
  void updateBaseNumParameters(std::shared_ptr<IntervalConstraint> interval);
  void updateMatrices();
  void updateQWithBaseNumParameters(size_t currChrNum, size_t minChrNum, size_t maxChrNum);
  void updateQWithGain(size_t currChrNum, size_t minChrNum);
  void updateQWithLoss(size_t currChrNum, size_t minChrNum);
  void updateQWithDupl(size_t currChrNum, size_t minChrNum, size_t maxChrNum = 0);
  void updateQWithDemiDupl(size_t currChrNum, size_t minChrNum, size_t maxChrNum);
  void updateEigenMatrices();
  void getParametersValues();
  void calculateExp_Qt(size_t pow, double s, size_t m, double v) const;
  void calculateExp_Qt(size_t pow, double* s, double v, bool ef = false) const;
  double getFirstNorm() const;
  double get_epsilon() const{ return 0.0001;};
  void updateConstRateParameter(double paramValueConst, double paramValueChange, string parameterName, std::shared_ptr<IntervalConstraint> interval);
  void updateLinearChangeParameter(double paramValueConst, double paramValueChange, string parameterName);
  void setNewBoundsForLinearParameters(double &constRate, double &changeRate, string paramNameConst, string paramNameLinear);
};
} // end of namespace bpp.

#endif  // _CHROMOSOMESUBSTITUTIONMODEL_H_