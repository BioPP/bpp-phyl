//
// File: COA.cpp
// Created by: Bastien Boussau
// Created on: Tue May 18 15:23:20 2010
//

/*
   Copyright or � or Copr. CNRS, (November 16, 2004)

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
#include "COA.h"

// From Utils:
#include <Utils/FileTools.h>
#include <Utils/TextTools.h>
#include <Utils/ApplicationTools.h>
#include <Utils/StringTokenizer.h>

// From NumCalc:
#include <NumCalc/VectorTools.h>
#include <NumCalc/MatrixTools.h>


// From the STL:
#include <iostream>
#include <fstream>
#include <string>

using namespace bpp;

using namespace std;


/******************************************************************************/

IncludingInterval COA::_WeightConstraint (-1.0, 1.0);


COA::COA(
  const ProteicAlphabet* alpha,
  string baseModel,
  vector<double> axisWeights,
  vector<double> equilibriumFrequencies) :
  /*,
     IncludingInterval WeightConstraint) :*/
  // AbstractSubstitutionModel(alpha),
  ProteinSubstitutionModel(),
  AbstractReversibleSubstitutionModel(alpha, "COA."),
  _baseModel(""),
  _P(0)
{
  _baseModel = baseModel;
  setNamespace(_baseModel+"+COA.");
  // _WeightConstraint = WeightConstraint;
  // _P = new RowMatrix<double>(20,20);
  int numberOfAxes = 19;
  // Setting the exchangeability matrix
  if (baseModel == "JC69")
  {
    for (unsigned int i = 0; i < 20; i++)
    {
      for (unsigned int j = 0; j < 20; j++)
      {
        generator_(i, j) = (i == j) ? -1. : 1. / 19.;
        exchangeability_(i, j) = generator_(i, j) * 20.;
      }
    }
  }
  else if (baseModel == "DSO78")
  {
    #include "__DSO78ExchangeabilityCode"
  }
  else if (baseModel == "JTT92")
  {
    #include "__JTT92ExchangeabilityCode"
  }
  else
  {
    readFromFile();
  }
  // Now we need to set the equilibrium frequencies

  // _WeightConstraint = new IncludingInterval (-1.0, 1.0);

  /* for (int i = 0; i < numberOfAxes; i++)
     {
     addParameter_(Parameter("COA.AxWei" + TextTools::toString(i), axisWeights[i])); // , &_WeightConstraint ));
     }*/
  // We set equilibrium frequencies to default empirical values (obtained on a dataset of 115 species with 3336 sites)

  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(0) + "_" + TextTools::toString(0), 1.08550953220348  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(0) + "_" + TextTools::toString(1), 1.34015833562751  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(0) + "_" + TextTools::toString(2), -1.34472986103380  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(0) + "_" + TextTools::toString(3), 0.531426663381176  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(0) + "_" + TextTools::toString(4), -3.00411856094694  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(0) + "_" + TextTools::toString(5), 0.5412256470285  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(0) + "_" + TextTools::toString(6), 0.647747317411624  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(0) + "_" + TextTools::toString(7), 0.480375873411076  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(0) + "_" + TextTools::toString(8), -0.207855736501920  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(0) + "_" + TextTools::toString(9), -1.17859613594355  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(0) + "_" + TextTools::toString(10), -0.241324993042156  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(0) + "_" + TextTools::toString(11), -1.58417814180993  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(0) + "_" + TextTools::toString(12), -0.0453770340560095  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(0) + "_" + TextTools::toString(13), -1.04721234149109  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(0) + "_" + TextTools::toString(14), 0.89922597116275  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(0) + "_" + TextTools::toString(15), -1.44534542029580  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(0) + "_" + TextTools::toString(16), 0.173221335465671  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(0) + "_" + TextTools::toString(17), 0.0840771707348402  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(0) + "_" + TextTools::toString(18), -0.81253607024061  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(0) + "_" + TextTools::toString(19), 0.710229305642463  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(1) + "_" + TextTools::toString(0), 0.0724807060062998  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(1) + "_" + TextTools::toString(1), 0.560472641691227  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(1) + "_" + TextTools::toString(2), -0.486180636302382  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(1) + "_" + TextTools::toString(3), -0.360854519138564  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(1) + "_" + TextTools::toString(4), 8.57476296883954  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(1) + "_" + TextTools::toString(5), 1.46542502082867  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(1) + "_" + TextTools::toString(6), -0.988341508075483  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(1) + "_" + TextTools::toString(7), -0.105119465762877  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(1) + "_" + TextTools::toString(8), 1.79941407042448  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(1) + "_" + TextTools::toString(9), -0.98670689011329  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(1) + "_" + TextTools::toString(10), 0.509379038112302  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(1) + "_" + TextTools::toString(11), -0.855446404296351  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(1) + "_" + TextTools::toString(12), 0.504462684745238  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(1) + "_" + TextTools::toString(13), 0.656990322583955  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(1) + "_" + TextTools::toString(14), -0.142887074289711  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(1) + "_" + TextTools::toString(15), 0.496329476572691  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(1) + "_" + TextTools::toString(16), 0.0199299798692343  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(1) + "_" + TextTools::toString(17), 0.389158484481491  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(1) + "_" + TextTools::toString(18), -0.939435997291933  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(1) + "_" + TextTools::toString(19), -0.128251656383470  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(2) + "_" + TextTools::toString(0), 0.732504450191948  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(2) + "_" + TextTools::toString(1), -0.888257979200525  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(2) + "_" + TextTools::toString(2), 1.19781360092311  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(2) + "_" + TextTools::toString(3), 0.994514860791484  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(2) + "_" + TextTools::toString(4), -3.05276055932937  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(2) + "_" + TextTools::toString(5), 1.26063936893623  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(2) + "_" + TextTools::toString(6), -1.44051062509273  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(2) + "_" + TextTools::toString(7), 0.429510919097404  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(2) + "_" + TextTools::toString(8), -0.542063098473156  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(2) + "_" + TextTools::toString(9), -0.696434635138525  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(2) + "_" + TextTools::toString(10), -0.244741979670493  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(2) + "_" + TextTools::toString(11), -0.392942468986463  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(2) + "_" + TextTools::toString(12), 0.300975050509667  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(2) + "_" + TextTools::toString(13), -0.224685368192731  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(2) + "_" + TextTools::toString(14), -1.08656707901677  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(2) + "_" + TextTools::toString(15), 2.03639305181187  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(2) + "_" + TextTools::toString(16), 1.25925646261773  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(2) + "_" + TextTools::toString(17), -2.76536721776282  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(2) + "_" + TextTools::toString(18), -2.2898075487343  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(2) + "_" + TextTools::toString(19), 0.0280562260106468  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(3) + "_" + TextTools::toString(0), -0.050791982456506  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(3) + "_" + TextTools::toString(1), -0.87200330836045  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(3) + "_" + TextTools::toString(2), -1.09490716978747  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(3) + "_" + TextTools::toString(3), 1.25321296302308  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(3) + "_" + TextTools::toString(4), 2.19571647563149  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(3) + "_" + TextTools::toString(5), -1.80462988553516  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(3) + "_" + TextTools::toString(6), 1.52198699771143  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(3) + "_" + TextTools::toString(7), 0.487260247075932  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(3) + "_" + TextTools::toString(8), 0.23575109936382  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(3) + "_" + TextTools::toString(9), -1.14782169408790  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(3) + "_" + TextTools::toString(10), -1.25666856034824  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(3) + "_" + TextTools::toString(11), 0.653881269853238  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(3) + "_" + TextTools::toString(12), 1.70792307661356  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(3) + "_" + TextTools::toString(13), 0.802917565579562  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(3) + "_" + TextTools::toString(14), -0.737707656781832  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(3) + "_" + TextTools::toString(15), 0.70085806740801  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(3) + "_" + TextTools::toString(16), -0.939489964261408  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(3) + "_" + TextTools::toString(17), -0.711875903852786  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(3) + "_" + TextTools::toString(18), -0.916022694855166  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(3) + "_" + TextTools::toString(19), 1.00231585038895  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(4) + "_" + TextTools::toString(0), -1.25720822176789  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(4) + "_" + TextTools::toString(1), 1.00407659957668  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(4) + "_" + TextTools::toString(2), 0.399251325864265  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(4) + "_" + TextTools::toString(3), 1.46845156800862  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(4) + "_" + TextTools::toString(4), 0.0262839131503755  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(4) + "_" + TextTools::toString(5), 0.050040648207608  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(4) + "_" + TextTools::toString(6), 1.26094785564366  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(4) + "_" + TextTools::toString(7), -0.285290843253106  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(4) + "_" + TextTools::toString(8), -0.50008664495708  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(4) + "_" + TextTools::toString(9), 0.727466063306291  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(4) + "_" + TextTools::toString(10), -0.55593078068898  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(4) + "_" + TextTools::toString(11), -1.61760307032272  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(4) + "_" + TextTools::toString(12), 0.144008909975884  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(4) + "_" + TextTools::toString(13), -0.395821711771923  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(4) + "_" + TextTools::toString(14), 0.0950634647927662  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(4) + "_" + TextTools::toString(15), 2.19578121255548  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(4) + "_" + TextTools::toString(16), -0.289288302317705  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(4) + "_" + TextTools::toString(17), -0.203835596438925  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(4) + "_" + TextTools::toString(18), 0.818354332878153  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(4) + "_" + TextTools::toString(19), -0.625685674406527  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(5) + "_" + TextTools::toString(0), 0.321786477604558  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(5) + "_" + TextTools::toString(1), 0.843636063216866  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(5) + "_" + TextTools::toString(2), -1.43668449329593  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(5) + "_" + TextTools::toString(3), -1.03810704899687  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(5) + "_" + TextTools::toString(4), -1.99685977745902  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(5) + "_" + TextTools::toString(5), -0.837944513532773  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(5) + "_" + TextTools::toString(6), -1.26087442916475  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(5) + "_" + TextTools::toString(7), 0.261999208527547  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(5) + "_" + TextTools::toString(8), -1.43253749515546  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(5) + "_" + TextTools::toString(9), -0.0921634206705537  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(5) + "_" + TextTools::toString(10), 0.698254101822732  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(5) + "_" + TextTools::toString(11), 0.237410095118684  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(5) + "_" + TextTools::toString(12), -0.968110842972956  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(5) + "_" + TextTools::toString(13), 0.653891044643201  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(5) + "_" + TextTools::toString(14), 0.236276257265647  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(5) + "_" + TextTools::toString(15), 2.19813897746200  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(5) + "_" + TextTools::toString(16), -2.15866429245791  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(5) + "_" + TextTools::toString(17), 2.98284821699641  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(5) + "_" + TextTools::toString(18), -0.226935265447352  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(5) + "_" + TextTools::toString(19), 0.53765226414757  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(6) + "_" + TextTools::toString(0), 0.681173229361595  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(6) + "_" + TextTools::toString(1), 0.81790483491224  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(6) + "_" + TextTools::toString(2), 0.736296814524209  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(6) + "_" + TextTools::toString(3), -0.719692837060817  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(6) + "_" + TextTools::toString(4), 1.08017094429282  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(6) + "_" + TextTools::toString(5), -1.56927636989044  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(6) + "_" + TextTools::toString(6), -1.53507787529079  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(6) + "_" + TextTools::toString(7), 0.249221962262765  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(6) + "_" + TextTools::toString(8), -0.457600144224005  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(6) + "_" + TextTools::toString(9), 1.00571053301677  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(6) + "_" + TextTools::toString(10), -1.49903745402771  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(6) + "_" + TextTools::toString(11), -0.243081444820963  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(6) + "_" + TextTools::toString(12), 3.92871755642722  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(6) + "_" + TextTools::toString(13), -0.74432561973885  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(6) + "_" + TextTools::toString(14), 0.228650238524935  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(6) + "_" + TextTools::toString(15), -0.285108171359608  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(6) + "_" + TextTools::toString(16), -0.299799838473454  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(6) + "_" + TextTools::toString(17), -2.30712327500503  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(6) + "_" + TextTools::toString(18), 0.712276934076681  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(6) + "_" + TextTools::toString(19), -0.0352724380109204  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(7) + "_" + TextTools::toString(0), 0.281728793148175  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(7) + "_" + TextTools::toString(1), 0.793959252070366  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(7) + "_" + TextTools::toString(2), 1.76320200787985  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(7) + "_" + TextTools::toString(3), -0.1247477518605  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(7) + "_" + TextTools::toString(4), 0.810631766456864  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(7) + "_" + TextTools::toString(5), -3.34216329587717  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(7) + "_" + TextTools::toString(6), -0.304166324789613  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(7) + "_" + TextTools::toString(7), 0.0597861934549203  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(7) + "_" + TextTools::toString(8), 0.680629745896586  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(7) + "_" + TextTools::toString(9), -0.712410289781704  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(7) + "_" + TextTools::toString(10), -0.0954732644354652  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(7) + "_" + TextTools::toString(11), -0.431316161629253  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(7) + "_" + TextTools::toString(12), -3.18423950294211  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(7) + "_" + TextTools::toString(13), 1.43317269648283  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(7) + "_" + TextTools::toString(14), -0.851267846192951  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(7) + "_" + TextTools::toString(15), -0.00114139827703831  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(7) + "_" + TextTools::toString(16), 0.785389130695595  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(7) + "_" + TextTools::toString(17), -2.65132828516872  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(7) + "_" + TextTools::toString(18), 1.51585575520904  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(7) + "_" + TextTools::toString(19), 0.432382137638077  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(8) + "_" + TextTools::toString(0), -0.560223097346833  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(8) + "_" + TextTools::toString(1), -0.458377116356005  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(8) + "_" + TextTools::toString(2), -0.804854776823822  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(8) + "_" + TextTools::toString(3), 2.12458679645293  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(8) + "_" + TextTools::toString(4), -2.20617989565429  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(8) + "_" + TextTools::toString(5), 0.493896166931432  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(8) + "_" + TextTools::toString(6), -1.88674173494950  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(8) + "_" + TextTools::toString(7), -0.330488301962051  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(8) + "_" + TextTools::toString(8), 2.53812498782889  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(8) + "_" + TextTools::toString(9), -0.100358666643086  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(8) + "_" + TextTools::toString(10), -0.456252042885703  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(8) + "_" + TextTools::toString(11), -0.219722540089655  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(8) + "_" + TextTools::toString(12), 0.286149874390503  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(8) + "_" + TextTools::toString(13), 1.74859760877058  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(8) + "_" + TextTools::toString(14), 0.76679791730802  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(8) + "_" + TextTools::toString(15), -0.463542432507303  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(8) + "_" + TextTools::toString(16), -0.125921353062139  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(8) + "_" + TextTools::toString(17), 2.58250347959997  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(8) + "_" + TextTools::toString(18), 2.65945358920765  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(8) + "_" + TextTools::toString(19), 0.512202114693986  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(9) + "_" + TextTools::toString(0), 0.194314051727008  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(9) + "_" + TextTools::toString(1), -0.826041505116499  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(9) + "_" + TextTools::toString(2), 0.943729826386134  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(9) + "_" + TextTools::toString(3), -0.407748199645375  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(9) + "_" + TextTools::toString(4), 1.84339501258909  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(9) + "_" + TextTools::toString(5), 1.51649406399173  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(9) + "_" + TextTools::toString(6), 0.119689094154147  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(9) + "_" + TextTools::toString(7), -0.480664768495801  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(9) + "_" + TextTools::toString(8), -4.01902257600276  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(9) + "_" + TextTools::toString(9), -0.631406796190173  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(9) + "_" + TextTools::toString(10), 0.0596156760118508  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(9) + "_" + TextTools::toString(11), -0.526231243209235  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(9) + "_" + TextTools::toString(12), 0.379898326154561  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(9) + "_" + TextTools::toString(13), -0.359308914077076  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(9) + "_" + TextTools::toString(14), -0.191311459287480  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(9) + "_" + TextTools::toString(15), 0.0389859061443903  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(9) + "_" + TextTools::toString(16), -0.140096055631894  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(9) + "_" + TextTools::toString(17), 0.542189863635839  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(9) + "_" + TextTools::toString(18), 3.59356938927178  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(9) + "_" + TextTools::toString(19), 1.19847442719103  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(10) + "_" + TextTools::toString(0), 1.89741860143244  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(10) + "_" + TextTools::toString(1), -0.837870955139807  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(10) + "_" + TextTools::toString(2), -0.809371066361505  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(10) + "_" + TextTools::toString(3), 1.82162281453909  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(10) + "_" + TextTools::toString(4), 2.42382269382627  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(10) + "_" + TextTools::toString(5), -0.354306125122515  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(10) + "_" + TextTools::toString(6), 0.0241819882405925  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(10) + "_" + TextTools::toString(7), 0.235227170538017  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(10) + "_" + TextTools::toString(8), -0.730076065220892  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(10) + "_" + TextTools::toString(9), 0.853226172450979  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(10) + "_" + TextTools::toString(10), 0.07666826093445  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(10) + "_" + TextTools::toString(11), -0.277902923594346  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(10) + "_" + TextTools::toString(12), -1.78972509877230  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(10) + "_" + TextTools::toString(13), -0.358579313414140  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(10) + "_" + TextTools::toString(14), -0.183356188966359  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(10) + "_" + TextTools::toString(15), -0.186260191262936  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(10) + "_" + TextTools::toString(16), -1.10452039631363  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(10) + "_" + TextTools::toString(17), -0.988400159124436  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(10) + "_" + TextTools::toString(18), 1.09660592489814  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(10) + "_" + TextTools::toString(19), -1.25781813751169  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(11) + "_" + TextTools::toString(0), -0.437257628318621  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(11) + "_" + TextTools::toString(1), 0.0195401160513908  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(11) + "_" + TextTools::toString(2), 0.690882190759627  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(11) + "_" + TextTools::toString(3), 1.43057961058004  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(11) + "_" + TextTools::toString(4), 0.5333351900158  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(11) + "_" + TextTools::toString(5), -2.03009455747631  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(11) + "_" + TextTools::toString(6), -0.574297077038491  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(11) + "_" + TextTools::toString(7), 0.233375109309212  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(11) + "_" + TextTools::toString(8), 0.101338559417307  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(11) + "_" + TextTools::toString(9), -0.909469381435248  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(11) + "_" + TextTools::toString(10), 1.61712149416552  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(11) + "_" + TextTools::toString(11), 0.305400397394848  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(11) + "_" + TextTools::toString(12), 1.33937956542265  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(11) + "_" + TextTools::toString(13), -3.70947659584285  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(11) + "_" + TextTools::toString(14), -0.0158639033956643  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(11) + "_" + TextTools::toString(15), -0.0586195853956851  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(11) + "_" + TextTools::toString(16), -0.00519331068773546  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(11) + "_" + TextTools::toString(17), 3.09652864065114  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(11) + "_" + TextTools::toString(18), 0.543921974295614  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(11) + "_" + TextTools::toString(19), -0.189360512042165  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(12) + "_" + TextTools::toString(0), -1.06142371674827  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(12) + "_" + TextTools::toString(1), 0.860272697621056  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(12) + "_" + TextTools::toString(2), -1.51161543404565  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(12) + "_" + TextTools::toString(3), 1.16080691393398  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(12) + "_" + TextTools::toString(4), 2.00861462540275  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(12) + "_" + TextTools::toString(5), 0.553566611231842  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(12) + "_" + TextTools::toString(6), -1.00056120795515  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(12) + "_" + TextTools::toString(7), 0.781482210057666  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(12) + "_" + TextTools::toString(8), -2.69321071098326  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(12) + "_" + TextTools::toString(9), 0.100778523584181  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(12) + "_" + TextTools::toString(10), -0.371172999919078  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(12) + "_" + TextTools::toString(11), 1.04776140966312  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(12) + "_" + TextTools::toString(12), -1.63042478034389  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(12) + "_" + TextTools::toString(13), -0.673578034598588  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(12) + "_" + TextTools::toString(14), 0.525158686844206  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(12) + "_" + TextTools::toString(15), -0.397998776275377  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(12) + "_" + TextTools::toString(16), 0.735788016384533  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(12) + "_" + TextTools::toString(17), -5.95906429069484  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(12) + "_" + TextTools::toString(18), -0.104699431217637  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(12) + "_" + TextTools::toString(19), 0.167233850933525  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(13) + "_" + TextTools::toString(0), -0.403071696121283  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(13) + "_" + TextTools::toString(1), 0.177682978303177  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(13) + "_" + TextTools::toString(2), 1.98012839206104  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(13) + "_" + TextTools::toString(3), 1.09760995985273  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(13) + "_" + TextTools::toString(4), -0.148294060689771  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(13) + "_" + TextTools::toString(5), 0.743129907965774  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(13) + "_" + TextTools::toString(6), -0.209795703460702  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(13) + "_" + TextTools::toString(7), -0.160069053087494  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(13) + "_" + TextTools::toString(8), 0.101364661818872  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(13) + "_" + TextTools::toString(9), 0.608585511529575  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(13) + "_" + TextTools::toString(10), 0.855716314502404  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(13) + "_" + TextTools::toString(11), -0.320684009375495  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(13) + "_" + TextTools::toString(12), 0.0345450325512344  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(13) + "_" + TextTools::toString(13), 0.529982908408623  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(13) + "_" + TextTools::toString(14), -0.646683566062795  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(13) + "_" + TextTools::toString(15), -1.33296941810514  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(13) + "_" + TextTools::toString(16), -2.40088568401337  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(13) + "_" + TextTools::toString(17), -4.2151530451272  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(13) + "_" + TextTools::toString(18), -2.06195167161038  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(13) + "_" + TextTools::toString(19), 0.907817403407018  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(14) + "_" + TextTools::toString(0), 0.274601024798303  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(14) + "_" + TextTools::toString(1), 0.434333814213829  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(14) + "_" + TextTools::toString(2), -0.712050917069825  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(14) + "_" + TextTools::toString(3), -0.377840154291803  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(14) + "_" + TextTools::toString(4), -2.28277788436393  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(14) + "_" + TextTools::toString(5), 0.5923206045253  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(14) + "_" + TextTools::toString(6), 0.479451170975527  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(14) + "_" + TextTools::toString(7), -0.664870775729664  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(14) + "_" + TextTools::toString(8), 1.35411330212889  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(14) + "_" + TextTools::toString(9), -1.07066249914665  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(14) + "_" + TextTools::toString(10), 0.900755534356667  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(14) + "_" + TextTools::toString(11), 0.504911139112662  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(14) + "_" + TextTools::toString(12), 1.69451283284855  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(14) + "_" + TextTools::toString(13), 0.272732844539596  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(14) + "_" + TextTools::toString(14), -0.293332101597302  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(14) + "_" + TextTools::toString(15), 0.55631910753811  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(14) + "_" + TextTools::toString(16), -0.545258894928866  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(14) + "_" + TextTools::toString(17), -8.71351339993846  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(14) + "_" + TextTools::toString(18), 2.38628193979761  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(14) + "_" + TextTools::toString(19), -0.750654256013073  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(15) + "_" + TextTools::toString(0), -0.269095427328304  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(15) + "_" + TextTools::toString(1), -0.348279550042306  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(15) + "_" + TextTools::toString(2), 0.347455660643089  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(15) + "_" + TextTools::toString(3), 0.378708018571084  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(15) + "_" + TextTools::toString(4), -0.258396490320314  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(15) + "_" + TextTools::toString(5), -1.27565186438219  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(15) + "_" + TextTools::toString(6), -0.00842241285552317  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(15) + "_" + TextTools::toString(7), 0.466873322843914  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(15) + "_" + TextTools::toString(8), -2.79337846612105  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(15) + "_" + TextTools::toString(9), -0.556002114276989  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(15) + "_" + TextTools::toString(10), 1.00266719006890  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(15) + "_" + TextTools::toString(11), -0.481483573192847  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(15) + "_" + TextTools::toString(12), 1.88311760795327  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(15) + "_" + TextTools::toString(13), 3.34282542371850  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(15) + "_" + TextTools::toString(14), 1.73189553428812  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(15) + "_" + TextTools::toString(15), -0.540645094480224  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(15) + "_" + TextTools::toString(16), 0.363545298499678  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(15) + "_" + TextTools::toString(17), 0.249576552256159  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(15) + "_" + TextTools::toString(18), -0.916760485629085  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(15) + "_" + TextTools::toString(19), -1.01638189002350  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(16) + "_" + TextTools::toString(0), 0.369275549020256  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(16) + "_" + TextTools::toString(1), 1.52578281236519  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(16) + "_" + TextTools::toString(2), 1.29894720163174  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(16) + "_" + TextTools::toString(3), 1.00986397844077  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(16) + "_" + TextTools::toString(4), -0.360453040799594  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(16) + "_" + TextTools::toString(5), 1.19949141901544  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(16) + "_" + TextTools::toString(6), 0.0103908901301893  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(16) + "_" + TextTools::toString(7), -0.775882721716194  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(16) + "_" + TextTools::toString(8), -1.27592510415571  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(16) + "_" + TextTools::toString(9), -1.05069729613604  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(16) + "_" + TextTools::toString(10), -0.966117485198648  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(16) + "_" + TextTools::toString(11), 1.27233612856449  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(16) + "_" + TextTools::toString(12), 0.0557923003065529  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(16) + "_" + TextTools::toString(13), 0.766420359798032  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(16) + "_" + TextTools::toString(14), -1.03194139571713  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(16) + "_" + TextTools::toString(15), -0.453163744237857  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(16) + "_" + TextTools::toString(16), -0.466259946464549  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(16) + "_" + TextTools::toString(17), 4.71620969482211  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(16) + "_" + TextTools::toString(18), -0.219783906996106  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(16) + "_" + TextTools::toString(19), -0.954337710601457  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(17) + "_" + TextTools::toString(0), -0.325302295058482  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(17) + "_" + TextTools::toString(1), -0.423733862472779  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(17) + "_" + TextTools::toString(2), 1.77628826449096  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(17) + "_" + TextTools::toString(3), -0.784392806820304  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(17) + "_" + TextTools::toString(4), -0.102788384731683  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(17) + "_" + TextTools::toString(5), 1.08803981933188  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(17) + "_" + TextTools::toString(6), 0.092062400135042  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(17) + "_" + TextTools::toString(7), 1.75809107502665  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(17) + "_" + TextTools::toString(8), 1.26038615491896  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(17) + "_" + TextTools::toString(9), -0.9429085386281  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(17) + "_" + TextTools::toString(10), -0.893389016013572  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(17) + "_" + TextTools::toString(11), 0.149785938234667  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(17) + "_" + TextTools::toString(12), -1.15572720159456  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(17) + "_" + TextTools::toString(13), -0.737544902742  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(17) + "_" + TextTools::toString(14), 1.95150848205868  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(17) + "_" + TextTools::toString(15), 0.193745829990836  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(17) + "_" + TextTools::toString(16), -1.32551510572315  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(17) + "_" + TextTools::toString(17), -0.401452705607341  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(17) + "_" + TextTools::toString(18), 0.994549347682586  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(17) + "_" + TextTools::toString(19), -0.621739763075379  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(18) + "_" + TextTools::toString(0), 0.489856752317616  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(18) + "_" + TextTools::toString(1), -0.217302946607416  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(18) + "_" + TextTools::toString(2), 0.637506466994677  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(18) + "_" + TextTools::toString(3), 0.369007428217951  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(18) + "_" + TextTools::toString(4), 0.95513668028049  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(18) + "_" + TextTools::toString(5), -0.626580802243678  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(18) + "_" + TextTools::toString(6), -0.0169647885481034  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(18) + "_" + TextTools::toString(7), -1.91426726078049  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(18) + "_" + TextTools::toString(8), 0.0587152331543925  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(18) + "_" + TextTools::toString(9), -0.223724412037975  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(18) + "_" + TextTools::toString(10), -0.410592553673277  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(18) + "_" + TextTools::toString(11), 0.435330941103377  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(18) + "_" + TextTools::toString(12), -0.839596760837366  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(18) + "_" + TextTools::toString(13), -0.82088943675775  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(18) + "_" + TextTools::toString(14), 3.24016100558309  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(18) + "_" + TextTools::toString(15), 0.729915598755108  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(18) + "_" + TextTools::toString(16), 0.0126898077775198  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(18) + "_" + TextTools::toString(17), -1.82328397706611  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(18) + "_" + TextTools::toString(18), -1.47045014080454  ));
  addParameter_(Parameter(_baseModel+"+COA.EqFre" + TextTools::toString(18) + "_" + TextTools::toString(19), 0.398816100687077  ));


  RowMatrix<double> M; // (20,20);
  M.resize(20, 20);


  // We build the matrix from the factors
  for (int i = 0; i < 19; i++)
  {
    for (int j = 0; j < 20; j++)
    {
      M(i,j) = getParameter("EqFre" + TextTools::toString(i) + "_" + TextTools::toString(j)).getValue();
    }
  }

  // One more line composed only of 1s to get a square matrix
  for (int j = 0; j < 20; j++)
  {
    M(19,j) = 1.0;
  }

  _P = new RowMatrix<double>;
  (MatrixTools::inv(M, *_P));

  // For now, we constrain the equilibrium frequencies to remain fixed,
  // so we do not need these parameters anymore.
  // In the future, this might be changed, in which case not removing these parameters would be useful.

  resetParameters_(); // We remove all parameters
  // and we put back the axis weights.

  for (int i = 0; i < numberOfAxes; i++)
  {
    addParameter_(Parameter(_baseModel+"+COA.AxWei" + TextTools::toString(i), axisWeights[i])); // , &_WeightConstraint ));
  }

  updateMatrices();
}

/******************************************************************************/

COA::COA(const COA& coa) :
  ProteinSubstitutionModel(coa),
  AbstractReversibleSubstitutionModel(coa),
  _baseModel(""),
  _P(0)
{
  _P = coa._P->clone();
  _baseModel = coa._baseModel;
  _WeightConstraint = coa._WeightConstraint;
}

/******************************************************************************/

COA & COA::operator=(const COA& coa)
{
  ProteinSubstitutionModel::operator=(coa);
  AbstractReversibleSubstitutionModel::operator=(coa);
  if (_P)
    delete _P;
  _P = coa._P->clone();
  _baseModel = coa._baseModel;
  _WeightConstraint = coa._WeightConstraint;
  return *this;
}

/******************************************************************************/

void COA::readFromFile()
{
  ifstream in(_baseModel.c_str(), ios::in);
  // Read exchangeability matrix:
  for (unsigned int i = 1; i < 20; i++)
  {
   string line = FileTools::getNextLine(in);
    StringTokenizer st(line);
    for (unsigned int j = 0; j < i; j++)
    {
   double s = TextTools::toDouble(st.nextToken());
      exchangeability_(i,j) = exchangeability_(j,i) = s;
    }
  }
  // Read frequencies: should be useless, here.
  /*	unsigned int fCount = 0;
     while(in && fCount < 20)
     {
     string line = FileTools::getNextLine(in);
     StringTokenizer st(line);
     while(st.hasMoreToken() && fCount < 20)
     {
      _freq[fCount] = TextTools::toDouble(st.nextToken());
      fCount++;
     }
     }
     double sf = VectorTools::sum(_freq);
     if(sf - 1 > 0.000001)
     {
     ApplicationTools::displayMessage("WARNING!!! Frequencies sum to " + TextTools::toString(sf) + ", frequencies have been scaled.");
     sf *= 1./sf;
     }
   */
  // Now build diagonal of the exchangeability matrix:
  for (unsigned int i = 0; i < 20; i++)
  {
   double sum = 0;
    for (unsigned int j = 0; j < 20; j++)
    {
      if (j != i)
        sum += exchangeability_(i,j);
    }
    exchangeability_(i,i) = -sum;
  }

  // Closing stream:
  in.close();
}


/******************************************************************************/
void COA::computeEquilibriumFrequencies()
{
  // This M matrix could be rebuilt every time we need
  // to update the substitution matrix. This would make sense, if we
  // want to optimize the factor vectors.
  // However, optimizing the factor vectors requires adding a lot of free parameters,
  // so for now we avoid that.

  /*RowMatrix<double> M;// (20,20);
     M.resize(20, 20);

     //We build the matrix from the factors
     for (int i = 0 ; i < 19; i++ ) {
     for (int j = 0 ; j < 20 ; j++ ) {
       M(i,j)=parameters_.getParameter("EqFre"+TextTools::toString(i)+ "_" + TextTools::toString(j))->getValue();
     }
     }

     //One more line composed only of 1s to get a square matrix
     for (int j = 0 ; j < 20 ; j++ ) {
     M(19,j) = 1.0;
     }

     RowMatrix<double> *P =new RowMatrix<double>(MatrixTools::inv(M));
   */


  // Now we get the coordinates
  vector<double> coord;

  for (int i = 0; i < 19; i++)
  {
   coord.push_back(getParameter("AxWei" + TextTools::toString(i)).getValue());
  }
  // We add the final 1:
  coord.push_back(1.0);

  // We now have the P matrix and the coordinate vectors, a simple multiplication will get us the frequencies.

  freq_ = prodMatrixVector(*_P,coord);

  /* cout <<"coord: ";
     VectorTools::print (coord);
     cout <<"freqs: ";
     VectorTools::print (_freq);*/
  bool neg = false;
  for (int i = 0; i < 20; i++)
  {
    if (freq_[i] < 0)
    {
      neg = true;
      freq_[i] = 0.000001;
    }
  }
  if (neg == true)
  {
    cerr << "WARNING: There were negative amino-acid frequencies !!!!!" << endl;
    double s = VectorTools::sum(freq_);
    for (int i = 0; i < 20; i++)
    {
      freq_[i] = freq_[i] / s;
    }
    cerr << "freqs: ";
    VectorTools::print (freq_);
  }
}

/******************************************************************************/

void COA::updateMatrices()
{
  computeEquilibriumFrequencies();
  AbstractReversibleSubstitutionModel::updateMatrices();
}


/******************************************************************************/
/* Function that computes the product of a matrix P of size 20 with a vector V of size 20, and returns a vector.*/

vector<double> COA::prodMatrixVector(RowMatrix<double>& P, vector<double>& V)
{
   vector<double> E(20, 0.0);

  for (int i = 0; i < 20; i++)  // factors are in lines
  {
    for (int j = 0; j < 20; j++)  // amino acids are in columns
    {
      E[i] = E[i] + P(i,j) * V[j];
    }
  }
  return E;
}

/******************************************************************************/
