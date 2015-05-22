//
// File: BppORateDistributionFormat.cpp
// Created by: Laurent Guéguen and Julien Dutheil
// Created on: Fri 16 november 2012, at 13:44
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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

#include "BppORateDistributionFormat.h"
#include "../Model/RateDistribution/ConstantRateDistribution.h"
#include "../Model/RateDistribution/GammaDiscreteRateDistribution.h"
#include "../Model/RateDistribution/GaussianDiscreteRateDistribution.h"
#include "../Model/RateDistribution/ExponentialDiscreteRateDistribution.h"

//From bpp-core:
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Numeric/Prob/InvariantMixedDiscreteDistribution.h>
#include <Bpp/Numeric/Prob/SimpleDiscreteDistribution.h>
#include <Bpp/Numeric/Prob/MixtureOfDiscreteDistributions.h>
#include <Bpp/Io/BppOParametrizableFormat.h>

using namespace bpp;

// From the STL:
#include <iomanip>

using namespace std;


DiscreteDistribution* BppORateDistributionFormat::read(
    const std::string& distDescription,
    bool parseArguments)
{
  unparsedArguments_.clear();
  string distName;
  map<string, string> args;
  KeyvalTools::parseProcedure(distDescription, distName, args);
  auto_ptr<DiscreteDistribution> rDist;

  if (distName == "Uniform")
    throw Exception("BppO Warning: Uniform distribution is deprecated, use Constant instead.");

  if ((distName == "InvariantMixed") || (distName == "Invariant"))
  {
    // We have to parse the nested distribution first:
    string nestedDistDescription = args["dist"];
    if (TextTools::isEmpty(nestedDistDescription))
      throw Exception("BppORateDistributionFormat::read. Missing argument 'dist' for distribution 'Invariant'.");
    if (verbose_)
      ApplicationTools::displayResult("Invariant Mixed distribution", distName );
    BppORateDistributionFormat nestedReader(false);
    DiscreteDistribution* nestedDistribution = nestedReader.read(nestedDistDescription, false);
    map<string, string> unparsedArgumentsNested(nestedReader.getUnparsedArguments());

    // Now we create the Invariant rate distribution:
    rDist.reset(new InvariantMixedDiscreteDistribution(nestedDistribution, 0.1, 0.000001));

    // Then we update the parameter set:
    for (map<string, string>::iterator it = unparsedArgumentsNested.begin();
         it != unparsedArgumentsNested.end(); it++)
    {
      unparsedArguments_["Invariant." + it->first] = it->second;
    }

    if (args.find("p") != args.end())
      unparsedArguments_["Invariant.p"] = args["p"];
  }
  else if (distName == "Constant")
  {
    if (!allowConstant_)
      throw Exception("BppORateDistributionFormat::read(). Constant distribution not allowed.");

    if (args.find("value") != args.end())
      throw Exception("Found argument 'value' in Constant distribution. Constant distribution is defined to have an average of 1.");
    rDist.reset(new ConstantRateDistribution());
  }
  else if (distName == "Simple")
  {
    if (args.find("values") == args.end())
      throw Exception("Missing argument 'values' in Simple distribution");
    if (args.find("probas") == args.end())
      throw Exception("Missing argument 'probas' in Simple distribution");
    vector<double> probas, values;

    string rf = args["values"];
    StringTokenizer strtok(rf.substr(1, rf.length() - 2), ",");
    while (strtok.hasMoreToken())
      values.push_back(TextTools::toDouble(strtok.nextToken()));

    rf = args["probas"];
    StringTokenizer strtok2(rf.substr(1, rf.length() - 2), ",");
    while (strtok2.hasMoreToken())
      probas.push_back(TextTools::toDouble(strtok2.nextToken()));

    std::map<size_t, std::vector<double> > ranges;

    if (args.find("ranges") != args.end())
    {
      string rr = args["ranges"];
      StringTokenizer strtok3(rr.substr(1, rr.length() - 2), ",");
      string desc;
      double deb, fin;
      size_t num;
      size_t po, pf, ppv;
      while (strtok3.hasMoreToken())
      {
        desc = strtok3.nextToken();
        po = desc.find("[");
        ppv = desc.find(";");
        pf = desc.find("]");
        num = (size_t)(TextTools::toInt(desc.substr(1, po - 1)));
        deb = TextTools::toDouble(desc.substr(po + 1, ppv - po - 1));
        fin = TextTools::toDouble(desc.substr(ppv + 1, pf - ppv - 1));
        vector<double> vd;
        vd.push_back(deb);
        vd.push_back(fin);
        ranges[num] = vd;
      }
    }
    if (ranges.size() == 0)
      rDist.reset(new SimpleDiscreteDistribution(values, probas));
    else
      rDist.reset(new SimpleDiscreteDistribution(values, ranges, probas));

    vector<string> v = rDist->getParameters().getParameterNames();

    for (size_t i = 0; i < v.size(); i++)
    {
      unparsedArguments_[v[i]] = TextTools::toString(rDist->getParameterValue(rDist->getParameterNameWithoutNamespace(v[i])));
    }
  }
  else if (distName == "Mixture")
  {
    if (args.find("probas") == args.end())
      throw Exception("Missing argument 'probas' in Mixture distribution");
    vector<double> probas;
    vector<DiscreteDistribution*> v_pdd;
    string rf = args["probas"];
    StringTokenizer strtok2(rf.substr(1, rf.length() - 2), ",");
    while (strtok2.hasMoreToken())
      probas.push_back(TextTools::toDouble(strtok2.nextToken()));

    vector<string> v_nestedDistrDescr;

    size_t nbd = 0;
    while (args.find("dist" + TextTools::toString(++nbd)) != args.end())
      v_nestedDistrDescr.push_back(args["dist" + TextTools::toString(nbd)]);

    if (v_nestedDistrDescr.size() != probas.size())
      throw Exception("Number of distributions (keyword 'dist" + TextTools::toString(probas.size()) + "') do not fit the number of probabilities");

    for (unsigned i = 0; i < v_nestedDistrDescr.size(); i++)
    {
      BppORateDistributionFormat nestedReader(false);
      auto_ptr<DiscreteDistribution> pdd(nestedReader.read(v_nestedDistrDescr[i], false));
      map<string, string> unparsedArgumentsNested(nestedReader.getUnparsedArguments());

      for (map<string, string>::iterator it = unparsedArgumentsNested.begin(); it != unparsedArgumentsNested.end(); it++)
      {
        unparsedArguments_[distName + "." + TextTools::toString(i + 1) + "_" + it->first] = it->second;
      }
      v_pdd.push_back(pdd.release());
    }
    rDist.reset(new MixtureOfDiscreteDistributions(v_pdd, probas));
  }
  else
  {
    if (args.find("n") == args.end())
      throw Exception("Missing argument 'n' (number of classes) in " + distName
                      + " distribution");
    size_t nbClasses = TextTools::to<size_t>(args["n"]);

    if (distName == "Gamma")
    {
      rDist.reset(new GammaDiscreteRateDistribution(nbClasses, 1.));

      if (args.find("alpha") != args.end())
        unparsedArguments_["Gamma.alpha"] = args["alpha"];
      if (args.find("beta") != args.end())
        throw Exception("Found argument 'beta' in Gamma distribution. Gamma distribution is defined to have an average of 1, with beta=alpha.");
    }
    else if (distName == "Gaussian")
    {
      rDist.reset(new GaussianDiscreteRateDistribution(nbClasses, 1));
      if (args.find("mu") != args.end())
        throw Exception("Found argument 'mu' in Gaussian distribution. Gaussian distribution is defined to have an average of 1, with mu=1.");
      if (args.find("sigma") != args.end())
        unparsedArguments_["Gaussian.sigma"] = args["sigma"];
    }
    else if (distName == "Exponential")
    {
      rDist.reset(new ExponentialDiscreteRateDistribution(nbClasses));
      if (args.find("lambda") != args.end())
        throw Exception("Found argument 'lambda' in Exponential distribution. Exponential distribution is defined to have an average of 1, with lambda=1.");
    }
    else
    {
      throw Exception("Unsupported rate distribution: " + distName + ".");
    }
  }
  if (verbose_)
  {
    ApplicationTools::displayResult("Distribution", distName);
    ApplicationTools::displayResult("Number of classes", TextTools::toString((int)rDist->getNumberOfCategories()));
  }

  if (parseArguments)
    initialize_(*rDist);

  return rDist.release();
}


void BppORateDistributionFormat::write(
    const DiscreteDistribution& dist,
    OutputStream& out,
    std::map<std::string, std::string>& globalAliases,
    std::vector<std::string>& writtenNames) const
{
  bool comma = false;

  const DiscreteDistribution* pd;

  out << dist.getName() + "(";

  const InvariantMixedDiscreteDistribution* invar = dynamic_cast<const InvariantMixedDiscreteDistribution*>(&dist);
  if (invar)
  {
    pd = invar->getVariableSubDistribution();
    out << "dist=";
    write(*pd, out, globalAliases, writtenNames);
    comma = true;
  }
  else
  {
    const MixtureOfDiscreteDistributions* mix = dynamic_cast<const MixtureOfDiscreteDistributions*>(&dist);
    if (mix)
    {
      size_t nd = mix->getNumberOfDistributions();
      for (size_t i = 0; i < nd; i++)
      {
        if (comma)
          out << ",";
        out << "dist" + TextTools::toString(i + 1) + "=";
        write(*mix->getNDistribution(i), out, globalAliases, writtenNames);
        comma = true;
      }
      out << ",probas=(";
      for (size_t i = 0; i < nd; i++)
      {
        out << mix->getNProbability(i);
        if (i != nd - 1)
          out << ",";
      }
      out << ")";
      for (size_t i = 1; i < nd; i++)
      {
        writtenNames.push_back(mix->getNamespace() + "theta" + TextTools::toString(i));
      }
    }
  }
  if (dynamic_cast<const ExponentialDiscreteRateDistribution*>(&dist) ||
      dynamic_cast<const GammaDiscreteRateDistribution*>(&dist) ||
      dynamic_cast<const GaussianDiscreteRateDistribution*>(&dist))
  {
    if (comma)
      out << ",";
    out << "n="  << dist.getNumberOfCategories();
    comma = true;
  }

  const SimpleDiscreteDistribution* ps = dynamic_cast<const SimpleDiscreteDistribution*>(&dist);
  if (ps)
  {
    size_t nd = ps->getNumberOfCategories();
    if (comma)
      out << ",";
    out << "values=(";
    for (size_t i = 0; i < nd; i++)
    {
      out << ps->getCategory(i);
      if (i != nd - 1)
        out << ",";
    }
    out << "),probas=(";
    for (size_t i = 0; i < nd; i++)
    {
      out << ps->getProbability(i);
      if (i != nd - 1)
        out << ",";
    }
    out << ")";

    const std::map<size_t, std::vector<double> > range = ps->getRanges();
    if (range.size() != 0)
    {
      out << ", ranges=(";
      std::map<size_t, std::vector<double> >::const_iterator it(range.begin());
      while (it != range.end())
      {
        out << "V" << TextTools::toString(it->first);
        out << "[" << TextTools::toString(it->second[0]) << ";" << TextTools::toString(it->second[1]) << "]";
        it++;
        if (it != range.end())
          out << ",";
      }
    }

    out << ")";

    for (size_t i = 1; i < nd; i++)
    {
      writtenNames.push_back(ps->getNamespace() + "theta" + TextTools::toString(i));
    }
    for (size_t i = 1; i < nd + 1; i++)
    {
      writtenNames.push_back(ps->getNamespace() + "V" + TextTools::toString(i));
    }

    comma = true;
  }


  // Writing the parameters
  BppOParametrizableFormat bOP;
  bOP.write(dynamic_cast<const ParameterAliasable*>(&dist), out, globalAliases, dist.getIndependentParameters().getParameterNames(), writtenNames, true, comma);
  out << ")";
}


