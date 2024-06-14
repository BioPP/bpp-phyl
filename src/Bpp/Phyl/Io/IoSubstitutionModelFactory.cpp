// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "BppOSubstitutionModelFormat.h"

using namespace bpp;

const std::string IOSubstitutionModelFactory::BPPO_FORMAT = "Bpp0";

ISubstitutionModel* IOSubstitutionModelFactory::createReader(const std::string& format)
{
  if (format == BPPO_FORMAT)
    return new BppOSubstitutionModelFormat(BppOSubstitutionModelFormat::ALL, true, true, true, true, 0);
  else
    throw Exception("Format " + format + " is not supported for input.");
}

OSubstitutionModel* IOSubstitutionModelFactory::createWriter(const std::string& format)
{
  if (format == BPPO_FORMAT)
    return new BppOSubstitutionModelFormat(BppOSubstitutionModelFormat::ALL, true, true, true, true, 0);
  else
    throw Exception("Format " + format + " is not supported for output.");
}
