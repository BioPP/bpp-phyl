// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "BppOSubstitutionModelFormat.h"

using namespace bpp;

const std::string IOSubstitutionModelFactory::BPPO_FORMAT = "Bpp0";

std::unique_ptr<ISubstitutionModel> IOSubstitutionModelFactory::createReader(const std::string& format)
{
  if (format == BPPO_FORMAT)
    return std::make_unique<BppOSubstitutionModelFormat>(BppOSubstitutionModelFormat::ALL, true, true, true, true, 0);
  else
    throw Exception("Format " + format + " is not supported for input.");
}

std::unique_ptr<OSubstitutionModel> IOSubstitutionModelFactory::createWriter(const std::string& format)
{
  if (format == BPPO_FORMAT)
    return std::make_unique<BppOSubstitutionModelFormat>(BppOSubstitutionModelFormat::ALL, true, true, true, true, 0);
  else
    throw Exception("Format " + format + " is not supported for output.");
}
