// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "BppOFrequencySetFormat.h"

using namespace bpp;

const std::string IOFrequencySetFactory::BPPO_FORMAT = "Bpp0";

std::unique_ptr<IFrequencySet> IOFrequencySetFactory::createReader(const std::string& format)
{
  if (format == BPPO_FORMAT)
    return std::make_unique<BppOFrequencySetFormat>(BppOFrequencySetFormat::ALL, true, 1);
  else
    throw Exception("Format " + format + " is not supported for input.");
}

std::unique_ptr<OFrequencySet> IOFrequencySetFactory::createWriter(const std::string& format)
{
  if (format == BPPO_FORMAT)
    return std::make_unique<BppOFrequencySetFormat>(BppOFrequencySetFormat::ALL, true, 1);
  else
    throw Exception("Format " + format + " is not supported for output.");
}
