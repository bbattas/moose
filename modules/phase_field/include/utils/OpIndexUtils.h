#pragma once
#include <algorithm>
#include <cctype>
#include <string>
#include "MooseError.h" // or whichever header provides mooseError

namespace OpIndexUtils
{

inline unsigned
parseOpIndex(const std::string & full, const std::string & base)
{
  if (full.rfind(base, 0) != 0)
    mooseError("Variable '", full, "' does not start with base '", base, "'.");
  const std::string tail = full.substr(base.size());
  if (tail.empty() ||
      !std::all_of(tail.begin(), tail.end(), [](unsigned char c) { return std::isdigit(c); }))
    mooseError("Variable '", full, "' must end with a non-negative integer index.");
  return static_cast<unsigned>(std::stoul(tail));
}

} // namespace OpIndexUtils
