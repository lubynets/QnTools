#include "H1F.hpp"
#include <iostream>

namespace Qn {
  
H1F Merge(const H1F &lhs, const H1F &rhs) {
  H1F merged = *(H1F*)lhs.Clone();
  merged.Add(&rhs);

  return merged;
}
}