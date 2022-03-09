#ifndef QNTOOLS_STATDISCRIMINATOR_H_
#define QNTOOLS_STATDISCRIMINATOR_H_

#include <cmath>
#include <vector>
#include <iostream>

#include "Rtypes.h"
#include "Stat.hpp"
#include "StatCollect.hpp"
#include "StatCalculate.hpp"

namespace Qn {

class StatDiscriminator : public Stat {
 public:
   
  /// Default constructor.
  StatDiscriminator() = default;
  
  /// Construct StatDiscriminator from StatCalculate.
  explicit StatDiscriminator(StatCalculate &sc);
  
  /// Construct StatDiscriminator from StatCollect.
  explicit StatDiscriminator(StatCollect &sc, ErrorType type = ErrorType::BOOTSTRAP);
  
  virtual ~StatDiscriminator();
  
  [[nodiscard]] double Mean() const { return value_; } 
  [[nodiscard]] double StandardErrorOfMean() const { return error_; }
  [[nodiscard]] double SumWeights() const { return weight_; }
  
  friend StatDiscriminator Merge(const StatDiscriminator &, const StatDiscriminator &);
  
  friend StatDiscriminator operator+(const StatDiscriminator &, const StatDiscriminator &);
  friend StatDiscriminator operator-(const StatDiscriminator &, const StatDiscriminator &);
  friend StatDiscriminator operator*(const StatDiscriminator &, const StatDiscriminator &);
  friend StatDiscriminator operator/(const StatDiscriminator &, const StatDiscriminator &);
  friend StatDiscriminator operator/(const StatDiscriminator &, double);
  friend StatDiscriminator operator*(const StatDiscriminator &, double);
  friend StatDiscriminator operator*(double, const StatDiscriminator &);
  friend StatDiscriminator Pow(const StatDiscriminator &, double);
  friend StatDiscriminator Sqrt(const StatDiscriminator &);
  
 protected:
  double value_{0};
  double error_{0};
  double weight_{0};
  
  /// \cond CLASSIMP
 ClassDef(StatDiscriminator, 1);
  /// \endcond
};

StatDiscriminator Merge(const StatDiscriminator &lhs, const StatDiscriminator &rhs);

StatDiscriminator operator+(const StatDiscriminator &lhs, const StatDiscriminator &rhs);
StatDiscriminator operator-(const StatDiscriminator &lhs, const StatDiscriminator &rhs);

StatDiscriminator operator*(const StatDiscriminator &lhs, const StatDiscriminator &rhs);
StatDiscriminator operator/(const StatDiscriminator &num, const StatDiscriminator &den);

StatDiscriminator operator*(const StatDiscriminator &operand, double scale);
StatDiscriminator operator*(double operand, const StatDiscriminator &rhs);
StatDiscriminator operator/(const StatDiscriminator &operand, double scale);

StatDiscriminator Pow(const StatDiscriminator &base, double exp);
StatDiscriminator Sqrt(const StatDiscriminator &operand);

}
#endif  // QNTOOLS_STATDISCRIMINATOR_H_