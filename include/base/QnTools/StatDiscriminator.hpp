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
  explicit StatDiscriminator(StatCollect &sc);
  
  virtual ~StatDiscriminator();
  
  double Mean() const;

  double MeanFromPropagation() const { return value_; }

  double MeanFromBootstrap() const;

  double StandardErrorOfMean() const;

  double SumWeights() const { return weight_; }

  double StdDevOfMeanFromPropagation() const { return error_; }

  /// Returns the standard error of the mean from bootstrapping using the variance statistic.
  double StdDevOfMeanFromBootstrapVariance() const { return std::sqrt(VarianceOfMeanFromBootstrap()); }

  /// Returns variance of the sample mean from bootstrapping using the variance statistic.
  double VarianceOfMeanFromBootstrap() const;

  double Neff() const { return n_effective_; }

  /// Retunrs the vector of means of bootstrap samples
  const std::vector<double>& GetSampleMeans() const { return sample_means_; }

  /// Retunrs the vector of weights of bootstrap samples
  const std::vector<double>& GetSampleWeights() const { return sample_weights_; }

  /// Sets the error type.
  void SetMeanType(ErrorType type) { type_of_mean_ = type; }

  /// Returns the error type.
  ErrorType GetMeanType() const { return type_of_mean_; }

  void SetValue(double value) { value_ = value; }
  void SetError(double value) { error_ = value; }
  void SetWeight(double value) { weight_ = value; }
  void SetVEW(double value, double error=0, double weight=1);
  void AddSampleMean(double value) { sample_means_.push_back(value); }
  void AddSampleWeight(double value) { sample_weights_.push_back(value); }

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

  std::vector<double> sample_means_; /// means of bootstrap samples
  std::vector<double> sample_weights_; /// weights of bootstrap samples

  double n_effective_{0};

  ErrorType type_of_mean_ = ErrorType::PROPAGATION;
  
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
