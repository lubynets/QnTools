#include "StatDiscriminator.hpp"

namespace Qn {

/**
 * Default constructor
 */

  
StatDiscriminator::~StatDiscriminator() = default;

StatDiscriminator::StatDiscriminator(StatCalculate &scalc, ErrorType type) {
  ErrorType type_buffer = scalc.GetErrorType();
  scalc.SetErrorType(type);
  weight_ = scalc.SumWeights();
  value_ = scalc.Mean();
  error_ = scalc.StandardErrorOfMean();
  scalc.SetErrorType(type_buffer);
}

StatDiscriminator::StatDiscriminator(StatCollect &scol, ErrorType type) {
  StatCalculate scalc(scol);
  scalc.SetErrorType(type);
  
  weight_ = scalc.SumWeights();           // how to call upper constructor in the current one?
  value_ = scalc.Mean();
  error_ = scalc.StandardErrorOfMean();
}

void StatDiscriminator::SetVEW(double value, double error, double weight) {
  value_ = value;
  error_ = error;
  weight_ = weight;
}

double StatDiscriminator::StandardErrorOfMean() const {
  if (GetErrorType() == ErrorType::PROPAGATION) return StdDevOfMeanFromPropagation();
  else                                          return StdDevOfMeanFromBootstrapVariance();
}

double StatDiscriminator::Mean() const {
  if (GetMeanType() == ErrorType::PROPAGATION) return MeanFromPropagation();
  else                                         return MeanFromBootstrap();
}

double StatDiscriminator::VarianceOfMeanFromBootstrap() const {
  Statistics stats;
  for (std::size_t i = 0; i < sample_means_.size(); ++i) {
    if (std::isnan(sample_means_[i]) && sample_weights_[i] > 0) {
      std::cerr << "Skipping NAN-value with non-zero weight for " << i << "-th sample"  << std::endl;
      continue;
    }
    stats.Fill(sample_means_[i], sample_weights_[i]);
  }
  return stats.Variance();
}

double StatDiscriminator::MeanFromBootstrap() const {
  Statistics stats;
  for (std::size_t i = 0; i < sample_means_.size(); ++i) {
    if (std::isnan(sample_means_[i]) && sample_weights_[i] > 0) {
      std::cerr << "Skipping NAN-value with non-zero weight for " << i << "-th sample"  << std::endl;
      continue;
    }
    stats.Fill(sample_means_[i], sample_weights_[i]);
  }
  return stats.Mean();
}

StatDiscriminator Merge(const StatDiscriminator &lhs, const StatDiscriminator &rhs) {
  StatDiscriminator merged;
  
  merged.weight_ = lhs.SumWeights() + rhs.SumWeights();
  merged.value_ = (lhs.MeanFromPropagation()*lhs.SumWeights() + rhs.MeanFromPropagation()*rhs.SumWeights()) / merged.SumWeights();
  merged.error_ = std::sqrt(lhs.StdDevOfMeanFromPropagation()*lhs.StdDevOfMeanFromPropagation()*lhs.SumWeights()*lhs.SumWeights() +
                            rhs.StdDevOfMeanFromPropagation()*rhs.StdDevOfMeanFromPropagation()*rhs.SumWeights()*rhs.SumWeights()) /
                  merged.SumWeights();

  // Bootstrap samples
  if (lhs.weight_ == 0.) {
    merged.sample_means_ = rhs.sample_means_;
    merged.sample_weights_ = rhs.sample_weights_;
  } else if (rhs.weight_ == 0.) {
    merged.sample_means_ = lhs.sample_means_;
    merged.sample_weights_ = lhs.sample_weights_;
  } else {
    for (size_t i = 0; i < rhs.sample_means_.size(); ++i) {
      auto lhs_mean = lhs.sample_means_[i];
      auto lhs_weight = lhs.sample_weights_[i];
      auto rhs_mean = rhs.sample_means_[i];
      auto rhs_weight = rhs.sample_weights_[i];
      auto merged_weight = lhs_weight + rhs_weight;
      auto lhs_wm = lhs_weight > 0? lhs_mean * lhs_weight : 0.0;
      auto rhs_wm = rhs_weight > 0? rhs_mean * rhs_weight : 0.0;
      auto merged_mean = merged_weight > 0? (lhs_wm + rhs_wm) / merged_weight : 0.0;
      merged.sample_weights_.push_back(merged_weight);
      merged.sample_means_.push_back(merged_mean);
    }
  }

  return merged;
}

StatDiscriminator operator+(const StatDiscriminator &lhs, const StatDiscriminator &rhs) {
  StatDiscriminator sum;
  
  sum.weight_ = lhs.SumWeights();
  sum.value_ = lhs.MeanFromPropagation() + rhs.MeanFromPropagation();
  sum.error_ = std::sqrt(lhs.StdDevOfMeanFromPropagation()*lhs.StdDevOfMeanFromPropagation() +
                         rhs.StdDevOfMeanFromPropagation()*rhs.StdDevOfMeanFromPropagation());

  // Bootstrap samples
  for (size_t i = 0; i < lhs.sample_means_.size(); ++i) {
    auto lhs_mean = lhs.sample_means_[i];
    auto lhs_weight = lhs.sample_weights_[i];
    auto rhs_mean = rhs.sample_means_[i];
    auto rhs_weight = rhs.sample_weights_[i];
    /* if any of arguments is not determined, sum is also not determined */
    if (lhs_weight <= 0 || rhs_weight <= 0) {
      sum.sample_weights_[i] = 0.0;
      sum.sample_means_[i] = 0.0;
    } else {
      sum.sample_weights_[i] = lhs.sample_weights_[i];
      sum.sample_means_[i] = lhs_mean + rhs_mean;
    }
  }

  return sum;
}

StatDiscriminator operator-(const StatDiscriminator &lhs, const StatDiscriminator &rhs) {
  StatDiscriminator difference;
  
  difference.weight_ = lhs.SumWeights();
  difference.value_ = lhs.MeanFromPropagation() - rhs.MeanFromPropagation();
  difference.error_ = std::sqrt(lhs.StdDevOfMeanFromPropagation()*lhs.StdDevOfMeanFromPropagation() +
                                rhs.StdDevOfMeanFromPropagation()*rhs.StdDevOfMeanFromPropagation());

  // Bootstrap samples
  for (size_t i = 0; i < lhs.sample_means_.size(); ++i) {
    auto lhs_mean = lhs.sample_means_[i];
    auto lhs_weight = lhs.sample_weights_[i];
    auto rhs_mean = rhs.sample_means_[i];
    auto rhs_weight = rhs.sample_weights_[i];
    /* if any of arguments is not determined, difference is also not determined */
    if (lhs_weight <= 0 || rhs_weight <= 0) {
      difference.sample_weights_[i] = 0;
      difference.sample_means_[i] = 0;
    } else {
      difference.sample_weights_[i] = lhs_weight;
      difference.sample_means_[i] = lhs_mean - rhs_mean;
    }
  }

  return difference;
}

StatDiscriminator operator*(const StatDiscriminator &lhs, const StatDiscriminator &rhs) {
  StatDiscriminator product;
  
  product.weight_ = lhs.SumWeights();
  product.value_ = lhs.MeanFromPropagation() * rhs.MeanFromPropagation();
  product.error_ = std::sqrt(lhs.StdDevOfMeanFromPropagation()*lhs.StdDevOfMeanFromPropagation()*rhs.MeanFromPropagation()*rhs.MeanFromPropagation() +
                             rhs.StdDevOfMeanFromPropagation()*rhs.StdDevOfMeanFromPropagation()*lhs.MeanFromPropagation()*lhs.MeanFromPropagation());

  // Bootstrap samples
  for (size_t i = 0; i < lhs.sample_means_.size(); ++i) {
    auto lhs_mean = lhs.sample_means_[i];
    auto lhs_weight = lhs.sample_weights_[i];
    auto rhs_mean = rhs.sample_means_[i];
    auto rhs_weight = rhs.sample_weights_[i];
    if (lhs_weight <= 0 || rhs_weight <= 0) {
      product.sample_weights_[i] = 0.;
      product.sample_means_[i] = 0;
    } else {
      product.sample_weights_[i] = lhs_weight;
      product.sample_means_[i] = lhs_mean * rhs_mean;
    }
  }

  return product;
}

StatDiscriminator operator/(const StatDiscriminator &num, const StatDiscriminator &den) {
  StatDiscriminator ratio;
  
  ratio.weight_ = num.SumWeights();
  ratio.value_ = num.MeanFromPropagation() / den.MeanFromPropagation();
  ratio.error_ = std::sqrt(num.MeanFromPropagation()*num.MeanFromPropagation()*den.StdDevOfMeanFromPropagation()*den.StdDevOfMeanFromPropagation() +
                           den.MeanFromPropagation()*den.MeanFromPropagation()*num.StdDevOfMeanFromPropagation()*num.StdDevOfMeanFromPropagation()) /
                 den.MeanFromPropagation() / den.MeanFromPropagation();

  // Bootstrap samples
  for (size_t i = 0; i < num.sample_means_.size(); ++i) {
    auto lhs_mean = num.sample_means_[i];
    auto lhs_weight = num.sample_weights_[i];
    auto rhs_mean = den.sample_means_[i];
    auto rhs_weight = den.sample_weights_[i];
    if (lhs_weight <= 0 || rhs_weight <= 0) {
      ratio.sample_weights_[i] = 0.0;
      ratio.sample_means_[i] = 0;
    } else {
      ratio.sample_weights_[i] = lhs_weight;
      ratio.sample_means_[i] = lhs_mean / rhs_mean;
    }
  }

  return ratio;
}

StatDiscriminator operator*(const StatDiscriminator &operand, double scale) {
  StatDiscriminator scaled;
  
  scaled.weight_ = operand.SumWeights();
  scaled.value_ = operand.MeanFromPropagation() * scale;
  scaled.error_ = operand.StdDevOfMeanFromPropagation() * scale;

  // Bootstrap samples
  for (size_t i = 0; i < operand.sample_means_.size(); ++i) {
    scaled.sample_weights_[i] = operand.sample_weights_[i];
    scaled.sample_means_[i] = operand.sample_means_[i] * scale;
  }

  return scaled;
}

StatDiscriminator operator*(double operand, const StatDiscriminator &rhs) {
  return rhs * operand;
}

StatDiscriminator operator/(const StatDiscriminator &operand, double scale) {
  return operand * (1./scale);
}

StatDiscriminator Pow(const StatDiscriminator &base, double exp) {
  StatDiscriminator result;
  
  result.weight_ = base.SumWeights();
  result.value_ = std::pow(base.MeanFromPropagation(), exp);
  result.error_ = exp * result.MeanFromPropagation() * base.StdDevOfMeanFromPropagation() / base.MeanFromPropagation();

  // Bootstrap samples
  for (size_t i = 0; i < base.sample_means_.size(); ++i) {
    result.sample_weights_[i] = base.sample_weights_[i];
    result.sample_means_[i] = std::pow(base.sample_means_[i],exp);
  }

  return result;  
}

StatDiscriminator Sqrt(const StatDiscriminator &operand) {
  return Pow(operand, 0.5);
}

}
