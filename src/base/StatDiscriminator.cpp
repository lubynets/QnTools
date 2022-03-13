#include "StatDiscriminator.hpp"

namespace Qn {

/**
 * Default constructor
 */

  
StatDiscriminator::~StatDiscriminator() = default;

StatDiscriminator::StatDiscriminator(StatCalculate &scalc) {
  weight_ = scalc.SumWeights();
  value_ = scalc.Mean();
  error_ = scalc.StandardErrorOfMean();
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

StatDiscriminator Merge(const StatDiscriminator &lhs, const StatDiscriminator &rhs) {
  StatDiscriminator merged;
  
  merged.weight_ = lhs.SumWeights() + rhs.SumWeights();
  merged.value_ = (lhs.Mean()*lhs.SumWeights() + rhs.Mean()*rhs.SumWeights()) / merged.SumWeights();
  merged.error_ = std::sqrt(lhs.StandardErrorOfMean()*lhs.StandardErrorOfMean()*lhs.SumWeights()*lhs.SumWeights() +
                            rhs.StandardErrorOfMean()*rhs.StandardErrorOfMean()*rhs.SumWeights()*rhs.SumWeights()) /
                  merged.SumWeights();

  return merged;
}

StatDiscriminator operator+(const StatDiscriminator &lhs, const StatDiscriminator &rhs) {
  StatDiscriminator sum;
  
  sum.weight_ = lhs.SumWeights();
  sum.value_ = lhs.Mean() + rhs.Mean();
  sum.error_ = std::sqrt(lhs.StandardErrorOfMean()*lhs.StandardErrorOfMean() +
                         rhs.StandardErrorOfMean()*rhs.StandardErrorOfMean());

  return sum;
}

StatDiscriminator operator-(const StatDiscriminator &lhs, const StatDiscriminator &rhs) {
  StatDiscriminator difference;
  
  difference.weight_ = lhs.SumWeights();
  difference.value_ = lhs.Mean() - rhs.Mean();
  difference.error_ = std::sqrt(lhs.StandardErrorOfMean()*lhs.StandardErrorOfMean() +
                                rhs.StandardErrorOfMean()*rhs.StandardErrorOfMean());

  return difference;
}

StatDiscriminator operator*(const StatDiscriminator &lhs, const StatDiscriminator &rhs) {
  StatDiscriminator product;
  
  product.weight_ = lhs.SumWeights();
  product.value_ = lhs.Mean() * rhs.Mean();
  product.error_ = std::sqrt(lhs.StandardErrorOfMean()*lhs.StandardErrorOfMean()*rhs.Mean()*rhs.Mean() +
                             rhs.StandardErrorOfMean()*rhs.StandardErrorOfMean()*lhs.Mean()*lhs.Mean());

  return product;
}

StatDiscriminator operator/(const StatDiscriminator &num, const StatDiscriminator &den) {
  StatDiscriminator ratio;
  
  ratio.weight_ = num.SumWeights();
  ratio.value_ = num.Mean() / den.Mean();
  ratio.error_ = std::sqrt(num.Mean()*num.Mean()*den.StandardErrorOfMean()*den.StandardErrorOfMean() +
                           den.Mean()*den.Mean()*num.StandardErrorOfMean()*num.StandardErrorOfMean()) /
                 den.Mean() / den.Mean();

  return ratio;
}

StatDiscriminator operator*(const StatDiscriminator &operand, double scale) {
  StatDiscriminator scaled;
  
  scaled.weight_ = operand.SumWeights();
  scaled.value_ = operand.Mean() * scale;
  scaled.error_ = operand.StandardErrorOfMean() * scale;

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
  result.value_ = std::pow(base.Mean(), exp);
  result.error_ = exp * result.Mean() * base.StandardErrorOfMean() / base.Mean();

  return result;  
}

StatDiscriminator Sqrt(const StatDiscriminator &operand) {
  return Pow(operand, 0.5);
}

}