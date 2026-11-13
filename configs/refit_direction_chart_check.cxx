// Mathematical validation of the exact header used by the Calypso refit.
// Input consists of real native-writer phi/theta values; no tracking is run.
#include "RefitDirectionChart.h"
#include <algorithm>
#include <iomanip>
#include <iostream>

Acts::Vector3 direction(const Acts::BoundVector& p) {
  const double phi = p[Acts::eBoundPhi], theta = p[Acts::eBoundTheta];
  return {std::cos(phi)*std::sin(theta), std::sin(phi)*std::sin(theta), std::cos(theta)};
}

int main() {
  std::size_t count = 0;
  double phi, theta, meanError = 0., inverseError = 0., derivativeError = 0.;
  while (std::cin >> phi >> theta) {
    Acts::BoundVector input = Acts::BoundVector::Zero();
    input[Acts::eBoundPhi] = phi;
    input[Acts::eBoundTheta] = theta;
    const auto [rotated, jacobian] = refitDirectionChart(input);
    const auto [recovered, inverse] = refitDirectionChart(rotated, true);
    meanError = std::max(meanError, (direction(input)-direction(recovered)).cwiseAbs().maxCoeff());
    Acts::BoundSquareMatrix expected = Acts::BoundSquareMatrix::Identity();
    if (std::sin(theta) < 0.) expected(Acts::eBoundTheta, Acts::eBoundTheta) = -1.;
    inverseError = std::max(inverseError, (inverse*jacobian-expected).cwiseAbs().maxCoeff());
    for (int column : {Acts::eBoundPhi, Acts::eBoundTheta}) {
      constexpr double step = 1.e-7;
      Acts::BoundVector plus = input, minus = input;
      plus[column] += step;
      minus[column] -= step;
      const auto high = refitDirectionChart(plus).first;
      const auto low = refitDirectionChart(minus).first;
      Acts::BoundVector difference = high-low;
      difference[Acts::eBoundPhi] = std::remainder(difference[Acts::eBoundPhi], 2.*std::acos(-1.));
      for (int row : {Acts::eBoundPhi, Acts::eBoundTheta}) {
        const double numeric = difference[row]/(2.*step);
        derivativeError = std::max(derivativeError,
            std::abs(numeric-jacobian(row, column))/(1.+std::abs(jacobian(row, column))));
      }
    }
    ++count;
  }
  std::cout << std::setprecision(17)
      << "{\"real_writer_angle_pairs\":" << count
      << ",\"direction_roundtrip_max_absolute\":" << meanError
      << ",\"jacobian_inverse_max_absolute\":" << inverseError
      << ",\"scaled_central_difference_max_error\":" << derivativeError << "}\n";
  return count > 500 && meanError < 1.e-12 && inverseError < 1.e-7 && derivativeError < 1.e-7 ? 0 : 1;
}
