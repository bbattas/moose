#pragma once

#include "InputParameters.h"
#include "MooseError.h"
#include "EulerAngleTxtFileReader.h"
#include "GBPairPacking.h"

#include <Eigen/Geometry>
#include <vector>

class GBMisorientationHelper
{
public:
  static InputParameters validParams();

protected:
  struct MisorientationData
  {
    Real theta = 0.0;
    Real polar_ax = 0.0;
    Real azim_ax = 0.0;

    Eigen::Quaternion<Real> q = Eigen::Quaternion<Real>(1.0, 0.0, 0.0, 0.0);
    Real qnorm = 0.0;
  };

  GBMisorientationHelper(const InputParameters & parameters);

  bool misorientationEnabled() const { return _enable_misorientation; }

  std::size_t grainNum() const;
  std::size_t misorientationTableSize() const;

  const MisorientationData & getMisorientationData(std::size_t grain_i, std::size_t grain_j) const;

  Real getMisorientationAngle(std::size_t grain_i, std::size_t grain_j) const;
  Real getMisorientationPolarAxis(std::size_t grain_i, std::size_t grain_j) const;
  Real getMisorientationAzimuthAxis(std::size_t grain_i, std::size_t grain_j) const;

private:
  std::size_t packedPairIndex(std::size_t grain_i, std::size_t grain_j) const;

  void buildMisorientationTable();

  void rotationSymmetryToQuaternion(const double O[3][3], Eigen::Quaternion<Real> & q);

  void defineSymmetryOperator();

  MisorientationData getMisorientationFromQuaternion(const Eigen::Quaternion<Real> & qi,
                                                     const Eigen::Quaternion<Real> & qj) const;

protected:
  const bool _enable_misorientation;

  const EulerAngleProvider * _euler = nullptr;

  std::vector<EulerAngles> _euler_angle;
  std::vector<Eigen::Quaternion<Real>> _quat_angle;

  std::vector<Eigen::Quaternion<Real>> _sym_quat;
  static constexpr int _o_sym = 24;

  std::vector<MisorientationData> _miso_table;
};
