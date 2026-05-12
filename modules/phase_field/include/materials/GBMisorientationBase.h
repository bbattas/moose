#pragma once

#include "Material.h"
#include "EulerAngleTxtFileReader.h"
#include "GBPairPacking.h"

#include <Eigen/Geometry>
#include <vector>

/**
 * Base material that builds and stores a global all-pair grain misorientation table.
 *
 * This class intentionally does not inspect the active grains/features at a qp and
 * intentionally does not declare qp-specific misorientation output properties. A
 * derived parent material should handle GrainTracker/order-parameter logic and call
 * getMisorientationData(i, j) when it needs theta, polar_ax, and azim_ax for a pair.
 */
class GBMisorientationBase : public Material
{
public:
  static InputParameters validParams();

  GBMisorientationBase(const InputParameters & parameters);

protected:
  struct MisorientationData
  {
    Real theta = 0.0;    ///< Misorientation angle in radians
    Real polar_ax = 0.0; ///< Misorientation-axis polar angle in radians
    Real azim_ax = 0.0;  ///< Misorientation-axis azimuth angle in radians

    Eigen::Quaternion<Real> q = Eigen::Quaternion<Real>(1.0, 0.0, 0.0, 0.0);
    Real qnorm = 0.0; ///< Norm of the vector part of q before axis normalization
  };

  /**
   * No qp-specific work is done in this base class.
   *
   * Derived materials should override computeQpProperties(), gather their active
   * grains/features, and then call getMisorientationData(i, j).
   */
  virtual void computeQpProperties() override {}

  /// Number of grains/orientations used to construct the table.
  std::size_t grainNum() const;

  /// Number of unique pair entries stored in the misorientation table.
  std::size_t misorientationTableSize() const;

  /// Packed flat-vector index for the unordered grain pair (grain_i, grain_j).
  std::size_t packedPairIndex(std::size_t grain_i, std::size_t grain_j) const;

  /// Full misorientation record for the unordered grain pair (grain_i, grain_j).
  const MisorientationData & getMisorientationData(std::size_t grain_i, std::size_t grain_j) const;

  /// Convenience accessors.
  Real getMisorientationAngle(std::size_t grain_i, std::size_t grain_j) const;
  Real getMisorientationPolarAxis(std::size_t grain_i, std::size_t grain_j) const;
  Real getMisorientationAzimuthAxis(std::size_t grain_i, std::size_t grain_j) const;

  /**
   * Direct table access for derived classes that need to share the same packed
   * index across several flat pair arrays.
   */
  const std::vector<MisorientationData> & misorientationTable() const { return _miso_table; }

private:
  /// Build the global all-pair misorientation lookup table once.
  void buildMisorientationTable();

  /// Convert a 3x3 symmetry rotation matrix to a quaternion.
  void rotationSymmetryToQuaternion(const double O[3][3], Eigen::Quaternion<Real> & q);

  /// Define the cubic crystal symmetry operators as quaternions.
  void defineSymmetryOperator();

  /// Compute the minimum misorientation data between two orientation quaternions.
  MisorientationData getMisorientationFromQuaternion(const Eigen::Quaternion<Real> & qi,
                                                     const Eigen::Quaternion<Real> & qj) const;

protected:
  /// UserObject providing Euler angles for grain/orientation ids.
  const EulerAngleProvider & _euler;

  /// Per-grain Euler angles and quaternions used to build the pair table.
  std::vector<EulerAngles> _euler_angle;
  std::vector<Eigen::Quaternion<Real>> _quat_angle;

  /// Cubic symmetry operators.
  std::vector<Eigen::Quaternion<Real>> _sym_quat;
  static constexpr int _o_sym = 24;

  /// Flat packed table of all unique grain-pair misorientation records.
  std::vector<MisorientationData> _miso_table;
};
