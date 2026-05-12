#include "GBMisorientationBase.h"

#include <algorithm>
#include <cmath>

InputParameters
GBMisorientationBase::validParams()
{
  InputParameters params = Material::validParams();

  params.addClassDescription(
      "Base material that builds a global all-pair grain misorientation lookup table.");

  params.addRequiredParam<UserObjectName>("euler_angle_provider",
                                          "Name of Euler angle provider user object.");

  return params;
}

GBMisorientationBase::GBMisorientationBase(const InputParameters & parameters)
  : Material(parameters), _euler(getUserObject<EulerAngleProvider>("euler_angle_provider"))
{
  buildMisorientationTable();
}

std::size_t
GBMisorientationBase::grainNum() const
{
  return _euler.getGrainNum();
}

std::size_t
GBMisorientationBase::misorientationTableSize() const
{
  return _miso_table.size();
}

std::size_t
GBMisorientationBase::packedPairIndex(std::size_t grain_i, std::size_t grain_j) const
{
  if (grain_i == grain_j)
    mooseError("GBMisorientationBase::packedPairIndex: identical grain ids grain_i=",
               grain_i,
               " grain_j=",
               grain_j);

  const auto n_grains = grainNum();
  if (grain_i >= n_grains || grain_j >= n_grains)
    mooseError("GBMisorientationBase::packedPairIndex: grain id outside Euler table. grain_i=",
               grain_i,
               " grain_j=",
               grain_j,
               " grainNum=",
               n_grains);

  return GBPairPacking::pack_unordered(grain_i, grain_j);
}

const GBMisorientationBase::MisorientationData &
GBMisorientationBase::getMisorientationData(std::size_t grain_i, std::size_t grain_j) const
{
  const auto k = packedPairIndex(grain_i, grain_j);

  if (k >= _miso_table.size())
    mooseError("GBMisorientationBase::getMisorientationData: packed index k=",
               k,
               " exceeds table size=",
               _miso_table.size(),
               " for grain_i=",
               grain_i,
               " grain_j=",
               grain_j);

  return _miso_table[k];
}

Real
GBMisorientationBase::getMisorientationAngle(std::size_t grain_i, std::size_t grain_j) const
{
  return getMisorientationData(grain_i, grain_j).theta;
}

Real
GBMisorientationBase::getMisorientationPolarAxis(std::size_t grain_i, std::size_t grain_j) const
{
  return getMisorientationData(grain_i, grain_j).polar_ax;
}

Real
GBMisorientationBase::getMisorientationAzimuthAxis(std::size_t grain_i, std::size_t grain_j) const
{
  return getMisorientationData(grain_i, grain_j).azim_ax;
}

void
GBMisorientationBase::buildMisorientationTable()
{
  defineSymmetryOperator();

  const auto n_grains = grainNum();

  _euler_angle.resize(n_grains);
  _quat_angle.resize(n_grains);

  for (std::size_t i = 0; i < n_grains; ++i)
  {
    _euler_angle[i] = _euler.getEulerAngles(i);

    // In the current EulerAngleTxtFileReader workflow, getEulerAngles() supplies radians,
    // while EulerAngles::toQuaternion() expects degrees.
    EulerAngles e_deg;
    e_deg.phi1 = _euler_angle[i].phi1 * 180.0 / libMesh::pi;
    e_deg.Phi = _euler_angle[i].Phi * 180.0 / libMesh::pi;
    e_deg.phi2 = _euler_angle[i].phi2 * 180.0 / libMesh::pi;

    _quat_angle[i] = e_deg.toQuaternion();
  }

  const std::size_t n_pairs = n_grains * (n_grains - 1) / 2;
  _miso_table.clear();
  _miso_table.resize(n_pairs);

  for (std::size_t j = 1; j < n_grains; ++j)
    for (std::size_t i = 0; i < j; ++i)
    {
      const auto k = GBPairPacking::pack_upper(i, j);
      _miso_table[k] = getMisorientationFromQuaternion(_quat_angle[i], _quat_angle[j]);
    }
}

void
GBMisorientationBase::rotationSymmetryToQuaternion(const double O[3][3],
                                                   Eigen::Quaternion<Real> & q)
{
  double q4 = 0.0;

  if ((1.0 + O[0][0] + O[1][1] + O[2][2]) > 0.0)
  {
    q4 = std::sqrt(1.0 + O[0][0] + O[1][1] + O[2][2]) / 2.0;
    q.w() = q4;
    q.x() = (O[2][1] - O[1][2]) / (4.0 * q4);
    q.y() = (O[0][2] - O[2][0]) / (4.0 * q4);
    q.z() = (O[1][0] - O[0][1]) / (4.0 * q4);
  }
  else if ((1.0 + O[0][0] - O[1][1] - O[2][2]) > 0.0)
  {
    q4 = std::sqrt(1.0 + O[0][0] - O[1][1] - O[2][2]) / 2.0;
    q.w() = (O[2][1] - O[1][2]) / (4.0 * q4);
    q.x() = q4;
    q.y() = (O[1][0] + O[0][1]) / (4.0 * q4);
    q.z() = (O[0][2] + O[2][0]) / (4.0 * q4);
  }
  else if ((1.0 - O[0][0] + O[1][1] - O[2][2]) > 0.0)
  {
    q4 = std::sqrt(1.0 - O[0][0] + O[1][1] - O[2][2]) / 2.0;
    q.w() = (O[0][2] - O[2][0]) / (4.0 * q4);
    q.x() = (O[1][0] + O[0][1]) / (4.0 * q4);
    q.y() = q4;
    q.z() = (O[2][1] + O[1][2]) / (4.0 * q4);
  }
  else if ((1.0 - O[0][0] - O[1][1] + O[2][2]) > 0.0)
  {
    q4 = std::sqrt(1.0 - O[0][0] - O[1][1] + O[2][2]) / 2.0;
    q.w() = (O[1][0] - O[0][1]) / (4.0 * q4);
    q.x() = (O[0][2] + O[2][0]) / (4.0 * q4);
    q.y() = (O[2][1] + O[1][2]) / (4.0 * q4);
    q.z() = q4;
  }
  else
    mooseError("GBMisorientationBase::rotationSymmetryToQuaternion: invalid symmetry matrix.");
}

void
GBMisorientationBase::defineSymmetryOperator()
{
  _sym_quat.resize(_o_sym);

  // Cubic crystal symmetry operators.
  double sym_rotation[24][3][3] = {
      {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}},    {{1, 0, 0}, {0, -1, 0}, {0, 0, -1}},
      {{1, 0, 0}, {0, 0, -1}, {0, 1, 0}},   {{1, 0, 0}, {0, 0, 1}, {0, -1, 0}},
      {{-1, 0, 0}, {0, 1, 0}, {0, 0, -1}},  {{-1, 0, 0}, {0, -1, 0}, {0, 0, 1}},
      {{-1, 0, 0}, {0, 0, -1}, {0, -1, 0}}, {{-1, 0, 0}, {0, 0, 1}, {0, 1, 0}},
      {{0, 1, 0}, {-1, 0, 0}, {0, 0, 1}},   {{0, 1, 0}, {0, 0, -1}, {-1, 0, 0}},
      {{0, 1, 0}, {1, 0, 0}, {0, 0, -1}},   {{0, 1, 0}, {0, 0, 1}, {1, 0, 0}},
      {{0, -1, 0}, {1, 0, 0}, {0, 0, 1}},   {{0, -1, 0}, {0, 0, -1}, {1, 0, 0}},
      {{0, -1, 0}, {-1, 0, 0}, {0, 0, -1}}, {{0, -1, 0}, {0, 0, 1}, {-1, 0, 0}},
      {{0, 0, 1}, {0, 1, 0}, {-1, 0, 0}},   {{0, 0, 1}, {1, 0, 0}, {0, 1, 0}},
      {{0, 0, 1}, {0, -1, 0}, {1, 0, 0}},   {{0, 0, 1}, {-1, 0, 0}, {0, -1, 0}},
      {{0, 0, -1}, {0, 1, 0}, {1, 0, 0}},   {{0, 0, -1}, {-1, 0, 0}, {0, 1, 0}},
      {{0, 0, -1}, {0, -1, 0}, {-1, 0, 0}}, {{0, 0, -1}, {1, 0, 0}, {0, -1, 0}}};

  for (int o = 0; o < _o_sym; ++o)
    rotationSymmetryToQuaternion(sym_rotation[o], _sym_quat[o]);
}

GBMisorientationBase::MisorientationData
GBMisorientationBase::getMisorientationFromQuaternion(const Eigen::Quaternion<Real> & qi_in,
                                                      const Eigen::Quaternion<Real> & qj_in) const
{
  Eigen::Quaternion<Real> qi = qi_in;
  Eigen::Quaternion<Real> qj = qj_in;

  qi.normalize();
  qj.normalize();

  constexpr Real eps = 1e-8;

  Real best_theta = 2.0 * libMesh::pi;
  Eigen::Quaternion<Real> qmin(1.0, 0.0, 0.0, 0.0);

  for (int o1 = 0; o1 < _o_sym; ++o1)
  {
    for (int o2 = 0; o2 < _o_sym; ++o2)
    {
      Eigen::Quaternion<Real> qib = _sym_quat[o1] * qi;
      Eigen::Quaternion<Real> qjb = _sym_quat[o2] * qj;

      qib.normalize();
      qjb.normalize();

      Eigen::Quaternion<Real> q = qib * qjb.conjugate();
      q.normalize();

      const bool in_fz =
          (0.0 <= q.z()) && (q.z() <= q.y()) && (q.y() <= q.x()) && (q.x() <= 1.0);

      const bool inv_in_fz =
          (0.0 <= -q.z()) && (-q.z() <= -q.y()) && (-q.y() <= -q.x()) && (-q.x() <= 1.0);

      if (!(in_fz || inv_in_fz))
        continue;

      const Real w = std::clamp(std::abs(q.w()), Real(-1.0), Real(1.0));
      const Real theta = 2.0 * std::acos(w);

      if (theta < best_theta - eps)
      {
        best_theta = theta;
        qmin = in_fz ? q : q.conjugate();
      }
    }
  }

  RealGradient axis(qmin.x(), qmin.y(), qmin.z());

  Real polar_ax = 0.0;
  Real azim_ax = 0.0;

  const Real vnorm = axis.norm();

  if (vnorm > 1e-12)
  {
    axis /= vnorm;

    polar_ax = std::acos(std::clamp(axis(2), Real(-1.0), Real(1.0)));

    azim_ax = std::atan2(axis(1), axis(0));
    if (azim_ax < 0.0)
      azim_ax += 2.0 * libMesh::pi;
  }

  MisorientationData out;
  out.theta = best_theta;
  out.polar_ax = polar_ax;
  out.azim_ax = azim_ax;
  out.q = qmin;
  out.qnorm = vnorm;

  return out;
}
