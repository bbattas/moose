#include "GBMisorientation.h"
#include "SolutionUserObject.h"

registerMooseObject("PhaseFieldApp", GBMisorientation);

InputParameters
GBMisorientation::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Calculate types of grain boundaries in a polycrystalline sample");
  params.addRequiredParam<UserObjectName>("grain_tracker",
                                          "The GrainTracker UserObject to get values from.");
  params.addRequiredParam<UserObjectName>("ebsd_reader", "The EBSDReader GeneralUserObject");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  params.addParam<Real>("miso_weight", 1, "0-1 Weight for misorientation in calcs.");
  params.addParam<MaterialPropertyName>(
      "gb_energy_iso_name", "sigma", "Isotropic GB energy before inclination dependence.");
  params.addParam<MaterialPropertyName>(
      "kappa", "kappa", "Gradient energy constant kappa material name.");
  params.addParam<MaterialPropertyName>(
      "mu", "mu", "Free energy thermodynamic parameter mu (or m in some formulations).");
  return params;
}

GBMisorientation::GBMisorientation(const InputParameters & parameters)
  : Material(parameters),
    _grain_tracker(getUserObject<GrainTracker>("grain_tracker")),
    _ebsd_reader(getUserObject<EBSDReader>("ebsd_reader")),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    // INPUTS
    _m_weight(getParam<Real>("miso_weight")),
    _gbe_iso(getMaterialProperty<Real>("gb_energy_iso_name")),
    _kappa(getMaterialProperty<Real>(getParam<MaterialPropertyName>("kappa"))),
    _mu(getMaterialProperty<Real>(getParam<MaterialPropertyName>("mu"))),
    _int_width(declareProperty<Real>("int_width")),
    _gamma(declareProperty<Real>("gamma_asymm")),
    // TESTING
    _gtnum(declareProperty<Real>("gt_num")),
    _other_out(declareProperty<Real>("other_out")),
    // EULER ANGLES
    _eul_a(declareProperty<Real>("eul_a")),
    _eul_b(declareProperty<Real>("eul_b")),
    _eul_c(declareProperty<Real>("eul_c")),
    // QUATERNIONS
    _quat_a(declareProperty<Real>("quat_a")),
    _quat_b(declareProperty<Real>("quat_b")),
    _quat_c(declareProperty<Real>("quat_c")),
    _quat_d(declareProperty<Real>("quat_d")),
    _quat_mag(declareProperty<Real>("quat_mag")),
    // Misorientation out
    _miso(declareProperty<Real>("misorientation")),
    _miso_polar(declareProperty<Real>("miso_axis_polar")),
    _miso_azim(declareProperty<Real>("miso_axis_azimuth")),
    // Energies
    _miso_ang_en(declareProperty<Real>("miso_ang_energy")),
    _miso_ax_en(declareProperty<Real>("miso_ax_energy")),
    _twist(declareProperty<Real>("twist_energy")),
    _tilt(declareProperty<Real>("tilt_energy")),
    _f_mis(declareProperty<Real>("f_miso"))
{
  getMisorientationAngles();
}

void
GBMisorientation::computeQpProperties()
{
  // Find out the number of boundary unique_id and save them
  _gb_pairs.clear();
  _gb_op_pairs.clear();
  //
  // _gb_i_vals.clear();

  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());
  for (auto i : index_range(op_to_grains))
  {
    if (op_to_grains[i] == FeatureFloodCount::invalid_id)
      continue;

    _gb_pairs.push_back(_ebsd_reader.getFeatureID(op_to_grains[i]));
    _gb_op_pairs.push_back((*_vals[i])[_qp]);
    // _gb_i_vals.push_back(i);
  }

  // Zero testing outputs
  _eul_a[_qp] = 0;
  _eul_b[_qp] = 0;
  _eul_c[_qp] = 0;
  _quat_a[_qp] = 0;
  _quat_b[_qp] = 0;
  _quat_c[_qp] = 0;
  _quat_d[_qp] = 0;
  _quat_mag[_qp] = 0;

  _gtnum[_qp] = _gb_pairs.size();
  _other_out[_qp] = 0;

  _miso_ang_en[_qp] = 0;
  _miso_ax_en[_qp] = 0;
  Real ang_cut = 62 * libMesh::pi / 180;
  _twist[_qp] = 0;
  _tilt[_qp] = 0;
  _f_mis[_qp] = 1;

  // Compute GB type by the number of id
  switch (_gb_pairs.size())
  {
    case 0:
    {
      _miso[_qp] = 0;
      _miso_polar[_qp] = 0;
      _miso_azim[_qp] = 0;
      _twist[_qp] = 0.7; // * _m_weight + (1 - _m_weight);
      _tilt[_qp] = 0.3;  // * _m_weight + (1 - _m_weight);
      _f_mis[_qp] = 1.0; // 0.3 + _twist[_qp];
      break;
    }

    case 1:
    {
      // _gb_pairs[0]
      auto e_angle = _ebsd_reader.getEulerAngles(_gb_pairs[0]);
      _eul_a[_qp] = e_angle.phi1;
      _eul_b[_qp] = e_angle.Phi;
      _eul_c[_qp] = e_angle.phi2;
      auto quat = e_angle.toQuaternion();
      _quat_a[_qp] = quat.w();
      _quat_b[_qp] = quat.x();
      _quat_c[_qp] = quat.y();
      _quat_d[_qp] = quat.z();
      _quat_mag[_qp] = quat.norm();
      //
      _miso[_qp] = 0;
      _miso_polar[_qp] = 0;
      _miso_azim[_qp] = 0;
      _twist[_qp] = 0.7; // * _m_weight + (1 - _m_weight);
      _tilt[_qp] = 0.3;  // * _m_weight + (1 - _m_weight);
      _f_mis[_qp] = 1.0; // 0.3 + _twist[_qp];
      break;
    }

    case 2:
    {
      auto idx = getLineNum(_gb_pairs[0], _gb_pairs[1]);
      _miso[_qp] = _misorientation_angles[idx];
      _miso_polar[_qp] = _miso_ax_polar[idx];
      _miso_azim[_qp] = _miso_ax_azimuth[idx];
      // Energies
      _miso_ang_en[_qp] = (_miso[_qp] / ang_cut) * (1 - std::log(_miso[_qp] / ang_cut));
      _miso_ax_en[_qp] = std::pow(std::abs(std::cos(_miso_polar[_qp])), 0.4) +
                         std::pow(std::abs(std::cos(_miso_azim[_qp] / 2)), 0.4);
      if (_miso_ang_en[_qp] > 1.0)
        _miso_ang_en[_qp] = 1.0;
      if (_miso_ax_en[_qp] > 1.0)
        _miso_ax_en[_qp] = 1.0;
      _twist[_qp] = 0.7 * _miso_ax_en[_qp] * _miso_ang_en[_qp] * _m_weight + (1 - _m_weight);
      _tilt[_qp] = 0.3 * _miso_ax_en[_qp] * _miso_ang_en[_qp] * _m_weight + (1 - _m_weight);
      // Normalized GBE
      _f_mis[_qp] = 0.3 + _twist[_qp];
      // OTHER
      auto quat = _qmin[idx];
      _quat_a[_qp] = quat.w();
      _quat_b[_qp] = quat.x();
      _quat_c[_qp] = quat.y();
      _quat_d[_qp] = quat.z();
      _quat_mag[_qp] = quat.norm();
      _other_out[_qp] = _other[idx];
      break;
    }

    default:
    {
      // _miso[_qp] = 0;
      _miso_polar[_qp] = 0;
      _miso_azim[_qp] = 0;
      // All grain pairs- Using average approach as per Yang 2025
      Real tj_miso = 0.0;
      Real tj_cross = 0.0;
      for (std::size_t idx1 = 0; idx1 < _gb_pairs.size(); ++idx1)
        for (std::size_t idx2 = idx1 + 1; idx2 < _gb_pairs.size(); ++idx2)
        {
          auto idx = getLineNum(_gb_pairs[idx1], _gb_pairs[idx2]);
          Real ang_en = (_misorientation_angles[idx] / ang_cut) *
                        (1 - std::log(_misorientation_angles[idx] / ang_cut));
          Real ax_en = std::pow(std::abs(std::cos(_miso_ax_polar[idx])), 0.4) +
                       std::pow(std::abs(std::cos(_miso_ax_azimuth[idx] / 2)), 0.4);
          if (ang_en > 1.0)
            ang_en = 1.0;
          if (ax_en > 1.0)
            ax_en = 1.0;
          Real cross = ax_en * ang_en; // * _m_weight + (1 - _m_weight);
          // Outputs
          tj_miso += _misorientation_angles[idx];
          tj_cross += cross;
        }
      // Average
      _miso[_qp] = tj_miso / _gb_pairs.size();
      _twist[_qp] = 0.7 * (tj_cross / _gb_pairs.size()) * _m_weight + (1 - _m_weight);
      _tilt[_qp] = 0.3 * (tj_cross / _gb_pairs.size()) * _m_weight + (1 - _m_weight);
      _f_mis[_qp] = 0.3 + _twist[_qp];
    }
  }

  // Calculate outputs from misorientation using Moelans approach
  // Hard-coded coefficients (for poly gamma)
  constexpr Real a1 = -3.0944; // coefficient for g2^4
  constexpr Real a2 = -1.8169; // coefficient for g2^3
  constexpr Real a3 = 10.323;  // coefficient for g2^2
  constexpr Real a4 = -8.1819; // coefficient for g2
  constexpr Real a5 = 2.0033;  // constant term
  // Gamma
  Real g = _gbe_iso[_qp] * _f_mis[_qp] / (std::sqrt(_kappa[_qp] * _mu[_qp]));
  Real g2 = g * g;
  Real pg = (((a1 * g2 + a2) * g2 + a3) * g2 + a4) * g2 + a5;
  _gamma[_qp] = 1 / pg; // 1.5;
  // IW
  Real f0_int =
      (((((0.0788 * pg - 0.4955) * pg + 1.2244) * pg - 1.5281) * pg + 1.0686) * pg - 0.5563) * pg +
      0.2907;
  _int_width[_qp] = (std::sqrt(_kappa[_qp] / _mu[_qp])) * (std::sqrt(1 / f0_int));
}
// Function to output total line number of Misorientation angle file
unsigned int
GBMisorientation::getTotalLineNum() const
{
  return _misorientation_angles.size();
}

// Function to output specific line number in Misorientation angle file
unsigned int
GBMisorientation::getLineNum(unsigned int grain_i, unsigned int grain_j)
{
  if (grain_i > grain_j)
    return grain_j + (grain_i - 1) * grain_i / 2;
  else
    return grain_i + (grain_j - 1) * grain_j / 2;
}

// Function to convert symmetry matrix to quaternion form
void
GBMisorientation::rotationSymmetryToQuaternion(const double O[3][3], Eigen::Quaternion<Real> & q)
{
  double q4 = 0;
  if ((1 + O[0][0] + O[1][1] + O[2][2]) > 0)
  {
    q4 = sqrt(1 + O[0][0] + O[1][1] + O[2][2]) / 2;
    q.w() = q4;
    q.x() = (O[2][1] - O[1][2]) / (4 * q4);
    q.y() = (O[0][2] - O[2][0]) / (4 * q4);
    q.z() = (O[1][0] - O[0][1]) / (4 * q4);
  }
  else if ((1 + O[0][0] - O[1][1] - O[2][2]) > 0)
  {
    q4 = sqrt(1 + O[0][0] - O[1][1] - O[2][2]) / 2;
    q.w() = (O[2][1] - O[1][2]) / (4 * q4);
    q.x() = q4;
    q.y() = (O[1][0] + O[0][1]) / (4 * q4);
    q.z() = (O[0][2] + O[2][0]) / (4 * q4);
  }
  else if ((1 - O[0][0] + O[1][1] - O[2][2]) > 0)
  {
    q4 = sqrt(1 - O[0][0] + O[1][1] - O[2][2]) / 2;
    q.w() = (O[0][2] - O[2][0]) / (4 * q4);
    q.x() = (O[1][0] + O[0][1]) / (4 * q4);
    q.y() = q4;
    q.z() = (O[2][1] + O[1][2]) / (4 * q4);
  }
  else if ((1 - O[0][0] - O[1][1] + O[2][2]) > 0)
  {
    q4 = sqrt(1 - O[0][0] - O[1][1] + O[2][2]) / 2;
    q.w() = (O[1][0] - O[0][1]) / (4 * q4);
    q.x() = (O[0][2] + O[2][0]) / (4 * q4);
    q.y() = (O[2][1] + O[1][2]) / (4 * q4);
    q.z() = q4;
  }
}

// Function to define the symmetry operator
void
GBMisorientation::defineSymmetryOperator()
{
  // grow by number of symmetric operators
  _sym_quat.resize(_o_sym);

  // cubic symmetry
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

  // initialize global operator
  for (int o = 0; o < _o_sym; o++)
    rotationSymmetryToQuaternion(sym_rotation[o], _sym_quat[o]);
}

GBMisorientation::MisorientationResult
GBMisorientation::getMisorientationFromQuaternion(const Eigen::Quaternion<Real> & qi_in,
                                                  const Eigen::Quaternion<Real> & qj_in)
{
  Eigen::Quaternion<Real> qi = qi_in;
  Eigen::Quaternion<Real> qj = qj_in;
  qi.normalize();
  qj.normalize();
  constexpr double eps = 1e-8;

  Real best_theta = 2 * libMesh::pi; // std::numeric_limits<Real>::max();
  Eigen::Quaternion<Real> qmin(1, 0, 0, 0);

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

      // Check if in functional zone
      bool in_fz = (0.0 <= q.z()) && (q.z() <= q.y()) && (q.y() <= q.x()) && (q.x() <= 1.0);
      bool inv_in_fz =
          (0.0 <= -q.z()) && (-q.z() <= -q.y()) && (-q.y() <= -q.x()) && (-q.x() <= 1.0);
      if (!(in_fz || inv_in_fz))
        continue;

      Real w = std::clamp(std::abs(q.w()), Real(-1), Real(1));
      Real theta = 2 * std::acos(w); // in [0, pi], pi/2 with the abs(q.w()) above

      if ((in_fz || inv_in_fz) && (theta < (best_theta - eps)))
      {
        best_theta = theta;
        qmin = in_fz ? q : q.conjugate();
      }
    }
  }

  // Shouldnt need this if forcing + and - test for in fz
  // // Fix sign to stabilize axis direction
  // if (qmin.w() < 0)
  //   qmin.coeffs() *= -1;

  // Axis from vector part
  RealGradient axis(qmin.x(), qmin.y(), qmin.z());
  Real polar_ax = 0.0;
  Real azim_ax = 0.0;

  Real vnorm = axis.norm();
  if (vnorm > 1e-12)
  {
    axis /= vnorm;
    // // theta_ax = std::acos(std::clamp(axis(2), Real(-1), Real(1))); // full 0-pi
    // theta_ax = std::acos(std::clamp(std::abs(axis(2)), 0.0, 1.0)); // half 0-pi/2
    // phi_ax = std::atan2(axis(1), axis(0));
    // if (phi_ax < 0)
    //   phi_ax += 2 * libMesh::pi;
    // New Way from VECTOR
    polar_ax = std::acos(axis(2)); // half 0-pi/2
    // azim_ax = std::atan2(axis(1), axis(0)) + libMesh::pi; // Lins method in VECTOR?
    azim_ax = std::atan2(axis(1), axis(0));
    if (azim_ax < 0)
      azim_ax += 2 * libMesh::pi;
  }

  MisorientationResult out;
  out.theta = best_theta;
  out.polar_ax = polar_ax;
  out.azim_ax = azim_ax;
  out.q = qmin;
  out.qnorm = vnorm; // will be ~1 because we normalized
  return out;
}

// Get the Misorientation angle
void
GBMisorientation::getMisorientationAngles()
{
  // Initialize symmetry operator as quaternion vectors
  defineSymmetryOperator();
  // Initialize parameters to calculate misorientation
  const auto grain_num = _ebsd_reader.getGrainNum();
  _euler_angle.resize(grain_num);
  _quat_angle.resize(grain_num);

  // Get Euler Angle of Orientation
  for (const auto i : make_range(grain_num))
  {
    auto grain_id = _ebsd_reader.getFeatureID(i);
    mooseAssert(grain_id < grain_num, "Feature ids cannot exceed max grain number");
    _euler_angle[grain_id] = _ebsd_reader.getEulerAngles(i);
    _quat_angle[grain_id] = _euler_angle[grain_id].toQuaternion();
  }

  for (const auto j : make_range(std::make_unsigned_t<int>(1), grain_num))
  {
    for (const auto i : make_range(j))
    {
      auto miso = getMisorientationFromQuaternion(_quat_angle[i], _quat_angle[j]);
      _misorientation_angles.push_back(miso.theta); // / libMesh::pi * 180);
      _miso_ax_polar.push_back(miso.polar_ax);
      _miso_ax_azimuth.push_back(miso.azim_ax);
      _qmin.push_back(miso.q);
      _other.push_back(miso.qnorm);
    }
  }
}
