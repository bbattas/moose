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
  params.addParam<Real>("angle_threshold", 15, "Max LAGB Misorientation angle");
  params.addParam<Real>("miso_weight", 1, "0-1 Weight for misorientation in calcs.");
  return params;
}

GBMisorientation::GBMisorientation(const InputParameters & parameters)
  : Material(parameters),
    _grain_tracker(getUserObject<GrainTracker>("grain_tracker")),
    _ebsd_reader(getUserObject<EBSDReader>("ebsd_reader")),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _angle_threshold(getParam<Real>("angle_threshold")),
    // _gb_type(declareADProperty<Real>("gb_type")),
    // INPUTS
    _m_weight(getParam<Real>("miso_weight")),
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
    _miso_theta(declareProperty<Real>("miso_axis_polar")),
    _miso_phi(declareProperty<Real>("miso_axis_azimuth")),
    // Energies
    _miso_ang_en(declareProperty<Real>("miso_ang_energy")),
    _miso_ax_en(declareProperty<Real>("miso_ax_energy")),
    _twist(declareProperty<Real>("twist_energy")),
    _tilt(declareProperty<Real>("tilt_energy"))
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

  // Compute GB type by the number of id
  switch (_gb_pairs.size())
  {
    case 0:
    {
      _miso[_qp] = 0;
      _miso_theta[_qp] = 0;
      _miso_phi[_qp] = 0;
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
      _miso[_qp] = 0;
      _miso_theta[_qp] = 0;
      _miso_phi[_qp] = 0;
      break;
    }

    case 2:
    {
      auto idx = getLineNum(_gb_pairs[0], _gb_pairs[1]);
      _miso[_qp] = _misorientation_angles[idx];
      _miso_theta[_qp] = _miso_ax_polar[idx];
      _miso_phi[_qp] = _miso_ax_azimuth[idx];
      // Energies
      _miso_ang_en[_qp] = (_miso[_qp] / ang_cut) * (1 - std::log(_miso[_qp] / ang_cut));
      _miso_ax_en[_qp] = std::pow(std::abs(std::cos(_miso_theta[_qp])), 0.4) +
                         std::pow(std::abs(std::cos(_miso_phi[_qp] / 2)), 0.4);
      _twist[_qp] = 0.7 * _miso_ax_en[_qp] * _miso_ang_en[_qp] * _m_weight + (1 - _m_weight);
      _tilt[_qp] = 0.3 * _miso_ax_en[_qp] * _miso_ang_en[_qp] * _m_weight + (1 - _m_weight);
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
      _miso[_qp] = 0;
      _miso_theta[_qp] = 0;
      _miso_phi[_qp] = 0;
      // _gb_type[_qp] = getTripleJunctionType();
      // _quat_a[_qp] = 2;
      // _quat_b[_qp] = 0;
      // _quat_c[_qp] = 0;
      // _quat_d[_qp] = 0;
    }
  }
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

// // Function to calculate the GB type in Triple junction
// Real
// GBMisorientation::getTripleJunctionType()
// {
//   unsigned int lagb_num = 0;
//   unsigned int hagb_num = 0;
//   Real ratio_base = 0.0;
//   Real ratio_lagb = 0.0;
//   for (unsigned int i = 1; i < _gb_pairs.size(); ++i)
//   {
//     for (unsigned int j = 0; j < i; ++j)
//     {
//       ratio_base += (_gb_op_pairs[i] * _gb_op_pairs[i] * _gb_op_pairs[j] * _gb_op_pairs[j]);
//       if (_misorientation_angles[getLineNum(_gb_pairs[j], _gb_pairs[i])] < _angle_threshold)
//       {
//         lagb_num += 1;
//         ratio_lagb += (_gb_op_pairs[i] * _gb_op_pairs[i] * _gb_op_pairs[j] * _gb_op_pairs[j]);
//       }
//       else
//         hagb_num += 1;
//     }
//   }
//   if (lagb_num == 0)
//     return 2;
//   else if (hagb_num == 0)
//     return 1;
//   else
//     return 2 - ratio_lagb / ratio_base;
// }

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

// // Function to return the misorientation of two quaternions
// GBMisorientation::MisorientationResult
// GBMisorientation::getMisorientationFromQuaternion(const Eigen::Quaternion<Real> & qi,
//                                                   const Eigen::Quaternion<Real> & qj)
// {
//   Real miso0, misom = 2.0 * libMesh::pi;
//   Eigen::Quaternion<Real> q, qib, qjb, qmin;
//   qmin.w() = 0;
//   qmin.x() = 0;
//   qmin.y() = 0;
//   qmin.z() = 0;

//   for (int o1 = 0; o1 < _o_sym; o1++)
//   {
//     for (int o2 = 0; o2 < _o_sym; o2++)
//     {
//       qib = _sym_quat[o1] * qi;
//       qjb = _sym_quat[o2] * qj;

//       // j-grain conjugate quaternion
//       qjb.x() = -qjb.x();
//       qjb.y() = -qjb.y();
//       qjb.z() = -qjb.z();
//       q = qib * qjb;
//       miso0 = round(2 * std::acos(q.w()) * 1e5) / 1e5;

//       if (miso0 > libMesh::pi)
//         miso0 = miso0 - 2 * libMesh::pi;
//       if (std::abs(miso0) < misom)
//       {
//         misom = std::abs(miso0);
//         qmin = q;
//       }
//     }
//   }

//   miso0 = 2 * std::acos(qmin.w());
//   if (miso0 > libMesh::pi)
//     miso0 = miso0 - 2 * libMesh::pi;
//   double theta = std::abs(miso0);

//   // return std::abs(miso0);
//   // Angles for polar/azimuth on q axis
//   // Fix sign ambiguity so w >= 0 (optional but recommended)
//   if (qmin.w() < 0.0)
//     qmin.coeffs() *= -1.0; // coeffs() = (x,y,z,w) in Eigen
//   RealGradient axis(qmin.x(), qmin.y(), qmin.z());
//   double s = std::sin(theta / 2.0);
//   double theta_ax = 0.0;
//   double phi_ax = 0.0;
//   if (s > 1e-12)
//   {
//     // remove sin(theta/2) from quat vector pieces
//     axis /= s;
//     axis /= axis.norm();
//     // Spherical angles of the axis relative to +z
//     theta_ax = std::acos(std::clamp(axis.z(), -1.0, 1.0)); // polar
//     phi_ax = std::atan2(axis.y(), axis.x());               // azimuth
//     if (phi_ax < 0.0)
//       phi_ax += 2 * libMesh::pi;
//   }

//   // OUTPUT
//   MisorientationResult m_out;
//   m_out.theta = theta;
//   m_out.theta_ax = theta_ax;
//   m_out.phi_ax = phi_ax;
//   m_out.q = qmin;
//   m_out.qnorm = qmin.norm();

//   return m_out;
// }

GBMisorientation::MisorientationResult
GBMisorientation::getMisorientationFromQuaternion(const Eigen::Quaternion<Real> & qi_in,
                                                  const Eigen::Quaternion<Real> & qj_in)
{
  Eigen::Quaternion<Real> qi = qi_in;
  Eigen::Quaternion<Real> qj = qj_in;
  qi.normalize();
  qj.normalize();

  Real best_theta = std::numeric_limits<Real>::max();
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

      Real w = std::clamp(std::abs(q.w()), Real(-1), Real(1));
      Real theta = 2 * std::acos(w); // in [0, pi]

      if (theta < best_theta)
      {
        best_theta = theta;
        qmin = q;
      }
    }
  }

  // Fix sign to stabilize axis direction
  if (qmin.w() < 0)
    qmin.coeffs() *= -1;

  // Axis from vector part
  RealGradient axis(qmin.x(), qmin.y(), qmin.z());
  Real theta_ax = 0.0;
  Real phi_ax = 0.0;

  Real vnorm = axis.norm();
  if (vnorm > 1e-12)
  {
    axis /= vnorm;
    // theta_ax = std::acos(std::clamp(axis(2), Real(-1), Real(1))); // full 0-pi
    theta_ax = std::acos(std::clamp(std::abs(axis(2)), 0.0, 1.0)); // half 0-pi/2
    phi_ax = std::atan2(axis(1), axis(0));
    if (phi_ax < 0)
      phi_ax += 2 * libMesh::pi;
  }

  MisorientationResult out;
  out.theta = best_theta;
  out.theta_ax = theta_ax;
  out.phi_ax = phi_ax;
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
      _miso_ax_polar.push_back(miso.theta_ax);
      _miso_ax_azimuth.push_back(miso.phi_ax);
      _qmin.push_back(miso.q);
      _other.push_back(miso.qnorm);
    }
  }
}
