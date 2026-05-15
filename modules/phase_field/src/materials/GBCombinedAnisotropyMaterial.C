#include "GBCombinedAnisotropyMaterial.h"

registerMooseObject("PhaseFieldApp", GBCombinedAnisotropyMaterial);

InputParameters
GBCombinedAnisotropyMaterial::validParams()
{
  InputParameters params = GBInclinationBase::validParams();

  params += GBMisorientationHelper::validParams();

  // Settings/Modes
  MooseEnum gb_mode("cos=0 inc=1 miso=2 full=3 iso=4", "full");
  params.addParam<MooseEnum>("gb_mode",
                             gb_mode,
                             "Which grain-boundary property model to use. "
                             "'cos' and 'inc' use only inclination data. "
                             "'miso' and 'full' also use misorientation data.");

  MooseEnum combine_form("avg=0 weighted=1", "avg");
  params.addParam<MooseEnum>(
      "tj_mode", combine_form, "Which TJ combination approach, weighted average or just average.");
  params.addParam<bool>(
      "L_2nd", false, "Is L 2nd order with full derivatives, requires 2nd order mesh and OPs.");
  params.addParam<bool>(
      "stiffness", true, "Include the 1st/2nd derivative of gbe wrt inc/gradeta, else = 0.");
  params.addParam<bool>("aniso_gbmob",
                        false,
                        "Apply gbe anisotropy also to the gb mobility (inferred: L*f -> L*f*f).");
  // Isotropic Input Properties
  params.addParam<MaterialPropertyName>(
      "kappa_name", "kappa", "Gradient energy constant kappa material name.");
  params.addParam<MaterialPropertyName>(
      "gbe_iso_name", "sigma_iso", "Isotropic GB energy before anisotropy is applied.");
  // VALUES
  params.addParam<Real>(
      "bulk_scalar",
      1.0,
      "Normalized GBE multiplier to use in non-gb regions for the interconnected calculations.");
  params.addParam<Real>("w_inc", 1, "Weight for inclination contribution in combined form [0-1].");
  params.addParam<Real>(
      "w_miso", 1, "Weight for misorientation contribution in combined form [0-1].");
  // Constants for COS inc function
  params.addParam<Real>("ifunc_a", 0.05, "Inclination function constant a.");
  // params.addParam<Real>("ifunc_b", 2, "Inclination function constant b.");
  // params.addParam<Real>("ifunc_c", 0.0, "Inclination function constant c.");

  params.addClassDescription(
      "Parent material combining inclination and optional misorientation data "
      "to compute pairwise grain-boundary properties.");

  return params;
}

GBCombinedAnisotropyMaterial::GBCombinedAnisotropyMaterial(const InputParameters & parameters)
  : GBInclinationBase(parameters),
    GBMisorientationHelper(parameters,
                           parameters.isParamValid("euler_angle_provider")
                               ? &getUserObject<EulerAngleProvider>("euler_angle_provider")
                               : nullptr),
    _gb_mode(getParam<MooseEnum>("gb_mode")),
    _tj_mode(getParam<MooseEnum>("tj_mode")),
    _aniso_L(getParam<bool>("L_2nd")),
    _stiffness(getParam<bool>("stiffness")),
    _aniso_mob(getParam<bool>("aniso_gbmob")),
    // Isotropic Inputs
    _kappa(getMaterialProperty<Real>("kappa_name")),
    _gbe_iso(getMaterialProperty<Real>("gbe_iso_name")),
    _mu(getMaterialProperty<Real>("mu")),
    _L0(getMaterialProperty<Real>("L0")),
    // Anisotropic Outputs
    _fgbe(declareProperty<std::vector<Real>>("fgbe")),
    _gbe_norm(declareProperty<Real>("gbe_norm")),
    // _L_ij(declareProperty<std::vector<Real>>("L_ij")),
    _L(declareProperty<Real>("L")),
    _gamma_ij(declareProperty<std::vector<Real>>("gamma_ij")),
    _dgamma_dgradeta(declareProperty<std::vector<RealGradient>>("dgamma_dgradeta")),
    _d2gamma_dgradeta2(declareProperty<std::vector<RealTensorValue>>("d2gamma_dgradeta2")),
    _gamma_asymm(declareProperty<Real>("gamma_asymm")),
    _int_width(declareProperty<Real>("int_width")),
    // CONSTANTS
    _bulk_mult(getParam<Real>("bulk_scalar")),
    _w_inc(getParam<Real>("w_inc")),
    _w_miso(getParam<Real>("w_miso")),
    // COS FUNCTION CONSTANTS
    _if_a(getParam<Real>("ifunc_a")),
    // _if_b(getParam<Real>("ifunc_b")),
    // _if_c(getParam<Real>("ifunc_c")),
    // HARDCODED skip param needed in kernel but should change that later
    _elem_noij(declareProperty<bool>("elem_no_ij")),
    _testout1(declareProperty<Real>("testout_1")),
    _thetaout(declareProperty<Real>("theta_out")),
    _noij_out(declareProperty<Real>("noij_out")),
    // Moelans L Derivatives
    _dL_deta(_op_num),
    _d2L_deta2(_op_num),
    _dL_dgradeta_name(_op_num),
    _dL_dgradeta(_op_num),
    _d2L_dgradeta2_name(_op_num),
    _d2L_dgradeta2(_op_num),
    _d2L_dgradetadeta(_op_num)

{
  if ((_gb_mode == MISO || _gb_mode == FULL) && !misorientationEnabled())
    mooseError(name(),
               ": gb_mode='miso' or gb_mode='full' requires "
               "enable_misorientation=true.");
  // L derivative dependence for full model
  if (_aniso_L)
  {
    for (unsigned int i = 0; i < _op_num; ++i)
    {
      _dL_deta[i] = &declarePropertyDerivative<Real>("L", coupledName("v", i)); //_vals_name[i]);
      _d2L_deta2[i].resize(_op_num);
      // First-derivative wrt grad eta_i
      _dL_dgradeta_name[i] = "dL_dgradeta_" + Moose::stringify(i);
      _dL_dgradeta[i] = &declareProperty<RealGradient>(_dL_dgradeta_name[i]);
      // Second-derivative wrt (grad eta_i, grad eta_*)
      _d2L_dgradeta2_name[i] = "d2L_dgradeta2_" + Moose::stringify(i);
      _d2L_dgradeta2[i] = &declareProperty<std::vector<RealTensorValue>>(_d2L_dgradeta2_name[i]);
      // deta parts
      _d2L_dgradetadeta[i].resize(_op_num);
      for (unsigned int j = 0; j < _op_num; ++j)
      {
        _d2L_deta2[i][j] =
            &declarePropertyDerivative<Real>("L", coupledName("v", i), coupledName("v", j));
        _d2L_dgradetadeta[i][j] =
            &declarePropertyDerivative<RealGradient>(_dL_dgradeta_name[i], coupledName("v", j));
      }
    }
  }
}

void
GBCombinedAnisotropyMaterial::computeQpProperties()
{
  // Builds theta_ij, polar_ij, ij_i, ij_j, no_ij_pairs, and derivatives.
  GBInclinationBase::computeQpProperties();

  // Angles and dangle/dgradeta
  const auto & theta = _theta_ij[_qp];
  const auto & dtheta = _dtheta_dgradeta[_qp];
  const auto & d2theta = _d2theta_dgradeta2[_qp];
  const auto & polar = _polar_ij[_qp];
  const auto & dpolar = _dpolar_dgradeta[_qp];
  const auto & d2polar = _d2polar_dgradeta2[_qp];

  const auto & ij_i = _ij_i[_qp];
  const auto & ij_j = _ij_j[_qp];

  // New Parameters
  // f(angle) and derivs wrt gradeta
  auto & finc = _fgbe[_qp];
  std::vector<RealGradient> dfinc;
  std::vector<RealTensorValue> d2finc;
  finc.assign(theta.size(), 0.0);
  dfinc.assign(theta.size(), RealGradient(0.0));
  d2finc.assign(theta.size(), RealTensorValue(0.0));
  // gamma_ij and derivs wrt gradeta
  auto & gamma = _gamma_ij[_qp];
  auto & dgamma = _dgamma_dgradeta[_qp];
  auto & d2gamma = _d2gamma_dgradeta2[_qp];
  gamma.assign(theta.size(), 0.0);
  dgamma.assign(theta.size(), RealGradient(0.0));
  d2gamma.assign(theta.size(), RealTensorValue(0.0));

  // If using full derivative dependence on L, zero all derivatives
  if (_aniso_L)
  {
    for (unsigned int i = 0; i < _op_num; ++i)
    {
      // First-derivative wrt grad eta_i: zero it at this qp
      (*_dL_dgradeta[i])[_qp] = RealGradient(0.0);
      auto & d2L_i = (*_d2L_dgradeta2[i])[_qp];
      if (d2L_i.size() != _op_num)
        d2L_i.resize(_op_num);
      std::fill(d2L_i.begin(), d2L_i.end(), RealTensorValue(0.0));
      (*_dL_deta[i])[_qp] = 0.0;
      for (unsigned int j = 0; j < _op_num; ++j)
      {
        (*_d2L_dgradetadeta[i][j])[_qp] = RealGradient(0.0);
        (*_d2L_deta2[j][i])[_qp] = 0.0;
      }
    }
  }

  // Hard-coded coefficients (for poly gamma)
  constexpr Real a1 = -3.0944; // coefficient for g2^4
  constexpr Real a2 = -1.8169; // coefficient for g2^3
  constexpr Real a3 = 10.323;  // coefficient for g2^2
  constexpr Real a4 = -8.1819; // coefficient for g2
  constexpr Real a5 = 2.0033;  // constant term

  // TEMPORARY
  _elem_noij[_qp] = false;
  _noij_out[_qp] = 0;
  _testout1[_qp] = -2;
  _thetaout[_qp] = -2;

  // ISOTROPIC BULK VALUES
  if (_no_ij_pairs[_qp])
  {
    _noij_out[_qp] = 1;
    // Gamma
    Real g = _gbe_iso[_qp] * _bulk_mult / (std::sqrt(_kappa[_qp] * _mu[_qp]));
    Real g2 = g * g;
    Real pg = (((a1 * g2 + a2) * g2 + a3) * g2 + a4) * g2 + a5;
    _gamma_asymm[_qp] = 1 / pg;
    // Calculate IW
    Real f0_int =
        (((((0.0788 * pg - 0.4955) * pg + 1.2244) * pg - 1.5281) * pg + 1.0686) * pg - 0.5563) *
            pg +
        0.2907;
    _int_width[_qp] = (std::sqrt(_kappa[_qp] / _mu[_qp])) * (std::sqrt(1 / f0_int));
    // Calculate L
    _L[_qp] = _L0[_qp] * _bulk_mult;
    _gbe_norm[_qp] = _bulk_mult;
    return;
  }

  // Looping sums
  Real hgb_tot = 0.0;
  Real iw_sum = 0.0;
  Real gamma_sum = 0.0;
  Real Lij_sum = 0.0;
  Real gbe_sum = 0.0;

  for (std::size_t k = 0; k < theta.size(); ++k)
  {
    // Possibly unnecessary extra check- noij comes from this
    if (theta[k] == -1.0)
      continue;

    const std::size_t i = ij_i[k];
    const std::size_t j = ij_j[k];

    switch (_gb_mode)
    {
      case COS:
      {
        const auto out = computeCosineOnlyGBE(theta[k], polar[k], _if_a, 2, 0);
        finc[k] = out.f;
        dfinc[k] = out.df_dtheta * dtheta[k];
        d2finc[k] = out.d2f_dtheta2 * libMesh::outer_product(dtheta[k], dtheta[k]) +
                    out.df_dtheta * d2theta[k];
        break;
      }

      case INC:
      {
        const auto out = computeInclinationGBE(theta[k], polar[k]);
        finc[k] = out.f;
        dfinc[k] = out.df_dtheta * dtheta[k] + out.df_dpolar * dpolar[k];
        d2finc[k] = out.d2f_dtheta2 * libMesh::outer_product(dtheta[k], dtheta[k]) +
                    out.df_dtheta * d2theta[k] + out.df_dpolar * d2polar[k] +
                    out.d2f_dthetadpolar * (libMesh::outer_product(dtheta[k], dpolar[k]) +
                                            libMesh::outer_product(dpolar[k], dtheta[k]));
        // + d2fdpolar2 * libMesh::outer_product(dpolar[k],dpolar[k]) // = 0
        break;
      }

      case MISO:
      {
        const auto & miso = getMisorientationData(i, j);
        const auto out = computeMisorientationGBE(miso);
        finc[k] = out.f;
        dfinc[k] = RealGradient(0.0); // out.df_dtheta * dtheta[k];
        d2finc[k] = RealTensorValue(0.0);
        // out.d2f_dtheta2 * libMesh::outer_product(dtheta[k], dtheta[k]) +
        // out.df_dtheta * d2theta[k];
        break;
      }

      case FULL:
      {
        const auto & miso = getMisorientationData(i, j);
        const auto out = computeFullGBE(theta[k], polar[k], miso, _w_inc, _w_miso);
        finc[k] = out.f;
        dfinc[k] = out.df_dtheta * dtheta[k] + out.df_dpolar * dpolar[k];
        d2finc[k] = out.d2f_dtheta2 * libMesh::outer_product(dtheta[k], dtheta[k]) +
                    out.df_dtheta * d2theta[k] + out.df_dpolar * d2polar[k] +
                    out.d2f_dthetadpolar * (libMesh::outer_product(dtheta[k], dpolar[k]) +
                                            libMesh::outer_product(dpolar[k], dtheta[k]));
        // + d2fdpolar2 * libMesh::outer_product(dpolar[k],dpolar[k]) // = 0
        break;
      }

      case ISO:
      {
        finc[k] = _bulk_mult;
        dfinc[k] = RealGradient(0.0);
        d2finc[k] = RealTensorValue(0.0);
        break;
      }

      default:
        mooseError(name(), ": unknown gb_mode = ", _gb_mode);
    }

    // #########################################################
    // Condense this k
    // Calculate gamma
    Real g = _gbe_iso[_qp] * finc[k] / (std::sqrt(_kappa[_qp] * _mu[_qp]));
    Real dg_df = _gbe_iso[_qp] / (std::sqrt(_kappa[_qp] * _mu[_qp]));
    Real g2 = g * g;
    // Polynomial (inverse gamma for 0.5 < gamma < 40)
    Real pg = (((a1 * g2 + a2) * g2 + a3) * g2 + a4) * g2 + a5;
    Real dpg = 4 * a1 * g2 * g2 * g2 + 3 * a2 * g2 * g2 + 2 * a3 * g2 + a4;
    Real d2pg = 12 * a1 * g2 * g2 + 6 * a2 * g2 + 2 * a3;
    gamma[k] = 1 / pg;
    if (_stiffness)
    {
      dgamma[k] = -2 * g * gamma[k] * dpg * dg_df * dfinc[k];
      d2gamma[k] = 2 * gamma[k] * gamma[k] * dg_df * dg_df *
                       (4 * g2 * gamma[k] * dpg * dpg - 2 * g2 * d2pg - dpg) *
                       libMesh::outer_product(dfinc[k], dfinc[k]) -
                   2 * g * gamma[k] * dpg * dg_df * d2finc[k];
    }
    // Calculate IW
    Real f0_int =
        (((((0.0788 * pg - 0.4955) * pg + 1.2244) * pg - 1.5281) * pg + 1.0686) * pg - 0.5563) *
            pg +
        0.2907;
    Real iw_ij = (std::sqrt(_kappa[_qp] / _mu[_qp])) * (std::sqrt(1 / f0_int));

    // Summation
    Real gb_ij_weight = 1;
    if (_tj_mode == WEIGHTED) // Weighted summation
      gb_ij_weight = (*_vals[i])[_qp] * (*_vals[i])[_qp] * (*_vals[j])[_qp] * (*_vals[j])[_qp];
    hgb_tot += gb_ij_weight;
    iw_sum += iw_ij * gb_ij_weight;
    gamma_sum += gamma[k] * gb_ij_weight;
    gbe_sum += finc[k];

    // L
    // Lij not used in kernel, so dont need to store full vector
    Real Lij = 0.0;
    Real dLdf = 0.0;
    Real d2Ldf2 = 0.0;
    if (_aniso_mob) // if gbmob is aniso same as gbe then L has the aniso f applied twice
    {
      Lij = _L0[_qp] * finc[k] * finc[k];
      dLdf = _L0[_qp] * finc[k]; // should it be 2x? there was a reason but i dont remember why
      d2Ldf2 = _L0[_qp];         // should it be 2x? there was a reason but i dont remember why
    }
    else
    {
      Lij = _L0[_qp] * finc[k];
      dLdf = _L0[_qp];
      d2Ldf2 = 0.0;
    }
    Lij_sum += Lij * gb_ij_weight;
    if (_aniso_L && _stiffness)
    {
      (*_dL_dgradeta[i])[_qp] += dLdf * dfinc[k] * gb_ij_weight;
      (*_dL_dgradeta[j])[_qp] -= dLdf * dfinc[k] * gb_ij_weight;
      auto & d2L_i = (*_d2L_dgradeta2[i])[_qp];
      auto & d2L_j = (*_d2L_dgradeta2[j])[_qp];

      d2L_i[i] += (d2Ldf2 * libMesh::outer_product(dfinc[k], dfinc[k]) + dLdf * d2finc[k]) *
                  gb_ij_weight; // (i,i)
      d2L_j[j] += (d2Ldf2 * libMesh::outer_product(dfinc[k], dfinc[k]) + dLdf * d2finc[k]) *
                  gb_ij_weight; // (j,j)
      d2L_i[j] -= (d2Ldf2 * libMesh::outer_product(dfinc[k], dfinc[k]) + dLdf * d2finc[k]) *
                  gb_ij_weight; // (i,j)
      d2L_j[i] -= (d2Ldf2 * libMesh::outer_product(dfinc[k], dfinc[k]) + dLdf * d2finc[k]) *
                  gb_ij_weight; // (j,i)
    }

    if (k == 0) // DEBUG
    {
      _testout1[_qp] = finc[k];
      _thetaout[_qp] = theta[k];
    }
  }
  // #########################################################
  // Combine all k
  _gamma_asymm[_qp] = gamma_sum / hgb_tot;
  _int_width[_qp] = iw_sum / hgb_tot;
  _L[_qp] = Lij_sum / hgb_tot;
  _gbe_norm[_qp] = gbe_sum / hgb_tot;
  if (_aniso_L && _stiffness)
  {
    for (unsigned int i = 0; i < _op_num; ++i)
    {
      (*_dL_dgradeta[i])[_qp] /= hgb_tot;
      for (unsigned int j = 0; j < _op_num; ++j)
      {
        (*_d2L_dgradeta2[i])[_qp][j] /= hgb_tot;
      }
    }
  }
}

GBCombinedAnisotropyMaterial::AngleFunctionResult
GBCombinedAnisotropyMaterial::computeCosineOnlyGBE(
    const Real theta_inc, const Real polar_inc, Real a, Real b, Real c) const
{
  AngleFunctionResult out;

  out.f = 1.0 + a * std::cos(b * (theta_inc + c));
  out.df_dtheta = -a * b * std::sin(b * (theta_inc + c));
  out.d2f_dtheta2 = -a * b * b * std::cos(b * (theta_inc + c));
  out.df_dpolar = 0.0;
  out.d2f_dpolar2 = 0.0;
  out.d2f_dthetadpolar = 0.0;

  return out;
}

GBCombinedAnisotropyMaterial::AngleFunctionResult
GBCombinedAnisotropyMaterial::computeInclinationGBE(const Real theta_inc,
                                                    const Real polar_inc) const
{
  // Smooth version of Lins combined inc portion, [0.5-1]
  // 0.5 + 0.5 * (1 + cos(2*theta)) * B
  AngleFunctionResult out;
  const Real c2 = std::cos(2.0 * theta_inc);
  const Real s2 = std::sin(2.0 * theta_inc);
  const Real B = 0.5 + (0.2 - 0.5) * polar_inc / (libMesh::pi / 2);
  const Real dB_dpolar = (0.2 - 0.5) / (libMesh::pi / 2);
  // Using a 50 percent range (0.5-1)
  out.f = 0.5 + 0.5 * (1.0 + c2) * B;
  out.df_dtheta = -s2 * B;
  out.d2f_dtheta2 = -2.0 * c2 * B;
  out.df_dpolar = 0.5 * (1.0 + c2) * dB_dpolar;
  out.d2f_dpolar2 = 0.0;
  out.d2f_dthetadpolar = -s2 * dB_dpolar;

  return out;
}

GBCombinedAnisotropyMaterial::AngleFunctionResult
GBCombinedAnisotropyMaterial::computeMisorientationGBE(const MisorientationData & miso) const
{
  // Available:
  //   miso.theta    = misorientation angle
  //   miso.polar_ax = misorientation-axis polar angle
  //   miso.azim_ax  = misorientation-axis azimuth angle
  //   miso.q        = quaternion
  //   miso.qnorm    = vector-part norm before axis normalization
  AngleFunctionResult out;
  Real ang_max = 62 * libMesh::pi / 180;
  const Real miso_rat = std::min(miso.theta / ang_max, 1.0);
  Real miso_ang_en = (miso_rat > 0.0) ? (miso_rat) * (1 - std::log(miso_rat)) : 0.0;
  Real miso_ax_en = std::pow(std::abs(std::cos(miso.polar_ax)), 0.4) +
                    std::pow(std::abs(std::cos(miso.azim_ax / 2)), 0.4);
  // Set upper limit to 1
  if (miso_ang_en > 1.0)
    miso_ang_en = 1.0;
  if (miso_ax_en > 1.0)
    miso_ax_en = 1.0;
  // Using a 50 percent range (0.5-1)
  // out.f = 0.3 + 0.7 * miso_ax_en * miso_ang_en;
  out.f = 0.5 + 0.5 * miso_ax_en * miso_ang_en;
  out.df_dtheta = 0.0;
  out.d2f_dtheta2 = 0.0;
  out.df_dpolar = 0.0;
  out.d2f_dpolar2 = 0.0;
  out.d2f_dthetadpolar = 0.0;

  return out;
}

GBCombinedAnisotropyMaterial::AngleFunctionResult
GBCombinedAnisotropyMaterial::computeFullGBE(const Real theta_inc,
                                             const Real polar_inc,
                                             const MisorientationData & miso,
                                             const Real w_inc,
                                             const Real w_miso) const
{
  AngleFunctionResult out;
  // MISORIENTATION
  Real ang_max = 62 * libMesh::pi / 180;
  const Real miso_rat = std::min(miso.theta / ang_max, 1.0);
  Real miso_ang_en = (miso_rat > 0.0) ? (miso_rat) * (1 - std::log(miso_rat)) : 0.0;
  Real miso_ax_en = std::pow(std::abs(std::cos(miso.polar_ax)), 0.4) +
                    std::pow(std::abs(std::cos(miso.azim_ax / 2)), 0.4);
  // Set upper limit to 1
  if (miso_ang_en > 1.0)
    miso_ang_en = 1.0;
  if (miso_ax_en > 1.0)
    miso_ax_en = 1.0;
  const Real twist = 0.5 * (miso_ax_en * miso_ang_en * w_miso + (1 - w_miso));
  const Real tilt = 0.2 * (miso_ax_en * miso_ang_en * w_miso + (1 - w_miso));
  // COMBINED SECTION
  // 0.5 + A * B * w_inc + twist * (1 - w_inc)
  // A (smooth) = 0.5 * (1 + cos (2 theta))
  // B = twist + polar/(pi/2) * (tilt - twist)
  // A
  const Real A = 0.5 * (1.0 + std::cos(2.0 * theta_inc));
  const Real dA_dtheta = -std::sin(2.0 * theta_inc);
  const Real d2A_dtheta2 = -2.0 * std::cos(2.0 * theta_inc);
  // B
  const Real B = twist + (tilt - twist) * polar_inc / (libMesh::pi / 2);
  const Real dB_dpolar = (tilt - twist) / (libMesh::pi / 2);
  // Condese full output
  out.f = 0.5 + A * B * w_inc + twist * (1.0 - w_inc);

  out.df_dtheta = dA_dtheta * B * w_inc;
  out.d2f_dtheta2 = d2A_dtheta2 * B * w_inc;

  out.df_dpolar = A * dB_dpolar * w_inc;
  out.d2f_dpolar2 = 0.0;

  out.d2f_dthetadpolar = dA_dtheta * dB_dpolar * w_inc;

  return out;
}
