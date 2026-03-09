#include "GBInclination.h"

registerMooseObject("PhaseFieldApp", GBInclination);

InputParameters
GBInclination::validParams()
{
  InputParameters params = GBInclinationBase::validParams();
  params.addClassDescription(
      "Child material to determine inclination dependent properties for AGG.");
  // params.addRequiredCoupledVarWithAutoBuild(
  //     "v", "var_name_base", "op_num", "Array of coupled variables");
  MooseEnum inc_func("cos=0 well=1 man_well=2 smooth_well=3 smooth_half=4", "cos");
  params.addParam<MooseEnum>(
      "inc_func",
      inc_func,
      "Which inclination function to use.  "
      "cos: 1 + a * cos(b(theta + c)); "
      "well: b wells starting from angle c of (1 - a) within angle d of any well, else (1 + a).");
  params.addParam<Real>("ifunc_a", 0.05, "Inclination function constant a.");
  params.addParam<Real>("ifunc_b", 2, "Inclination function constant b.");
  params.addParam<Real>("ifunc_c", 0.0, "Inclination function constant c.");
  params.addParam<Real>("ifunc_d", 10.0, "Inclination function constant d.");
  // GB ENERGY
  params.addParam<MaterialPropertyName>(
      "gb_energy_iso_name", "sigma_iso", "Isotropic GB energy before inclination dependence.");
  params.addParam<MaterialPropertyName>(
      "gb_energy_aniso_name", "sigma", "Inclination dependent GB energy output.");
  // FREE ENERGY PROPS
  params.addParam<MaterialPropertyName>(
      "kappa", "kappa", "Gradient energy constant kappa material name.");
  params.addParam<Real>(
      "free_energy_m", 1, "Free energy function constant m (or mu in PF kernels).");
  // OPTIONAL Anisotropic L
  params.addParam<bool>("aniso_L", false, "Is L anisotropic, else L=L0.");
  params.addParam<bool>(
      "noDeriv_L",
      false,
      "Use anisotropic L but without the derivatives (for kernel with variable_L=false)");
  params.addParam<bool>(
      "aniso_gbmob", false, "Apply gbe anisotropy also to the gb mobility (inferred).");
  params.addParam<bool>("stiffness", true, "Include the stiffness (d2gamma/dgradeta2), esle = 0.");
  params.addParam<MaterialPropertyName>(
      "L0", "L0", "AC mobility prefactor/reference value material.");
  MooseEnum combine_form("weighted=0 avg=1", "weighted");
  params.addParam<MooseEnum>("combine_gb_form",
                             combine_form,
                             "Which GB combination approach, weighted average or just average.");
  // params.addParam<UserObjectName>("grain_tracker",
  //                                 "The GrainTracker UserObject to get values from.");
  // params.addParam<UserObjectName>("ffc", "The FFC UserObject to get values from.");
  // MooseEnum angular_func("atan_2D=0 atan_3D=1 acos=2 atan_half=3", "atan_2D");
  // params.addParam<MooseEnum>(
  //     "angular_func",
  //     angular_func,
  //     "Which angular distance function to use. "
  //     "atan_2D: oriented angle atan2(y,x) in [0,2pi); "
  //     "atan_3D: oriented azimuth around +x (atan2(z,y)) in [0,2pi) *UNTESTED*;"
  //     "acos: acos(x) in [0,pi]; "
  //     "atan_half: atan2(sqrt(y^2+z^2), x) in [0,pi].");
  // params.addParam<Real>("intol", 100, "hgbalpha tolerance");
  // params.addParam<Real>("altol", 100, "alpha tolerance");
  // params.addParam<Real>("gt_tol", 0.001, "alpha tolerance");
  // params.addParam<MaterialPropertyName>("hgb", "hgb", "Name of gb switching function.");
  return params;
}

GBInclination::GBInclination(const InputParameters & parameters)
  : GBInclinationBase(parameters),
    _testout1(declareProperty<Real>("testout1")),
    _testout2(declareProperty<Real>("testout2")),
    _testoutgrad(declareProperty<RealGradient>("testoutgrad")),
    _testouttens(declareProperty<RealTensorValue>("testouttens")),
    // Which inclination function to use
    _inc_func(getParam<MooseEnum>("inc_func")),
    _if_a(getParam<Real>("ifunc_a")),
    _if_b(getParam<Real>("ifunc_b")),
    _if_c(getParam<Real>("ifunc_c")),
    _if_d(getParam<Real>("ifunc_d")),
    // Inclination (1 + cos) output
    _inclination(declareProperty<std::vector<Real>>("inclination")),
    // gamma_ij
    _gamma_ij(declareProperty<std::vector<Real>>("gamma_ij")),
    _dgamma_dgradeta(declareProperty<std::vector<RealGradient>>("dgamma_dgradeta")),
    _d2gamma_dgradeta2(declareProperty<std::vector<RealTensorValue>>("d2gamma_dgradeta2")),
    // GB Energy
    _gbe_iso(getMaterialProperty<Real>("gb_energy_iso_name")),
    _gbe_aniso(declareProperty<Real>(getParam<MaterialPropertyName>("gb_energy_aniso_name"))),
    // Other Free Energy Props
    _kappa(getMaterialProperty<Real>(getParam<MaterialPropertyName>("kappa"))),
    _const_m(getParam<Real>("free_energy_m")),
    _mu(declareProperty<Real>("mu")),
    // Elemental rounding needed or not here
    _int_width(_limit_umag ? declareProperty<Real>("int_width")
                           : declareProperty<Real>("int_width_qp")),
    _gamma_qp(_limit_umag ? declareProperty<Real>("gamma_asymm")
                          : declareProperty<Real>("gamma_qp")),
    _elem_noij(_limit_umag ? &declareProperty<bool>("elem_no_ij") : nullptr),
    // Optional Anisotropic L
    _aniso_L(getParam<bool>("aniso_L")),
    _no_deriv_L(getParam<bool>("noDeriv_L")),
    _aniso_mob(getParam<bool>("aniso_gbmob")),
    _stiffness(getParam<bool>("stiffness")),
    _gb_combo(getParam<MooseEnum>("combine_gb_form")),
    _L0(getMaterialProperty<Real>("L0")),
    _L_ij(declareProperty<std::vector<Real>>("L_ij")),
    // _dL_dgradeta(declareProperty<std::vector<RealGradient>>("dL_dgradeta")),
    // _d2L_dgradeta2(declareProperty<std::vector<RealTensorValue>>("d2L_dgradeta2")),
    _L(declareProperty<Real>("L")),
    // _aniso_L ? declareProperty<Real>("L_qp") : declareProperty<Real>("L")),
    // Moelans L Derivatives
    _dL_deta(_op_num),
    _d2L_deta2(_op_num),
    _dL_dgradeta_name(_op_num),
    _dL_dgradeta(_op_num),
    _d2L_dgradeta2_name(_op_num),
    _d2L_dgradeta2(_op_num),
    _d2L_dgradetadeta(_op_num)

{
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
GBInclination::computeQpProperties()
{
  GBInclinationBase::computeQpProperties(); // force it to actually execute the parent

  // Flatpack vector matrices definitions
  const unsigned num_pairs = GBPairPacking::count_upper(_op_num);
  auto & theta = _theta_ij[_qp];
  auto & dtheta = _dtheta_dgradeta[_qp];
  auto & d2theta = _d2theta_dgradeta2[_qp];

  // Declare inclination flatpacked vector
  auto & finc = _inclination[_qp];
  if (finc.size() != num_pairs)
    finc.resize(num_pairs);
  std::fill(finc.begin(), finc.end(), 1.0);

  // Declare derivative pieces we need for inclination function
  std::vector<RealGradient> dfinc_dgradeta;
  std::vector<RealTensorValue> d2finc_dgradeta2;
  dfinc_dgradeta.clear();
  d2finc_dgradeta2.clear();
  if (dfinc_dgradeta.size() != num_pairs)
    dfinc_dgradeta.resize(num_pairs);
  if (d2finc_dgradeta2.size() != num_pairs)
    d2finc_dgradeta2.resize(num_pairs);
  std::fill(dfinc_dgradeta.begin(), dfinc_dgradeta.end(), RealGradient(0.0));
  std::fill(d2finc_dgradeta2.begin(), d2finc_dgradeta2.end(), RealTensorValue(0.0));

  // Declare/build gamma
  auto & gamma = _gamma_ij[_qp];
  auto & dgamma = _dgamma_dgradeta[_qp];
  auto & d2gamma = _d2gamma_dgradeta2[_qp];
  // Size once per QP (no-op if already correct); then initialize values
  if (gamma.size() != num_pairs)
    gamma.resize(num_pairs);
  if (dgamma.size() != num_pairs)
    dgamma.resize(num_pairs);
  if (d2gamma.size() != num_pairs)
    d2gamma.resize(num_pairs);
  std::fill(gamma.begin(), gamma.end(), 1.0);
  std::fill(dgamma.begin(), dgamma.end(), RealGradient(0.0));
  std::fill(d2gamma.begin(), d2gamma.end(), RealTensorValue(0.0));

  // Declare/build L_ij
  auto & Lij = _L_ij[_qp];
  // auto & dL = _dL_dgradeta[_qp];
  // auto & d2L = _d2L_dgradeta2[_qp];
  // Size once per QP (no-op if already correct); then initialize values
  if (Lij.size() != num_pairs)
    Lij.resize(num_pairs);
  // if (dL.size() != num_pairs)
  //   dL.resize(num_pairs);
  // if (d2L.size() != num_pairs)
  //   d2L.resize(num_pairs);
  std::fill(Lij.begin(), Lij.end(), 0.0);
  // std::fill(dL.begin(), dL.end(), RealGradient(0.0));
  // std::fill(d2L.begin(), d2L.end(), RealTensorValue(0.0));

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

  // IW initialization
  std::vector<Real> iw_ij;
  iw_ij.clear();
  if (iw_ij.size() != num_pairs)
    iw_ij.resize(num_pairs);
  std::fill(iw_ij.begin(), iw_ij.end(), 0.0);

  // hgb storage
  Real hgb_tot = 0.0;
  Real iw_sum = 0.0;
  Real gamma_sum = 0.0;
  Real Lij_sum = 0.0;

  _testout1[_qp] = -1; //(*_vals[i])[_qp];
  _testout2[_qp] = -1; //(*_vals[j])[_qp];
  _testoutgrad[_qp] = RealGradient(0.0);

  if (!_no_ij_pairs[_qp])
  {
    for (std::size_t k = 0; k < theta.size(); ++k) // theta.size()
    {
      if (theta[k] == -1.0)
        continue; // skip this k=ij
      // Maybe also add an OP value check for i and j to skip? should be handled by ffc/gt though

      // Actually compute for this k=ij
      const std::size_t i = _ij_i[_qp][k];
      const std::size_t j = _ij_j[_qp][k];

      // choose which inclination function
      switch (_inc_func)
      {
        case 0:
        {
          // cos: 1 + a * cos(b(theta + c))
          const Real fcos = 1 + _if_a * std::cos(_if_b * (theta[k] + _if_c));
          const Real dfcos = -1 * _if_a * _if_b * std::sin(_if_b * (theta[k] + _if_c));
          const Real d2fcos = -1 * _if_a * _if_b * _if_b * std::cos(_if_b * (theta[k] + _if_c));
          if (fcos >= 0.0)
          {
            finc[k] = fcos;
            dfinc_dgradeta[k] = dfcos * dtheta[k];
            d2finc_dgradeta2[k] =
                d2fcos * libMesh::outer_product(dtheta[k], dtheta[k]) + dfcos * d2theta[k];
          }
          break;
        }

        case 1:
        {
          // well with b minimums of width +/-d with the first offset by angle c
          const unsigned int nb = static_cast<unsigned int>(_if_b);
          if (nb == 0)
            break;
          bool in_well = false;
          const Real twopi = 2.0 * libMesh::pi;
          const Real rad_d = _if_d * twopi / 360;
          const Real spacing = twopi / _if_b;
          for (unsigned int m = 0; m < nb && !in_well; ++m)
          {
            const Real center = (_if_c * twopi / 360) + m * spacing;
            double d = fmod(theta[k] - center, twopi);
            if (d < 0.0)
              d += twopi;
            d = std::min(d, twopi - d);
            if (d < rad_d)
              in_well = true;
          }
          const Real fwell = in_well ? (1.0 - _if_a) : (1.0 + _if_a);
          finc[k] = fwell;
          break;
        }

        case 2:
        {
          // hard coded well at 2 pole 90/270 with 10 deg
          bool in_well = false;
          const Real twopi = 2.0 * libMesh::pi;
          const Real rad_d = 10 * twopi / 360;
          const Real rad_w1 = 90 * twopi / 360;
          const Real rad_w2 = 270 * twopi / 360;
          const bool well_1 = (std::abs(theta[k] - rad_w1) < rad_d);
          const bool well_2 = (std::abs(theta[k] - rad_w2) < rad_d);
          in_well = (well_1 || well_2);
          const Real fwell = in_well ? (1.0 - _if_a) : (1.0 + _if_a);
          finc[k] = fwell;
          break;
        }

        case 3:
        {
          // smoothed out well at 2 pole 90/270 with 10 deg
          const Real twopi = 2.0 * libMesh::pi;
          const Real td = theta[k] * 360 / twopi;
          Real w = 0.0;
          if (td > 78 && td < 82)
          {
            Real s = (td - 78) / 4;
            w = s * s * s * (6.0 * s * s - 15.0 * s + 10.0);
          }
          else if (td > 258 && td < 262)
          {
            Real s = (td - 258) / 4;
            w = s * s * s * (6.0 * s * s - 15.0 * s + 10.0);
          }
          else if (td > 98 && td < 102)
          {
            Real s = (td - 98) / 4;
            w = 1 - (s * s * s * (6.0 * s * s - 15.0 * s + 10.0));
          }
          else if (td > 278 && td < 282)
          {
            Real s = (td - 278) / 4;
            w = 1 - (s * s * s * (6.0 * s * s - 15.0 * s + 10.0));
          }
          else if (td >= 82 && td <= 98)
          {
            w = 1.0;
          }
          else if (td >= 262 && td <= 278)
          {
            w = 1.0;
          }
          const Real fwell = 1.0 + _if_a - (2 * w * _if_a);
          finc[k] = fwell;
          break;
        }

        case 4:
        {
          // smoothed out well at 2 pole 90/270 with 10 deg
          const Real twopi = 2.0 * libMesh::pi;
          const Real td = theta[k] * 360 / twopi;
          Real w = 0.0;
          if (td > 79 && td < 81)
          {
            Real s = (td - 79) / 2;
            w = s * s * s * (6.0 * s * s - 15.0 * s + 10.0);
          }
          else if (td > 259 && td < 261)
          {
            Real s = (td - 259) / 2;
            w = s * s * s * (6.0 * s * s - 15.0 * s + 10.0);
          }
          else if (td > 99 && td < 101)
          {
            Real s = (td - 99) / 2;
            w = 1 - (s * s * s * (6.0 * s * s - 15.0 * s + 10.0));
          }
          else if (td > 279 && td < 281)
          {
            Real s = (td - 279) / 2;
            w = 1 - (s * s * s * (6.0 * s * s - 15.0 * s + 10.0));
          }
          else if (td >= 81 && td <= 99)
          {
            w = 1.0;
          }
          else if (td >= 261 && td <= 279)
          {
            w = 1.0;
          }
          const Real fwell = 1.0 + _if_a - (2 * w * _if_a);
          finc[k] = fwell;
          break;
        }

        default:
          mooseError("Unknown inc_func = ", _inc_func);
          break;
      }
      // Calculate gamma
      Real g = _gbe_iso[_qp] * finc[k] / (std::sqrt(_kappa[_qp] * _const_m));
      Real dg_df = _gbe_iso[_qp] / (std::sqrt(_kappa[_qp] * _const_m));
      Real g2 = g * g;
      // Polynomial (inverse gamma for 0.5 < gamma < 40)
      Real pg = (((a1 * g2 + a2) * g2 + a3) * g2 + a4) * g2 + a5;
      Real dpg = 4 * a1 * g2 * g2 * g2 + 3 * a2 * g2 * g2 + 2 * a3 * g2 + a4;
      Real d2pg = 12 * a1 * g2 * g2 + 6 * a2 * g2 + 2 * a3;
      // Save gamma_ij
      gamma[k] = 1 / pg;
      if (_stiffness)
      {
        dgamma[k] = -2 * g * gamma[k] * dpg * dg_df * dfinc_dgradeta[k];
        d2gamma[k] = 2 * gamma[k] * gamma[k] * dg_df * dg_df *
                         (4 * g2 * gamma[k] * dpg * dpg - 2 * g2 * d2pg - dpg) *
                         libMesh::outer_product(dfinc_dgradeta[k], dfinc_dgradeta[k]) -
                     2 * g * gamma[k] * dpg * dg_df * d2finc_dgradeta2[k];
      }

      // Calculate IW
      Real f0_int =
          (((((0.0788 * pg - 0.4955) * pg + 1.2244) * pg - 1.5281) * pg + 1.0686) * pg - 0.5563) *
              pg +
          0.2907;
      iw_ij[k] = (std::sqrt(_kappa[_qp] / _const_m)) * (std::sqrt(1 / f0_int));

      Real gb_ij_weight = 0.0;
      if (_gb_combo == 0) // Weighted summation
        gb_ij_weight = (*_vals[i])[_qp] * (*_vals[i])[_qp] * (*_vals[j])[_qp] * (*_vals[j])[_qp];
      else
        gb_ij_weight = 1;

      // sum based on weight/avg multiplier
      hgb_tot += gb_ij_weight;
      iw_sum += iw_ij[k] * gb_ij_weight;
      gamma_sum += gamma[k] * gb_ij_weight;

      // Optionally anisotropic L
      if (_aniso_L || _no_deriv_L)
      {
        // Real gb2 = (*_vals[i])[_qp] * (*_vals[i])[_qp] * (*_vals[j])[_qp] * (*_vals[j])[_qp];
        // could do the Lij and derivativs in the kernel instead?
        Real dLdf = 0.0;
        Real d2Ldf2 = 0.0;
        if (_aniso_mob) // if gbmob is aniso same as gbe then L has the aniso f applied twice
        {
          Lij[k] = _L0[_qp] * finc[k] * finc[k];
          dLdf = _L0[_qp] * finc[k];
          d2Ldf2 = _L0[_qp];
        }
        else
        {
          Lij[k] = _L0[_qp] * finc[k];
          dLdf = _L0[_qp];
          d2Ldf2 = 0.0;
        }
        Lij_sum += Lij[k] * gb_ij_weight;
        if (_aniso_L && _stiffness)
        {
          // dL[k] = _L0[_qp] * dfinc_dgradeta[k];
          // d2L[k] = _L0[_qp] * d2finc_dgradeta2[k];
          (*_dL_dgradeta[i])[_qp] += dLdf * dfinc_dgradeta[k] * gb_ij_weight;
          (*_dL_dgradeta[j])[_qp] -= dLdf * dfinc_dgradeta[k] * gb_ij_weight;
          auto & d2L_i = (*_d2L_dgradeta2[i])[_qp];
          auto & d2L_j = (*_d2L_dgradeta2[j])[_qp];

          d2L_i[i] += (d2Ldf2 * libMesh::outer_product(dfinc_dgradeta[k], dfinc_dgradeta[k]) +
                       dLdf * d2finc_dgradeta2[k]) *
                      gb_ij_weight; // (i,i)
          d2L_j[j] += (d2Ldf2 * libMesh::outer_product(dfinc_dgradeta[k], dfinc_dgradeta[k]) +
                       dLdf * d2finc_dgradeta2[k]) *
                      gb_ij_weight; // (j,j)
          d2L_i[j] -= (d2Ldf2 * libMesh::outer_product(dfinc_dgradeta[k], dfinc_dgradeta[k]) +
                       dLdf * d2finc_dgradeta2[k]) *
                      gb_ij_weight; // (i,j)
          d2L_j[i] -= (d2Ldf2 * libMesh::outer_product(dfinc_dgradeta[k], dfinc_dgradeta[k]) +
                       dLdf * d2finc_dgradeta2[k]) *
                      gb_ij_weight; // (j,i)
        }
      }
    }
  }
  // Summation and whole outputs
  // Check if no ij pairs at qp use finc = 1 for calculation of condensed output
  if ((_no_ij_pairs[_qp]) || (hgb_tot < 1e-6))
  {
    Real g = _gbe_iso[_qp] / (std::sqrt(_kappa[_qp] * _const_m));
    Real g2 = g * g;
    Real pg = (((a1 * g2 + a2) * g2 + a3) * g2 + a4) * g2 + a5;
    _gamma_qp[_qp] = 1 / pg; // 1.5;
    // IW
    Real f0_int =
        (((((0.0788 * pg - 0.4955) * pg + 1.2244) * pg - 1.5281) * pg + 1.0686) * pg - 0.5563) *
            pg +
        0.2907;
    _int_width[_qp] = (std::sqrt(_kappa[_qp] / _const_m)) * (std::sqrt(1 / f0_int));
  }
  else
  {
    _int_width[_qp] = iw_sum / hgb_tot;
    _gamma_qp[_qp] = gamma_sum / hgb_tot;
  }
  _mu[_qp] = _const_m;

  // AC Mobility
  if (!_aniso_L)
  {
    if (_no_deriv_L)
    {
      if ((_no_ij_pairs[_qp]) || (hgb_tot < 1e-6))
        _L[_qp] = _L0[_qp];
      else
        _L[_qp] = Lij_sum / hgb_tot;
    }
    else
      _L[_qp] = _L0[_qp];
  }
  else if (_no_ij_pairs[_qp] || (hgb_tot < 1e-6))
  {
    // aniso L but skipping this qp
    _L[_qp] = _L0[_qp];
  }
  else
  {
    // eta summation of moelans L_ij
    _L[_qp] = Lij_sum / hgb_tot;
    if (_stiffness)
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

  if (_limit_umag)
    (*_elem_noij)[_qp] = false;

  if (_no_ij_pairs[_qp])
    _testout1[_qp] = 0;
  else
    _testout1[_qp] = 1;
  _testout2[_qp] = hgb_tot;
  // if ((_ij_i[_qp]).size() > 0)
  // {
  //   _testout1[_qp] = _ij_i[_qp][0];
  //   _testout2[_qp] = _ij_j[_qp][0];
  // }

  _testoutgrad[_qp] = dgamma[1];
  _testouttens[_qp] = d2gamma[1];
}
