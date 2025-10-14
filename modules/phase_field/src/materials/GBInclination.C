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
  MooseEnum inc_func("cos=0", "cos");
  params.addParam<MooseEnum>("inc_func",
                             inc_func,
                             "Which inclination function to use.  "
                             "cos: 1 + a * cos(b(theta + c));.");
  params.addParam<Real>("ifunc_a", 0.05, "Inclination function constant a.");
  params.addParam<Real>("ifunc_b", 2, "Inclination function constant b.");
  params.addParam<Real>("ifunc_c", 0.0, "Inclination function constant c.");
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
    _mu(declareProperty<Real>("mu"))
// _op_num(coupledComponents("v")),
// _vals(coupledValues("v")),
// _vals_name(_op_num),
// // Angular distance to the x axis
// _theta_ij(declareProperty<std::vector<Real>>("theta_ij")),
// _dtheta_dgradeta(declareProperty<std::vector<RealGradient>>("dtheta_dgradeta")),
// _d2theta_dgradeta2(declareProperty<std::vector<RealTensorValue>>("d2theta_dgradeta2")),
// // _dtheta_dgradeta(declareProperty<RealGradient>("dincdgrad_eta")),
// // Grain Tracker/FFC for GB identification
// _gb_case(getParam<MooseEnum>("gb_id_method")),
// _grain_tracker(isParamValid("grain_tracker") ? &getUserObject<GrainTracker>("grain_tracker")
//                                              : nullptr),
// _ffc_tracker(isParamValid("ffc") ? &getUserObject<FeatureFloodCount>("ffc") : nullptr),
// _gt_tol(getParam<Real>("gt_tol")),
// _gtnum(declareProperty<Real>("gtnum")),
// _angular_func(getParam<MooseEnum>("angular_func")),
// _intol(getParam<Real>("intol")),
// _altol(getParam<Real>("altol")),
// _hgb(getMaterialProperty<Real>(getParam<MaterialPropertyName>("hgb")))

{
  // if (_op_num < 2)
  //   mooseError("Inclination properties requires op_num >= 2");

  // _vals.resize(_op_num);
  // _grad_vals.resize(_op_num);

  // for (unsigned int i = 0; i < _op_num; ++i)
  // {
  //   _vals[i] = &coupledValue("v", i);
  //   _vals_name[i] = coupledName("v", i);
  //   _grad_vals[i] = &coupledGradient("v", i);
  // }
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
  // mooseAssert(theta.size() == _k2ij.size(),
  //             "theta/_k2ij size mismatch"); // move to error checking and makr sure both = num
  //             pairs

  // Declare inclination flatpacked vector
  auto & finc = _inclination[_qp];
  if (finc.size() != num_pairs)
    finc.resize(num_pairs);
  std::fill(finc.begin(), finc.end(), 0.0);

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
  std::fill(gamma.begin(), gamma.end(), -1.0);
  std::fill(dgamma.begin(), dgamma.end(), RealGradient(0.0));
  std::fill(d2gamma.begin(), d2gamma.end(), RealTensorValue(0.0));

  // Hard-coded coefficients (for poly gamma)
  constexpr Real a1 = -3.0944; // coefficient for g2^4
  constexpr Real a2 = -1.8169; // coefficient for g2^3
  constexpr Real a3 = 10.323;  // coefficient for g2^2
  constexpr Real a4 = -8.1819; // coefficient for g2
  constexpr Real a5 = 2.0033;  // constant term

  _testout1[_qp] = -1; //(*_vals[i])[_qp];
  _testout2[_qp] = -1; //(*_vals[j])[_qp];
  _testoutgrad[_qp] = RealGradient(0.0);

  for (std::size_t k = 0; k < theta.size(); ++k) // theta.size()
  {
    if (theta[k] == -1.0)
      continue; // skip this k=ij
    // Maybe also add an OP value check for i and j to skip? should be handled by ffc/gt though

    // Actually compute for this k=ij
    const std::size_t i = _k2i[k];
    const std::size_t j = _k2j[k];
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

      default:
        mooseError("Unknown inc_func = ", _inc_func);
        break;
    }
    // Calculate gamma
    Real g = _gbe_iso[_qp] * finc[k] / (std::sqrt(_kappa * _const_m));
    Real dg_df = _gbe_iso[_qp] / (std::sqrt(_kappa * _const_m));
    Real g2 = g * g;
    // Polynomial (inverse gamma for 0.5 < gamma < 40)
    Real poly_g = (((a1 * g2 + a2) * g2 + a3) * g2 + a4) * g2 + a5;
    Real dpoly_g_dg2 = 4 * a1 * g2 * g2 * g2 + 3 * a2 * g2 * g2 + 2 * a3 * g2 + a4;
    Real d2poly_g_dg22 = 12 * a1 * g2 * g2 + 6 * a2 * g2 + 2 * a3;
    // Save gamma_ij
    gamma[k] = 1 / poly_g;
    // DERIVS and IW and L?
  }
  _testout1[_qp] = finc[0];
  _testoutgrad[_qp] = dfinc_dgradeta[0];
  _testouttens[_qp] = d2finc_dgradeta2[0];
}
