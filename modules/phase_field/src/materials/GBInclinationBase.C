//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GBInclinationBase.h"

registerMooseObject("PhaseFieldApp", GBInclinationBase);

InputParameters
GBInclinationBase::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Inclination dependent properties for AGG.");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  MooseEnum gb_id_method("graintracker=0 ffc=1 all=2", "graintracker");
  params.addParam<MooseEnum>(
      "gb_id_method", gb_id_method, "Which GB OP identification/selection option.");
  params.addParam<UserObjectName>("grain_tracker",
                                  "The GrainTracker UserObject to get values from.");
  params.addParam<UserObjectName>("ffc", "The FFC UserObject to get values from.");
  // params.addParam<UserObjectName>("ebsd_reader", "The EBSDReader GeneralUserObject");
  // params.addParam<Real>("delta_ij", 0.05, "Anisotropy weight in cos function");
  // params.addParam<Real>(
  //     "theta_prefactor", 4.0, "Multiplier in cos function (cos(n(theta + theta0)))");
  // params.addParam<Real>("inc_ij_0", 0.0, "Inclination function offset in cos function");
  // params.addParam<std::vector<MaterialPropertyName>>(
  //     "dgamma_grad_eta_names",
  //     std::vector<MaterialPropertyName>(),
  //     "Interfacial / grain boundary gamma parameter names (leave empty for gamma0... gammaN)");
  // params.addParam<std::vector<MaterialPropertyName>>(
  //     "d2gamma_grad_eta2_names",
  //     std::vector<MaterialPropertyName>(),
  //     "Interfacial / grain boundary gamma parameter names (leave empty for gamma0... gammaN)");
  // params.addParam<MaterialPropertyName>("gb_energy_input",
  //                                       "GB energy before inclination dependence");
  // params.addParam<MaterialPropertyName>(
  //     "gb_energy", "gb_energy", "Inclination dependent GB energy output.");
  // params.addParam<Real>("kappa", 1, "Gradient energy constant kappa value");
  // params.addParam<Real>("free_energy_m", 1, "Free energy function constant m");
  // params.addParam<bool>(
  //     "L_of_eta", false, "Is L a function of eta, requiring those derivatives to be defined.");
  // params.addParam<bool>("aniso_L", false, "Is L anisotropic, else L=L0.");
  // params.addParam<MaterialPropertyName>("L0", "AC mobility prefactor/reference value material.");
  // params.addParam<MaterialPropertyName>("gamma0",
  //                                       "gamma prefactor/reference value material (for mu
  //                                       calc).");
  // params.addParam<MaterialPropertyName>("L_name", "L_aniso", "Name of anisotropic L output.");
  // params.addParam<MaterialPropertyName>(
  //     "gamma_name", "gamma_aniso", "Name of anisotropic gamma output.");
  // params.addParam<MaterialPropertyName>("mu_name", "mu_aniso", "Name of anisotropic mu output.");
  // params.addParam<bool>(
  //     "continuous", false, "Disregard GT and calculate for all variables everywhere.");
  MooseEnum angular_func("atan_2D=0 atan_3D=1 acos=2 atan_half=3", "atan_2D");
  params.addParam<MooseEnum>(
      "angular_func",
      angular_func,
      "Which angular distance function to use. "
      "atan_2D: oriented angle atan2(y,x) in [0,2pi); "
      "atan_3D: oriented azimuth around +x (atan2(z,y)) in [0,2pi) *UNTESTED*;"
      "acos: acos(x) in [0,pi]; "
      "atan_half: atan2(sqrt(y^2+z^2), x) in [0,pi].");
  // MooseEnum alphacase("base=0 ngb=1 thresh=2 exclude=3 zero=4 subzero=5 hgb=6 both=7", "base");
  // params.addParam<MooseEnum>("alphacase", alphacase, "Which alpha option to use.");
  params.addParam<Real>("intol", 1e-6, "hgbalpha tolerance");
  params.addParam<Real>("altol", 1e-6, "alpha tolerance");
  params.addParam<Real>("gt_tol", 0.001, "alpha tolerance");
  params.addParam<MaterialPropertyName>("hgb", "hgb", "Name of gb switching function.");
  // params.addParam<bool>("moelans_mu", true, "Is mu defined per moelans aniso, else = sigma.");

  return params;
}

GBInclinationBase::GBInclinationBase(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    // : Material(parameters),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _vals_name(_op_num),
    // Inclination angular
    _inclination(declareProperty<Real>("inclination_function")),
    // Angular distance to the x axis
    _theta_ij(declareProperty<std::vector<Real>>("theta_ij")),
    _dtheta_dgradeta(declareProperty<std::vector<RealGradient>>("dtheta_dgradeta")),
    _d2theta_dgradeta2(declareProperty<std::vector<RealTensorValue>>("d2theta_dgradeta2")),
    // _dtheta_dgradeta(declareProperty<RealGradient>("dincdgrad_eta")),
    // Grain Tracker/FFC for GB identification
    _gb_case(getParam<MooseEnum>("gb_id_method")),
    _grain_tracker(isParamValid("grain_tracker") ? &getUserObject<GrainTracker>("grain_tracker")
                                                 : nullptr),
    _ffc_tracker(isParamValid("ffc") ? &getUserObject<FeatureFloodCount>("ffc") : nullptr),
    _gt_tol(getParam<Real>("gt_tol")),
    // _ebsd_reader(getUserObject<EBSDReader>("ebsd_reader")),
    // _delta_ij(getParam<Real>("delta_ij")),
    // _theta_pre(getParam<Real>("theta_prefactor")),
    // _inc_ij_0(getParam<Real>("inc_ij_0")),
    // _gamma(declareProperty<Real>(getParam<MaterialPropertyName>("gamma_name"))),
    // _dgammadgrad_eta_name(getParam<std::vector<MaterialPropertyName>>("dgamma_grad_eta_names")),
    // _dgammadgrad_eta(_op_num),
    // _d2gammadgrad_eta2_name(getParam<std::vector<MaterialPropertyName>>("d2gamma_grad_eta2_names")),
    // _d2gammadgrad_eta2(_op_num),
    // _moelans_mu(getParam<bool>("moelans_mu")),
    // _mu0(getParam<Real>("mu0")),
    // Inclincation dependent L
    // _aniso_L(getParam<bool>("aniso_L")),
    // _L_of_eta(getParam<bool>("L_of_eta")),
    // _L_name(getParam<MaterialPropertyName>("L_name")),
    // _L(declareProperty<Real>(_L_name)),
    // _dLdeta(_op_num),
    // _d2Ldetadeta(_op_num),
    // _dLdgrad_eta_name(_op_num),
    // _dLdgrad_eta(_op_num),
    // _d2Ldgrad_eta2_name(_op_num),
    // _d2Ldgrad_eta2(_op_num),
    // _d2Ldgrad_etadeta(_op_num),
    // TEMP TEST OUTPUTS
    // _opout(declareProperty<Real>("opout")),
    // _opout2(declareProperty<Real>("opout2")),
    // _testout(declareProperty<Real>("testout")),
    // _testout2(declareProperty<Real>("testout2")),
    // _testout3(declareProperty<Real>("testout3")),
    // _alpha_out(declareProperty<Real>("alpha_out")),
    _gtnum(declareProperty<Real>("gtnum")),
    // _num_pairs(GBPairPacking::count_upper_incl_diag(_op_num)),
    // _altnum(declareProperty<Real>("altnum")),
    // _atens(declareProperty<RealTensorValue>("atens")),
    // _t2tens(declareProperty<RealTensorValue>("t2tens")),
    // _ngbtens(declareProperty<RealTensorValue>("ngbtens")),
    // _testoutgrad(declareProperty<RealGradient>("testoutgrad")),
    // _testoutgrad2(declareProperty<RealTensorValue>("testoutgrad2")),
    // _inclin(declareProperty<RealGradient>("inclin_vec")),
    // _dadb(declareProperty<RealGradient>("dadb")),
    // _d2adb2(declareProperty<RealTensorValue>("d2adb2")),
    // _d3adb3(declareProperty<RankTwoTensor>("d3adb3")),
    // Other material properties
    // _gbe(getMaterialProperty<Real>("gb_energy_input")),
    // _gbe_inc(declareProperty<Real>(getParam<MaterialPropertyName>("gb_energy"))),
    // _kappa(getParam<Real>("kappa")),
    // _const_m(getParam<Real>("free_energy_m")),
    // _L0(getMaterialProperty<Real>("L0")),
    // _gamma0(getMaterialProperty<Real>("gamma0")),
    // _int_width(declareProperty<Real>("int_width")),
    // _mu(declareProperty<Real>(getParam<MaterialPropertyName>("mu_name"))),
    // _gb_test_grad(declareProperty<RealGradient>("gb_test_grad")),
    // _alt_vec(declareProperty<RealGradient>("alt_vec")),
    // _continuous(getParam<bool>("continuous")),
    _angular_func(getParam<MooseEnum>("angular_func")),
    // _alpha_case(getParam<MooseEnum>("alphacase")),
    _intol(getParam<Real>("intol")),
    _altol(getParam<Real>("altol")),
    _hgb(getMaterialProperty<Real>(getParam<MaterialPropertyName>("hgb")))

{
  if (_op_num == 0)
    mooseError("Model requires op_num > 0");

  // if (_L_of_eta)
  //   mooseError("Inclination L currently only coded for grad u dependent, not u dependent.");

  // if (_dgammadgrad_eta_name.size() != 0 && _dgammadgrad_eta_name.size() != _op_num)
  //   paramError("dgamma_grad_eta_names",
  //              "Specify either as many entries as op_num values or none at all for auto-naming
  //              the " "gamma gradients with respect to gradient of grain OPs.");
  // if (_d2gammadgrad_eta2_name.size() != 0 && _d2gammadgrad_eta2_name.size() != _op_num)
  //   paramError("d2gamma_grad_eta2_names",
  //              "Specify either as many entries as op_num values or none at all for auto-naming
  //              the " "gamma gradients with respect to gradient of grain OPs.");

  // ADD SOME ERROR FOR CASE AND GT VS FFC INPUTS SINCE THEYRE OPTIONAL

  // // automatic names for the gamma properties
  // if (_dgammadgrad_eta_name.size() == 0)
  // {
  //   _dgammadgrad_eta_name.resize(_op_num);
  //   for (unsigned int i = 0; i < _op_num; i++)
  //     _dgammadgrad_eta_name[i] = "dgammadgrad_eta_" + Moose::stringify(i);
  // }
  // if (_d2gammadgrad_eta2_name.size() == 0)
  // {
  //   _d2gammadgrad_eta2_name.resize(_op_num);
  //   for (unsigned int i = 0; i < _op_num; i++)
  //     _d2gammadgrad_eta2_name[i] = "d2gammadgrad_eta2_" + Moose::stringify(i);
  // }

  _vals.resize(_op_num);
  _grad_vals.resize(_op_num);
  // _incl_tens.resize(_op_num);
  // _ang_dist.resize(_op_num);

  // _num_pairs = GBPairPacking::count_upper_incl_diag(_op_num);
  // _theta_ij.resize(_num_pairs);

  for (unsigned int i = 0; i < _op_num; ++i)
  {
    _vals[i] = &coupledValue("v", i);
    _vals_name[i] = coupledName("v", i);
    _grad_vals[i] = &coupledGradient("v", i);
    // // Build the ij tensor/matrix of inclinations
    // _incl_tens[i].resize(_op_num);
    // _ang_dist[i].resize(_op_num);
    // _dgammadgrad_eta[i] = &declareProperty<RealGradient>(_dgammadgrad_eta_name[i]);
    // _d2gammadgrad_eta2[i] = &declareProperty<RealTensorValue>(_d2gammadgrad_eta2_name[i]);
    // // L derivatives
    // _dLdeta[i] = &declarePropertyDerivative<Real>(_L_name, _vals_name[i]);
    // _d2Ldetadeta[i].resize(_op_num);
    // _dLdgrad_eta_name[i] = "dLdgrad_eta_" + Moose::stringify(i);
    // _dLdgrad_eta[i] = &declareProperty<RealGradient>(_dLdgrad_eta_name[i]);
    // _d2Ldgrad_eta2_name[i] = "d2Ldgrad_eta2_" + Moose::stringify(i);
    // _d2Ldgrad_eta2[i] = &declareProperty<RealTensorValue>(_d2Ldgrad_eta2_name[i]);
    // _d2Ldgrad_etadeta[i].resize(_op_num);
    // // _d2Ldetadgrad_eta[i] =
    // //     &declarePropertyDerivative<RealGradient>(_dLdgrad_eta_name[i], _vals_name[i]);
    // // _d2Ldetadgrad_eta[i].resize(_op_num); // only using same eta for eta and grad eta here
    // for (unsigned int j = 0; j <= i; ++j)
    // {
    //   _d2Ldetadeta[j][i] = &declarePropertyDerivative<Real>(_L_name, _vals_name[j],
    //   _vals_name[i]); _d2Ldgrad_etadeta[i][j] =
    //       &declarePropertyDerivative<RealGradient>(_dLdgrad_eta_name[i], _vals_name[j]);
    // }
  }
}

void
GBInclinationBase::computeQpProperties()
{
  const unsigned num_pairs = GBPairPacking::count_upper(_op_num);
  auto & theta = _theta_ij[_qp];
  auto & dtheta = _dtheta_dgradeta[_qp];
  auto & d2theta = _d2theta_dgradeta2[_qp];
  // Size once per QP (no-op if already correct); then initialize values
  if (theta.size() != num_pairs)
    theta.resize(num_pairs);
  if (dtheta.size() != num_pairs)
    dtheta.resize(num_pairs);
  if (d2theta.size() != num_pairs)
    d2theta.resize(num_pairs);
  std::fill(theta.begin(), theta.end(), -1.0);
  std::fill(dtheta.begin(), dtheta.end(), RealGradient(0.0));
  std::fill(d2theta.begin(), d2theta.end(), RealTensorValue(0.0));
  // // Fill
  // for (unsigned j = 0; j < _op_num; ++j)
  //   for (unsigned i = 0; i <= j; ++i)
  //   {
  //     const unsigned k = GBPairPacking::pack_upper_incl_diag(i, j);
  //     theta[k] = k;
  //     // _theta_ij[_qp][k] = k;
  //   }

  // // Zero out the derivatives for gamma wrt grad_eta
  // for (unsigned int i = 0; i < _op_num; ++i)
  //   (*_dgammadgrad_eta[i])[_qp] = RealGradient(0.0, 0.0, 0.0);

  // Find out the number of boundary unique_id and save them
  _gb_ij_list.clear();
  switch (_gb_case)
  {
    case 0:
    {
      // Grain Tracker Case
      if (!_grain_tracker)
        mooseError(name(), ": gb_id_method=graintracker requires 'grain_tracker' user object.");
      const auto & op_to_grains = (*_grain_tracker).getVarToFeatureVector(_current_elem->id());
      for (auto i : index_range(op_to_grains))
      {
        if (op_to_grains[i] == FeatureFloodCount::invalid_id)
          continue;

        _gb_ij_list.push_back(i);
      }
      break;
    }

    case 1:
    {
      // FFC Case
      if (!_ffc_tracker)
        mooseError(name(), ": gb_id_method=ffc requires 'ffc' user object.");
      const auto & op_to_grains = (*_ffc_tracker).getVarToFeatureVector(_current_elem->id());
      for (auto i : index_range(op_to_grains))
      {
        if (op_to_grains[i] == FeatureFloodCount::invalid_id)
          continue;

        _gb_ij_list.push_back(i);
      }
      break;
    }

    case 2:
    {
      // Variable threshold cutoffs- BROKEN/FAILS
      for (unsigned int i = 0; i < _op_num; ++i)
      {
        if ((*_vals[i])[_qp] >= _gt_tol)
        {
          _gb_ij_list.push_back(i);
        }
      }
      break;
    }

    default:
      mooseError("Unknown gb_case = ", _gb_case);
      break;
  }

  _gtnum[_qp] = theta.size(); // _gb_ij_list.size();
  std::sort(_gb_ij_list.begin(), _gb_ij_list.end());
  _hgb_pairs.clear();
  _inc_pairs.clear();

  // Inclination theta (angle) derivatives wrt gradient of eps
  // std::vector<RealGradient> dtheta_dgeta_list;
  // std::vector<RealTensorValue> d2theta_dgeta2_list;
  // dtheta_dgeta_list.clear();
  // d2theta_dgeta2_list.clear();
  std::vector<std::array<unsigned int, 2>> active_ij_pairs;
  active_ij_pairs.clear();

  switch (_gb_ij_list.size())
  {
    case 0:
      // _inclination_distance[_qp] = -1.0; // angular distance out of 0-2pi range (not counted)
      _inclination[_qp] = 1.0; // f = 1 + cos() = 1
      break;
    case 1:
      // _inclination_distance[_qp] = -1.0; // angular distance out of 0-2pi range (not counted)
      _inclination[_qp] = 1.0; // f = 1 + cos() = 1
      break;

    default:
      // do all ij pairs (i>j) if more than 2 vars/features
      for (std::size_t idx1 = 0; idx1 < _gb_ij_list.size(); ++idx1)
        for (std::size_t idx2 = idx1 + 1; idx2 < _gb_ij_list.size(); ++idx2)
        {
          unsigned int i = _gb_ij_list[idx1];
          unsigned int j = _gb_ij_list[idx2];
          const unsigned k = GBPairPacking::pack_upper(i, j);
          // if j-i points outward from lower grain number
          // if i-j like in paper points inward on lower grain number
          // Can also invert it by doing ngb = -ngb
          RealGradient ngb = (*_grad_vals[i])[_qp] - (*_grad_vals[j])[_qp];
          Real R = 0.0;
          Real a_dist = 0.0;
          RealGradient dtheta_dgetai(0.0, 0.0, 0.0);
          RealTensorValue d2theta_dgetai2(0.0);

          if (ngb.norm() > 1.0e-10)
          {
            RealGradient uxyz(0.0, 0.0, 0.0);
            Real alpha = 0.0;
            bool alpha_skip = false;
            Real hovtol = _hgb[_qp] / _altol;
            Real invtol = 1 / _intol;
            // BOTH: hgb and alpha- use hgb*alpha as well as straight alpha cutoff
            uxyz = ngb;
            ngb /= ngb.norm(); // normalize ngb now that we have uxyz for the un-normalized
            if ((uxyz.norm() < hovtol) || (uxyz.norm() < invtol)) // hgb (altol)
            {
              alpha = 0.0;
              alpha_skip = true;
            }
            else
            {
              alpha = 1 / uxyz.norm();
            }

            // Now calculate the inclination using arc-trig functions
            if (alpha_skip)
            {
              // exit this part of the for loop to the next in the loop
            }
            else
            {
              // "atan_2D=0 atan_3D=1 acos=2 atan_half=3"
              switch (_angular_func)
              {
                case 0: // atan_2D
                {
                  // ATAN returning full 360 degrees instead of just 0-180 for +/-180
                  // --- Oriented azimuth around +x, referenced from +y, increasing toward +z ---
                  const Real x = ngb(0);
                  const Real y = ngb(1);

                  // Value ([-pi, pi] -> [0, 2pi))
                  Real angle = std::atan2(y, x);
                  if (angle < 0.0)
                    angle += 2.0 * libMesh::pi;
                  a_dist = angle;

                  // --- First derivatives wrt phi (= components of normalized ngb) ---
                  // angle = atan2(z, y)
                  // dθ/dx = -y/(x^2+y^2), dθ/dy = x/(x^2+y^2), dθ/dz = 0
                  const Real r2 = x * x + y * y;
                  const Real r2_eps = 1e-20; // robust near axis
                  const Real rdenom = std::max(r2, r2_eps);

                  RealGradient dthetadphi(0.0, 0.0, 0.0);
                  dthetadphi(0) = -y / rdenom; // dθ/dx
                  dthetadphi(1) = x / rdenom;  // dθ/dy
                  dthetadphi(2) = 0.0;         // dθ/dz

                  // --- Second derivatives (Hessian) wrt phi ---
                  // ∂²θ/∂x² =  2xy / r2^2
                  // ∂²θ/∂y² = -2xy / r2^2
                  // ∂²θ/∂x∂y = ∂²θ/∂y∂x = (y^2 - x^2) / r2^2
                  const Real inv_r2 = 1.0 / rdenom;
                  const Real inv_r2_2 = inv_r2 * inv_r2;

                  RealTensorValue d2thetadphi2(0.0);
                  d2thetadphi2(0, 0) = 2.0 * x * y * inv_r2_2;
                  d2thetadphi2(1, 1) = -2.0 * x * y * inv_r2_2;
                  const Real dxy = (y * y - x * x) * inv_r2_2;
                  d2thetadphi2(0, 1) = dxy;
                  d2thetadphi2(1, 0) = dxy;
                  // z-rows/cols remain zero

                  // --- Chain rule to grad(eta_i) and its second derivatives ---
                  // COPIED FROM PREVIOUS PARTS
                  RealTensorValue dij(1, 0, 0, 0, 1, 0, 0, 0, 1);
                  RankTwoTensor dphi_dgradetai;
                  RankThreeTensor d2phi_dgradetai2;
                  for (unsigned int p = 0; p < 3; ++p)
                  {
                    for (unsigned int q = 0; q < 3; ++q)
                    {
                      dphi_dgradetai(p, q) =
                          alpha * dij(p, q) - alpha * alpha * alpha * uxyz(p) * uxyz(q);
                      for (unsigned int r = 0; r < 3; ++r)
                      {
                        d2phi_dgradetai2(p, q, r) =
                            alpha * alpha * alpha *
                            (3 * alpha * alpha * uxyz(p) * uxyz(q) * uxyz(r) - uxyz(p) * dij(q, r) -
                             uxyz(q) * dij(p, r) - uxyz(r) * dij(p, q));
                      }
                    }
                    for (unsigned int q = 0; q < 3; ++q)
                    {
                      dtheta_dgetai(p) += dthetadphi(q) * dphi_dgradetai(q, p);
                      for (unsigned int r = 0; r < 3; ++r)
                        for (unsigned int s = 0; s < 3; ++s)
                        {
                          d2theta_dgetai2(q, r) +=
                              d2thetadphi2(p, s) * dphi_dgradetai(p, r) * dphi_dgradetai(s, q) +
                              dthetadphi(p) * d2phi_dgradetai2(p, q, r);
                        }
                    }
                  }
                  // Finally, after all the checks and calcs we save this pair to ij -> k
                  theta[k] = a_dist;
                  dtheta[k] = dtheta_dgetai;    // this is d/d \nabla \eta_i = -d/dj
                  d2theta[k] = d2theta_dgetai2; // this is d2/d \nabla \eta_i^2 = d2/djj = -d2/dij

                  _inclination[_qp] = 1.0;
                  break;
                }
                default:
                  mooseError("Unknown angular_func = ", _angular_func);
                  break;
              }
            }

            // if (!alpha_skip)
            // {
            // _inc_pairs.push_back(a_dist); // Radians // * 180 / M_PI for degrees
            // // _dincdgradetai_pairs.push_back(dincij_dgradetai);
            // dtheta_dgeta_list.push_back(dinc_dgetai);
            // d2theta_dgeta2_list.push_back(d2inc_dgetai2);
            // _hgb_pairs.push_back((*_vals[i])[_qp] * (*_vals[i])[_qp] * (*_vals[j])[_qp] *
            //                      (*_vals[j])[_qp]);
            // }
          }
          // else
          // {
          //   // IF ngb < cutoff we dont count this pair
          //   ngb = 0.0;
          //   // If norms of gradients too small/flat
          //   // a_dist = (libMesh::pi / (2 * _theta_pre)) - _inc_ij_0;
          //   // a_dist = -1;
          // }
        }
  }
}
