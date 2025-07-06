//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GGInclinationMaterial.h"

registerMooseObject("PhaseFieldApp", GGInclinationMaterial);

InputParameters
GGInclinationMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Inclination dependent properties for AGG.");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  params.addParam<MaterialPropertyName>("inclination_name",
                                        "Name of inclination cos function material output");
  params.addParam<UserObjectName>("grain_tracker",
                                  "The GrainTracker UserObject to get values from.");
  // params.addParam<UserObjectName>("ebsd_reader", "The EBSDReader GeneralUserObject");
  params.addParam<Real>("delta_ij", 0.05, "Anisotropy weight in cos function");
  params.addParam<Real>(
      "theta_prefactor", 4.0, "Multiplier in cos function (cos(n(theta + theta0)))");
  params.addParam<Real>("inc_ij_0", 0.0, "Inclination function offset in cos function");
  params.addParam<std::vector<MaterialPropertyName>>(
      "dgamma_grad_eta_names",
      std::vector<MaterialPropertyName>(),
      "Interfacial / grain boundary gamma parameter names (leave empty for gamma0... gammaN)");
  params.addParam<std::vector<MaterialPropertyName>>(
      "d2gamma_grad_eta2_names",
      std::vector<MaterialPropertyName>(),
      "Interfacial / grain boundary gamma parameter names (leave empty for gamma0... gammaN)");
  params.addParam<MaterialPropertyName>("gb_energy_input",
                                        "GB energy before inclination dependence");
  params.addParam<MaterialPropertyName>(
      "gb_energy", "gb_energy", "Inclination dependent GB energy output.");
  params.addParam<Real>("kappa", 1, "Gradient energy constant kappa value");
  params.addParam<Real>("free_energy_m", 1, "Free energy function constant m");
  params.addParam<bool>(
      "L_of_eta", false, "Is L a function of eta, requiring those derivatives to be defined.");
  params.addParam<MaterialPropertyName>("L0", "AC mobility prefactor/reference value material.");
  params.addParam<MaterialPropertyName>("gamma0",
                                        "gamma prefactor/reference value material (for mu calc).");
  params.addParam<MaterialPropertyName>("L_name", "L_aniso", "Name of anisotropic L output.");
  params.addParam<MaterialPropertyName>(
      "gamma_name", "gamma_aniso", "Name of anisotropic gamma output.");
  params.addParam<MaterialPropertyName>("mu_name", "mu_aniso", "Name of anisotropic mu output.");
  params.addParam<bool>(
      "continuous", false, "Disregard GT and calculate for all variables everywhere.");
  MooseEnum angular_func("atan=0 acos=1", "atan");
  params.addParam<MooseEnum>(
      "angular_func", angular_func, "Which angular distance function to use.");
  params.addParam<Real>("gradtol", 1e-6, "Gradient magnitude tolerance");
  return params;
}

GGInclinationMaterial::GGInclinationMaterial(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    // : Material(parameters),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _vals_name(_op_num),
    // Inclination cos function
    _inclination(declareProperty<Real>(getParam<MaterialPropertyName>("inclination_name"))),
    // Angular distance to the x axis
    _inclination_distance(declareProperty<Real>("inclination_distance")),
    // Inclination vector for polar plots
    // _inclination_vector(declareProperty<Real>("inclination_vector")),
    // Grain Tracker/EBSD for GB identification
    _grain_tracker(getUserObject<GrainTracker>("grain_tracker")),
    // _ebsd_reader(getUserObject<EBSDReader>("ebsd_reader")),
    _delta_ij(getParam<Real>("delta_ij")),
    _theta_pre(getParam<Real>("theta_prefactor")),
    _inc_ij_0(getParam<Real>("inc_ij_0")),
    _gamma(declareProperty<Real>(getParam<MaterialPropertyName>("gamma_name"))),
    _dgammadgrad_eta_name(getParam<std::vector<MaterialPropertyName>>("dgamma_grad_eta_names")),
    _dgammadgrad_eta(_op_num),
    _d2gammadgrad_eta2_name(getParam<std::vector<MaterialPropertyName>>("d2gamma_grad_eta2_names")),
    _d2gammadgrad_eta2(_op_num),
    // Inclincation dependent L
    _L_of_eta(getParam<bool>("L_of_eta")),
    _L_name(getParam<MaterialPropertyName>("L_name")),
    _L(declareProperty<Real>(_L_name)),
    _dLdeta(_op_num),
    _d2Ldetadeta(_op_num),
    _dLdgrad_eta_name(_op_num),
    _dLdgrad_eta(_op_num),
    _d2Ldgrad_eta2_name(_op_num),
    _d2Ldgrad_eta2(_op_num),
    _d2Ldgrad_etadeta(_op_num),
    // TEMP TEST OUTPUTS
    _testout(declareProperty<Real>("testout")),
    _testout2(declareProperty<Real>("testout2")),
    _testoutgrad(declareProperty<RealGradient>("testoutgrad")),
    _testoutgrad2(declareProperty<RealTensorValue>("testoutgrad2")),
    _inclin(declareProperty<RealGradient>("inclin_vec")),
    _dadb(declareProperty<RealGradient>("dadb")),
    _d2adb2(declareProperty<RealTensorValue>("d2adb2")),
    _d3adb3(declareProperty<RankTwoTensor>("d3adb3")),
    // Other material properties
    _gbe(getMaterialProperty<Real>("gb_energy_input")),
    _gbe_inc(declareProperty<Real>(getParam<MaterialPropertyName>("gb_energy"))),
    _kappa(getParam<Real>("kappa")),
    _const_m(getParam<Real>("free_energy_m")),
    _L0(getMaterialProperty<Real>("L0")),
    _gamma0(getMaterialProperty<Real>("gamma0")),
    _int_width(declareProperty<Real>("int_width")),
    _mu(declareProperty<Real>(getParam<MaterialPropertyName>("mu_name"))),
    _continuous(getParam<bool>("continuous")),
    _angular_func(getParam<MooseEnum>("angular_func")),
    _intol(getParam<Real>("gradtol"))

{
  if (_op_num == 0)
    mooseError("Model requires op_num > 0");

  if (_L_of_eta)
    mooseError("Inclination L currently only coded for grad u dependent, not u dependent.");

  if (_dgammadgrad_eta_name.size() != 0 && _dgammadgrad_eta_name.size() != _op_num)
    paramError("dgamma_grad_eta_names",
               "Specify either as many entries as op_num values or none at all for auto-naming the "
               "gamma gradients with respect to gradient of grain OPs.");
  if (_d2gammadgrad_eta2_name.size() != 0 && _d2gammadgrad_eta2_name.size() != _op_num)
    paramError("d2gamma_grad_eta2_names",
               "Specify either as many entries as op_num values or none at all for auto-naming the "
               "gamma gradients with respect to gradient of grain OPs.");

  // automatic names for the gamma properties
  if (_dgammadgrad_eta_name.size() == 0)
  {
    _dgammadgrad_eta_name.resize(_op_num);
    for (unsigned int i = 0; i < _op_num; i++)
      _dgammadgrad_eta_name[i] = "dgammadgrad_eta_" + Moose::stringify(i);
  }
  if (_d2gammadgrad_eta2_name.size() == 0)
  {
    _d2gammadgrad_eta2_name.resize(_op_num);
    for (unsigned int i = 0; i < _op_num; i++)
      _d2gammadgrad_eta2_name[i] = "d2gammadgrad_eta2_" + Moose::stringify(i);
  }

  _vals.resize(_op_num);
  _grad_vals.resize(_op_num);
  _incl_tens.resize(_op_num);
  _ang_dist.resize(_op_num);

  for (unsigned int i = 0; i < _op_num; ++i)
  {
    _vals[i] = &coupledValue("v", i);
    _vals_name[i] = coupledName("v", i);
    _grad_vals[i] = &coupledGradient("v", i);
    // Build the ij tensor/matrix of inclinations
    _incl_tens[i].resize(_op_num);
    _ang_dist[i].resize(_op_num);
    _dgammadgrad_eta[i] = &declareProperty<RealGradient>(_dgammadgrad_eta_name[i]);
    _d2gammadgrad_eta2[i] = &declareProperty<RealTensorValue>(_d2gammadgrad_eta2_name[i]);
    // L derivatives
    _dLdeta[i] = &declarePropertyDerivative<Real>(_L_name, _vals_name[i]);
    _d2Ldetadeta[i].resize(_op_num);
    _dLdgrad_eta_name[i] = "dLdgrad_eta_" + Moose::stringify(i);
    _dLdgrad_eta[i] = &declareProperty<RealGradient>(_dLdgrad_eta_name[i]);
    _d2Ldgrad_eta2_name[i] = "d2Ldgrad_eta2_" + Moose::stringify(i);
    _d2Ldgrad_eta2[i] = &declareProperty<RealTensorValue>(_d2Ldgrad_eta2_name[i]);
    _d2Ldgrad_etadeta[i].resize(_op_num);
    // _d2Ldetadgrad_eta[i] =
    //     &declarePropertyDerivative<RealGradient>(_dLdgrad_eta_name[i], _vals_name[i]);
    // _d2Ldetadgrad_eta[i].resize(_op_num); // only using same eta for eta and grad eta here
    for (unsigned int j = 0; j <= i; ++j)
    {
      _d2Ldetadeta[j][i] = &declarePropertyDerivative<Real>(_L_name, _vals_name[j], _vals_name[i]);
      _d2Ldgrad_etadeta[i][j] =
          &declarePropertyDerivative<RealGradient>(_dLdgrad_eta_name[i], _vals_name[j]);
    }
  }
}

void
GGInclinationMaterial::computeQpProperties()
{
  // Zero out the derivatives for gamma wrt grad_eta
  for (unsigned int i = 0; i < _op_num; ++i)
    (*_dgammadgrad_eta[i])[_qp] = RealGradient(0.0, 0.0, 0.0);
  // From the ComputeGBMisorientationType:
  // Find out the number of boundary unique_id and save them
  // _gb_pairs.clear();
  // _gb_op_pairs.clear();
  _gb_ij_pairs.clear();
  if (_continuous)
  {
    // Give it every OP index once
    _gb_ij_pairs.resize(_op_num);
    std::iota(_gb_ij_pairs.begin(), _gb_ij_pairs.end(), 0);
  }
  else
  {
    const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());
    for (auto i : index_range(op_to_grains))
    {
      if (op_to_grains[i] == FeatureFloodCount::invalid_id)
        continue;

      // Real i_norm = (*_grad_vals[i])[_qp].norm();
      // if (i_norm >= 0.01)
      //   _gb_ij_pairs.push_back(i);
      // continue;

      // _gb_pairs.push_back(_ebsd_reader.getFeatureID(op_to_grains[i]));
      // _gb_op_pairs.push_back((*_vals[i])[_qp]);
      _gb_ij_pairs.push_back(i);
    }
  }

  // Make a copy and sort
  _gb_ij_sorted = _gb_ij_pairs;
  std::sort(_gb_ij_sorted.begin(), _gb_ij_sorted.end());
  _hgb_pairs.clear();
  _inc_pairs.clear();
  // _dincdgradetai_pairs.clear();
  // std::vector<Real> gb_i_list;
  // gb_i_list.clear();

  std::vector<RealGradient> dinc_dgeta_list;
  std::vector<RealTensorValue> d2inc_dgeta2_list;
  dinc_dgeta_list.clear();
  d2inc_dgeta2_list.clear();
  _testout2[_qp] = 0.0;
  _dadb[_qp] = RealGradient(0.0, 0.0, 0.0);
  _d2adb2[_qp] = RealTensorValue(0.0);
  _d3adb3[_qp] = RealTensorValue(0.0);

  switch (_gb_ij_pairs.size())
  {
    case 0:
      break;
    case 1:
      // _inclination[_qp] = _pre_inc[_qp];
      // _inclination_distance[_qp] = 0.0;
      _inclination_distance[_qp] = (libMesh::pi / (2 * _theta_pre)) - _inc_ij_0;
      _inclination[_qp] = 1.0;
      break;
    // case 2:
    //   unsigned int i = _gb_ij_sorted[0];
    //   unsigned int j = _gb_ij_sorted[1];
    //   RealGradient ngb = (*_grad_vals[i])[_qp] - (*_grad_vals[j])[_qp];
    //   Real R = 0.0;
    //   Real a_dist = 0.0;
    //   _hgb_pairs.push_back((*_vals[i])[_qp] * (*_vals[i])[_qp] * (*_vals[j])[_qp] *
    //                        (*_vals[j])[_qp]);
    //   if (ngb.norm() > 1.0e-10)
    //   {
    //     ngb /= ngb.norm();
    //     R = std::sqrt((ngb(1) * ngb(1)) + (ngb(2) * ngb(2)));
    //     a_dist = std::atan2(R, ngb(0));
    //   }
    //   else
    //   {
    //     ngb = 0.0;
    //   }
    //   _inc_pairs.push_back(a_dist * 180 / M_PI); // Degrees
    //   break;
    default:
      // do all ij pairs if more than 2 vars/features
      for (std::size_t idx1 = 0; idx1 < _gb_ij_sorted.size(); ++idx1)
        for (std::size_t idx2 = idx1 + 1; idx2 < _gb_ij_sorted.size(); ++idx2)
        {
          unsigned int i = _gb_ij_sorted[idx1];
          unsigned int j = _gb_ij_sorted[idx2];
          RealGradient ngb = (*_grad_vals[i])[_qp] - (*_grad_vals[j])[_qp];
          Real R = 0.0;
          Real a_dist = 0.0;
          RealGradient dincij_dgradetai(0.0, 0.0, 0.0);
          RealGradient dinc_dgetai(0.0, 0.0, 0.0);
          RealTensorValue d2inc_dgetai2(0.0);
          _hgb_pairs.push_back((*_vals[i])[_qp] * (*_vals[i])[_qp] * (*_vals[j])[_qp] *
                               (*_vals[j])[_qp]);
          // _testout2[_qp] = ngb.norm();
          // Real i_norm = (*_grad_vals[i])[_qp].norm();
          // Real j_norm = (*_grad_vals[j])[_qp].norm();
          // if (i_norm >= _intol)
          //   _testout2[_qp] = (*_grad_vals[j])[_qp].norm();
          // else
          //   _testout2[_qp] = 2.0;
          // Real j_norm = (*_grad_vals[j])[_qp].norm();
          if (ngb.norm() >
              1.0e-10) // &&
                       // (i_norm >= _intol ||
                       //  j_norm >= _intol)) // && i_norm >= 0.01 && j_norm >= 0.01) // 1.0e-10
          {
            ngb /= ngb.norm();       // Really dont think this should be here?
            RealGradient uxyz = ngb; // thought the math needed uxyz in derivs not normalized?
            // Real alpha = 1 / uxyz.norm();
            const Real alpha_tol = 0.1;
            Real alpha = 1 / std::max(uxyz.norm(), alpha_tol);
            // ngb /= ngb.norm();
            if (_angular_func == 0)
            {
              R = std::sqrt((ngb(1) * ngb(1)) + (ngb(2) * ngb(2)));
              a_dist = std::atan2(R, ngb(0));
              // Derivs- Only if there is a non-zero tilt (R exists & nonzero)
              // else we get NaN in dy and dz
              if (R > 1e-10)
              {
                // Derivatives
                RealTensorValue dij(1, 0, 0, 0, 1, 0, 0, 0, 1);
                RankTwoTensor dphi_dgradetai;
                RankThreeTensor d2phi_dgradetai2;
                RealGradient dAdphi(1.0, 0.0, 0.0);
                RealGradient dBdphi(0.0, ngb(1) / R, ngb(2) / R);
                RealTensorValue d2Adphiji(0, 0, 0, 0, 0, 0, 0, 0, 0);
                RealTensorValue d2Bdphiji(0,
                                          0,
                                          0,
                                          0,
                                          ngb(2) * ngb(2) / (R * R * R),
                                          -ngb(1) * ngb(2) / (R * R * R),
                                          0,
                                          -ngb(1) * ngb(2) / (R * R * R),
                                          ngb(1) * ngb(1) / (R * R * R));
                Real dincdA = -R;
                Real dincdB = ngb(0);
                Real d2incdA2 = 2 * R * ngb(0);
                Real d2incdB2 = -2 * R * ngb(0);
                Real d2incdAB = R * R - ngb(0) * ngb(0);
                RealGradient dincdphi(0.0, 0.0, 0.0);
                RealTensorValue d2incdphi2(0.0);
                for (unsigned int p = 0; p < 3; ++p)
                {
                  dincdphi(p) = dincdA * dAdphi(p) + dincdB * dBdphi(p);
                  for (unsigned int q = 0; q < 3; ++q)
                  {
                    dphi_dgradetai(p, q) =
                        alpha * dij(p, q) - alpha * alpha * alpha * uxyz(p) * uxyz(q);
                    d2incdphi2(p, q) = (d2incdA2 * dAdphi(q) * dAdphi(p)) +
                                       d2incdAB * (dAdphi(q) * dBdphi(p) + dAdphi(p) * dBdphi(q)) +
                                       (d2incdB2 * dBdphi(q) * dBdphi(p)) +
                                       (dincdA * d2Adphiji(p, q)) + (dincdB * d2Bdphiji(p, q));
                    for (unsigned int r = 0; r < 3; ++r)
                    {
                      d2phi_dgradetai2(p, q, r) =
                          alpha * alpha * alpha *
                          (3 * alpha * alpha * uxyz(p) * uxyz(q) * uxyz(r) - uxyz(p) * dij(q, r) -
                           uxyz(q) * dij(p, r) - uxyz(r) * dij(p, q));
                    }
                  }
                  // }
                  // for (unsigned int p = 0; p < 3; ++p)
                  // {
                  for (unsigned int q = 0; q < 3; ++q)
                  {
                    dinc_dgetai(p) += dincdphi(q) * dphi_dgradetai(q, p);
                    for (unsigned int r = 0; r < 3; ++r)
                      for (unsigned int s = 0; s < 3; ++s)
                      {
                        d2inc_dgetai2(q, r) +=
                            d2incdphi2(p, s) * dphi_dgradetai(p, r) * dphi_dgradetai(s, q) +
                            dincdphi(p) * d2phi_dgradetai2(p, q, r);
                      }
                  }
                }
                if (idx1 == 0 && idx2 == 1)
                {
                  _testout2[_qp] = alpha; // uxyz.norm();
                  _dadb[_qp] = uxyz;      // dincdphi;
                  _d2adb2[_qp] = d2incdphi2;
                  _d3adb3[_qp] = dphi_dgradetai;
                }
                // Real temp_dincldgradetaix =
                //     (-alpha / R) * ((1 - (uxyz(0) * uxyz(0) * alpha * alpha)) +
                //                     (ngb(1) / ngb(0)) * (uxyz(0) * uxyz(1) * alpha * alpha) +
                //                     (ngb(2) / ngb(0)) * (uxyz(0) * uxyz(2) * alpha * alpha));
                // Real temp_dincldgradetaiy =
                //     (-alpha / R) * ((uxyz(0) * uxyz(1) * alpha * alpha) +
                //                     (ngb(1) / ngb(0)) * (1 - (uxyz(1) * uxyz(1) * alpha * alpha))
                //                     + (ngb(2) / ngb(0)) * (uxyz(1) * uxyz(2) * alpha * alpha));
                // Real temp_dincldgradetaiz =
                //     (-alpha / R) * ((uxyz(0) * uxyz(2) * alpha * alpha) +
                //                     (ngb(1) / ngb(0)) * (uxyz(2) * uxyz(1) * alpha * alpha) +
                //                     (ngb(2) / ngb(0)) * (1 - (uxyz(2) * uxyz(2) * alpha *
                //                     alpha)));
                // dincij_dgradetai =
                //     RealGradient(temp_dincldgradetaix, temp_dincldgradetaiy,
                //     temp_dincldgradetaiz);
              }
            }
            else if (_angular_func == 1)
            {
              // Trim values outside allowed range (which shouldnt exist anyway in ngb)
              Real cx = std::max(-1.0, std::min(1.0, ngb(0)));
              a_dist = std::acos(cx);
              // Derivatives
              // acos (or inc) to phi
              // d/dx acos(x) = -1/sqrt(1-x^2)
              Real eps = 1e-3;
              Real den_root = std::sqrt(1 - cx * cx + eps * eps);
              RealGradient dincdphi(0.0, 0.0, 0.0);
              Real d2incdphi2_xx = 0.0;
              if (den_root > 1e-10)
              {
                dincdphi(0) = -1 / den_root;
                d2incdphi2_xx = -cx / (den_root * den_root * den_root);
              }
              RealTensorValue d2incdphi2(0.0);
              d2incdphi2(0, 0) = d2incdphi2_xx;
              // phi to gradeta
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
                // }
                // for (unsigned int p = 0; p < 3; ++p)
                // {
                for (unsigned int q = 0; q < 3; ++q)
                {
                  dinc_dgetai(p) += dincdphi(q) * dphi_dgradetai(q, p);
                  for (unsigned int r = 0; r < 3; ++r)
                    for (unsigned int s = 0; s < 3; ++s)
                    {
                      d2inc_dgetai2(q, r) +=
                          d2incdphi2(p, s) * dphi_dgradetai(p, r) * dphi_dgradetai(s, q) +
                          dincdphi(p) * d2phi_dgradetai2(p, q, r);
                    }
                }
              }
              if (idx1 == 0 && idx2 == 1)
              {
                _dadb[_qp] = dincdphi;
                _d2adb2[_qp] = d2incdphi2;
                _d3adb3[_qp] = dphi_dgradetai;
              }
            }
            else
            {
              mooseError("Angular function must be 0 or 1");
            }

            // _inc_pairs.push_back(a_dist); // Radians // * 180 / M_PI for degrees
            // // _dincdgradetai_pairs.push_back(dincij_dgradetai);
            // dinc_dgeta_list.push_back(dinc_dgetai);
            // d2inc_dgeta2_list.push_back(d2inc_dgetai2);
          }
          else
          {
            ngb = 0.0;
            // If norms of gradients too small/flat
            // a_dist = (libMesh::pi / (2 * _theta_pre)) - _inc_ij_0;
            // a_dist = -1;
          }
          _inc_pairs.push_back(a_dist); // Radians // * 180 / M_PI for degrees
          // _dincdgradetai_pairs.push_back(dincij_dgradetai);
          dinc_dgeta_list.push_back(dinc_dgetai);
          d2inc_dgeta2_list.push_back(d2inc_dgetai2);
        }
  }
  // Combine the inclination now!
  // temp weighted combination of just the angle
  std::vector<RealGradient> dINC_dGradEta(_op_num, RealGradient(0.0));
  std::vector<RealGradient> dfinc_dgeta(_op_num, RealGradient(0.0));
  std::vector<RealTensorValue> d2finc_dgeta2(_op_num, RealTensorValue(0.0));
  if (_inc_pairs.size() > 0)
  {
    Real numer = 0.0;
    Real denom = 0.0;
    Real ang_numer = 0.0;
    for (std::size_t n = 0; n < _inc_pairs.size(); ++n)
    {
      // Real inc_func = 1 + _delta_ij * std::cos(4 * (_inc_pairs[n] + _inc_ij_0));
      numer += (1 + _delta_ij * std::cos(_theta_pre * (_inc_pairs[n] + _inc_ij_0))) * _hgb_pairs[n];
      denom += _hgb_pairs[n];
      ang_numer += _inc_pairs[n] * _hgb_pairs[n];
      // Derivatives
      // dinc_dgradeta_i += (-4 * _delta_ij * std::sin(4.0 * (_inc_pairs[n] + _inc_ij_0))) *
      // _hgb_pairs[n] * _dincdgradetai_pairs[n];
    }
    if (denom > 1.0e-10)
    {
      _inclination_distance[_qp] = ang_numer / denom;
      _inclination[_qp] = numer / denom;

      // Derivatives
      std::size_t pair_index = 0;
      for (std::size_t idx1 = 0; idx1 < _gb_ij_sorted.size(); ++idx1)
      {
        for (std::size_t idx2 = idx1 + 1; idx2 < _gb_ij_sorted.size(); ++idx2)
        {
          unsigned int p = _gb_ij_sorted[idx1];
          unsigned int q = _gb_ij_sorted[idx2];

          // 1 + delta * cos(...):
          //   derivative wrt inc is   dF_dinc = -4 * delta * sin(4*(inc + inc0))
          //   (if that is your function)
          // Real inc_n = _inc_pairs[pair_index];
          Real dF_dinc = -1 * _theta_pre * _delta_ij *
                         std::sin(_theta_pre * (_inc_pairs[pair_index] + _inc_ij_0));
          Real d2F_dinc2 = -1 * _theta_pre * _theta_pre * _delta_ij *
                           std::cos(_theta_pre * (_inc_pairs[pair_index] + _inc_ij_0));

          // Weighted partial for p:
          //   w_n * dF_dinc * dinc/d(grad eta_p)
          // const RealGradient & dinc_dGradEta_p = _dincdgradetai_pairs[pair_index];
          // RealGradient partial_p =
          //     _hgb_pairs[pair_index] * dF_dinc * _dincdgradetai_pairs[pair_index];
          RealGradient df_dgeta_p = _hgb_pairs[pair_index] * dF_dinc * dinc_dgeta_list[pair_index];
          RealTensorValue d2f_dgeta2_p =
              _hgb_pairs[pair_index] *
              (d2F_dinc2 * libMesh::outer_product(dinc_dgeta_list[pair_index],
                                                  dinc_dgeta_list[pair_index]) +
               dF_dinc * d2inc_dgeta2_list[pair_index]);

          // Weighted partial for q:
          //   w_n * dF_dinc * dinc/d(grad eta_q)
          //   But inc = grad p - grad q => derivative wrt grad q is just - partial_p
          // RealGradient partial_q = -partial_p;
          // RealGradient df_dgeta_q = -df_dgeta_p;

          // Add to the correct OP i’s derivative
          // dINC_dGradEta[p] += partial_p;
          // dINC_dGradEta[q] += partial_q;
          // Alternate
          dfinc_dgeta[p] += df_dgeta_p;
          dfinc_dgeta[q] -= df_dgeta_p;
          d2finc_dgeta2[p] += d2f_dgeta2_p;
          d2finc_dgeta2[q] += d2f_dgeta2_p;

          pair_index++;
        }
      }
      for (unsigned int i = 0; i < _op_num; ++i)
      {
        // dINC_dGradEta[i] /= denom;
        dfinc_dgeta[i] /= denom;
        d2finc_dgeta2[i] /= denom;
      }
    }
  }
  else
  {
    // _inclination_distance[_qp] = 0.0;
    _inclination_distance[_qp] = (libMesh::pi / (2 * _theta_pre)) - _inc_ij_0;
    _inclination[_qp] = 1.0;
  }

  _gbe_inc[_qp] = _inclination[_qp] * _gbe[_qp];
  Real g = _gbe_inc[_qp] / (std::sqrt(_kappa * _const_m));
  Real dg_dfinc = _gbe[_qp] / (std::sqrt(_kappa * _const_m));
  Real g2 = g * g;
  // Hard-coded coefficients (for example purposes)
  constexpr Real a1 = -3.0944; // coefficient for g2^4
  constexpr Real a2 = -1.8169; // coefficient for g2^3
  constexpr Real a3 = 10.323;  // coefficient for g2^2
  constexpr Real a4 = -8.1819; // coefficient for g2
  constexpr Real a5 = 2.0033;  // constant term

  Real poly_g = (((a1 * g2 + a2) * g2 + a3) * g2 + a4) * g2 + a5;
  Real dpoly_g_dg2 = 4 * a1 * g2 * g2 * g2 + 3 * a2 * g2 * g2 + 2 * a3 * g2 + a4;
  Real d2poly_g_dg22 = 12 * a1 * g2 * g2 + 6 * a2 * g2 + 2 * a3;

  _gamma[_qp] = 1 / poly_g;
  // Int width
  Real f0_int =
      (((((0.0788 * poly_g - 0.4955) * poly_g + 1.2244) * poly_g - 1.5281) * poly_g + 1.0686) *
           poly_g -
       0.5563) *
          poly_g +
      0.2907;
  _int_width[_qp] = (std::sqrt(_kappa / _const_m)) * (std::sqrt(1 / f0_int));
  // Build Anisotropic L
  // If gb mobility were a different f we would use that instead of 2 of the gb energy ones
  // These are basically placeholders so that the L derivations will still work if we change these
  Real & finc_m = _inclination[_qp];
  std::vector<RealGradient> & dfinc_m_dgeta = dfinc_dgeta;
  std::vector<RealTensorValue> & d2finc_m_dgeta2 = d2finc_dgeta2;
  // Now for L
  _L[_qp] = _L0[_qp] * _inclination[_qp] * finc_m;

  for (unsigned int i = 0; i < _op_num; ++i)
  {
    // gamma derivatives
    (*_dgammadgrad_eta[i])[_qp] =
        -2 * g * _gamma[_qp] * _gamma[_qp] * dpoly_g_dg2 * dg_dfinc * dfinc_dgeta[i];
    (*_d2gammadgrad_eta2[i])[_qp] =
        8 * _gamma[_qp] * _gamma[_qp] * _gamma[_qp] * dpoly_g_dg2 * dpoly_g_dg2 * g2 * dg_dfinc *
            dg_dfinc * libMesh::outer_product(dfinc_dgeta[i], dfinc_dgeta[i]) -
        _gamma[_qp] * (d2poly_g_dg22 * 4 * g2 * dg_dfinc * dg_dfinc *
                           libMesh::outer_product(dfinc_dgeta[i], dfinc_dgeta[i]) +
                       dpoly_g_dg2 * (2 * dg_dfinc * dg_dfinc *
                                          libMesh::outer_product(dfinc_dgeta[i], dfinc_dgeta[i]) +
                                      2 * g * dg_dfinc * d2finc_dgeta2[i]));
    if (_L_of_eta)
    {
      // Put the actual derivative here if L = f(u)
      // (*_dLdeta[i])[_qp] = 0.0;
      // for (unsigned int j = 0; j < i; ++j)
      //   (*_d2Ldetadeta[j][i])[_qp] = 0.0;
      // (*_d2Ldetadgrad_eta[i])[_qp] = RealGradient(0.0);
      // (*_dLdgrad_eta[i])[_qp] =
      //     _L0[_qp] * (_inclination[_qp] * dfinc_m_dgeta[i] + dfinc_dgeta[i] * finc_m);
      // (*_d2Ldgrad_eta2[i])[_qp] =
      //     _L0[_qp] * (_inclination[_qp] * d2finc_m_dgeta2[i] +
      //                 2 * libMesh::outer_product(dfinc_m_dgeta[i], dfinc_dgeta[i]) +
      //                 d2finc_dgeta2[i] * finc_m);
    }
    else
    {
      // Not a function of u
      (*_dLdeta[i])[_qp] = 0.0;
      for (unsigned int j = 0; j < i; ++j)
      {
        (*_d2Ldetadeta[j][i])[_qp] = 0.0;
        (*_d2Ldgrad_etadeta[i][j])[_qp] = RealGradient(0.0);
      }
      (*_dLdgrad_eta[i])[_qp] =
          _L0[_qp] * (_inclination[_qp] * dfinc_m_dgeta[i] + dfinc_dgeta[i] * finc_m);
      (*_d2Ldgrad_eta2[i])[_qp] =
          _L0[_qp] * (_inclination[_qp] * d2finc_m_dgeta2[i] +
                      2 * libMesh::outer_product(dfinc_m_dgeta[i], dfinc_dgeta[i]) +
                      d2finc_dgeta2[i] * finc_m);
    }
  }

  // mu calc
  _mu[_qp] = _L0[_qp] * _gamma0[_qp] * std::sqrt(_kappa / _const_m) * finc_m;

  // dINC_dGradEta[i];

  // REMEMBER INC here is the f = cos (phi^{\prime})
  _testout[_qp] = 0.0;
  _testoutgrad[_qp] = RealGradient(0.0);
  _testoutgrad2[_qp] = RealTensorValue(0.0);
  if (_inc_pairs.size() > 0)
  {
    _testout[_qp] = _inc_pairs[0];
    _testoutgrad[_qp] = dinc_dgeta_list[0];    // dfinc_dgeta[0];
    _testoutgrad2[_qp] = d2inc_dgeta2_list[0]; // d2finc_dgeta2[0];
  }
  RealGradient inclin = (*_grad_vals[0])[_qp] - (*_grad_vals[1])[_qp];
  Real mag = inclin.norm();
  if (mag > 1e-4)
    _inclin[_qp] = inclin / mag;
  else
    _inclin[_qp] = RealGradient(0.0);
  // const Real inc_tol = 1e-8;
  // if (std::abs(_inclin[_qp].norm() - 1.0) <= inc_tol)
  //   _testout2[_qp] = 1.0; //_gamma[_qp]; //_inc_pairs[1];
  // else
  //   _testout2[_qp] = 0.0;
  // _testout2[_qp] = _inclin[_qp].norm();

  // dinc_dgeta_list[0];
  // _testout[_qp] = (*_grad_vals[0])[_qp].norm();
  // // _testout2[_qp] = testout_x;
  // _testout2[_qp] = testout_den;

  // _incder_temp[_qp] = dfinc_dgeta[0];
}
