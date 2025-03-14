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
  params.addParam<UserObjectName>("ebsd_reader", "The EBSDReader GeneralUserObject");
  params.addParam<Real>("delta_ij", 0.05, "Anisotropy weight in cos function");
  params.addParam<Real>("inc_ij_0", 0.0, "Inclination function offset in cos function");
  params.addParam<std::vector<MaterialPropertyName>>(
      "gamma_grad_eta_names",
      std::vector<MaterialPropertyName>(),
      "Interfacial / grain boundary gamma parameter names (leave empty for gamma0... gammaN)");
  params.addParam<MaterialPropertyName>("gb_energy_input",
                                        "GB energy before inclination dependence");
  params.addParam<MaterialPropertyName>(
      "gb_energy", "gb_energy", "Inclination dependent GB energy output.");
  params.addParam<Real>("kappa", 1, "Gradient energy constant kappa value");
  params.addParam<Real>("free_energy_m", 1, "Free energy function constant m");
  return params;
}

GGInclinationMaterial::GGInclinationMaterial(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    // : Material(parameters),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    // Inclination cos function
    _inclination(declareProperty<Real>(getParam<MaterialPropertyName>("inclination_name"))),
    // Angular distance to the x axis
    _inclination_distance(declareProperty<Real>("inclination_distance")),
    // Grain Tracker/EBSD for GB identification
    _grain_tracker(getUserObject<GrainTracker>("grain_tracker")),
    _ebsd_reader(getUserObject<EBSDReader>("ebsd_reader")),
    _delta_ij(getParam<Real>("delta_ij")),
    _inc_ij_0(getParam<Real>("inc_ij_0")),
    _gamma(declareProperty<Real>("gamma_inc")),
    _dgammadgrad_eta_name(getParam<std::vector<MaterialPropertyName>>("gamma_grad_eta_names")),
    _dgammadgrad_eta(_op_num),
    // TEMP TEST OUTPUTS
    _testout(declareProperty<Real>("testout")),
    _testout2(declareProperty<Real>("testout2")),
    _gbe(getMaterialProperty<Real>("gb_energy_input")),
    _gbe_inc(declareProperty<Real>(getParam<MaterialPropertyName>("gb_energy"))),
    _kappa(getParam<Real>("kappa")),
    _const_m(getParam<Real>("free_energy_m"))

{
  if (_op_num == 0)
    mooseError("Model requires op_num > 0");

  if (_dgammadgrad_eta_name.size() != 0 && _dgammadgrad_eta_name.size() != _op_num)
    paramError("gamma_grad_eta_names",
               "Specify either as many entries as op_num values or none at all for auto-naming the "
               "gamma gradients with respect to gradient of grain OPs.");

  // automatic names for the gamma properties
  if (_dgammadgrad_eta_name.size() == 0)
  {
    _dgammadgrad_eta_name.resize(_op_num);
    for (unsigned int i = 0; i < _op_num; i++)
      _dgammadgrad_eta_name[i] = "dgammadgrad_eta" + Moose::stringify(i);
  }

  _vals.resize(_op_num);
  _grad_vals.resize(_op_num);
  _incl_tens.resize(_op_num);
  _ang_dist.resize(_op_num);

  for (unsigned int i = 0; i < _op_num; ++i)
  {
    _vals[i] = &coupledValue("v", i);
    _grad_vals[i] = &coupledGradient("v", i);
    // Build the ij tensor/matrix of inclinations
    _incl_tens[i].resize(_op_num);
    _ang_dist[i].resize(_op_num);
    _dgammadgrad_eta[i] = &declareProperty<RealGradient>(_dgammadgrad_eta_name[i]);
  }
}

void
GGInclinationMaterial::computeQpProperties()
{
  // From the ComputeGBMisorientationType:
  // Find out the number of boundary unique_id and save them
  _gb_pairs.clear();
  _gb_op_pairs.clear();
  _gb_ij_pairs.clear();

  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());
  for (auto i : index_range(op_to_grains))
  {
    if (op_to_grains[i] == FeatureFloodCount::invalid_id)
      continue;

    _gb_pairs.push_back(_ebsd_reader.getFeatureID(op_to_grains[i]));
    _gb_op_pairs.push_back((*_vals[i])[_qp]);
    _gb_ij_pairs.push_back(i);
  }

  // Make a copy and sort
  _gb_ij_sorted = _gb_ij_pairs;
  std::sort(_gb_ij_sorted.begin(), _gb_ij_sorted.end());
  _hgb_pairs.clear();
  _inc_pairs.clear();

  switch (_gb_pairs.size())
  {
    case 0:
      break;
    case 1:
      // _inclination[_qp] = _pre_inc[_qp];
      _inclination_distance[_qp] = 0.0;
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
          _hgb_pairs.push_back((*_vals[i])[_qp] * (*_vals[i])[_qp] * (*_vals[j])[_qp] *
                               (*_vals[j])[_qp]);
          if (ngb.norm() > 1.0e-10)
          {
            ngb /= ngb.norm();
            R = std::sqrt((ngb(1) * ngb(1)) + (ngb(2) * ngb(2)));
            a_dist = std::atan2(R, ngb(0));
          }
          else
          {
            ngb = 0.0;
          }
          _inc_pairs.push_back(a_dist); // Radians // * 180 / M_PI for degrees
        }
  }
  // Combine the inclination now!
  // temp weighted combination of just the angle
  if (_inc_pairs.size() > 0)
  {
    Real numer = 0.0;
    Real denom = 0.0;
    Real ang_numer = 0.0;
    for (std::size_t n = 0; n < _inc_pairs.size(); ++n)
    {
      // Real inc_func = 1 + _delta_ij * std::cos(4 * (_inc_pairs[n] + _inc_ij_0));
      numer += (1 + _delta_ij * std::cos(4 * (_inc_pairs[n] + _inc_ij_0))) * _hgb_pairs[n];
      denom += _hgb_pairs[n];
      ang_numer += _inc_pairs[n] * _hgb_pairs[n];
    }
    _inclination_distance[_qp] = ang_numer / denom;
    _inclination[_qp] = numer / denom;
  }
  else
  {
    _inclination_distance[_qp] = 0.0;
    _inclination[_qp] = 1.0;
  }

  _gbe_inc[_qp] = _inclination[_qp] * _gbe[_qp];
  Real g = _gbe_inc[_qp] / (std::sqrt(_kappa * _const_m));
  Real g2 = g * g;
  // Hard-coded coefficients (for example purposes)
  constexpr Real a1 = -3.0944; // coefficient for g2^4
  constexpr Real a2 = -1.8169; // coefficient for g2^3
  constexpr Real a3 = 10.323;  // coefficient for g2^2
  constexpr Real a4 = -8.1819; // coefficient for g2
  constexpr Real a5 = 2.0033;  // constant term

  Real poly_g = (((a1 * g2 + a2) * g2 + a3) * g2 + a4) * g2 + a5;

  _gamma[_qp] = 1 / poly_g;
  // // Test gamma function
  // _gamma[_qp] = (((*_grad_vals[0])[_qp]).norm()) * (((*_grad_vals[0])[_qp]).norm());
  for (unsigned int i = 0; i < _op_num; ++i)
    (*_dgammadgrad_eta[i])[_qp] = 2 * ((*_grad_vals[i])[_qp]);

  // _testout[_qp] = (*_grad_vals[0])[_qp](0);
  // _testout2[_qp] = (*_dgammadgrad_eta[0])[_qp](0);
  _testout[_qp] = _gbe_inc[_qp];
  _testout2[_qp] = g;
  // // Calcluate for
  // // Compute GB type by the number of id
  // // _grains_on_gb[_qp] = _gb_pairs.size();
  // _grains_on_gb[_qp] = 0;
  // _gb_id[_qp] = 0;
  // switch (_gb_pairs.size())
  // {
  //   case 0:
  //     break;
  //   case 1:
  //     _gb_id[_qp] = _gb_ij_pairs[0];
  //     break;
  //   case 2:
  //     // get type by Misorientation angle
  //     _gb_id[_qp] = _gb_ij_pairs[0];
  //     _grains_on_gb[_qp] = _gb_ij_pairs[1];
  //     break;
  //   default:
  //     // get continuous type at triple junction
  //     _gb_id[_qp] = -1;
  // }

  // // Vector of vectors to 1 ij pair output
  // for (unsigned int i = 0; i < _op_num - 1; ++i)
  //   for (unsigned int j = i + 1; j < _op_num; ++j)
  //   {
  //     RealGradient ngb = (*_grad_vals[i])[_qp] - (*_grad_vals[j])[_qp];
  //     Real R = 0.0;
  //     if (ngb.norm() > 1.0e-10)
  //     {
  //       ngb /= ngb.norm();
  //       R = std::sqrt((ngb(1) * ngb(1)) + (ngb(2) * ngb(2)));
  //       _ang_dist[i][j] = std::atan2(R, ngb(0));
  //     }
  //     else
  //     {
  //       ngb = 0.0;
  //       _ang_dist[i][j] = 0.0;
  //     }
  //     // Store to inclination matrix
  //     _incl_tens[i][j] = ngb;
  //   }

  // _inclination[_qp] = _ang_dist[_i_value][_j_value]; // Radians
  // _inclination[_qp] = _ang_dist[_i_value][_j_value] * 180 / M_PI; // Degrees
}
