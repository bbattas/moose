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
  params.addClassDescription(
      "Phase field parameters for polynomial free energy for single component systems");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  params.addParam<MaterialPropertyName>("inclination_name", "Name of inclination material");
  params.addParam<UserObjectName>("grain_tracker",
                                  "The GrainTracker UserObject to get values from.");
  params.addParam<UserObjectName>("ebsd_reader", "The EBSDReader GeneralUserObject");
  params.addParam<unsigned int>("i_value", 0, "i of inclination_{ij} for output");
  params.addParam<unsigned int>("j_value", 0, "j of inclination_{ij} for output");
  // params.addCoupledVar("T", "Temperature variable in Kelvin");
  // params.addRequiredCoupledVar("c", "Concentration");
  // params.addRequiredParam<Real>(
  //     "int_width", "The interfacial width of void surface in the length scale of the problem");
  // params.addParam<Real>(
  //     "length_scale", 1.0e-9, "defines the base length scale of the problem in m");
  // params.addParam<Real>("time_scale", 1.0e-9, "defines the base time scale of the problem");
  // MooseEnum poly_order("4 6 8");
  // params.addRequiredParam<MooseEnum>(
  //     "polynomial_order", poly_order, "Order of polynomial free energy");
  // params.addRequiredParam<Real>("D0", "Diffusivity prefactor for vacancies in m^2/s");
  // params.addRequiredParam<Real>("Em", "Migration energy in eV");
  // params.addRequiredParam<Real>("Ef", "Formation energy in eV");
  // params.addRequiredParam<Real>("surface_energy", "Surface energy in J/m2");
  return params;
}

GGInclinationMaterial::GGInclinationMaterial(const InputParameters & parameters)
  : Material(parameters),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _inclination(declareProperty<Real>(getParam<MaterialPropertyName>("inclination_name"))),
    // Grain Tracker/EBSD for GB identification
    _grain_tracker(getUserObject<GrainTracker>("grain_tracker")),
    _ebsd_reader(getUserObject<EBSDReader>("ebsd_reader")),
    _grains_on_gb(declareProperty<Real>("grains_on_gb")),
    _gb_id(declareProperty<Real>("gb_id")),
    // RealGradient
    _i_value(getParam<unsigned int>("i_value")),
    _j_value(getParam<unsigned int>("j_value"))
// RealTensorValue
// _c(coupledValue("c")),
// _T(coupledValue("T")),
// _M(declareProperty<Real>("M")),
// _grad_M(declareProperty<RealGradient>("grad_M")),
// _kappa(declareProperty<Real>("kappa")),
// _c_eq(declareProperty<Real>("c_eq")),
// _W(declareProperty<Real>("barr_height")),
// _Qstar(declareProperty<Real>("Qstar")),
// _D(declareProperty<Real>("D")),
// _int_width(getParam<Real>("int_width")),
// _length_scale(getParam<Real>("length_scale")),
// _time_scale(getParam<Real>("time_scale")),
// _order(getParam<MooseEnum>("polynomial_order")),
// _D0(getParam<Real>("D0")),
// _Em(getParam<Real>("Em")),
// _Ef(getParam<Real>("Ef")),
// _surface_energy(getParam<Real>("surface_energy")),
// _JtoeV(6.24150974e18), // joule to eV conversion
// _kb(8.617343e-5)       // Boltzmann constant in eV/K
{
  if (_op_num == 0)
    mooseError("Model requires op_num > 0");

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

  // Compute GB type by the number of id
  // _grains_on_gb[_qp] = _gb_pairs.size();
  _grains_on_gb[_qp] = 0;
  _gb_id[_qp] = 0;
  switch (_gb_pairs.size())
  {
    case 0:
      break;
    case 1:
      _gb_id[_qp] = _gb_ij_pairs[0];
      break;
    case 2:
      // get type by Misorientation angle
      _gb_id[_qp] = _gb_ij_pairs[0];
      _grains_on_gb[_qp] = _gb_ij_pairs[1];
      break;
    default:
      // get continuous type at triple junction
      _gb_id[_qp] = -1;
  }

  // Vector of vectors to 1 ij pair output
  for (unsigned int i = 0; i < _op_num - 1; ++i)
    for (unsigned int j = i + 1; j < _op_num; ++j)
    {
      RealGradient ngb = (*_grad_vals[i])[_qp] - (*_grad_vals[j])[_qp];
      Real R = 0.0;
      if (ngb.norm() > 1.0e-10)
      {
        ngb /= ngb.norm();
        R = std::sqrt((ngb(1) * ngb(1)) + (ngb(2) * ngb(2)));
        _ang_dist[i][j] = std::atan2(R, ngb(0));
      }
      else
      {
        ngb = 0.0;
        _ang_dist[i][j] = 0.0;
      }
      // Store to inclination matrix
      _incl_tens[i][j] = ngb;
    }

  // _inclination[_qp] = _ang_dist[_i_value][_j_value]; // Radians
  _inclination[_qp] = _ang_dist[_i_value][_j_value] * 180 / M_PI; // Degrees
}
