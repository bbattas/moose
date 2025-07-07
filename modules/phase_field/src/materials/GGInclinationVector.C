//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GGInclinationVector.h"

registerMooseObject("PhaseFieldApp", GGInclinationVector);

InputParameters
GGInclinationVector::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Inclination dependent properties for AGG.");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  params.addParam<UserObjectName>("grain_tracker",
                                  "The GrainTracker UserObject to get values from.");
  return params;
}

GGInclinationVector::GGInclinationVector(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    // : Material(parameters),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _vals_name(_op_num),
    // Inclination cos function
    // _inclination(declareProperty<Real>(getParam<MaterialPropertyName>("inclination_name"))),
    // // Angular distance to the x axis
    // _inclination_distance(declareProperty<Real>("inclination_distance")),
    // Inclination vector for polar plots
    _inclination_vector(declareProperty<RealGradient>("inclination_vector")),
    _ang_dist(declareProperty<Real>("ang_dist")),
    // Grain Tracker/EBSD for GB identification
    _grain_tracker(getUserObject<GrainTracker>("grain_tracker"))
{
  if (_op_num == 0)
    mooseError("Model requires op_num > 0");

  _vals.resize(_op_num);
  _grad_vals.resize(_op_num);

  for (unsigned int i = 0; i < _op_num; ++i)
  {
    _vals[i] = &coupledValue("v", i);
    _vals_name[i] = coupledName("v", i);
    _grad_vals[i] = &coupledGradient("v", i);
  }
}

void
GGInclinationVector::computeQpProperties()
{
  _gb_ij_pairs.clear();
  _hgb_pairs.clear();
  _inc_vec_pairs.clear();
  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());
  for (auto i : index_range(op_to_grains))
  {
    if (op_to_grains[i] == FeatureFloodCount::invalid_id)
      continue;

    _gb_ij_pairs.push_back(i);
  }

  // Make a copy and sort
  _gb_ij_sorted = _gb_ij_pairs;
  std::sort(_gb_ij_sorted.begin(), _gb_ij_sorted.end());

  switch (_gb_ij_pairs.size())
  {
    case 0:
      _inclination_vector[_qp] = RealGradient(0.0, 0.0, 0.0);
      _ang_dist[_qp] = -1;
      break;
    case 1:
      _inclination_vector[_qp] = RealGradient(0.0, 0.0, 0.0);
      _ang_dist[_qp] = -1;
      break;
    default:
      // do all ij pairs if more than 2 vars/features
      for (std::size_t idx1 = 0; idx1 < _gb_ij_sorted.size(); ++idx1)
        for (std::size_t idx2 = idx1 + 1; idx2 < _gb_ij_sorted.size(); ++idx2)
        {
          unsigned int i = _gb_ij_sorted[idx1];
          unsigned int j = _gb_ij_sorted[idx2];
          RealGradient ngb = (*_grad_vals[i])[_qp] - (*_grad_vals[j])[_qp];
          // _hgb_pairs.push_back((*_vals[i])[_qp] * (*_vals[i])[_qp] * (*_vals[j])[_qp] *
          //                      (*_vals[j])[_qp]);
          if (ngb.norm() > 1.0e-10)
          {
            ngb /= ngb.norm(); // Really dont think this should be here?
            _inc_vec_pairs.push_back(ngb);
            _hgb_pairs.push_back((*_vals[i])[_qp] * (*_vals[i])[_qp] * (*_vals[j])[_qp] *
                                 (*_vals[j])[_qp]);
          }
          else
          {
            ngb = 0.0;
          }
        }

      // Now take the assumed multiple and convert
      if (_hgb_pairs.size() > 0)
      {
        RealGradient numer(0.0, 0.0, 0.0);
        Real denom = 0.0;
        Real ang_numer = 0.0;
        for (std::size_t n = 0; n < _hgb_pairs.size(); ++n)
        {
          RealGradient ngb_n = _inc_vec_pairs[n];
          numer += ngb_n * _hgb_pairs[n];
          denom += _hgb_pairs[n];
          Real R = std::sqrt((ngb_n(1) * ngb_n(1)) + (ngb_n(2) * ngb_n(2)));
          Real a_dist = std::atan2(R, ngb_n(0));
          ang_numer += a_dist * _hgb_pairs[n];
        }
        if (denom > 1.0e-10)
        {
          _inclination_vector[_qp] = numer / denom;
          _ang_dist[_qp] = ang_numer / denom;
        }
        else
        {
          _inclination_vector[_qp] = RealGradient(0.0, 0.0, 0.0);
          _ang_dist[_qp] = -1;
        }
      }
      else
      {
        _inclination_vector[_qp] = RealGradient(0.0, 0.0, 0.0);
        _ang_dist[_qp] = -1;
      }
  }
}
