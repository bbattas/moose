//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElementalGammaWellMaterial.h"
#include "libmesh/quadrature.h"

registerMooseObject("PhaseFieldApp", ElementalGammaWellMaterial);

InputParameters
ElementalGammaWellMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("New version of elemental gamma rounding for the well funciton only, "
                             "with no dL but constant GBMob so L varies.");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  params.addParam<bool>(
      "skip", false, "Skip all the averaging calculations (for testing purposes).");
  params.addParam<bool>("round_down", false, "Use low inc function values over high.");
  return params;
}

ElementalGammaWellMaterial::ElementalGammaWellMaterial(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _gamma_in(getMaterialProperty<Real>("gamma_qp")),
    _gamma_out(declareProperty<Real>("gamma_asymm")),
    _int_width_in(getMaterialProperty<Real>("int_width_qp")),
    _int_width_out(declareProperty<Real>("int_width")),
    _L_in(getMaterialProperty<Real>("L_qp")),
    _L_out(declareProperty<Real>("L")),
    _skip(getParam<bool>("skip")),
    _round_down(getParam<bool>("round_down"))
{
}

void
ElementalGammaWellMaterial::computeProperties()
{
  // initialize placeholders
  bool gamma_change = false;
  Real min_gamma = 0.0;
  Real max_gamma = 0.0;
  Real min_iw = 0.0;
  Real max_iw = 0.0;
  Real min_L = 0.0;
  Real max_L = 0.0;

  // If skipping, just copy input gamma/int_width to outputs on all qps
  if (_skip)
  {
    for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
    {
      _gamma_out[_qp] = _gamma_in[_qp];
      _int_width_out[_qp] = _int_width_in[_qp];
      _L_out[_qp] = _L_in[_qp];
    }

    return; // Skip all remaining computations
  }
  else
  {
    // initialize min/max from first qp
    min_gamma = _gamma_in[0];
    max_gamma = _gamma_in[0];

    min_iw = _int_width_in[0];
    max_iw = _int_width_in[0];

    min_L = _L_in[0];
    max_L = _L_in[0];

    // Check for uneven values in the qps
    for (_qp = 1; _qp < _qrule->n_points(); ++_qp)
    {
      // check if any qp is different from the first one
      if (_gamma_in[_qp] != _gamma_in[0])
        gamma_change = true;

      // update min / max
      if (_gamma_in[_qp] < min_gamma)
      {
        min_gamma = _gamma_in[_qp];
        max_iw = _int_width_in[_qp];
        min_L = _L_in[_qp];
      }
      if (_gamma_in[_qp] > max_gamma)
      {
        max_gamma = _gamma_in[_qp];
        min_iw = _int_width_in[_qp];
        max_L = _L_in[_qp];
      }
    }

    // Make change
    for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
    {
      if (gamma_change)
      {
        if (_round_down)
        {
          _gamma_out[_qp] = min_gamma;
          _int_width_out[_qp] = max_iw;
          _L_out[_qp] = min_L;
        }
        else
        {
          _gamma_out[_qp] = max_gamma;
          _int_width_out[_qp] = min_iw;
          _L_out[_qp] = max_L;
        }
      }
      else
      {
        _gamma_out[_qp] = _gamma_in[_qp];
        _int_width_out[_qp] = _int_width_in[_qp];
        _L_out[_qp] = _L_in[_qp];
      }
    }
  }
}
