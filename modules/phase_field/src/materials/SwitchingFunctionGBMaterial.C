//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SwitchingFunctionGBMaterial.h"
#include "MooseException.h"

registerMooseObject("PhaseFieldApp", SwitchingFunctionGBMaterial);
registerMooseObject("PhaseFieldApp", ADSwitchingFunctionGBMaterial);

template <bool is_ad>
InputParameters
SwitchingFunctionGBMaterialTempl<is_ad>::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<MaterialPropertyName>(
      "h_name", "Name of the switching function material property for the grain boundaries");
  params.addRequiredCoupledVar("grain_ops", "Vector of order parameters for the given phase");
  // params.addRequiredCoupledVar("all_ops", "Vector of all order parameters for all phases that you
  // want derivatives wrt");
  params.addClassDescription("Calculates the switching function for a grain boundary in a "
                             "multi-phase, multi-order parameter model");
  return params;
}

template <bool is_ad>
SwitchingFunctionGBMaterialTempl<is_ad>::SwitchingFunctionGBMaterialTempl(
    const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _h_name(this->getParam<MaterialPropertyName>("h_name")),
    _num_eta_gb(coupledComponents("grain_ops")),
    _eta_gb(coupledGenericValues<is_ad>("grain_ops")),
    _eta_gb_names(coupledNames("grain_ops")),
    // _num_eta(coupledComponents("all_ops")),
    // _eta(coupledGenericValues<is_ad>("all_ops")),
    // _eta_names(coupledNames("all_ops")),
    // _is_p(_num_eta),
    _prop_h(declareGenericProperty<Real, is_ad>(_h_name)),
    _prop_dh(_num_eta_gb),
    _prop_d2h(_num_eta_gb)
{
  // Declare h derivative properties
  for (unsigned int i = 0; i < _num_eta_gb; ++i)
    _prop_d2h[i].resize(_num_eta_gb, NULL);

  for (unsigned int i = 0; i < _num_eta_gb; ++i)
  {
    _prop_dh[i] = &this->template declarePropertyDerivative<Real, is_ad>(_h_name, _eta_gb_names[i]);
    for (unsigned int j = i; j < _num_eta_gb; ++j)
    {
      _prop_d2h[i][j] = _prop_d2h[j][i] = &this->template declarePropertyDerivative<Real, is_ad>(
          _h_name, _eta_gb_names[i], _eta_gb_names[j]);
    }
  }

  // // Determine which order parameters in the list of all etas belong to phase p
  // for (unsigned int i = 0; i < _num_eta; ++i)
  // {
  //   _is_p[i] = false;
  //   for (unsigned int j = 0; j < _num_eta_p; ++j)
  //   {
  //     if (_eta_names[i] == _eta_p_names[j])
  //       _is_p[i] = true;
  //   }
  // }
}

template <bool is_ad>
void
SwitchingFunctionGBMaterialTempl<is_ad>::computeQpProperties()
{
  GenericReal<is_ad> hgb = 0.0;

  for (unsigned int i = 0; i < _num_eta_gb; ++i)
    for (unsigned int j = i + 1; j < _num_eta_gb; ++j)
    {
      hgb += (*_eta_gb[i])[_qp] * (*_eta_gb[i])[_qp] * (*_eta_gb[j])[_qp] * (*_eta_gb[j])[_qp];
    }

  _prop_h[_qp] = 16 * hgb;

  // For derivatives: sum of other squares
  // GenericReal<is_ad> sum_other = 0.0;
  std::vector<GenericReal<is_ad>> sum_other(_num_eta_gb, 0.0);
  for (unsigned int i = 0; i < _num_eta_gb; ++i)
  {
    GenericReal<is_ad> sum_tmp = 0.0;
    for (unsigned int j = 0; j < _num_eta_gb; ++j)
    {
      if (i != j)
      {
        sum_tmp += (*_eta_gb[j])[_qp] * (*_eta_gb[j])[_qp];
      }
    }
    sum_other[i] = sum_tmp;
  }

  for (unsigned int i = 0; i < _num_eta_gb; ++i)
  {
    // First derivatives
    (*_prop_dh[i])[_qp] = 32 * (*_eta_gb[i])[_qp] * sum_other[i];

    // Second derivatives
    for (unsigned int j = 0; j < _num_eta_gb; ++j)
    {
      if (i == j)
      {
        (*_prop_d2h[i][j])[_qp] = 32 * sum_other[i];
      }
      else
      {
        (*_prop_d2h[i][j])[_qp] = 64 * (*_eta_gb[i])[_qp] * (*_eta_gb[j])[_qp];
      }
    }
  }
}

// explicit instantiation
template class SwitchingFunctionGBMaterialTempl<true>;
template class SwitchingFunctionGBMaterialTempl<false>;
