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
  params.addClassDescription("Calculates the switching function for a grain boundary in a "
                             "multi-phase, multi-order parameter model");
  params.addRequiredParam<MaterialPropertyName>(
      "h_name", "Name of the switching function material property for the grain boundaries");
  // params.addRequiredCoupledVar("grain_ops", "Vector of order parameters for the given phase");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  params.addParam<Real>("hgb_threshold", 0.0, "Lower limit cutoff for hgb.");
  MooseEnum func_type("base=0 lambda=1 combined=2", "base");
  params.addParam<MooseEnum>("func_type", func_type, "Which hgb option to use.");

  return params;
}

template <bool is_ad>
SwitchingFunctionGBMaterialTempl<is_ad>::SwitchingFunctionGBMaterialTempl(
    const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _h_name(this->getParam<MaterialPropertyName>("h_name")),
    // Old Manual op input
    // _num_eta_gb(coupledComponents("grain_ops")),
    // _eta_gb(coupledGenericValues<is_ad>("grain_ops")),
    // _eta_gb_names(coupledNames("grain_ops")),
    // Automatic op input
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _vals_names(coupledNames("v")),
    // Details
    _hgb_threshold(getParam<Real>("hgb_threshold")),
    _func_type(getParam<MooseEnum>("func_type")),
    // Output
    _prop_h(declareGenericProperty<Real, is_ad>(_h_name)),
    _prop_dh(_op_num),
    _prop_d2h(_op_num)
{
  // Declare h derivative properties
  for (unsigned int i = 0; i < _op_num; ++i)
    _prop_d2h[i].resize(_op_num, NULL);

  for (unsigned int i = 0; i < _op_num; ++i)
  {
    _prop_dh[i] = &this->template declarePropertyDerivative<Real, is_ad>(_h_name, _vals_names[i]);
    for (unsigned int j = i; j < _op_num; ++j)
    {
      _prop_d2h[i][j] = _prop_d2h[j][i] = &this->template declarePropertyDerivative<Real, is_ad>(
          _h_name, _vals_names[i], _vals_names[j]);
    }
  }
}

template <bool is_ad>
void
SwitchingFunctionGBMaterialTempl<is_ad>::computeQpProperties()
{
  // OTHERWISE we actually calculate
  GenericReal<is_ad> hgb = 0.0;
  // GenericReal<is_ad> hgb_val = 0.0;

  switch (_func_type)
  {
    case 0:
    {
      // Previous use: $16 \sum \eta_i^2 \eta_j^2$
      // HGB
      for (unsigned int i = 0; i < _op_num; ++i)
        for (unsigned int j = i + 1; j < _op_num; ++j)
        {
          hgb += (*_vals[i])[_qp] * (*_vals[i])[_qp] * (*_vals[j])[_qp] * (*_vals[j])[_qp];
        }
      hgb *= 16;

      // DERIVATIVES
      for (unsigned int i = 0; i < _op_num; ++i)
      {
        // Calculate sum of non i OPs squares for derivatives
        GenericReal<is_ad> sum_other = 0.0;
        for (unsigned int j = 0; j < _op_num; ++j)
        {
          if (i != j)
          {
            sum_other += (*_vals[j])[_qp] * (*_vals[j])[_qp];
          }
        }
        // First derivatives
        (*_prop_dh[i])[_qp] = 32 * (*_vals[i])[_qp] * sum_other;
        // Second derivatives only for same var twice (rest assuming 0)
        (*_prop_d2h[i][i])[_qp] = 32 * sum_other;
      }
      break;
    }

    case 1:
    {
      // Previous use: $4 (1 - \sum \eta_i^2)^2$
      // HGB
      GenericReal<is_ad> sum_gr = 0.0;
      for (unsigned int i = 0; i < _op_num; ++i)
        sum_gr += (*_vals[i])[_qp] * (*_vals[i])[_qp];

      hgb = 4 * (1 - sum_gr) * (1 - sum_gr);

      // DERIVATIVES
      for (unsigned int i = 0; i < _op_num; ++i)
      {
        // First derivatives
        (*_prop_dh[i])[_qp] = 16 * (*_vals[i])[_qp] * (sum_gr - 1);
        // Second derivatives only for same var twice (rest assuming 0)
        (*_prop_d2h[i][i])[_qp] = 16 * (sum_gr - 1 + (2 * (*_vals[i])[_qp] * (*_vals[i])[_qp]));
      }
      break;
    }

    case 2:
    {
      // HGB 0: $16 \sum \eta_i^2 \eta_j^2$
      GenericReal<is_ad> hgb_cross = 0.0;
      for (unsigned int i = 0; i < _op_num; ++i)
        for (unsigned int j = i + 1; j < _op_num; ++j)
        {
          hgb_cross += (*_vals[i])[_qp] * (*_vals[i])[_qp] * (*_vals[j])[_qp] * (*_vals[j])[_qp];
        }
      hgb_cross *= 16;

      // Combined version for varied Iw
      // HGB 1: $4 (1 - \sum \eta_i^2)^2$
      GenericReal<is_ad> hgb_ssq = 0.0;
      GenericReal<is_ad> sum_gr = 0.0;
      for (unsigned int i = 0; i < _op_num; ++i)
        sum_gr += (*_vals[i])[_qp] * (*_vals[i])[_qp];

      hgb_ssq = 4 * (1 - sum_gr) * (1 - sum_gr);

      // Combination
      // hgb_cross = std::clamp(hgb_cross, 0.0, 1.0);
      // hgb_ssq = std::clamp(hgb_ssq, 0.0, 1.0);
      hgb = (hgb_cross + hgb_ssq) / 2;
      hgb = std::clamp(hgb, GenericReal<is_ad>(0), GenericReal<is_ad>(1));

      // // DERIVATIVES- 0
      // for (unsigned int i = 0; i < _op_num; ++i)
      // {
      //   // Calculate sum of non i OPs squares for derivatives
      //   GenericReal<is_ad> sum_other = 0.0;
      //   for (unsigned int j = 0; j < _op_num; ++j)
      //   {
      //     if (i != j)
      //     {
      //       sum_other += (*_vals[j])[_qp] * (*_vals[j])[_qp];
      //     }
      //   }
      //   // First derivatives
      //   (*_prop_dh[i])[_qp] = 32 * (*_vals[i])[_qp] * sum_other;
      //   // Second derivatives only for same var twice (rest assuming 0)
      //   (*_prop_d2h[i][i])[_qp] = 32 * sum_other;
      // }
      //
      // // DERIVATIVES - 1
      // for (unsigned int i = 0; i < _op_num; ++i)
      // {
      //   // First derivatives
      //   (*_prop_dh[i])[_qp] = 16 * (*_vals[i])[_qp] * (sum_gr - 1);
      //   // Second derivatives only for same var twice (rest assuming 0)
      //   (*_prop_d2h[i][i])[_qp] = 16 * (sum_gr - 1 + (2 * (*_vals[i])[_qp] * (*_vals[i])[_qp]));
      // }
      break;
    }

    default:
      mooseError("Unknown func_type = ", _func_type);
      break;
  }

  // If h_val is below the threshold, set _prop_h to 0 and return early
  if (hgb < _hgb_threshold)
  {
    _prop_h[_qp] = 0.0;

    // Skip all derivative calculations as they're assumed to be zero
    // Unsure but it seems like it actually needs to explicitly be set to 0?
    for (unsigned int i = 0; i < _op_num; ++i)
    {
      (*_prop_dh[i])[_qp] = 0.0;
      (*_prop_d2h[i][i])[_qp] = 0.0;
    }
    return;
  }
  // ELSE: Define hgb and its derivatives
  _prop_h[_qp] = hgb;
}

// explicit instantiation
template class SwitchingFunctionGBMaterialTempl<true>;
template class SwitchingFunctionGBMaterialTempl<false>;
