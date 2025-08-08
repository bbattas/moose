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
  params.addParam<Real>("hgb_threshold", 0.0, "Lower limit cutoff for hgb.");
  // params.addRequiredCoupledVar("all_ops", "Vector of all order parameters for all phases that you
  // want derivatives wrt");
  MooseEnum func_type("base=0 moelans=1", "base");
  params.addParam<MooseEnum>("func_type", func_type, "Which hgb option to use.");
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
    _hgb_threshold(getParam<Real>("hgb_threshold")),
    _func_type(getParam<MooseEnum>("func_type")),
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
  // OTHERWISE we actually calculate
  GenericReal<is_ad> hgb = 0.0;
  // GenericReal<is_ad> hgb_val = 0.0;

  switch (_func_type)
  {
    case 0:
    {
      // Previous use: $16 \sum \eta_i^2 \eta_j^2$
      // HGB
      for (unsigned int i = 0; i < _num_eta_gb; ++i)
        for (unsigned int j = i + 1; j < _num_eta_gb; ++j)
        {
          hgb += (*_eta_gb[i])[_qp] * (*_eta_gb[i])[_qp] * (*_eta_gb[j])[_qp] * (*_eta_gb[j])[_qp];
        }
      hgb *= 16;

      // DERIVATIVES
      for (unsigned int i = 0; i < _num_eta_gb; ++i)
      {
        // Calculate sum of non i OPs squares for derivatives
        GenericReal<is_ad> sum_other = 0.0;
        for (unsigned int j = 0; j < _num_eta_gb; ++j)
        {
          if (i != j)
          {
            sum_other += (*_eta_gb[j])[_qp] * (*_eta_gb[j])[_qp];
          }
        }
        // First derivatives
        (*_prop_dh[i])[_qp] = 32 * (*_eta_gb[i])[_qp] * sum_other;
        // Second derivatives only for same var twice (rest assuming 0)
        (*_prop_d2h[i][i])[_qp] = 32 * sum_other;
      }
      break;
    }

    case 1:
    {
      // Using the MOELANS version $\eta_i^2 \eta_j^2 / \sum \eta_i^2 \eta_j^2$
      GenericReal<is_ad> numer = 0.0;
      GenericReal<is_ad> denom = 0.0;
      for (unsigned int i = 0; i < _num_eta_gb; ++i)
        for (unsigned int j = i + 1; j < _num_eta_gb; ++j)
        {
          numer +=
              (*_eta_gb[i])[_qp] * (*_eta_gb[i])[_qp] * (*_eta_gb[j])[_qp] * (*_eta_gb[j])[_qp];
          denom +=
              (*_eta_gb[i])[_qp] * (*_eta_gb[i])[_qp] * (*_eta_gb[j])[_qp] * (*_eta_gb[j])[_qp];
        }
      // numer = (*_eta_gb[0])[_qp] * (*_eta_gb[0])[_qp] * (*_eta_gb[1])[_qp] * (*_eta_gb[1])[_qp];
      if (denom < 1e-8)
        hgb = 0.0;
      else
        hgb = numer / denom;
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
    for (unsigned int i = 0; i < _num_eta_gb; ++i)
    {
      (*_prop_dh[i])[_qp] = 0.0;
      (*_prop_d2h[i][i])[_qp] = 0.0;
      // for (unsigned int j = 0; j < _num_eta_gb; ++j)
      // {
      //   (*_prop_d2h[i][j])[_qp] = 0.0;
      // }
    }
    return;
  }
  // ELSE: Define hgb and its derivatives
  _prop_h[_qp] = hgb;

  // // // For derivatives: sum of other squares
  // // // GenericReal<is_ad> sum_other = 0.0;
  // // std::vector<GenericReal<is_ad>> sum_other(_num_eta_gb, 0.0);
  // // for (unsigned int i = 0; i < _num_eta_gb; ++i)
  // // {
  // //   GenericReal<is_ad> sum_tmp = 0.0;
  // //   for (unsigned int j = 0; j < _num_eta_gb; ++j)
  // //   {
  // //     if (i != j)
  // //     {
  // //       sum_tmp += (*_eta_gb[j])[_qp] * (*_eta_gb[j])[_qp];
  // //     }
  // //   }
  // //   sum_other[i] = sum_tmp;
  // // }

  // for (unsigned int i = 0; i < _num_eta_gb; ++i)
  // {
  //   // Calculate sum of non i OPs squares for derivatives
  //   GenericReal<is_ad> sum_other = 0.0;
  //   for (unsigned int j = 0; j < _num_eta_gb; ++j)
  //   {
  //     if (i != j)
  //     {
  //       sum_other += (*_eta_gb[j])[_qp] * (*_eta_gb[j])[_qp];
  //     }
  //   }
  //   // First derivatives
  //   (*_prop_dh[i])[_qp] = 32 * (*_eta_gb[i])[_qp] * sum_other;
  //   // Second derivatives only for same var twice (rest assuming 0)
  //   (*_prop_d2h[i][i])[_qp] = 32 * sum_other;

  //   // // Manually set second derivatives to zero for i != j and assign the diagonal term
  //   // for (unsigned int j = 0; j < _num_eta_gb; ++j)
  //   // {
  //   //   if (i == j)
  //   //   {
  //   //     (*_prop_d2h[i][j])[_qp] = 32 * sum_other;
  //   //   }
  //   //   else
  //   //   {
  //   //     (*_prop_d2h[i][j])[_qp] = 0.0; // Previously assumed it defaulted to 0 if undefined?
  //   //   }
  //   // }
  // }

  // // for (unsigned int i = 0; i < _num_eta_gb; ++i)
  // // {
  // //   // First derivatives
  // //   (*_prop_dh[i])[_qp] = 32 * (*_eta_gb[i])[_qp] * sum_other[i];

  // //   // Second derivatives
  // //   for (unsigned int j = 0; j < _num_eta_gb; ++j)
  // //   {
  // //     if (i == j)
  // //     {
  // //       (*_prop_d2h[i][j])[_qp] = 32 * sum_other[i];
  // //     }
  // //     else // Removed the cross term second derivatives
  // //     {
  // //       (*_prop_d2h[i][j])[_qp] = 0.0;//64 * (*_eta_gb[i])[_qp] * (*_eta_gb[j])[_qp];
  // //     }
  // //   }
  // // }
}

// explicit instantiation
template class SwitchingFunctionGBMaterialTempl<true>;
template class SwitchingFunctionGBMaterialTempl<false>;
