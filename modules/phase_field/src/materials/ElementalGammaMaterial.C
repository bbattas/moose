//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElementalGammaMaterial.h"
#include "libmesh/quadrature.h"

registerMooseObject("PhaseFieldApp", ElementalGammaMaterial);

InputParameters
ElementalGammaMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Generates a diffusion tensor to distinguish between the bulk, grain "
                             "boundaries, and surfaces");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  params.addParam<MaterialPropertyName>(
      "gb_energy_iso_name", "sigma_iso", "Isotropic GB energy before inclination dependence.");
  params.addParam<MaterialPropertyName>(
      "kappa", "kappa", "Gradient energy constant kappa material name.");
  params.addParam<Real>(
      "free_energy_m", 1, "Free energy function constant m (or mu in PF kernels).");
  return params;
}

ElementalGammaMaterial::ElementalGammaMaterial(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _gamma_in(getMaterialProperty<Real>("gamma_a")),
    _gamma_out(declareProperty<Real>("gamma_asymm")),
    _no_ij_pairs(getMaterialProperty<bool>("no_ij_pairs")),
    _elem_no_ij(declareProperty<bool>("elem_no_ij")),
    _int_noij(declareProperty<Real>("int_noij")),
    _gbe_iso(getMaterialProperty<Real>("gb_energy_iso_name")),
    _kappa(getMaterialProperty<Real>(getParam<MaterialPropertyName>("kappa"))),
    _const_m(getParam<Real>("free_energy_m"))
{
  // if (_op_num == 0)
  //   mooseError("Model requires op_num > 0");

  // _vals.resize(_op_num);
  // _grad_vals.resize(_op_num);
  // for (unsigned int i = 0; i < _op_num; ++i)
  // {
  //   _vals[i] = &coupledValue("v", i);
  //   _grad_vals[i] = &coupledGradient("v", i);
  //   _vals_name[i] = coupledName("v", i);
  //   if (!isCoupledConstant(_vals_name[i]))
  //     _dDdeta[i] = &declarePropertyDerivative<RealTensorValue>(_diffusivity_name, _vals_name[i]);
  // }
}

void
ElementalGammaMaterial::computeProperties()
{
  // Hard-coded coefficients (for poly gamma)
  constexpr Real a1 = -3.0944; // coefficient for g2^4
  constexpr Real a2 = -1.8169; // coefficient for g2^3
  constexpr Real a3 = 10.323;  // coefficient for g2^2
  constexpr Real a4 = -8.1819; // coefficient for g2
  constexpr Real a5 = 2.0033;  // constant term

  bool elem_skip = false;
  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    const bool qp_noij = _no_ij_pairs[_qp];
    if (qp_noij)
    {
      elem_skip = true;
      break;
    }
  }

  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    _elem_no_ij[_qp] = elem_skip;
    if (elem_skip)
    {
      _int_noij[_qp] = 1;
      // Gamma
      Real g = _gbe_iso[_qp] / (std::sqrt(_kappa[_qp] * _const_m));
      Real g2 = g * g;
      Real pg = (((a1 * g2 + a2) * g2 + a3) * g2 + a4) * g2 + a5;
      _gamma_out[_qp] = 1 / pg;
    }
    else
    {
      _int_noij[_qp] = 0;
      _gamma_out[_qp] = _gamma_in[_qp];
    }
  }
}
