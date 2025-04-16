//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ACInterfaceAnisoGamma.h"

registerMooseObject("PhaseFieldApp", ACInterfaceAnisoGamma);

InputParameters
ACInterfaceAnisoGamma::validParams()
{
  InputParameters params = JvarMapKernelInterface<Kernel>::validParams();
  params.addClassDescription("Gradient energy Allen-Cahn Kernel");
  params.addParam<MaterialPropertyName>("mob_name", "L", "The mobility used with the kernel");
  params.addParam<MaterialPropertyName>(
      "dgamma_dgradop_name",
      "The derivative of gamma with respect to the gradient of the variable being operated on");
  params.addParam<MaterialPropertyName>(
      "d2gamma_dgradop2_name",
      "The 2nd derivative of gamma with respect to the gradient of the variable being operated on");
  // params.addParam<MaterialPropertyName>("kappa_name", "kappa_op", "The kappa used with the
  // kernel");
  params.addParam<bool>("variable_L",
                        true,
                        "The mobility is a function of any MOOSE variable (if "
                        "this is set to false L must be constant over the "
                        "entire domain!)");
  params.addRequiredCoupledVar("v",
                               "Array of coupled order parameter names for OTHER order parameters");
  return params;
}

ACInterfaceAnisoGamma::ACInterfaceAnisoGamma(const InputParameters & parameters)
  : DerivativeMaterialInterface<JvarMapKernelInterface<Kernel>>(parameters),
    _L(getMaterialProperty<Real>("mob_name")),
    _gamma(getMaterialProperty<Real>("gamma_aniso")), // gamma_asymm is the normal name change later
    _dgammadgrad_op(getMaterialProperty<RealGradient>("dgamma_dgradop_name")),
    _d2gammadgrad_op2(getMaterialProperty<RealTensorValue>("d2gamma_dgradop2_name")),
    // Grain OP list as per ACGrGrPoly
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _vals_var(coupledIndices("v")),
    _variable_L(getParam<bool>("variable_L")),
    _dLdop(getMaterialPropertyDerivative<Real>("mob_name", _var.name())),
    _d2Ldop2(getMaterialPropertyDerivative<Real>("mob_name", _var.name(), _var.name())),
    // _dkappadop(getMaterialPropertyDerivative<Real>("kappa_name", _var.name())),
    _dLdarg(_n_args),
    _d2Ldargdop(_n_args),
    _d2Ldarg2(_n_args),
    // _dkappadarg(_n_args),
    // _arg(_n_args),
    _gradarg(_n_args)
{
  // Get mobility and kappa derivatives and coupled variable gradients
  for (unsigned int i = 0; i < _n_args; ++i)
  {
    MooseVariable * ivar = _coupled_standard_moose_vars[i];
    const VariableName iname = ivar->name();
    if (iname == _var.name())
    {
      if (isCoupled("args"))
        paramError("args",
                   "The kernel variable should not be specified in the coupled `args` parameter.");
      else
        paramError("coupled_variables",
                   "The kernel variable should not be specified in the coupled `coupled_variables` "
                   "parameter.");
    }

    _dLdarg[i] = &getMaterialPropertyDerivative<Real>("mob_name", i);
    // _dkappadarg[i] = &getMaterialPropertyDerivative<Real>("kappa_name", i);
    _d2Ldargdop[i] = &getMaterialPropertyDerivative<Real>("mob_name", iname, _var.name());

    // _arg[i] = ivar;
    _gradarg[i] = &(ivar->gradSln());

    _d2Ldarg2[i].resize(_n_args);
    for (unsigned int j = 0; j < _n_args; ++j)
      _d2Ldarg2[i][j] = &getMaterialPropertyDerivative<Real>("mob_name", i, j);
  }
}

void
ACInterfaceAnisoGamma::initialSetup()
{
  validateCoupling<Real>("mob_name");
  // validateCoupling<Real>("kappa_name");
}

RealGradient
ACInterfaceAnisoGamma::gradL() // no changes yet
{
  RealGradient g = _grad_u[_qp] * _dLdop[_qp];
  for (unsigned int i = 0; i < _n_args; ++i)
    g += (*_gradarg[i])[_qp] * (*_dLdarg[i])[_qp];
  return g;
}

RealGradient
ACInterfaceAnisoGamma::nablaLPsi() // RH $L \nabla \psi$
{
  // sum is the product rule gradient \f$ \nabla (L\psi) \f$
  RealGradient sum = _L[_qp] * _grad_test[_i][_qp];

  if (_variable_L)
    sum += gradL() * _test[_i][_qp];

  return sum;
}

// RealGradient
// ACInterfaceAnisoGamma::kappaNablaLPsi()
// {
//   return _kappa[_qp] * nablaLPsi();
// }

Real
ACInterfaceAnisoGamma::sumSqEtaj()
{
  // Sum all other (grain?) order parameters
  // Only works assuming the grains are the only args/coupled variables
  // SumOPj += (*_coupled_standard_moose_vars[i])[_qp] * (*_coupled_standard_moose_vars[i])[_qp];
  Real SumOPj = 0.0;
  for (unsigned int i = 0; i < _op_num; ++i)
    SumOPj += (*_vals[i])[_qp] * (*_vals[i])[_qp];
  return SumOPj;
}

Real
ACInterfaceAnisoGamma::computeQpResidual()
{
  return _u[_qp] * _u[_qp] * _dgammadgrad_op[_qp] * sumSqEtaj() * nablaLPsi();
}

Real
ACInterfaceAnisoGamma::computeQpJacobian()
{
  // dsum is the derivative of R wrt u (without L variable dependence for now)
  RealGradient ddir = 2 * _u[_qp] * _L[_qp] * _dgammadgrad_op[_qp] * sumSqEtaj() * _phi[_j][_qp];

  // dind is the derivative of R wrt \nabla u (without L variable dependence for now)
  RealTensorValue dind =
      _u[_qp] * _u[_qp] * _L[_qp] * _d2gammadgrad_op2[_qp] * sumSqEtaj(); // * _grad_phi[_j][_qp];

  RealGradient dind_cont = dind * _grad_phi[_j][_qp]; // Tensor*vector = contracted to vector

  // return dot(ddir, _grad_test[_i][_qp]) + dot(dind, _grad_test[_i][_qp]);
  Real jac1 = ddir * _grad_test[_i][_qp]; // Vector*Vector = Real
  Real jac2 = dind_cont * _grad_test[_i][_qp];
  return jac1 + jac2;

  // // dsum is the derivative \f$ \frac\partial{\partial \eta} \left( \nabla (L\psi) \right) \f$
  // RealGradient dsum =
  //     (_dkappadop[_qp] * _L[_qp] + _kappa[_qp] * _dLdop[_qp]) * _phi[_j][_qp] *
  //     _grad_test[_i][_qp];

  // // compute the derivative of the gradient of the mobility
  // if (_variable_L)
  // {
  //   RealGradient dgradL =
  //       _grad_phi[_j][_qp] * _dLdop[_qp] + _grad_u[_qp] * _phi[_j][_qp] * _d2Ldop2[_qp];

  //   for (unsigned int i = 0; i < _n_args; ++i)
  //     dgradL += (*_gradarg[i])[_qp] * _phi[_j][_qp] * (*_d2Ldargdop[i])[_qp];

  //   dsum += (_kappa[_qp] * dgradL + _dkappadop[_qp] * _phi[_j][_qp] * gradL()) * _test[_i][_qp];
  // }

  // return _grad_phi[_j][_qp] * kappaNablaLPsi() + _grad_u[_qp] * dsum;
}

Real
ACInterfaceAnisoGamma::computeQpOffDiagJacobian(unsigned int jvar)
{
  // // get the coupled variable jvar is referring to
  // const unsigned int cvar = mapJvarToCvar(jvar);

  // // dsum is the derivative \f$ \frac\partial{\partial \eta} \left( \nabla (L\psi) \right) \f$
  // RealGradient dsum = ((*_dkappadarg[cvar])[_qp] * _L[_qp] + _kappa[_qp] * (*_dLdarg[cvar])[_qp])
  // * _phi[_j][_qp] * _grad_test[_i][_qp];

  // // compute the derivative of the gradient of the mobility
  // if (_variable_L)
  // {
  //   RealGradient dgradL = _grad_phi[_j][_qp] * (*_dLdarg[cvar])[_qp] +
  //                         _grad_u[_qp] * _phi[_j][_qp] * (*_d2Ldargdop[cvar])[_qp];

  //   for (unsigned int i = 0; i < _n_args; ++i)
  //     dgradL += (*_gradarg[i])[_qp] * _phi[_j][_qp] * (*_d2Ldarg2[cvar][i])[_qp];

  //   dsum += (_kappa[_qp] * dgradL + _dkappadop[_qp] * _phi[_j][_qp] * gradL()) * _test[_i][_qp];
  // }

  // return _grad_u[_qp] * dsum;

  // get the coupled variable jvar is referring to
  // const unsigned int cvar = mapJvarToCvar(jvar);

  // Temp version with only using etas
  for (unsigned int i = 0; i < _op_num; ++i)
    if (jvar == _vals_var[i])
    {

      // dsum is the derivative of R wrt eta_j (without L variable dependence for now)
      RealGradient ddir =
          2 * (*_vals[i])[_qp] * _u[_qp] * _u[_qp] * _L[_qp] * _dgammadgrad_op[_qp] * _phi[_j][_qp];

      // dind is the derivative of R wrt \nabla eta_j (without L variable dependence for now)
      RealTensorValue dind = _u[_qp] * _u[_qp] * _L[_qp] * (-_d2gammadgrad_op2[_qp]) * sumSqEtaj();

      RealGradient dind_cont = dind * _grad_phi[_j][_qp]; // Tensor*vector = contracted to vector

      Real jac1 = ddir * _grad_test[_i][_qp]; // Vector*Vector = Real
      Real jac2 = dind_cont * _grad_test[_i][_qp];
      return jac1 + jac2;
    }

  return 0.0;
}
