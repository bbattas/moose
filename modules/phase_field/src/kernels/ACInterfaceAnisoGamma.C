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
      "gamma_name", "gamma_aniso", "The gamma used with the kernel");
  params.addParam<MaterialPropertyName>(
      "dgamma_dgradop_name",
      "The derivative of gamma with respect to the gradient of the variable being operated on");
  params.addParam<MaterialPropertyName>(
      "d2gamma_dgradop2_name",
      "The 2nd derivative of gamma with respect to the gradient of the variable being operated on");
  params.addParam<MaterialPropertyName>(
      "dL_dgradop_name",
      "The derivative of L with respect to the gradient of the variable being operated on");
  params.addParam<MaterialPropertyName>(
      "d2L_dgradop2_name",
      "The 2nd derivative of L with respect to the gradient of the variable being operated on");
  // params.addParam<MaterialPropertyName>("kappa_name", "kappa_op", "The kappa used with the
  // kernel");
  params.addParam<bool>("variable_L",
                        true,
                        "The mobility is a function of any MOOSE variable (if "
                        "this is set to false L must be constant over the "
                        "entire domain!)");
  params.addRequiredCoupledVar("v",
                               "Array of coupled order parameter names for OTHER order parameters");
  // params.addRequiredCoupledVar("var_name_base", "Array of coupled variable names");
  params.addRequiredParam<std::string>("var_name_base", "specifies the base name of the variables");
  return params;
}

ACInterfaceAnisoGamma::ACInterfaceAnisoGamma(const InputParameters & parameters)
  : DerivativeMaterialInterface<JvarMapKernelInterface<Kernel>>(parameters),
    // Second order values needed for $L=f(\nabla u)$
    _second_u(second()),
    _second_test(secondTest()),
    _second_phi(secondPhi()),
    _L(getMaterialProperty<Real>("mob_name")),
    _gamma(getMaterialProperty<Real>(getParam<MaterialPropertyName>("gamma_name"))),
    _dgammadgrad_op(getMaterialProperty<RealGradient>("dgamma_dgradop_name")),
    _d2gammadgrad_op2(getMaterialProperty<RealTensorValue>("d2gamma_dgradop2_name")),
    // Grain OP list as per ACGrGrPoly
    _grain_ids(),
    // _grain_set(),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _vals_var(coupledIndices("v")),
    _var_name_base(getParam<std::string>("var_name_base")),
    _variable_L(getParam<bool>("variable_L")),
    _dLdop(getMaterialPropertyDerivative<Real>("mob_name", _var.name())),
    _d2Ldop2(getMaterialPropertyDerivative<Real>("mob_name", _var.name(), _var.name())),
    _dLdgrad_op(getMaterialProperty<RealGradient>("dL_dgradop_name")),
    _d2Ldgrad_op2(getMaterialProperty<RealTensorValue>("d2L_dgradop2_name")),
    _dLdopdgrad_op(getMaterialPropertyDerivative<RealGradient>("dL_dgradop_name", _var.name())),
    _dLdarg(_n_args),
    _d2Ldargdop(_n_args),
    _d2Ldarg2(_n_args),
    _dLdgradarg(_n_args),
    _d2Ldgradarg2(_n_args),
    _d2Ldargdgradarg(_n_args),
    // _dkappadarg(_n_args),
    // _arg(_n_args),
    _gradarg(_n_args),
    _second_arg(_n_args)
{
  // Get mobility and kappa derivatives and coupled variable gradients
  mooseWarning("NAMEBASE: ", _var_name_base);
  // std::vector<unsigned int> _grain_ids;
  for (unsigned int i = 0; i < _n_args; ++i)
  {
    mooseWarning("Into initial");
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

    _gradarg[i] = &(ivar->gradSln());
    // _gradarg[i] = &coupledGradient(ivar->name());
    _second_arg[i] = &(ivar->secondSln());

    _d2Ldarg2[i].resize(_n_args);
    for (unsigned int j = 0; j < _n_args; ++j)
      _d2Ldarg2[i][j] = &getMaterialPropertyDerivative<Real>("mob_name", i, j);

    // ID the number at the end of the args to determine which dLdgrad to use
    // Check if coupled var starts with the var_name_base (grain op)
    if (MooseUtils::beginsWith(iname, _var_name_base) && iname != _var.name())
    {
      _grain_ids.push_back(i);
      // _grain_set.insert(i);
      // Pull the number at the end of the grain op
      const std::size_t pos = iname.find_last_not_of("0123456789");
      if (pos == std::string::npos || pos + 1 == iname.size())
        mooseError("Variable '",
                   iname,
                   "' does not end in a numeric suffix, despite starting with var_name_base.");
      const std::string digits = iname.substr(pos + 1);
      // mooseWarning("digits = ", digits);
      _dLdgradarg[i] = &getMaterialProperty<RealGradient>("dLdgrad_eta_" + digits);
      _d2Ldgradarg2[i] = &getMaterialProperty<RealTensorValue>("d2Ldgrad_eta2_" + digits);
      _d2Ldargdgradarg[i] =
          &getMaterialPropertyDerivative<RealGradient>("dLdgrad_eta_" + digits, i);
    }
    // else
    // {
    //   // Placeholder for now
    //   _dLdgradarg[i] = 0.0;
    //   _d2Ldgradarg2[i] = 0.0;
    // }

    // mooseWarning("debug: _n_args = ", _coupled_standard_moose_vars[i]->name());
    // mooseWarning("debug: number = ", "dLdgrad_eta_" + digits);
  }
  // mooseWarning("debug: _n_args = ", _n_args);
}

void
ACInterfaceAnisoGamma::initialSetup()
{
  validateCoupling<Real>("mob_name");
  // validateCoupling<Real>("kappa_name");
}

RealGradient
ACInterfaceAnisoGamma::gradL() // Includes grad op dependence
{
  RealGradient g = _grad_u[_qp] * _dLdop[_qp];
  g += _second_u[_qp] * _dLdgrad_op[_qp];
  for (unsigned int i = 0; i < _n_args; ++i)
  {
    g += (*_gradarg[i])[_qp] * (*_dLdarg[i])[_qp];
    g += (*_second_arg[i])[_qp] * (*_dLdgradarg[i])[_qp];
  }
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

Real
ACInterfaceAnisoGamma::sumSqEtaj()
{
  // Sum all other (grain?) order parameters
  Real SumOPj = 0.0;
  for (unsigned int j : _grain_ids)
    SumOPj += (*_coupled_standard_moose_vars[j]).sln()[_qp] *
              (*_coupled_standard_moose_vars[j]).sln()[_qp];
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
  // dsum is the derivative of R wrt u
  Real ddir = 2 * _u[_qp] * _dgammadgrad_op[_qp] * sumSqEtaj() * _phi[_j][_qp] * nablaLPsi();

  // dind is the derivative of R wrt \nabla u
  Real dind =
      _u[_qp] * _u[_qp] * _d2gammadgrad_op2[_qp] * sumSqEtaj() * _grad_phi[_j][_qp] * nablaLPsi();

  if (_variable_L)
  {
    // Grad L partials
    static const RealTensorValue I(1, 0, 0, 0, 1, 0, 0, 0, 1);
    RealGradient dgradLdu = _d2Ldop2[_qp] * _grad_u[_qp] + _dLdopdgrad_op[_qp] * _second_u[_qp];
    RealTensorValue dgradLdgradu = libMesh::outer_product(_dLdopdgrad_op[_qp], _grad_u[_qp]) +
                                   I * _dLdop[_qp] + _d2Ldgrad_op2[_qp] * _second_u[_qp];
    // Direct L dependence
    ddir += _dLdop[_qp] * _u[_qp] * _u[_qp] * _dgammadgrad_op[_qp] * sumSqEtaj() * _phi[_j][_qp] *
            _grad_test[_i][_qp];
    // Direct grad L dependence
    ddir += dgradLdu * _u[_qp] * _u[_qp] * _dgammadgrad_op[_qp] * sumSqEtaj() * _phi[_j][_qp] *
            _test[_i][_qp];
    // Indirect L dependence (of grad u)
    dind += _dLdgrad_op[_qp] * _u[_qp] * _u[_qp] * _dgammadgrad_op[_qp] * sumSqEtaj() *
            _grad_phi[_j][_qp] * _grad_test[_i][_qp];
    // Indirect grad L dependence (of grad u)
    dind += dgradLdgradu * _u[_qp] * _u[_qp] * _dgammadgrad_op[_qp] * sumSqEtaj() *
            _grad_phi[_j][_qp] * _test[_i][_qp];
  }

  return ddir + dind;
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
  const unsigned int cvar = mapJvarToCvar(jvar);

  // Check if the offdiag wrt variable is one of the grain ops
  bool is_grain = std::find(_grain_ids.begin(), _grain_ids.end(), cvar) != _grain_ids.end();
  if (is_grain)
  {
    // Direct dependence on arg
    Real ddir = 2 * (*_coupled_standard_moose_vars[cvar]).sln()[_qp] * _u[_qp] * _u[_qp] *
                _dgammadgrad_op[_qp] * _phi[_j][_qp] * nablaLPsi();
    // Indirect dependence on grad_arg
    Real dind = _u[_qp] * _u[_qp] * (-_d2gammadgrad_op2[_qp]) * sumSqEtaj() * _grad_phi[_j][_qp] *
                nablaLPsi();

    // if (_variable_L)
    // {
    //   // Need for loops like in ACInt for example
    //   // Grad L partials
    //   static const RealTensorValue I(1, 0, 0, 0, 1, 0, 0, 0, 1);
    //   RealGradient dgradLdeta =
    //       (*_d2Ldarg2[cvar])[_qp] * _grad_u[_qp] + _dLdopdgrad_op[_qp] * _second_u[_qp];
    //   RealTensorValue dgradLdgradeta = libMesh::outer_product(_dLdopdgrad_op[_qp], _grad_u[_qp])
    //   +
    //                                    I * _dLdop[_qp] + _d2Ldgrad_op2[_qp] * _second_u[_qp];
    //   // Direct L dependence
    //   ddir += _dLdop[_qp] * _u[_qp] * _u[_qp] * _dgammadgrad_op[_qp] * sumSqEtaj() *
    //   _phi[_j][_qp] *
    //           _grad_test[_i][_qp];
    //   // Direct grad L dependence
    //   ddir += dgradLdu * _u[_qp] * _u[_qp] * _dgammadgrad_op[_qp] * sumSqEtaj() * _phi[_j][_qp] *
    //           _test[_i][_qp];
    //   // Indirect L dependence (of grad u)
    //   dind += _dLdgrad_op[_qp] * _u[_qp] * _u[_qp] * _dgammadgrad_op[_qp] * sumSqEtaj() *
    //           _grad_phi[_j][_qp] * _grad_test[_i][_qp];
    //   // Indirect grad L dependence (of grad u)
    //   dind += dgradLdgradu * _u[_qp] * _u[_qp] * _dgammadgrad_op[_qp] * sumSqEtaj() *
    //           _grad_phi[_j][_qp] * _test[_i][_qp];
    // }
    return ddir + dind;
  }

  // // Temp version with only using etas
  // for (unsigned int i = 0; i < _op_num; ++i)
  //   if (jvar == _vals_var[i])
  //   {
  //     // dsum is the derivative of R wrt eta_j
  //     Real ddir = 2 * (*_vals[i])[_qp] * _u[_qp] * _u[_qp] * _dgammadgrad_op[_qp] * _phi[_j][_qp]
  //     *
  //                 nablaLPsi();

  //     // dind is the derivative of R wrt \nabla eta_j (without L variable dependence for now)
  //     Real dind = _u[_qp] * _u[_qp] * (-_d2gammadgrad_op2[_qp]) * sumSqEtaj() *
  //     _grad_phi[_j][_qp] *
  //                 nablaLPsi();

  //     // RealGradient dind_cont = dind * _grad_phi[_j][_qp]; // Tensor*vector = contracted to vector

  //     // Real jac1 = ddir * _grad_test[_i][_qp]; // Vector*Vector = Real
  //     // Real jac2 = dind_cont * _grad_test[_i][_qp];
  //     return ddir + dind;
  //     // // dsum is the derivative of R wrt eta_j
  //     // RealGradient ddir =
  //     //     2 * (*_vals[i])[_qp] * _u[_qp] * _u[_qp] * _L[_qp] * _dgammadgrad_op[_qp] *
  //     //     _phi[_j][_qp];

  //     // // dind is the derivative of R wrt \nabla eta_j (without L variable dependence for now)
  //     // RealTensorValue dind = _u[_qp] * _u[_qp] * _L[_qp] * (-_d2gammadgrad_op2[_qp]) *
  //     // sumSqEtaj();

  //     // RealGradient dind_cont = dind * _grad_phi[_j][_qp]; // Tensor*vector = contracted to vector

  //     // Real jac1 = ddir * _grad_test[_i][_qp]; // Vector*Vector = Real
  //     // Real jac2 = dind_cont * _grad_test[_i][_qp];
  //     // return jac1 + jac2;
  //   }

  return 0.0;
}
