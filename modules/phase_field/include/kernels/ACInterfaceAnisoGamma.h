#pragma once

#include "Kernel.h"
#include "JvarMapInterface.h"
#include "DerivativeMaterialInterface.h"

/**
 * Compute the Allen-Cahn interface term with the weak form residual
 * \f$ \left( \kappa_i \nabla\eta_i, \nabla (L_i \psi) \right) \f$
 */
class ACInterfaceAnisoGamma : public DerivativeMaterialInterface<JvarMapKernelInterface<Kernel>>
{
public:
  static InputParameters validParams();

  ACInterfaceAnisoGamma(const InputParameters & parameters);
  virtual void initialSetup();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  RealGradient gradL();
  // RealGradient gradKappa();

  /// the \f$ \nabla(L\psi) \f$ term
  RealGradient nablaLPsi();

  Real sumSqEtaj();

  // /// the \f$ \kappa\nabla(L\psi) \f$ term
  // RealGradient kappaNablaLPsi();

  ///@{ Variables for second order derivatives
  const VariableSecond & _second_u;
  const VariableTestSecond & _second_test;
  const VariablePhiSecond & _second_phi;
  ///@}

  /// Mobility
  const MaterialProperty<Real> & _L;

  /// FE Constant m/mu
  const MaterialProperty<Real> & _mu;

  // Gamma
  const MaterialProperty<Real> & _gamma;
  const MaterialProperty<RealGradient> & _dgammadgrad_op;
  const MaterialProperty<RealTensorValue> & _d2gammadgrad_op2;

  const bool _skip_off;
  // Adding a mask for GB only and to smooth it out a bit
  const MaterialProperty<Real> * _mask;
  const bool _mask_tf;

  std::vector<unsigned int> _grain_ids;
  // std::unordered_set<unsigned int> grain_set;

  /// Grain op values from input v = (all other grain ops)
  // const unsigned int _op_num;
  // const std::vector<const VariableValue *> _vals;
  // const std::vector<unsigned int> _vals_var;
  const std::string _var_name_base;

  // /// Interfacial parameter
  // const MaterialProperty<Real> & _kappa;

  /// flag set if L is a function of non-linear variables in args
  const bool _variable_L;

  /// @{ Mobility derivatives w.r.t. order parameter
  const MaterialProperty<Real> * _dLdop;
  const MaterialProperty<Real> * _d2Ldop2;
  /// @}

  const MaterialProperty<RealGradient> * _dLdgrad_op;
  const MaterialProperty<RealTensorValue> * _d2Ldgrad_op2;
  const MaterialProperty<RealGradient> * _d2Ldopdgrad_op;

  // /// kappa derivative w.r.t. order parameter
  // const MaterialProperty<Real> & _dkappadop;

  /// @{ Mobility derivative w.r.t. other coupled variables
  std::vector<const MaterialProperty<Real> *> _dLdarg;
  std::vector<const MaterialProperty<Real> *> _d2Ldargdop;
  std::vector<std::vector<const MaterialProperty<Real> *>> _d2Ldarg2;
  /// @}
  std::vector<const MaterialProperty<RealGradient> *> _dLdgradarg;
  std::vector<const MaterialProperty<RealGradient> *> _d2Ldgradargdop;
  std::vector<const MaterialProperty<RealGradient> *> _d2Ldargdgradop;
  std::vector<const MaterialProperty<RealTensorValue> *> _d2Ldgradarg2;
  std::vector<std::vector<const MaterialProperty<RealGradient> *>> _d2Ldgradargdarg;

  // /// kappa derivative w.r.t. other coupled variables
  // std::vector<const MaterialProperty<Real> *> _dkappadarg;
  /// Arg values for etsj^2
  // std::vector<const VariableValue *> _arg;
  // const std::vector<const VariableValue *> _arg;

  std::vector<const VariableValue *> _eta_vals;
  /// Gradients for all coupled variables
  std::vector<const VariableGradient *> _gradarg;
  std::vector<const VariableSecond *> _second_arg;
};
