#pragma once

#include "Kernel.h"
#include "JvarMapInterface.h"
#include "DerivativeMaterialInterface.h"
#include "OpIndexUtils.h"

/**
 * Compute the Allen-Cahn interface term with the weak form residual
 * \f$ \left( \kappa_i \nabla\eta_i, \nabla (L_i \psi) \right) \f$
 */
class ACInterfaceInclinationGamma
  : public DerivativeMaterialInterface<JvarMapKernelInterface<Kernel>>
{
public:
  static InputParameters validParams();

  ACInterfaceInclinationGamma(const InputParameters & parameters);
  virtual void initialSetup();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  RealGradient gradL();

  /// the \f$ \nabla(L\psi) \f$ term
  RealGradient nablaLPsi();

  ///@{ Variables for second order derivatives
  const VariableSecond & _second_u;
  const VariableTestSecond & _second_test;
  const VariablePhiSecond & _second_phi;
  ///@}

  const unsigned int _op_num;

  const std::string _var_name_base;

  const MaterialProperty<bool> & _no_ij_pairs;
  const MaterialProperty<std::vector<unsigned int>> & _ij_i;
  const MaterialProperty<std::vector<unsigned int>> & _ij_j;

  /// FE Constant m/mu
  const MaterialProperty<Real> & _mu;

  // Gamma
  const MaterialProperty<std::vector<RealGradient>> & _dgammaij_dgradeta;
  const MaterialProperty<std::vector<RealTensorValue>> & _d2gammaij_dgradeta2;

  /// Mobility
  const MaterialProperty<Real> & _L;
  // /// Interfacial parameter
  // const MaterialProperty<Real> & _kappa;

  /// flag set if L is a function of non-linear variables in args
  const bool _variable_L;

  const bool _skip_off;

  const bool _debug_kernel;

  // Adding a mask for GB only and to smooth it out a bit
  const MaterialProperty<Real> * _mask;
  const bool _mask_tf;

  /// Variable L
  /// wrt u
  const MaterialProperty<Real> * _dLdu;
  const MaterialProperty<Real> * _d2Ldu2;
  // gradient dependent
  std::vector<const MaterialProperty<RealGradient> *> _dL_dgradeta;
  std::vector<const MaterialProperty<RealGradient> *> _d2L_dgradudarg;
  const MaterialProperty<std::vector<RealTensorValue>> * _d2L_dgradudgradeta;
  std::vector<std::vector<const MaterialProperty<RealGradient> *>> _d2L_dgradetadarg;
  // // wrt grad etas (just need the u componenets though)
  // const MaterialProperty<std::vector<RealGradient>> * _dLdgrad_eta;
  // const MaterialProperty<std::vector<RealTensorValue>> * _d2Ldgrad_eta2;
  // const MaterialProperty<std::vector<RealGradient>> * _d2Ldudgrad_eta;
  // /// wrt n_args
  std::vector<const MaterialProperty<Real> *> _dLdarg;
  std::vector<const MaterialProperty<Real> *> _d2Ldargdu;
  std::vector<std::vector<const MaterialProperty<Real> *>> _d2Ldarg2;
  // std::vector<const MaterialProperty<std::vector<RealGradient>> *> _d2Ldargdgrad_eta;

  /// Gradients for all coupled variables
  std::vector<const VariableGradient *> _gradarg;
  std::vector<const VariableSecond *> _second_arg;

  // New adds
  // My order-parameter index for u = eta_i
  unsigned _my_op = std::numeric_limits<unsigned>::max();

  // For each OP index k, pointer to that variable's value array at qps
  std::vector<const VariableValue *> _eta_by_op;

  // Map OP index -> jvar id (so OffDiagJacobians can match quickly)
  std::vector<int> _op_to_jvar;

  // Convenience: get eta_k and grad phi_j at this qp
  inline Real eta_at(unsigned k) const { return (*_eta_by_op[k])[_qp]; }

  // std::unordered_map<unsigned int, const VariableValue *> _grain_val_by_k;    // key = k
  // std::unordered_map<unsigned int, const VariableValue *> _grain_val_by_argi; // key = i

  // // Inline, const, and noexcept (optional).
  // inline const VariableValue * get_val_by_k(unsigned int j) const noexcept
  // {
  //   if (auto it = _grain_val_by_k.find(j); it != _grain_val_by_k.end())
  //     return it->second;
  //   return nullptr;
  // }

  // inline const VariableValue * get_val_by_argi(unsigned int arg_i) const noexcept
  // {
  //   if (auto it = _grain_val_by_argi.find(arg_i); it != _grain_val_by_argi.end())
  //     return it->second;
  //   return nullptr;
  // }
};
