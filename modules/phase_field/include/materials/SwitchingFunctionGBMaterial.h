#pragma once

#include "Material.h"
#include "DerivativeMaterialInterface.h"

// Forward Declarations

/**
 * SwitchingFunctionGBMaterial is a grain boundary switching function for a multi-phase,
 * multi-order parameter system.
 */
template <bool is_ad>
class SwitchingFunctionGBMaterialTempl : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  SwitchingFunctionGBMaterialTempl(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  /// Name of the function
  MaterialPropertyName _h_name;

  // /// Order parameters for phase alpha
  // const unsigned int _num_eta_gb;
  // const std::vector<const GenericVariableValue<is_ad> *> _eta_gb;
  // const std::vector<VariableName> _eta_gb_names;

  /// order parameters
  const unsigned int _op_num;
  const std::vector<const VariableValue *> _vals;
  const std::vector<VariableName> _vals_names;

  Real _hgb_threshold;

  int _func_type;

  // /// Order parameters for all phases (including alpha)
  // const unsigned int _num_eta;
  // const std::vector<const GenericVariableValue<is_ad> *> _eta;
  // const std::vector<VariableName> _eta_names;

  // /// List of which order parameters in the full list of all etas belong to phase p
  // std::vector<bool> _is_p;

  /// Switching function and derivatives
  GenericMaterialProperty<Real, is_ad> & _prop_h;
  std::vector<GenericMaterialProperty<Real, is_ad> *> _prop_dh;
  std::vector<std::vector<GenericMaterialProperty<Real, is_ad> *>> _prop_d2h;
};

typedef SwitchingFunctionGBMaterialTempl<false> SwitchingFunctionGBMaterial;
typedef SwitchingFunctionGBMaterialTempl<true> ADSwitchingFunctionGBMaterial;
