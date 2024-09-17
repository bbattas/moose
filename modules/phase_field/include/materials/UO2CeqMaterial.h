#pragma once

#include "Material.h"
#include "DerivativeMaterialInterface.h"

/**
 * Calculates UO2 equilibrium concentration with or without irradiation for
 * vacancies or interstitials. Currently hard coded for 1600K.
 */

class UO2CeqMaterial : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  UO2CeqMaterial(const InputParameters & parameters);

protected:
  // From PCDTB
  virtual void computeQpProperties();

  /// number of solid phase order paramters
  const unsigned int _neta;

  /// solid phase order parameters
  std::vector<const VariableValue *> _eta;
  std::vector<VariableName> _eta_name;

  // const VariableValue & _T;

  // Void phase
  // const VariableValue & _c;
  // VariableName _c_name;

  // New ceq property
  MaterialPropertyName _ceq_name;
  MaterialProperty<Real> & _ceqs;
  std::vector<MaterialProperty<Real> *> _dceqsdeta;
  std::vector<std::vector<MaterialProperty<Real> *>> _d2ceqsdeta2;

  /// hgb
  const MaterialPropertyName _hgb_name;
  const MaterialProperty<Real> & _hgb;
  std::vector<const MaterialProperty<Real> *> _dhgbdeta;
  std::vector<std::vector<const MaterialProperty<Real> *>> _d2hgbdeta2;

  // Vacancy or interstitial, thermal or irradiation
  MooseEnum _vac_int;
};
