#pragma once

#include "Material.h"
#include "DerivativeMaterialInterface.h"

/**
 * Calculates UO2 ballistic mixing materials for vacancies
 * and interstitials.
 */

class UO2MixingMaterial : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  UO2MixingMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  /// number of solid phase order paramters
  const unsigned int _neta;

  /// solid phase order parameters
  std::vector<const VariableValue *> _eta;
  std::vector<VariableName> _eta_name;

  // CHEMICAL POTENTIALS
  /// chemical potential vacancies
  const VariableValue & _wu;
  const NonlinearVariableName _wu_name;
  /// chemical potential interstitials
  const VariableValue & _wi;
  const NonlinearVariableName _wi_name;

  // CONCENTRATIONS
  /// concentration vacancies
  const VariableValue & _cu;
  const NonlinearVariableName _cu_name;
  /// concentration interstitials
  const VariableValue & _ci;
  const NonlinearVariableName _ci_name;

  /// void phase order parameter
  const VariableValue & _phi;
  const NonlinearVariableName _phi_name;

  // const VariableValue & _T;

  // Susceptability*Va and Derivatives
  const MaterialProperty<Real> & _chiu;
  const MaterialProperty<Real> & _dchiudphi;
  const MaterialProperty<Real> & _dchiudwu;
  const MaterialProperty<Real> & _chii;
  const MaterialProperty<Real> & _dchiidphi;
  const MaterialProperty<Real> & _dchiidwi;

  // Output Mixing Vacancy value
  MaterialPropertyName _bmv_name;
  MaterialProperty<Real> & _bmv;
  MaterialProperty<Real> & _dbmvdphi;
  MaterialProperty<Real> & _dbmvdwu;
  MaterialProperty<Real> & _dbmvdcu;

  // Output Mixing Interstitial value
  MaterialPropertyName _bmi_name;
  MaterialProperty<Real> & _bmi;
  MaterialProperty<Real> & _dbmidphi;
  MaterialProperty<Real> & _dbmidwi;
  MaterialProperty<Real> & _dbmidci;

  // Constants
  const Real _fdot;
  const Real _Va;

  // Cascade Constants
  const Real _Nc;
  const Real _Vc;
  const Real _tc;
  const Real _Dc;
};
