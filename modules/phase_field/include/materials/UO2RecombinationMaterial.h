#pragma once

#include "Material.h"
#include "DerivativeMaterialInterface.h"

/**
 * Calculates UO2 vacancy and interstitial recombination based on
 * the concentration variables. Hardcoded for 1600K for now.
 */

class UO2RecombinationMaterial : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  UO2RecombinationMaterial(const InputParameters & parameters);

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

  // Susceptability*Va
  const MaterialProperty<Real> & _chiu;
  const MaterialProperty<Real> & _chii;

  // Output Recombination
  MaterialPropertyName _rec_name;
  MaterialProperty<Real> & _rec;
  MaterialProperty<Real> & _drecdphi;
  MaterialProperty<Real> & _drecdwu;
  MaterialProperty<Real> & _drecdwi;
  MaterialProperty<Real> & _drecdcu;
  MaterialProperty<Real> & _drecdci;

  // Constants
  const Real _hv_thresh;
  const Real _Va;
  const Real _kB;

  // Cascade Constants
  const Real _znum;
  const Real _a0;
  const Real _Di0;
  const Real _EiB;

  // Switching functions
  const MaterialProperty<Real> & _hv;
  const MaterialProperty<Real> & _hs;
  const MaterialProperty<Real> & _dhs;

  // ENUM to select which version of recombination ifs to use
  MooseEnum _if_case;
  // hs interface narrowing divisor if applicable
  const Real _switch;
};
