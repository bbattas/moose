//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"
#include "DerivativeMaterialInterface.h"

/**
 * This material calculates necessary parameters for the grand potential sintering model.
 * Especially those related to switching functions and thermodynamics.
 */
class GrandPotentialMultiSinteringMaterial : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  GrandPotentialMultiSinteringMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  /// number of solid phase order paramters
  const unsigned int _neta;

  /// solid phase order parameters
  std::vector<const VariableValue *> _eta;
  std::vector<VariableName> _eta_name;

  /// chemical potential vacancies
  const VariableValue & _wu;
  const NonlinearVariableName _wu_name;
  /// chemical potential interstitials
  const VariableValue & _wi;
  const NonlinearVariableName _wi_name;

  /// void phase order parameter
  const VariableValue & _phi;
  const NonlinearVariableName _phi_name;

  /// equilibrium vacancy concentration
  const MaterialPropertyName _csu_eq_name;
  const MaterialProperty<Real> & _csu_eq;
  std::vector<const MaterialProperty<Real> *> _dcsu_eq;
  std::vector<std::vector<const MaterialProperty<Real> *>> _d2csu_eq;
  /// equilibrium interstitial concentration
  const MaterialPropertyName _csi_eq_name;
  const MaterialProperty<Real> & _csi_eq;
  std::vector<const MaterialProperty<Real> *> _dcsi_eq;
  std::vector<std::vector<const MaterialProperty<Real> *>> _d2csi_eq;

  /// temperature
  const VariableValue & _T;

  /// vacancy void energy coefficient
  const MaterialProperty<Real> & _kvu;
  /// vacancy solid energy coefficient
  const MaterialProperty<Real> & _ksu;
  /// interstitial void energy coefficient
  const MaterialProperty<Real> & _kvi;
  /// interstitial solid energy coefficient
  const MaterialProperty<Real> & _ksi;

  /// void phase switching function
  MaterialProperty<Real> & _hv;
  MaterialProperty<Real> & _dhv;
  MaterialProperty<Real> & _d2hv;

  /// solid phase switching function
  MaterialProperty<Real> & _hs;
  MaterialProperty<Real> & _dhs;
  MaterialProperty<Real> & _d2hs;

  /// vacancy susceptibility
  MaterialProperty<Real> & _chiu;
  MaterialProperty<Real> & _dchiudphi;
  MaterialProperty<Real> & _dchiudw;
  MaterialProperty<Real> & _d2chiudphi2;
  MaterialProperty<Real> & _d2chiudw2;
  MaterialProperty<Real> & _d2chiudphidw;
  /// interstitial susceptibility
  MaterialProperty<Real> & _chii;
  MaterialProperty<Real> & _dchiidphi;
  MaterialProperty<Real> & _dchiidw;
  MaterialProperty<Real> & _d2chiidphi2;
  MaterialProperty<Real> & _d2chiidw2;
  MaterialProperty<Real> & _d2chiidphidw;

  /// void phase vacancy density
  MaterialProperty<Real> & _rhovu;
  MaterialProperty<Real> & _drhovudw;
  /// void phase interstial density
  MaterialProperty<Real> & _rhovi;
  MaterialProperty<Real> & _drhovidw;

  /// solid phase vacancy density
  MaterialProperty<Real> & _rhosu;
  MaterialProperty<Real> & _drhosudw;
  MaterialProperty<Real> & _d2rhosudw2;
  std::vector<MaterialProperty<Real> *> _drhosu;
  std::vector<MaterialProperty<Real> *> _d2rhosudwdeta;
  std::vector<std::vector<MaterialProperty<Real> *>> _d2rhosu;
  /// solid phase interstial density
  MaterialProperty<Real> & _rhosi;
  MaterialProperty<Real> & _drhosidw;
  MaterialProperty<Real> & _d2rhosidw2;
  std::vector<MaterialProperty<Real> *> _drhosi;
  std::vector<MaterialProperty<Real> *> _d2rhosidwdeta;
  std::vector<std::vector<MaterialProperty<Real> *>> _d2rhosi;

  /// void phase potential density
  MaterialProperty<Real> & _omegav;
  MaterialProperty<Real> & _domegavdwu;
  MaterialProperty<Real> & _d2omegavdwu2;
  MaterialProperty<Real> & _domegavdwi;
  MaterialProperty<Real> & _d2omegavdwi2;

  /// solid phase potential density
  MaterialProperty<Real> & _omegas;
  MaterialProperty<Real> & _domegasdwu;
  MaterialProperty<Real> & _d2omegasdwu2;
  MaterialProperty<Real> & _domegasdwi;
  MaterialProperty<Real> & _d2omegasdwi2;
  std::vector<MaterialProperty<Real> *> _domegasdeta;
  std::vector<MaterialProperty<Real> *> _d2omegasdwudeta;
  std::vector<MaterialProperty<Real> *> _d2omegasdwideta;
  std::vector<std::vector<MaterialProperty<Real> *>> _d2omegasdetadeta;

  /// energy barrier coefficient
  MaterialProperty<Real> & _mu;
  MaterialProperty<Real> & _dmu;
  MaterialProperty<Real> & _d2mu;

  /// gradient energy coefficient
  MaterialProperty<Real> & _kappa;
  MaterialProperty<Real> & _dkappa;
  MaterialProperty<Real> & _d2kappa;

  /// interface profile coefficient
  MaterialProperty<Real> & _gamma;

  /// Body Force coefficient for mass conservation in conc and chempot coupling
  MaterialProperty<Real> & _hv_c_min;
  MaterialProperty<Real> & _hs_c_min;

  /// MatReaction Force coefficient for mass conservation in conc and chempot coupling
  MaterialProperty<Real> & _hv_over_kVa;
  MaterialProperty<Real> & _hs_over_kVa;

  /// surface energy
  // const Real _sigma_s;
  const MaterialProperty<Real> & _sigma_s;

  /// grain boundary energy
  // const Real _sigma_gb;
  const MaterialProperty<Real> & _sigma_gb;

  /// interface width
  const Real _int_width;

  /// Parameter to determine accuracy of surface/GB phase switching function
  const Real _switch;

  /// Atomic volume of species
  const Real _Va;

  /// Type of energy function to use for the solid phase
  const MooseEnum _solid_energy;

  // Removing these into the actual material due to changing sigma gb and s to materials/variables
  // /// mu value on surfaces
  // const Real _mu_s;

  // /// mu value on grain boundaries
  // const Real _mu_gb;

  // /// kappa value on surfaces
  // const Real _kappa_s;

  // /// kappa value on grain boundaries
  // const Real _kappa_gb;

  /// Boltzmann constant
  const Real _kB;

  /// strict mass conservation flag
  bool _mass_conservation;
};
