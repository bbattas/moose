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

// Forward Declarations

/**
 * Inclination dependent properties for AGG
 */
class GGInclinationMaterial : public DerivativeMaterialInterface<Material> // Material
{
public:
  static InputParameters validParams();

  GGInclinationMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  /// total number of grains
  const unsigned int _op_num;

  std::vector<const VariableValue *> _vals;
  std::vector<VariableName> _vals_name;
  std::vector<const VariableGradient *> _grad_vals;
  std::vector<std::vector<RealGradient>> _incl_tens;
  std::vector<std::vector<Real>> _ang_dist;

  // MaterialProperty<RealGradient> & _inclination;
  MaterialProperty<Real> & _inclination;
  MaterialProperty<Real> & _inclination_distance;
  // MaterialProperty<Real> & _ang_dist;
  // std::vector<std::vector<Real>> & _inclination;
  // RealTensorValue

  /// Grain tracker object
  const GrainTracker & _grain_tracker;

  /// EBSD reader user object
  // const EBSDReader & _ebsd_reader;

  // Inclination function constants
  Real _delta_ij;
  Real _theta_pre;
  Real _inc_ij_0;

  // gamma testing
  MaterialProperty<Real> & _gamma;

  /// gamma gradient names
  std::vector<MaterialPropertyName> _dgammadgrad_eta_name;
  /// All the actual gradients of gamma with respect to eta
  std::vector<MaterialProperty<RealGradient> *> _dgammadgrad_eta;
  // Second derivatives
  std::vector<MaterialPropertyName> _d2gammadgrad_eta2_name;
  std::vector<MaterialProperty<RealTensorValue> *> _d2gammadgrad_eta2;

  // L mobility parameters
  const bool _L_of_eta;
  const MaterialPropertyName _L_name;
  MaterialProperty<Real> & _L;
  std::vector<MaterialProperty<Real> *> _dLdeta;
  std::vector<std::vector<MaterialProperty<Real> *>> _d2Ldetadeta;
  // Gradient derivatives
  std::vector<MaterialPropertyName> _dLdgrad_eta_name;
  std::vector<MaterialProperty<RealGradient> *> _dLdgrad_eta;
  std::vector<MaterialPropertyName> _d2Ldgrad_eta2_name;
  std::vector<MaterialProperty<RealTensorValue> *> _d2Ldgrad_eta2;
  std::vector<std::vector<MaterialProperty<RealGradient> *>> _d2Ldgrad_etadeta;

  // MaterialProperty<Real> & _testout;
  // MaterialProperty<Real> & _testout2;
  // MaterialProperty<RealGradient> & _incder_temp;

  const MaterialProperty<Real> & _gbe;
  MaterialProperty<Real> & _gbe_inc;

  Real _kappa;
  Real _const_m;

  const MaterialProperty<Real> & _L0;
  const MaterialProperty<Real> & _gamma0;

  MaterialProperty<Real> & _int_width;
  MaterialProperty<Real> & _mu;

  // other stuff

  /// parameters to store the EBSD id and corresponding value on GB
  // std::vector<unsigned int> _gb_pairs;
  std::vector<Real> _gb_op_pairs;
  std::vector<unsigned int> _gb_ij_pairs;
  std::vector<unsigned int> _gb_ij_sorted;

  // For storing in inclination calc for combination at the end
  std::vector<Real> _hgb_pairs;
  std::vector<Real> _inc_pairs;
  // std::vector<RealGradient> _dincdgradetai_pairs;
};
