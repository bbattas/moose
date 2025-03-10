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

// Forward Declarations

/**
 * Calculated properties for a single component phase field model using polynomial free energies
 */
class GGInclinationMaterial : public Material
{
public:
  static InputParameters validParams();

  GGInclinationMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  /// total number of grains
  const unsigned int _op_num;

  std::vector<const VariableValue *> _vals;
  std::vector<const VariableGradient *> _grad_vals;
  std::vector<std::vector<RealGradient>> _incl_tens;
  std::vector<std::vector<Real>> _ang_dist;

  // MaterialProperty<RealGradient> & _inclination;
  MaterialProperty<Real> & _inclination;
  // MaterialProperty<Real> & _ang_dist;
  // std::vector<std::vector<Real>> & _inclination;
  // RealTensorValue

  /// Grain tracker object
  const GrainTracker & _grain_tracker;

  /// EBSD reader user object
  const EBSDReader & _ebsd_reader;

  MaterialProperty<Real> & _grains_on_gb;
  MaterialProperty<Real> & _gb_id;

  /// parameters to store the EBSD id and corresponding value on GB
  std::vector<unsigned int> _gb_pairs;
  std::vector<Real> _gb_op_pairs;
  std::vector<unsigned int> _gb_ij_pairs;

  // ij for temp output
  const unsigned int _i_value;
  const unsigned int _j_value;

  // ///Variable values
  // const VariableValue & _c;
  // const VariableValue & _T;

  // ///Mateiral property declarations
  // MaterialProperty<Real> & _M;
  // MaterialProperty<RealGradient> & _grad_M;

  // MaterialProperty<Real> & _kappa;
  // MaterialProperty<Real> & _c_eq;
  // MaterialProperty<Real> & _W;
  // MaterialProperty<Real> & _Qstar;
  // MaterialProperty<Real> & _D;

  // ///Input parameters
  // Real _int_width;
  // Real _length_scale;
  // Real _time_scale;
  // MooseEnum _order;
  // Real _D0;
  // Real _Em;
  // Real _Ef;
  // Real _surface_energy;

  // const Real _JtoeV;
  // const Real _kb;
};
