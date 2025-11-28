#pragma once

#include "GBInclinationBase.h"

// Forward Declarations

/**
 * Inclination definition for
 */
class GBInclination : public GBInclinationBase
{
public:
  static InputParameters validParams();

  GBInclination(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  MaterialProperty<Real> & _testout1;
  MaterialProperty<Real> & _testout2;
  MaterialProperty<RealGradient> & _testoutgrad;
  MaterialProperty<RealTensorValue> & _testouttens;

  // Enum for which inclination function to use
  int _inc_func;

  // Inclination function constants
  const Real _if_a;
  const Real _if_b;
  const Real _if_c;
  const Real _if_d;

  // Inclination output (1+cos)
  MaterialProperty<std::vector<Real>> & _inclination;

  // gamma_ij and derivatives flatpacked to a vector with GBPairPacking
  MaterialProperty<std::vector<Real>> & _gamma_ij;
  MaterialProperty<std::vector<RealGradient>> & _dgamma_dgradeta;
  MaterialProperty<std::vector<RealTensorValue>> & _d2gamma_dgradeta2;

  // GBE
  const MaterialProperty<Real> & _gbe_iso; // input
  MaterialProperty<Real> & _gbe_aniso;     // output

  // Other Free Energy terms
  const MaterialProperty<Real> & _kappa;
  Real _const_m;
  MaterialProperty<Real> & _mu;

  MaterialProperty<Real> & _int_width;

  MaterialProperty<Real> & _gamma_qp;

  // AC Mobility
  const bool _aniso_L;
  const MaterialProperty<Real> & _L0;
  MaterialProperty<std::vector<Real>> & _L_ij;
  // MaterialProperty<std::vector<RealGradient>> & _dL_dgradeta;
  // MaterialProperty<std::vector<RealTensorValue>> & _d2L_dgradeta2;
  MaterialProperty<Real> & _L;

  /// aniso L Derivatives
  std::vector<MaterialProperty<Real> *> _dL_deta;
  std::vector<std::vector<MaterialProperty<Real> *>> _d2L_deta2;
  // d grad eta
  std::vector<MaterialPropertyName> _dL_dgradeta_name;
  std::vector<MaterialProperty<RealGradient> *> _dL_dgradeta;
  std::vector<MaterialPropertyName> _d2L_dgradeta2_name;
  std::vector<MaterialProperty<std::vector<RealTensorValue>> *> _d2L_dgradeta2;
  std::vector<std::vector<MaterialProperty<RealGradient> *>> _d2L_dgradetadeta;

  /// total number of grains
  // const unsigned int _op_num;

  // std::vector<const VariableValue *> _vals;
  // std::vector<VariableName> _vals_name;
  // std::vector<const VariableGradient *> _grad_vals;

  // // Inclination angular distance and derivatives flatpacked to a vector with GBPairPacking
  // MaterialProperty<std::vector<Real>> & _theta_ij;
  // MaterialProperty<std::vector<RealGradient>> & _dtheta_dgradeta;
  // MaterialProperty<std::vector<RealTensorValue>> & _d2theta_dgradeta2;

  // /// Enum for grain op identification
  // int _gb_case;
  // /// Grain tracker object
  // const GrainTracker * _grain_tracker;
  // const FeatureFloodCount * _ffc_tracker;
  // const Real _gt_tol;

  // // Number of OPs counted at each _qp
  // MaterialProperty<Real> & _gtnum;

  // // Vector of opi/opj pairs for iteration
  // std::vector<unsigned int> _gb_ij_list;

  // // Which trig approach to use for angular distance calculation
  // int _angular_func;

  // // Thresholds for the inverse magnitude of the gradient difference calc
  // const Real _intol;
  // const Real _altol;

  // // GB switching function for use with altol
  // const MaterialProperty<Real> & _hgb;
};
