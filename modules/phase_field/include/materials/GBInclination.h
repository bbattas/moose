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

  // Enum for which inclination function to use
  int _inc_func;

  // Inclination function constants
  const Real _if_a;
  const Real _if_b;
  const Real _if_c;

  // Inclination output (1+cos)
  MaterialProperty<std::vector<Real>> & _inclination;

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
