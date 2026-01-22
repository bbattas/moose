#pragma once

#include "Material.h"
#include "DerivativeMaterialInterface.h"
#include "GBPairPacking.h"

// Forward Declarations

/**
 * Inclination definition for
 */
class GBInclinationBase : public DerivativeMaterialInterface<Material> // Material
{
public:
  static InputParameters validParams();

  GBInclinationBase(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  /// total number of grains
  const unsigned int _op_num;

  std::vector<const VariableValue *> _vals;
  std::vector<VariableName> _vals_name;
  std::vector<const VariableGradient *> _grad_vals;

  // Inclination angular distance and derivatives flatpacked to a vector with GBPairPacking
  MaterialProperty<std::vector<Real>> & _theta_ij;
  MaterialProperty<std::vector<RealGradient>> & _dtheta_dgradeta;
  MaterialProperty<std::vector<RealTensorValue>> & _d2theta_dgradeta2;

  // Save the ij thats associated with each one so we can skip a sqrt calc later to unpack
  MaterialProperty<std::vector<unsigned int>> & _ij_i; // use unsigned int if dont need -1
  MaterialProperty<std::vector<unsigned int>> & _ij_j;

  /// Enum for grain op identification
  int _gb_case;
  /// Grain tracker object
  const GrainTracker * _grain_tracker;
  const FeatureFloodCount * _ffc_tracker;
  const Real _gt_tol;

  // Number of OPs counted at each _qp
  MaterialProperty<Real> & _gtnum;

  // Vector of opi/opj pairs for iteration
  std::vector<unsigned int> _gb_ij_list;

  // Which trig approach to use for angular distance calculation
  int _angular_func;

  // Thresholds for the inverse magnitude of the gradient difference calc
  const Real _intol;
  const Real _altol;

  // GB switching function for use with altol
  const MaterialProperty<Real> & _hgb;

  // Check if there are no actual ij pairs at this qp
  MaterialProperty<bool> & _no_ij_pairs;

  MaterialProperty<Real> & _testout3;
  MaterialProperty<RealGradient> & _aval;
  MaterialProperty<RealGradient> & _ival;
  MaterialProperty<RealGradient> & _acut;
  MaterialProperty<RealGradient> & _icut;

  // Use limit instead of zero
  const bool _limit_umag;
};
