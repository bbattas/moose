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
  // Azimuth angle [0,2pi)
  MaterialProperty<std::vector<Real>> & _theta_ij;
  MaterialProperty<std::vector<RealGradient>> & _dtheta_dgradeta;
  MaterialProperty<std::vector<RealTensorValue>> & _d2theta_dgradeta2;
  // POLAR angle (from z axis [0,pi/2))
  MaterialProperty<std::vector<Real>> & _polar_ij;
  MaterialProperty<std::vector<RealGradient>> & _dpolar_dgradeta;
  MaterialProperty<std::vector<RealTensorValue>> & _d2polar_dgradeta2;

  // compact active pair OP indices
  MaterialProperty<std::vector<unsigned int>> & _ij_i;
  MaterialProperty<std::vector<unsigned int>> & _ij_j;
  // OP activity flags at this qp
  // op_is_present[op] = 1 if GrainTracker/FFC says this OP has a valid feature on this elem
  // op_has_active_pair[op] = 1 if this OP appears in at least one saved ij pair
  MaterialProperty<std::vector<unsigned char>> & _op_is_present;
  MaterialProperty<std::vector<unsigned char>> & _op_has_active_pair;
  // GT feature id for misorientation
  MaterialProperty<std::vector<unsigned int>> & _ug_i;
  MaterialProperty<std::vector<unsigned int>> & _ug_j;

  /// Enum for grain op identification
  int _gb_case;

  /// Grain tracker object
  const GrainTracker * _grain_tracker;
  const FeatureFloodCount * _ffc_tracker;

  // Vector of opi/opj pairs for iteration
  std::vector<std::pair<unsigned int, unsigned int>> _gb_pairs;

  // Which trig approach to use for angular distance calculation
  int _angular_func;

  // Thresholds for the inverse magnitude of the gradient difference calc
  const Real _alphatol;
  const Real _hgbatol;

  // GB switching function for use with altol
  const MaterialProperty<Real> & _hgb;

  // Check if there are no actual ij pairs at this qp
  MaterialProperty<bool> & _no_ij_pairs;

  // GT Count Checker
  MaterialProperty<Real> & _gtnum;

  // Use limit instead of zero
  const bool _limit_umag;
};
