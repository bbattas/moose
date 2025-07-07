#pragma once

#include "Material.h"
#include "DerivativeMaterialInterface.h"

// Forward Declarations

/**
 * Inclination vector temporary output until i restructure ggincmat to better function
 */
class GGInclinationVector : public DerivativeMaterialInterface<Material> // Material
{
public:
  static InputParameters validParams();

  GGInclinationVector(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  /// total number of grains
  const unsigned int _op_num;

  std::vector<const VariableValue *> _vals;
  std::vector<VariableName> _vals_name;
  std::vector<const VariableGradient *> _grad_vals;

  MaterialProperty<RealGradient> & _inclination_vector;
  MaterialProperty<Real> & _ang_dist;

  /// Grain tracker object
  const GrainTracker & _grain_tracker;

  /// parameters to store the EBSD id and corresponding value on GB
  std::vector<unsigned int> _gb_ij_pairs;
  std::vector<unsigned int> _gb_ij_sorted;

  // For storing in inclination calc for combination at the end
  std::vector<Real> _hgb_pairs;
  std::vector<RealGradient> _inc_vec_pairs;
  std::vector<Real> _inc_pairs;
};
