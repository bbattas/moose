#pragma once

// MOOSE includes
#include "AuxKernel.h"

/**
 * Extract a component from the gradient of a variable
 */
class VariableGradUnitNorm : public VectorAuxKernel
{
public:
  static InputParameters validParams();

  /**
   * Class constructor
   * @param parameters Input parameters for the object
   */
  VariableGradUnitNorm(const InputParameters & parameters);

protected:
  virtual RealVectorValue computeValue() override;

private:
  /// Reference to the gradient of the coupled variable
  const VariableGradient & _gradient;

  // const VariableGradient & _grad_curve;

  /// Desired component
  int _component;
};
