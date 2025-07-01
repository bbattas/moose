#pragma once

// MOOSE includes
#include "AuxKernel.h"

/**
 * ECalculate the OP gradient based inclination
 */
class InclinationAux : public AuxKernel
{
public:
  static InputParameters validParams();

  /**
   * Class constructor
   * @param parameters Input parameters for the object
   */
  InclinationAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

private:
  /// Reference to the gradient of the coupled variable
  const VariableGradient & _grad1;
  const VariableGradient & _grad2;
  /// Desired component
  int _component;
  // Other
  // RealGradient raw;
  // Real mag;
};
