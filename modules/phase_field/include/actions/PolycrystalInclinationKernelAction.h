#pragma once

#include "Action.h"

/**
 * Action that sets up ACInterfaceInclinationGamma kernel for all ops.
 * Requires GBInclination material for properties.
 */
class PolycrystalInclinationKernelAction : public Action
{
public:
  static InputParameters validParams();

  PolycrystalInclinationKernelAction(const InputParameters & params);

  virtual void act();

protected:
  /// number of grains to create
  const unsigned int _op_num;

  /// base name for the order parameter variables
  const std::string _var_name_base;
};
