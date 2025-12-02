#pragma once
#include "AuxKernel.h"
#include "GBPairPacking.h"

class MatVectorijSum : public AuxKernel
{
public:
  static InputParameters validParams();
  MatVectorijSum(const InputParameters &);

protected:
  virtual Real computeValue() override;

  const unsigned int _op_num;
  const std::vector<const VariableValue *> _vals;

  const MaterialProperty<std::vector<Real>> & _vec_prop;
};
