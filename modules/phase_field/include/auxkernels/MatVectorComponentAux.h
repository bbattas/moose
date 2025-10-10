#pragma once
#include "AuxKernel.h"
#include "GBPairPacking.h"

class MatVectorComponentAux : public AuxKernel
{
public:
  static InputParameters validParams();
  MatVectorComponentAux(const InputParameters &);

protected:
  virtual Real computeValue() override;

  const unsigned _i;
  const unsigned _j;

  const bool _gradient;      // true => read from _vec_grad_prop, false => _vec_prop
  const unsigned _component; // only used when _gradient=true (0=x,1=y,2=z)

  // One of these will be set based on the 'gradient' flag; the other stays nullptr
  const MaterialProperty<std::vector<Real>> * _vec_prop;              // vector<Real>
  const MaterialProperty<std::vector<RealGradient>> * _vec_grad_prop; // vector<RealGradient>
};
