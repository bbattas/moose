#include "InclinationAux.h"

registerMooseObject("PhaseFieldApp", InclinationAux);

InputParameters
InclinationAux::validParams()
{
  MooseEnum component("x=0 y=1 z=2");
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Creates a field consisting of one component of the gradient of a coupled variable.");
  params.addRequiredCoupledVar("var1", "The variable from which to compute the gradient component");
  params.addRequiredCoupledVar("var2", "The variable from which to compute the gradient component");
  params.addParam<MooseEnum>("component", component, "The gradient component to compute");
  return params;
}

InclinationAux::InclinationAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _grad1(coupledGradient("var1")),
    _grad2(coupledGradient("var2")),
    _component(getParam<MooseEnum>("component"))
{
}

Real
InclinationAux::computeValue()
{
  RealGradient raw = _grad1[_qp] - _grad2[_qp];
  Real mag = raw.norm();
  raw /= mag;
  return raw(_component);
}
