#include "VariableGradUnitNorm.h"

registerMooseObject("PhaseFieldApp", VariableGradUnitNorm);

InputParameters
VariableGradUnitNorm::validParams()
{
  MooseEnum component("x=0 y=1 z=2");
  InputParameters params = VectorAuxKernel::validParams();
  params.addClassDescription(
      "Creates a field consisting of one component of the gradient of a coupled variable.");
  params.addRequiredCoupledVar("gradient_variable",
                               "The variable from which to compute the gradient component");
  params.addParam<MooseEnum>("component", component, "The gradient component to compute");
  return params;
}

VariableGradUnitNorm::VariableGradUnitNorm(const InputParameters & parameters)
  : VectorAuxKernel(parameters),
    _gradient(coupledGradient("gradient_variable")),
    _component(getParam<MooseEnum>("component"))
{
}

RealVectorValue
VariableGradUnitNorm::computeValue()
{
  RealVectorValue value;
  value = _gradient[_qp] / _gradient[_qp].norm();
  return value;
}
