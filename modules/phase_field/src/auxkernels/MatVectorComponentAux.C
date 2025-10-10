#include "MatVectorComponentAux.h"

registerMooseObject("PhaseFieldApp", MatVectorComponentAux);

InputParameters
MatVectorComponentAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredParam<MaterialPropertyName>("property", "Vector material property name");
  params.addRequiredParam<unsigned int>("i", "Component index i of (i,j) to extract");
  params.addRequiredParam<unsigned int>("j", "Component index j of (i,j) to extract");
  params.addParam<bool>("gradient",
                        false,
                        "If true, treat 'property' as vector<RealGradient>; "
                        "otherwise as vector<Real>.");
  params.addParam<unsigned int>("component",
                                0,
                                "For gradient=true: component to extract (0=x,1=y,2=z). "
                                "Ignored when gradient=false.");
  return params;
}

MatVectorComponentAux::MatVectorComponentAux(const InputParameters & params)
  : AuxKernel(params),
    // _vec_prop(getMaterialProperty<std::vector<Real>>(getParam<MaterialPropertyName>("property"))),
    _i(getParam<unsigned int>("i")),
    _j(getParam<unsigned int>("j")),
    _gradient(getParam<bool>("gradient")),
    _component(getParam<unsigned int>("component")),
    _vec_prop(_gradient ? nullptr
                        : &getMaterialProperty<std::vector<Real>>(
                              getParam<MaterialPropertyName>("property"))),
    _vec_grad_prop(_gradient ? &getMaterialProperty<std::vector<RealGradient>>(
                                   getParam<MaterialPropertyName>("property"))
                             : nullptr)
{
  if (_gradient && !_vec_grad_prop)
    mooseError(name(), ": gradient=true but gradient property lookup failed.");
  if (!_gradient && !_vec_prop)
    mooseError(name(), ": gradient=false but scalar-vector property lookup failed.");
}

Real
MatVectorComponentAux::computeValue()
{
  // const auto & v = _vec_prop[_qp];
  const unsigned k = GBPairPacking::pack_upper(_i, _j);

  if (_gradient)
  {
    // if (!_vec_grad_prop)
    //   mooseError(name(), ": gradient=true but gradient property pointer is null.");
    const auto & v = (*_vec_grad_prop)[_qp];
    if (k >= v.size())
      return std::numeric_limits<Real>::quiet_NaN();

    // RealGradient supports operator()(comp) or operator[]
    return v[k](_component);
  }
  else
  {
    // if (!_vec_prop)
    //   mooseError(name(), ": gradient=false but scalar-vector property pointer is null.");
    const auto & v = (*_vec_prop)[_qp];
    return (k < v.size()) ? v[k] : std::numeric_limits<Real>::quiet_NaN();
  }

  // return (k < v.size()) ? v[k] : std::numeric_limits<Real>::quiet_NaN();
}
