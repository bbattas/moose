#include "MatVectorijSum.h"

registerMooseObject("PhaseFieldApp", MatVectorijSum);

InputParameters
MatVectorijSum::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  params.addRequiredParam<MaterialPropertyName>("property", "Vector material property name");
  return params;
}

MatVectorijSum::MatVectorijSum(const InputParameters & params)
  : AuxKernel(params),
    // _vec_prop(getMaterialProperty<std::vector<Real>>(getParam<MaterialPropertyName>("property"))),
    // _i(getParam<unsigned int>("i")),
    // _j(getParam<unsigned int>("j")),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _vec_prop(getMaterialProperty<std::vector<Real>>(getParam<MaterialPropertyName>("property")))
{
}

Real
MatVectorijSum::computeValue()
{
  const auto & v = _vec_prop[_qp];
  // const unsigned k = GBPairPacking::pack_upper(_i, _j);
  Real num = 0.0;
  Real den = 0.0;
  for (std::size_t k = 0; k < v.size(); ++k)
  {
    auto [i, j] = GBPairPacking::unpack_upper(k);
    num += v[k] * (*_vals[i])[_qp] * (*_vals[i])[_qp] * (*_vals[j])[_qp] * (*_vals[j])[_qp];
    den += (*_vals[i])[_qp] * (*_vals[i])[_qp] * (*_vals[j])[_qp] * (*_vals[j])[_qp];
  }
  return (den > 1e-6) ? num / den : std::numeric_limits<Real>::quiet_NaN();
}
