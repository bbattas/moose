#include "PolycrystalInclinationKernelAction.h"
#include "Factory.h"
#include "Conversion.h"
#include "FEProblem.h"

registerMooseAction("PhaseFieldApp", PolycrystalInclinationKernelAction, "add_kernel");

InputParameters
PolycrystalInclinationKernelAction::validParams()
{
  InputParameters params = Action::validParams();
  params.addClassDescription("Set up ACInterfaceInclinationGamma kernels for all grains.");
  params.addRequiredParam<unsigned int>(
      "op_num", "specifies the total number of grains (deformed + recrystallized) to create");
  params.addRequiredParam<std::string>("var_name_base", "specifies the base name of the variables");
  // params.addParam<VariableName>("c", "Name of coupled concentration variable");
  // params.addParam<Real>("en_ratio", 1.0, "Ratio of surface to GB energy");
  // params.addParam<unsigned int>("ndef", 0, "specifies the number of deformed grains to create");
  // params.addParam<bool>("implicit", true, "Whether kernels are implicit or not");
  // params.addParam<bool>(
  //     "use_displaced_mesh", false, "Whether to use displaced mesh in the kernels");
  params.addParam<MaterialPropertyName>(
      "hgb_mask", "Name of GB switching function, if used. If not specified defaults to no mask.");
  params.addParam<bool>("variable_mobility",
                        false,
                        "The mobility is a function of any MOOSE variable (if "
                        "this is set to false, L must be constant over the "
                        "entire domain!)");
  // params.addCoupledVar("args", "Vector of nonlinear variable arguments that L depends on");
  // params.deprecateCoupledVar("args", "coupled_variables", "02/27/2024");

  return params;
}

PolycrystalInclinationKernelAction::PolycrystalInclinationKernelAction(
    const InputParameters & params)
  : Action(params),
    _op_num(getParam<unsigned int>("op_num")),
    _var_name_base(getParam<std::string>("var_name_base"))
{
}

void
PolycrystalInclinationKernelAction::act()
{
  auto variable_L = getParam<bool>("variable_mobility");
  //
  for (unsigned int op = 0; op < _op_num; ++op)
  {
    // Create variable names
    std::string var_name = _var_name_base + Moose::stringify(op);
    std::vector<VariableName> cv;
    cv.resize(_op_num - 1);

    // Other ops for coupled variables
    unsigned int ind = 0;
    for (unsigned int j = 0; j < _op_num; ++j)
      if (j != op)
        cv[ind++] = _var_name_base + Moose::stringify(j);

    // Set up ACInterfaceInclinationGamma kernel
    {
      InputParameters params = _factory.getValidParams("ACInterfaceInclinationGamma");
      params.set<NonlinearVariableName>("variable") = var_name;
      params.set<std::vector<VariableName>>("coupled_variables") = cv;
      params.set<bool>("variable_L") = variable_L;
      if (isParamValid("hgb_mask"))
        params.set<MaterialPropertyName>("mask_name") = getParam<MaterialPropertyName>("hgb_mask");
      params.applyParameters(parameters());

      std::string kernel_name = "ACInterfaceInclination_" + var_name;
      _problem->addKernel("ACInterfaceInclinationGamma", kernel_name, params);
    }
  }
}
