##############################################################################
# File: 16_2grain_withVoid_masscons.i
# File Location: /examples/sintering/paper3/19_massLoss_tests/16_2grain_withVoid_masscons
# Created Date: Friday June 7th 2024
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Friday June 7th 2024
# Modified By: Brandon Battas
# -----
# Description:
#  Other tests were wild, so im starting from the beginning with the simplest
#   input that has worked and im adding more, starting with void phase and
#   our various materials and whatnot
#
##############################################################################

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 60
  ny = 60
  xmin = -15
  xmax = 15
  ymin = -15
  ymax = 15
  uniform_refine = 0
  parallel_type = DISTRIBUTED
[]

[GlobalParams]
  profile = TANH
  int_width = 4
  op_num = 2
  var_name_base = gr
[]

[Variables]
  [w]
    initial_condition = 0
  []
  [c]
    # initial_condition = 0
  []
  [gr0]
  []
  [gr1]
  []
[]

[ICs]
  [IC_c]
    type = SolutionIC
    solution_uo = exo_soln
    variable = c
    from_variable = 'c_auxvar'
  []
  # [IC_c]
  #   type = SmoothCircleIC
  #   variable = c
  #   x1 = -30
  #   y1 = 0
  #   radius = 35
  #   invalue = 0.1
  #   outvalue = 0.9
  # []
  [IC_gr0]
    type = SmoothCircleIC
    variable = gr0
    x1 = -30
    y1 = 0
    radius = 35
    invalue = 1
    outvalue = 0
  []
  [IC_gr1]
    type = SmoothCircleIC
    variable = gr1
    x1 = -30
    y1 = 0
    radius = 35
    invalue = 0
    outvalue = 1
  []
[]


[UserObjects]
  [exo_soln]
    type = SolutionUserObject
    mesh = 'solutionIC_genC_out.e'
    system_variables = 'c_auxvar'
    timestep = LATEST
  []
[]


[Modules]
  [PhaseField]
    [GrandPotential]
      switching_function_names = 'ha hb'
      # Chempot
      chemical_potentials = 'w'
      mobilities = 'Dchi'
      anisotropic = 'false'
      susceptibilities = 'chi'
      free_energies_w = 'rhoa rhob'
      # Grains
      mobility_name_gr = L
      kappa_gr = kappa
      gamma_gr = gab
      free_energies_gr = 'omegaa omegab'
      # # Other OPs (Phi)
      # additional_ops = 'phi'
      # mobility_name_op = L
      # kappa_op = kappa
      # gamma_grxop = gab
      free_energies_op = '' #empty when no phi'omegaa omegab'
      # Mass Conservation
      mass_conservation = true
      concentrations = 'c'
      hj_over_kVa = 'hoverk_a hoverk_b'
      hj_c_min = 'caeq_mask cbeq_mask'
    []
  []
[]

[Materials]
  [ha_test]
    type = SwitchingFunctionMultiPhaseMaterial
    h_name = ha
    all_etas = 'gr0 gr1'
    phase_etas = 'gr0'
  []
  [hb_test]
    type = SwitchingFunctionMultiPhaseMaterial
    h_name = hb
    all_etas = 'gr0 gr1'
    phase_etas = 'gr1'
  []
  [hgb]
    type = DerivativeParsedMaterial
    property_name = hgb
    coupled_variables = 'gr0 gr1'
    expression = '16 * ( (gr0 * gr1)^2 )' #+ (gr0 * etad0)^2 + (gr1 * etad0)^2)
    # expression = '4*(1 - (gr0^2 + gr1^2 + etad0^2))^2'
  []
  [caeq]
    type = DerivativeParsedMaterial
    property_name = caeq
    constant_names = 'cb cgb'
    constant_expressions = '0.1 0.2'
    coupled_variables = 'gr0 gr1'
    material_property_names = 'hgb(gr0,gr1)'
    expression = 'hgb * cgb + (1 - hgb) * cb'
  []
  [cbeq]
    type = DerivativeParsedMaterial
    property_name = cbeq
    constant_names = 'cb cgb'
    constant_expressions = '0.9 0.95'
    coupled_variables = 'gr0 gr1'
    material_property_names = 'hgb(gr0,gr1)'
    expression = 'hgb * cgb + (1 - hgb) * cb'
  []
  [omegaa]
    type = DerivativeParsedMaterial
    coupled_variables = 'w gr0 gr1' #was just w now includes etas for derivatives
    property_name = omegaa
    material_property_names = 'Vm ka caeq(gr0,gr1)'
    expression = '-0.5*w^2/Vm^2/ka-w/Vm*caeq'
    derivative_order = 2
  []
  [omegab]
    type = DerivativeParsedMaterial
    coupled_variables = 'w gr0 gr1' #was just w now includes etas for derivatives
    property_name = omegab
    material_property_names = 'Vm kb cbeq(gr0,gr1)'
    expression = '-0.5*w^2/Vm^2/kb-w/Vm*cbeq'
    derivative_order = 2
  []
  [rhoa]
    type = DerivativeParsedMaterial
    coupled_variables = 'w gr0 gr1' #was just w now includes etas
    property_name = rhoa
    material_property_names = 'Vm ka caeq(gr0,gr1)'
    expression = 'w/Vm^2/ka + caeq/Vm'
    derivative_order = 2
  []
  [rhob]
    type = DerivativeParsedMaterial
    coupled_variables = 'w gr0 gr1' #was just w now includes etas
    property_name = rhob
    material_property_names = 'Vm kb cbeq(gr0,gr1)'
    expression = 'w/Vm^2/kb + cbeq/Vm'
    derivative_order = 2
  []
  [cmat]
    type = ParsedMaterial
    material_property_names = 'Vm rhoa(w,gr0,gr1) rhob(w,gr0,gr1) ha(gr0,gr1) hb(gr0,gr1)'
    expression = 'Vm * (ha * rhoa + hb * rhob)'
    property_name = cmat
  []
  [const]
    type = GenericConstantMaterial
    prop_names = 'kappa_c  kappa   L   D    Vm    ka  caeq_c  kb  cbeq_c  gab  mu  tgrad_corr_mult'
    prop_values = '0        1       1.0 1.0  1.0  10.0  0.1  10.0  0.9   1.5  1.0     0.0'
  []
  [Mobility]
    type = DerivativeParsedMaterial
    property_name = Dchi
    material_property_names = 'D chi(gr0,gr1) Vm'
    coupled_variables = 'gr0 gr1'
    expression = 'D*chi*Vm' #Vm needed because evolving c instead of rho
    derivative_order = 2
  []
  [chi]
    type = DerivativeParsedMaterial
    property_name = chi
    material_property_names = 'Vm ha(gr0,gr1) ka hb(gr0,gr1) kb'
    expression = '(ha/ka + hb/kb) / Vm^2'
    coupled_variables = 'gr0 gr1'
    derivative_order = 2
  []
  # For c evolution instead of rho for mass conservation
  [hoverk_a]
    type = DerivativeParsedMaterial
    material_property_names = 'ha(gr0,gr1) Vm ka'
    property_name = hoverk_a
    expression = 'ha / Vm / ka'
  []
  [hoverk_b]
    type = DerivativeParsedMaterial
    material_property_names = 'hb(gr0,gr1) Vm kb'
    property_name = hoverk_b
    expression = 'hb / Vm / kb'
  []
  [caeq_mask]
    type = DerivativeParsedMaterial
    property_name = caeq_mask
    material_property_names = 'ha(gr0,gr1) caeq(gr0,gr1)'
    coupled_variables = 'gr0 gr1'
    expression = 'ha * caeq'
  []
  [cbeq_mask]
    type = DerivativeParsedMaterial
    property_name = cbeq_mask
    material_property_names = 'hb(gr0,gr1) cbeq(gr0,gr1)'
    coupled_variables = 'gr0 gr1'
    expression = 'hb * cbeq'
  []
[]

[Postprocessors]
  [c_total]
    type = ElementIntegralVariablePostprocessor
    variable = c
  []
  [gr0_total]
    type = ElementIntegralVariablePostprocessor
    variable = gr0
  []
  [gr1_total]
    type = ElementIntegralVariablePostprocessor
    variable = gr1
  []
  [hgb_total]
    type = ElementIntegralMaterialProperty
    mat_prop = hgb
  []
  [cmat_total]
    type = ElementIntegralMaterialProperty
    mat_prop = cmat
  []
  [runtime]
    type = PerfGraphData
    section_name = "Root"
    data_type = TOTAL
  []
  [timestep]
    type = TimestepSize
    outputs = csv
  []
  [nonlinear]
    type = NumNonlinearIterations
    outputs = csv
  []
  [linear]
    type = NumLinearIterations
    outputs = csv
  []
  [residuals]
    type = NumResidualEvaluations
    outputs = csv
  []
  [n_DOFs]
    type = NumDOFs
    outputs = csv
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  nl_max_its = 15
  scheme = bdf2
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -sub_pc_type'
  petsc_options_value = 'asm lu'
  l_max_its = 20
  l_tol = 1.0e-6
  nl_rel_tol = 1.0e-8
  start_time = 0.0
  # end_time = 1000
  steady_state_detection = true
  num_steps = 5
  nl_abs_tol = 1e-10
  dt = 1
[]

[Outputs]
  csv = true
  exodus = true
  checkpoint = false
[]

