##############################################################################
# File: 07_2grain_noVoid_masscons.i
# File Location: /examples/sintering/paper3/19_massLoss_tests/07_2grain_noVoid_masscons
# Created Date: Thursday June 6th 2024
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Thursday June 6th 2024
# Modified By: Brandon Battas
# -----
# Description:
#  Same as 05 but using the mass conservation approach
#  Testing 2 grains with a curved boundary and no void phase
#  Using a constant ceq for each phase (no gb dependence)
#  Based on the GrandPotential3Phase and same_masscons but with plenty of changes
#  Might have overspecified the derivites in materials() but just playing safe
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
[]

[GlobalParams]
  profile = TANH
  int_width = 4
[]

[Variables]
  [w]
    initial_condition = 0
  []
  [c]
  []
  [etaa0]
  []
  [etab0]
  []
[]

[ICs]
  [IC_c]
    type = SmoothCircleIC
    variable = c
    x1 = -30
    y1 = 0
    radius = 35
    invalue = 0.1
    outvalue = 0.9
  []
  [IC_etaa0]
    type = SmoothCircleIC
    variable = etaa0
    x1 = -30
    y1 = 0
    radius = 35
    invalue = 1
    outvalue = 0
  []
  [IC_etab0]
    type = SmoothCircleIC
    variable = etab0
    x1 = -30
    y1 = 0
    radius = 35
    invalue = 0
    outvalue = 1
  []
[]

[Kernels]
  # Order parameter eta_alpha0
  [ACa0_bulk]
    type = ACGrGrMulti
    variable = etaa0
    v = 'etab0'
    gamma_names = 'gab'
  []
  [ACa0_sw]
    type = ACSwitching
    variable = etaa0
    Fj_names = 'omegaa omegab'
    hj_names = 'ha     hb'
    coupled_variables = 'etab0 w'
  []
  [ACa0_int]
    type = ACInterface
    variable = etaa0
    kappa_name = kappa
  []
  [ea0_dot]
    type = TimeDerivative
    variable = etaa0
  []
  # Order parameter eta_beta0
  [ACb0_bulk]
    type = ACGrGrMulti
    variable = etab0
    v = 'etaa0'
    gamma_names = 'gab'
  []
  [ACb0_sw]
    type = ACSwitching
    variable = etab0
    Fj_names = 'omegaa omegab'
    hj_names = 'ha     hb'
    coupled_variables = 'etaa0 w'
  []
  [ACb0_int]
    type = ACInterface
    variable = etab0
    kappa_name = kappa
  []
  [eb0_dot]
    type = TimeDerivative
    variable = etab0
  []
  #Concentration
  [c_dot]
    type = TimeDerivative
    variable = c
  []
  [Diffusion]
    type = MatDiffusion
    variable = c
    v = w
    diffusivity = Dchi
    args = 'etaa0 etab0'
  []
  #The following relate chemical potential to composition using Eq. (22)
  [w_rxn]
    type = MatReaction
    variable = w
    v = c
    mob_name = -1
  []
  [ca_rxn]
    type = MatReaction
    variable = w
    mob_name = 'hoverk_a'
    args = 'etaa0 etab0'
  []
  [ca_bodyforce]
    type = MaskedBodyForce
    variable = w
    mask = caeq_mask #ha
    coupled_variables = 'etaa0 etab0'
    value = 1 #0.1 #caeq
  []
  [cb_rxn]
    type = MatReaction
    variable = w
    mob_name = 'hoverk_b'
    args = 'etaa0 etab0'
  []
  [cb_bodyforce]
    type = MaskedBodyForce
    variable = w
    mask = cbeq_mask #hb
    coupled_variables = 'etaa0 etab0'
    value = 1 #0.9 #cbeq
  []
  # Not needed anymore if evolving based on C
  # [coupled_etaa0dot]
  #   type = CoupledSwitchingTimeDerivative
  #   variable = w
  #   v = etaa0
  #   Fj_names = 'rhoa rhob'
  #   hj_names = 'ha   hb'
  #   coupled_variables = 'etaa0 etab0 '
  # []
  # [coupled_etab0dot]
  #   type = CoupledSwitchingTimeDerivative
  #   variable = w
  #   v = etab0
  #   Fj_names = 'rhoa rhob'
  #   hj_names = 'ha   hb'
  #   coupled_variables = 'etaa0 etab0'
  # []
[]

[Materials]
  [ha_test]
    type = SwitchingFunctionMultiPhaseMaterial
    h_name = ha
    all_etas = 'etaa0 etab0'
    phase_etas = 'etaa0'
  []
  [hb_test]
    type = SwitchingFunctionMultiPhaseMaterial
    h_name = hb
    all_etas = 'etaa0 etab0'
    phase_etas = 'etab0'
  []
  [hgb]
    type = DerivativeParsedMaterial
    property_name = hgb
    coupled_variables = 'etaa0 etab0'
    expression = '16 * ( (etaa0 * etab0)^2 )' #+ (etaa0 * etad0)^2 + (etab0 * etad0)^2)
    # expression = '4*(1 - (etaa0^2 + etab0^2 + etad0^2))^2'
  []
  [caeq]
    type = DerivativeParsedMaterial
    property_name = caeq
    constant_names = 'cb cgb'
    constant_expressions = '0.1 0.2'
    coupled_variables = 'etaa0 etab0'
    material_property_names = 'hgb(etaa0,etab0)'
    expression = 'cb' #'hgb * cgb + (1 - hgb) * cb'
  []
  [cbeq]
    type = DerivativeParsedMaterial
    property_name = cbeq
    constant_names = 'cb cgb'
    constant_expressions = '0.9 0.95'
    coupled_variables = 'etaa0 etab0'
    material_property_names = 'hgb(etaa0,etab0)'
    expression = 'cb' #'hgb * cgb + (1 - hgb) * cb'
  []
  [omegaa]
    type = DerivativeParsedMaterial
    coupled_variables = 'w etaa0 etab0' #was just w now includes etas for derivatives
    property_name = omegaa
    material_property_names = 'Vm ka caeq(etaa0,etab0)'
    expression = '-0.5*w^2/Vm^2/ka-w/Vm*caeq'
    derivative_order = 2
  []
  [omegab]
    type = DerivativeParsedMaterial
    coupled_variables = 'w etaa0 etab0' #was just w now includes etas for derivatives
    property_name = omegab
    material_property_names = 'Vm kb cbeq(etaa0,etab0)'
    expression = '-0.5*w^2/Vm^2/kb-w/Vm*cbeq'
    derivative_order = 2
  []
  [rhoa]
    type = DerivativeParsedMaterial
    coupled_variables = 'w etaa0 etab0' #was just w now includes etas
    property_name = rhoa
    material_property_names = 'Vm ka caeq(etaa0,etab0)'
    expression = 'w/Vm^2/ka + caeq/Vm'
    derivative_order = 2
  []
  [rhob]
    type = DerivativeParsedMaterial
    coupled_variables = 'w etaa0 etab0' #was just w now includes etas
    property_name = rhob
    material_property_names = 'Vm kb cbeq(etaa0,etab0)'
    expression = 'w/Vm^2/kb + cbeq/Vm'
    derivative_order = 2
  []
  [cmat]
    type = ParsedMaterial
    material_property_names = 'Vm rhoa(w,etaa0,etab0) rhob(w,etaa0,etab0) ha(etaa0,etab0) hb(etaa0,etab0)'
    expression = 'Vm * (ha * rhoa + hb * rhob)'
    property_name = cmat
    outputs = 'exodus'
  []
  [const]
    type = GenericConstantMaterial
    prop_names = 'kappa_c  kappa   L   D    Vm    ka  caeq_c  kb  cbeq_c  gab  mu  tgrad_corr_mult'
    prop_values = '0        1       1.0 1.0  1.0  10.0  0.1  10.0  0.9   1.5  1.0     0.0'
  []
  [Mobility]
    type = DerivativeParsedMaterial
    property_name = Dchi
    material_property_names = 'D chi(etaa0,etab0) Vm'
    coupled_variables = 'etaa0 etab0'
    expression = 'D*chi*Vm' #Vm needed because evolving c instead of rho
    derivative_order = 2
  []
  [chi]
    type = DerivativeParsedMaterial
    property_name = chi
    material_property_names = 'Vm ha(etaa0,etab0) ka hb(etaa0,etab0) kb'
    expression = '(ha/ka + hb/kb) / Vm^2'
    coupled_variables = 'etaa0 etab0'
    derivative_order = 2
  []
  # For c evolution instead of rho for mass conservation
  [hoverk_a]
    type = DerivativeParsedMaterial
    material_property_names = 'ha(etaa0,etab0) Vm ka'
    property_name = hoverk_a
    expression = 'ha / Vm / ka'
  []
  [hoverk_b]
    type = DerivativeParsedMaterial
    material_property_names = 'hb(etaa0,etab0) Vm kb'
    property_name = hoverk_b
    expression = 'hb / Vm / kb'
  []
  [caeq_mask]
    type = DerivativeParsedMaterial
    property_name = caeq_mask
    material_property_names = 'ha(etaa0,etab0) caeq(etaa0,etab0)'
    coupled_variables = 'etaa0 etab0'
    expression = 'ha * caeq'
  []
  [cbeq_mask]
    type = DerivativeParsedMaterial
    property_name = cbeq_mask
    material_property_names = 'hb(etaa0,etab0) cbeq(etaa0,etab0)'
    coupled_variables = 'etaa0 etab0'
    expression = 'hb * cbeq'
  []
[]

[Postprocessors]
  [c_total]
    type = ElementIntegralVariablePostprocessor
    variable = c
  []
  [etaa0_total]
    type = ElementIntegralVariablePostprocessor
    variable = etaa0
  []
  [etab0_total]
    type = ElementIntegralVariablePostprocessor
    variable = etab0
  []
  [hgb_total]
    type = ElementIntegralMaterialProperty
    mat_prop = hgb
  []
  [cmat_total]
    type = ElementIntegralMaterialProperty
    mat_prop = cmat
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
  petsc_options_iname = -pc_type
  petsc_options_value = asm
  l_max_its = 15
  l_tol = 1.0e-3
  nl_rel_tol = 1.0e-8
  start_time = 0.0
  num_steps = 500
  nl_abs_tol = 1e-10
  dt = 1.0
[]

[Outputs]
  csv = true
  exodus = true
  checkpoint = false
[]
