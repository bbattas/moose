##############################################################################
# File: 01_bicr_base.i
# File Location: /examples/agg/01_initial_testing/12_bicrystal/01_bicr_base
# Created Date: Monday May 12th 2025
# Author: Battas,Brandon Scott (bbattas@ufl.edu)
# -----
# Last Modified: Friday May 16th 2025
# Modified By: Brandon Battas
# -----
# Description:
#  circular bicrystal input to test iso and aniso things
#  Base input (anisotropic)
#
#
##############################################################################

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 60
  ny = 60
  nz = 0
  xmin = 0
  xmax = 16 #1000
  ymin = 0
  ymax = 16 #1000
  zmin = 0
  zmax = 0
  elem_type = QUAD4
  parallel_type = DISTRIBUTED
  uniform_refine = 0
  second_order = true
[]

[GlobalParams]
  op_num = 2
  var_name_base = 'gr'
  # profile = TANH
  # int_width = int_width
[]

[UserObjects]
  [grain_tracker]
    type = GrainTracker
    # variable = eta
    threshold = 0.01 #0.2
    connecting_threshold = 0.08 #0.08
    compute_halo_maps = true
    compute_var_to_feature_map = true
    execute_on = 'initial timestep_end'
  []
[]

[ICs]
  # [PolycrystalICs]
  #   [BicrystalCircleGrainIC]
  #     radius = 4
  #     x = 8
  #     y = 8
  #   []
  # []
  [gr0_IC]
    type = SmoothCircleIC
    invalue = 1
    outvalue = 0
    radius = 4
    variable = gr0
    x1 = 8
    y1 = 8
    int_width = 3.2
  []
  [gr1_IC]
    type = SmoothCircleIC
    invalue = 0
    outvalue = 1
    radius = 4
    variable = gr1
    x1 = 8
    y1 = 8
    int_width = 3.2
  []
  # [gr0_IC]
  #   type = BoundingBoxIC
  #   inside = 1
  #   outside = 0
  #   variable = gr0
  #   x1 = 8
  #   y1 = 8
  #   x2 = 24
  #   y2 = 24
  #   int_width = 1.5
  # []
  # [gr1_IC]
  #   type = BoundingBoxIC
  #   inside = 0
  #   outside = 1
  #   variable = gr1
  #   x1 = 8
  #   y1 = 8
  #   x2 = 24
  #   y2 = 24
  #   int_width = 1.5
  # []
[]

[BCs]
  [Periodic]
    [All]
      auto_direction = 'x y'
    []
  []
[]

[AuxVariables]
  [bnds]
    order = FIRST
    family = LAGRANGE
  []
  [unique_grains]
    order = CONSTANT
    family = MONOMIAL
  []
  [var_indices]
    order = CONSTANT
    family = MONOMIAL
  []
  [ghost_regions]
    order = CONSTANT
    family = MONOMIAL
  []
  # Gradient components
  [gr0_x]
    order = CONSTANT
    family = MONOMIAL
  []
  [gr0_y]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Modules]
  [PhaseField]
    [GrainGrowth]
      mobility = L_aniso
      kappa = kappa
      order = SECOND
      family = LAGRANGE
    []
  []
[]

[Kernels]
  [gr0_ACIaniso]
    type = ACInterfaceAnisoGamma
    variable = gr0
    # v = 'gr1 gr2 gr3'
    coupled_variables = 'gr1'
    gamma_name = gamma_asymm
    dgamma_dgradop_name = dgammadgrad_eta_0
    d2gamma_dgradop2_name = d2gammadgrad_eta2_0
    mob_name = L_aniso
    dL_dgradop_name = dLdgrad_eta_0
    d2L_dgradop2_name = d2Ldgrad_eta2_0
    variable_L = false
    skip_off = true
  []
  [gr1_ACIaniso]
    type = ACInterfaceAnisoGamma
    variable = gr1
    # v = 'gr0 gr2 gr3'
    coupled_variables = 'gr0'
    gamma_name = gamma_asymm
    dgamma_dgradop_name = dgammadgrad_eta_1
    d2gamma_dgradop2_name = d2gammadgrad_eta2_1
    mob_name = L_aniso
    dL_dgradop_name = dLdgrad_eta_1
    d2L_dgradop2_name = d2Ldgrad_eta2_1
    variable_L = false
    skip_off = true
  []
[]

[AuxKernels]
  [BndsCalc]
    type = BndsCalcAux
    variable = bnds
    execute_on = 'initial timestep_end'
  []
  [unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
    execute_on = timestep_begin
  []
  [var_indices]
    type = FeatureFloodCountAux
    variable = var_indices
    flood_counter = grain_tracker
    field_display = VARIABLE_COLORING
    execute_on = timestep_begin
  []
  [ghosted_entities]
    type = FeatureFloodCountAux
    variable = ghost_regions
    flood_counter = grain_tracker
    field_display = GHOSTED_ENTITIES
    execute_on = 'initial timestep_end'
  []
  # gr0 gradient pieces
  [gr0_gradx]
    type = VariableGradientComponent
    variable = gr0_x
    component = x
    gradient_variable = gr0
  []
  [gr0_grady]
    type = VariableGradientComponent
    variable = gr0_y
    component = y
    gradient_variable = gr0
  []
[]

[Materials]
  [constants]
    type = GenericConstantMaterial
    prop_names = 'L0 kappa const_m gamma_iso iw_iso gbe_iso'
    prop_values = '1.0 0.3 0.9375   1.5       1.55    0.25'
  []
  # [iso_constants]
  #   type = GenericConstantMaterial
  #   prop_names = 'L0 kappa const_m gamma_asymm gamma_iso iw_iso gbe_iso mu'
  #   prop_values = '1.0 0.3 0.9375   1.5           1.5     1.55   0.25  0.85'
  # []
  [incl_test01]
    type = GGInclinationMaterial
    grain_tracker = grain_tracker
    gb_energy_input = gbe_iso
    kappa = 0.3 #kappa
    free_energy_m = 0.9375 #const_m
    L0 = L0
    gamma0 = gamma_iso
    theta_prefactor = 2
    continuous = false
    # Output Names
    inclination_name = inclination_mat
    L_name = L_aniso
    gamma_name = gamma_asymm
    mu_name = mu
    output_properties = 'gamma_asymm L_aniso mu gb_energy int_width inclination_mat inclination_distance'
    outputs = 'exodus'
  []
  [sumgr]
    type = ParsedMaterial
    property_name = sumgr
    coupled_variables = 'gr0 gr1'
    expression = 'gr0 + gr1'
    outputs = 'exodus'
  []
  # Variable Gradients
  [ng_gr0]
    type = VariableGradientMaterial
    prop = ng_gr0
    variable = gr0
    outputs = exodus
  []
  [ng_gr1]
    type = VariableGradientMaterial
    prop = ng_gr1
    variable = gr1
  []
  # Free energies
  [fe_bulk_manual]
    type = ParsedMaterial
    property_name = fe_bulk_manual
    coupled_variables = 'gr0 gr1'
    material_property_names = 'gamma_asymm'
    expression = 'gmeta:= gr0^2 * gr1^2;
    etaover:= 0.25*gr0^4 - 0.5*gr0^2 + 0.25*gr1^4 - 0.5*gr1^2;
    etaover + gamma_asymm * gmeta + 0.25'
    outputs = 'exodus'
  []
  [fe_grad_manual]
    type = ParsedMaterial
    property_name = fe_grad_manual
    coupled_variables = 'gr0 gr1'
    material_property_names = 'kappa ng_gr0 ng_gr1'
    expression = '0.5 * kappa * (ng_gr0^2 + ng_gr1^2)'
    outputs = 'exodus'
  []
  [fe_tot_manual]
    type = ParsedMaterial
    property_name = fe_tot_manual
    coupled_variables = 'gr0 gr1'
    material_property_names = 'fe_bulk_manual fe_grad_manual'
    expression = 'fe_bulk_manual + fe_grad_manual'
    outputs = 'exodus'
  []
  # [gr0_grad]
  #   type = VariableGradientMaterial
  #   property_name = gr0_grad


[]

[Postprocessors]
  [runtime]
    type = PerfGraphData
    section_name = "Root"
    data_type = TOTAL
  []
  [nl_its]
    type = NumNonlinearIterations
  []
  [l_its]
    type = NumLinearIterations
  []
  [max_mpi_memory]
    type = MemoryUsage
    value_type = max_process
    report_peak_value = True
    mem_units = MEGABYTES
    execute_on = 'NONLINEAR TIMESTEP_END'
  []
  [ctr_grain_area]
    type = ElementIntegralVariablePostprocessor
    variable = gr1
    execute_on = 'initial timestep_end'
  []
  [tot_grain_area]
    type = ElementIntegralMaterialProperty
    mat_prop = sumgr
  []
  [gbe_total]
    type = ElementIntegralMaterialProperty
    mat_prop = gb_energy
  []
  [mu_total]
    type = ElementIntegralMaterialProperty
    mat_prop = mu
  []
  [fe_bulk]
    type = ElementIntegralMaterialProperty
    mat_prop = fe_bulk_manual
  []
  [fe_grad]
    type = ElementIntegralMaterialProperty
    mat_prop = fe_grad_manual
  []
  [fe_tot]
    type = ElementIntegralMaterialProperty
    mat_prop = fe_tot_manual
  []
  # Variable Residuals
  [R_gr0]
    type = VariableResidual
    variable = gr0
  []
  [R_gr1]
    type = VariableResidual
    variable = gr1
  []
[]

# [VectorPostprocessors]
#   [0_line]
#     type = LineValueSampler
#     variable = gr0
#     start_point = '0 8 0'
#     end_point = '16 8 0'
#     sort_by = x
#     num_points = 61
#     outputs = csv
#   []
# []

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = 'PJFNK'

  # petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  # petsc_options_value = 'hypre boomeramg 31'
  petsc_options_iname = '-pc_type -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  petsc_options_value = 'asm        preonly       lu           2'
  # petsc_options_iname = '-pc_type -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  # petsc_options_value = 'asm        preonly       ilu           2'

  nl_max_its = 20
  l_max_its = 60
  l_tol = 1e-05
  nl_rel_tol = 1e-8 #default is 1e-8
  # nl_abs_tol = 1e-14 #only needed when near equilibrium or veeeery small dt

  start_time = 0.0
  end_time = 30
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 6
    linear_iteration_ratio = 1e5
    dt = 1
  []
[]

[Outputs]
  exodus = true
  console = true
  csv = true
[]
