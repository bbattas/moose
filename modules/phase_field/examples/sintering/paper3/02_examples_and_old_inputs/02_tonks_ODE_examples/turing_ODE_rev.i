[Mesh]
  type = GeneratedMesh
  nx = 450
  ny = 450
  xmax = 300
  ymax = 300
  dim = 2
[]

[Modules]
  [PhaseField]
    [Conserved]
      [cv]
        solve_type = FORWARD_SPLIT
        kappa = kappa
        free_energy = F
        mobility = Mv
        coupled_variables = ci
      []
    []
  []
[]

[Variables]
  [ci_scalar]
    family = SCALAR
    order = FIRST
    initial_condition = 0.01
  []
[]

[ICs]
  [cv_IC]
    type = RandomIC
    seed = 54932402
    max = 0.2525
    min = 0.2475
    variable = cv
  []
[]

[Kernels]
  [cv_gen]
    type = BodyForce
    value = 1.25e-3
    variable = cv
  []
  [cv_recombination]
    type = MatReaction
    variable = cv
    args = ci
    mob_name = R_for_cv
  []
[]

[Functions]
  [ci_func]
    type = ParsedFunction
    symbol_names = 'ci_scalar_reporter'
    symbol_values = 'ci_scalar_reporter'
    expression = 'ci_scalar_reporter'
    execute_on = 'INITIAL LINEAR'
  []
[]

[ScalarKernels]
  [ci_scalar_dot]
    type = ODETimeDerivative
    variable = ci_scalar
  []
  [ci_rest]
    type = ParsedODEKernel
    # function = '1.25e-3 - 0.5'
    expression = '-(c-cv_av*ci_scalar*a)'
    variable = ci_scalar
    postprocessors = cv_av
    constant_names = 'c a'
    constant_expressions = '1.25e-3 0.5'
  []
[]

[AuxVariables]
  [ci]
  []
[]

[AuxKernels]
  [ci]
    type = FunctionAux
    variable = ci
    function = ci_func
    execute_on = 'INITIAL LINEAR TIMESTEP_END'
  []
[]

[Materials]
  [constants]
    type = GenericConstantMaterial
    prop_names = 'kappa Mv R'
    prop_values = '1    1  -0.5'
  []
  [R_for_cv]
    type = DerivativeParsedMaterial
    material_property_names = 'R'
    coupled_variables = ci
    property_name = R_for_cv
    expression = 'R * ci'
  []
  [fbulk]
    type = DerivativeParsedMaterial
    property_name = F
    expression = 'ci^2 * (1 - ci)^2 + cv^2 * (1 - cv)^2'
    coupled_variables = 'ci cv'
  []
[]

[BCs]
  [Periodic]
    [All]
      auto_direction = 'x y'
    []
  []
[]

[Postprocessors]
  # [cv_total]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = cv
  #   execute_on = 'initial timestep_end'
  # []
  [cv_av]
    type = ElementAverageValue
    variable = cv
    execute_on = 'TIMESTEP_END'
  []
  [dt]
    type = TimestepSize
  []
  [num_voids]
    type = FeatureFloodCount
    variable = cv
    threshold = 0.5
    compute_var_to_feature_map = true
  []
  [mesh_volume]
    type = VolumePostprocessor
    execute_on = 'initial'
  []
  [porosity]
    type = FeatureVolumeFraction
    feature_volumes = void_volumes
    mesh_volume = mesh_volume
  []
  [run_time]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
  []
  [ci_scalar_reporter]
    type = ScalarVariable
    variable = ci_scalar
    execute_on = 'INITIAL LINEAR TIMESTEP_END'
  []
[]

[VectorPostprocessors]
  [void_volumes]
    type = FeatureVolumeVectorPostprocessor
    flood_counter = num_voids
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = NEWTON
  end_time = 40000
  # end_time = 1000
  petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap'
  petsc_options_value = 'asm lu 2'

  nl_max_its = 8
  l_tol = 1e-4
  l_max_its = 40
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-9

  automatic_scaling = true
  compute_scaling_once = false
  line_search = none

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
    optimal_iterations = 6
    linear_iteration_ratio = 300
  []
[]

[Outputs]
  file_base = ODE_rndIC
  perf_graph = true
  exodus = true
  csv = true
  # checkpoint = true
  interval = 3
[]
