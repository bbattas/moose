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
      []
    []
  []
[]

[AuxVariables]
  [cv_av]
  []
  [ci]
  []
[]

[AuxKernels]
  [cv_av_kernel]
    type = FunctionAux
    variable = cv_av
    function = cv_av_function
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [ci_kernel]
    type = ParsedAux
    variable = ci
    coupled_variables = cv_av
    expression = '1.25e-3/(0.5*cv_av)'
    # expression = '1.25e-3/(0.5*0.25)'
  []
[]

[Functions]
  [cv_av_function]
    type = ParsedFunction
    symbol_names = 'av_cv'
    symbol_values = 'av_cv'
    expression = 'av_cv'
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

[Materials]
  [constants]
    type = GenericConstantMaterial
    prop_names = 'kappa Mv R    c'
    prop_values = '1    1  -0.5 1.25e-3'
  []
  [R_for_cv]
    type = ParsedMaterial
    material_property_names = 'R'
    coupled_variables = ci
    property_name = R_for_cv
    function = 'R * ci'
  []
  [fbulk]
    type = DerivativeParsedMaterial
    f_name = F
    function = 'ci^2 * (1 - ci)^2 + cv^2 * (1 - cv)^2'
    coupled_variables = 'cv ci'
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
  # [max_cv]
  #   type = ElementExtremeValue
  #   variable = cv
  #   value_type = max
  # []
  # [min_cv]
  #   type = ElementExtremeValue
  #   variable = cv
  #   value_type = min
  # []
  [av_cv]
    type = ElementAverageValue
    variable = cv
    execute_on = 'initial TIMESTEP_BEGIN'
  []
  [run_time]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
  []
[]

[VectorPostprocessors]
  [void_volumes]
    type = FeatureVolumeVectorPostprocessor
    flood_counter = num_voids
    output_centroids = true
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = NEWTON

  #End time and time stepper
  # end_time = 1000
  end_time = 40000

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
    optimal_iterations = 6
    linear_iteration_ratio = 1000
  []

  #Numerical solution options
  petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap'
  petsc_options_value = 'asm lu 2'

  nl_max_its = 8
  l_tol = 1e-4
  l_max_its = 40
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-9
  automatic_scaling = true
  compute_scaling_once = false
[]

[Outputs]
  perf_graph = true
  exodus = true
  csv = true
  file_base = 1PDE_rndIC
  # checkpoint = true #Checkpoint files make it so you can restart simulations
  interval = 3 #Number of simulations between data output
[]
