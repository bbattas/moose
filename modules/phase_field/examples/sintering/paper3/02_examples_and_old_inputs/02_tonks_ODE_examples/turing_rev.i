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
      [ci]
        solve_type = FORWARD_SPLIT
        kappa = kappa
        free_energy = F
        mobility = Mi
        coupled_variables = cv
      []
    []
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
  [ci_IC]
    type = RandomIC
    seed = 694923043
    max = 0.0101
    min = 0.0099
    variable = ci
  []
[]

[Kernels]
  [ci_gen]
    type = MaskedBodyForce
    value = 1.25e-3
    variable = ci
    mask = gen_mask
    coupled_variables = cv
  []
  [ci_recombination]
    type = MatReaction
    variable = ci
    args = cv
    mob_name = R_for_ci
  []
  [cv_gen]
    type = MaskedBodyForce
    value = 1.25e-3
    variable = cv
    mask = gen_mask
    coupled_variables = ci
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
    prop_names = 'kappa Mi Mv R'
    prop_values = '1 5 1 -0.5'
  []
  [R_for_ci]
    type = DerivativeParsedMaterial
    material_property_names = 'R'
    coupled_variables = cv
    property_name = R_for_ci
    function = 'R * cv'
  []
  [R_for_cv]
    type = DerivativeParsedMaterial
    material_property_names = 'R'
    coupled_variables = ci
    property_name = R_for_cv
    function = 'R * ci'
  []
  [fbulk]
    type = DerivativeParsedMaterial
    f_name = F
    function = 'ci^2 * (1 - ci)^2 + cv^2 * (1 - cv)^2'
    coupled_variables = 'ci cv'
  []
  [gen_mask]
    type = DerivativeParsedMaterial
    f_name = gen_mask
    function = '(1 - cv)^2'
    coupled_variables = 'cv'
  []
[]

[BCs]
  [Periodic]
    [All]
      auto_direction = 'x y'
    []
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
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
  [max_cv]
    type = ElementExtremeValue
    variable = cv
    value_type = max
  []
  [min_cv]
    type = ElementExtremeValue
    variable = cv
    value_type = min
  []
  [max_ci]
    type = ElementExtremeValue
    variable = ci
    value_type = max
  []
  [min_ci]
    type = ElementExtremeValue
    variable = ci
    value_type = min
  []
  [av_cv]
    type = ElementAverageValue
    variable = cv
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
  file_base = Dratio_5_rndIC
  # checkpoint = true #Checkpoint files make it so you can restart simulations
  interval = 3 #Number of simulations between data output
[]
