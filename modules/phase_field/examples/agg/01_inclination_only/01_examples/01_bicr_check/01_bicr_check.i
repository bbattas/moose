##############################################################################
# File: 01_bicr_check.i
# File Location: /examples/agg/01_inclination_only/01_examples/01_bicr_check
# Created Date: Wednesday March 25th 2026
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Wednesday March 25th 2026
# Modified By: Brandon Battas
# -----
# Description:
#  Test input to make sure the AGG Inclination parts are working in the new branch
#
#
#
##############################################################################

amag = 0.0

[Mesh]
  [gmg]
    type = DistributedRectilinearMeshGenerator
    dim = 2
    nx = 200
    ny = 200
    xmin = 0
    xmax = 200
    ymin = 0
    ymax = 200
  []
  parallel_type = DISTRIBUTED
[]

[GlobalParams]
  op_num = 2
  var_name_base = gr
[]

[Variables]
  [PolycrystalVariables]
    # order = SECOND
  []
[]

[ICs]
  [gr0_IC]
    type = SmoothSuperellipsoidIC
    invalue = 1
    outvalue = 0
    variable = gr0
    a = 80
    b = 80
    n = 2
    x1 = 100
    y1 = 100
    int_width = 6
  []
  [gr1_IC]
    type = SmoothSuperellipsoidIC
    invalue = 0
    outvalue = 1
    variable = gr1
    a = 80
    b = 80
    n = 2
    x1 = 100
    y1 = 100
    int_width = 6
  []
[]

[UserObjects]
  [grain_tracker]
    type = GrainTracker
    threshold = 0.01
    connecting_threshold = 0.01
    compute_var_to_feature_map = true
    compute_halo_maps = false # Only necessary for displaying HALOS
    execute_on = 'initial timestep_end'
    halo_level = 6
    # use_less_than_threshold_comparison = true
    # polycrystal_ic_uo = ebsd
  []
  [term]
    type = Terminator
    expression = 'grain_tracker < 2'
  []
[]

[AuxVariables]
  [bnds]
  []
  [unique_grains]
    order = CONSTANT
    family = MONOMIAL
  []
  [contour]
  []
[]

[Kernels]
  [PolycrystalKernel]
    variable_mobility = false
  []
  [PolycrystalInclinationKernel]
    variable_mobility = false
    hgb_mask = hgb
  []
[]

[AuxKernels]
  [bnds_aux]
    type = BndsCalcAux
    variable = bnds
    execute_on = 'initial timestep_end'
  []
  [unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
    execute_on = 'initial timestep_end'
  []
  [contour]
    type = ParsedAux
    variable = contour
    coupled_variables = 'gr0 gr1'
    expression = 'gr0 * gr0 / ((gr0 * gr0) + (gr1 * gr1))'
    execute_on = 'initial timestep_end'
  []
[]

[BCs]
  [Periodic]
    [All]
      auto_direction = 'x y'
    []
  []
[]

[Materials]
  [constants]
    type = GenericConstantMaterial
    prop_names = 'L0         gamma_iso kappa_op  iw_iso gbe_iso'
    prop_values = '1.15382e-6   1.5    2.07337e7   6   4.60748e6'
  []
  [GBInc_cos]
    type = GBInclination
    gb_id_method = GRAINTRACKER
    grain_tracker = grain_tracker
    angular_func = ATAN_2D
    intol = 0
    altol = 5
    limit_umag = true
    # Inclination function
    inc_func = COS
    ifunc_a = ${amag}
    ifunc_b = 2
    ifunc_c = 0
    ifunc_d = 0
    # General
    combine_gb_form = avg
    aniso_L = false
    aniso_gbmob = false
    # Other Properties
    gb_energy_iso_name = gbe_iso
    kappa = kappa_op
    free_energy_m = 5.521269e6
    # output_properties = 'gtnum int_width gamma_asymm L'
    # outputs = 'exodus'
  []
  [hgb]
    type = SwitchingFunctionGBMaterial
    h_name = hgb
    hgb_threshold = 0
    func_type = COMBINED
    # output_properties = 'hgb'
    # outputs = 'exodus'
  []
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
  [timestep]
    type = TimestepSize
  []
  [n_elements]
    type = NumElements
    execute_on = 'initial timestep_end'
  []
  [n_nodes]
    type = NumNodes
    execute_on = 'initial timestep_end'
  []
  [DOFs]
    type = NumDOFs
    execute_on = 'initial timestep_end'
  []
  [avg_grain_volumes]
    type = AverageGrainVolume
    feature_counter = grain_tracker
    execute_on = 'initial timestep_end'
  []
  [max_mpi_memory]
    type = MemoryUsage
    value_type = max_process
    report_peak_value = True
    mem_units = MEGABYTES
    execute_on = 'NONLINEAR TIMESTEP_END'
  []
  [gr0_total]
    type = ElementIntegralVariablePostprocessor
    variable = gr0
  []
  [gr1_total]
    type = ElementIntegralVariablePostprocessor
    variable = gr1
  []

  [gamma_min]
    type = ElementExtremeMaterialProperty
    mat_prop = gamma_asymm
    value_type = MIN
  []
  [gamma_max]
    type = ElementExtremeMaterialProperty
    mat_prop = gamma_asymm
    value_type = MAX
  []
  [iw_min]
    type = ElementExtremeMaterialProperty
    mat_prop = int_width
    value_type = MIN
  []
  [iw_max]
    type = ElementExtremeMaterialProperty
    mat_prop = int_width
    value_type = MAX
  []
[]

# [VectorPostprocessors]
#   [horizontal]
#     type = LineValueSampler
#     end_point = '200 100 0'
#     num_points = 200
#     sort_by = X
#     start_point = '0 100 0'
#     variable = 'gr0 gr1 contour'
#     outputs = 'vpp'
#   []
#   [vertical]
#     type = LineValueSampler
#     end_point = '100 200 0'
#     num_points = 200
#     sort_by = Y
#     start_point = '100 0 0'
#     variable = 'gr0 gr1 contour'
#     outputs = 'vpp'
#   []
# []

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient # Type of executioner, here it is transient with an adaptive time step
  scheme = bdf2 # Type of time integration (2nd order backward euler), defaults to 1st order backward euler

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  # Uses newton iteration to solve the problem.
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -mat_mffd_type'
  petsc_options_value = 'hypre boomeramg 101 ds'
  # petsc_options_iname = '-pc_type -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  # petsc_options_value = 'asm        preonly       lu           2'

  l_max_its = 60 # Max number of linear iterations
  l_tol = 1e-4 # Relative tolerance for linear solves
  nl_max_its = 12 # Max number of nonlinear iterations
  nl_rel_tol = 1e-8 #1e-10 # Relative tolerance for nonlienar solves
  # nl_abs_tol = 1e-10

  start_time = 0.0
  end_time = 60
  num_steps = 3

  automatic_scaling = true
  compute_scaling_once = false
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.5
    optimal_iterations = 6
    linear_iteration_ratio = 1e5
  []
[]

[Outputs]
  exodus = true # Exodus file will be outputted
  #nemesis = true
  console = true
  csv = true
  checkpoint = false
  # [vpp]
  #   type = CSV
  # []
  # [checkpoint]
  #   type = Checkpoint
  #   num_files = 5
  # []
  # [console]
  #   type = Console
  #   max_rows = 20 # Will print the 20 most recent postprocessor values to the screen
  # []
  # [pgraph]
  #   type = PerfGraphOutput
  #   execute_on = 'initial final' # Default is "final"
  #   level = 2 # Default is 1
  #   heaviest_branch = true # Default is false
  #   heaviest_sections = 7 # Default is 0
  # []
  # file_base = 01_a${amag}_circle/01_a${amag}_circle_out
[]
