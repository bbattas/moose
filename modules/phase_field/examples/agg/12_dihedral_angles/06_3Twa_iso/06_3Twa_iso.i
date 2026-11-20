##############################################################################
# File: 06_3Twa_iso.i
# File Location: /examples/agg/12_dihedral_angles/06_3Twa_iso
# Created Date: Friday November 20th 2026
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Friday November 20th 2026
# Modified By: Brandon Battas
# -----
# Description:
#  weighted average version for tjs
#
#
#
##############################################################################

[Mesh]
  [gmg]
    type = DistributedRectilinearMeshGenerator
    dim = 2
    nx = 120
    ny = 200
    xmin = 0
    xmax = 60
    ymin = 0
    ymax = 100
  []
  parallel_type = DISTRIBUTED
  second_order = false
  uniform_refine = 0
[]

[GlobalParams]
  op_num = 3
  var_name_base = gr
[]

[Variables]
  [PolycrystalVariables]
  []
[]

[ICs]
  [gr0_IC]
    type = BoundingBoxIC
    variable = gr0
    x1 = -10
    y1 = -10
    x2 = 30
    y2 = 50
    inside = 1
    outside = 0
  []
  [gr1_IC]
    type = BoundingBoxIC
    variable = gr1
    x1 = 30
    y1 = -10
    x2 = 90
    y2 = 50
    inside = 1
    outside = 0
  []
  [gr2_IC]
    type = BoundingBoxIC
    variable = gr2
    x1 = -10
    y1 = 50
    x2 = 90
    y2 = 170
    inside = 1
    outside = 0
  []
[]

[UserObjects]
  # [euler_file]
  #   type = EulerAngleTxtFileReader
  #   file_name = '../00_sub/bicr_miso_low.txt'
  # []
  [grain_tracker]
    type = GrainTracker
    # variable = 'gr0 gr1 gr2 gr3 gr4'
    threshold = 0.01
    connecting_threshold = 0.01
    compute_var_to_feature_map = true
    compute_halo_maps = false
    execute_on = 'initial timestep_end'
    halo_level = 6
    # use_less_than_threshold_comparison = true
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
  # [halos]
  #   order = CONSTANT
  #   family = MONOMIAL
  # []
  [c0]
  []
  [c1]
  []
  [c2]
  []
[]

[Kernels]
  [PolycrystalKernel]
    variable_mobility = false
    # order = SECOND
  []
  # [PolycrystalInclinationKernel]
  #   variable_mobility = false
  #   hgb_mask = hgb
  # []
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
  # [halos]
  #   type = FeatureFloodCountAux
  #   variable = halos
  #   flood_counter = grain_tracker
  #   field_display = HALOS
  #   execute_on = 'initial timestep_end'
  # []
  [c0]
    type = ParsedAux
    variable = c0
    coupled_variables = 'gr0 gr1 gr2'
    expression = 'gr0 * gr0 / ((gr0 * gr0) + (gr1 * gr1) + (gr2 * gr2))'
    execute_on = 'initial timestep_end'
  []
  [c1]
    type = ParsedAux
    variable = c1
    coupled_variables = 'gr0 gr1 gr2'
    expression = 'gr1 * gr1 / ((gr0 * gr0) + (gr1 * gr1) + (gr2 * gr2))'
    execute_on = 'initial timestep_end'
  []
  [c2]
    type = ParsedAux
    variable = c2
    coupled_variables = 'gr0 gr1 gr2'
    expression = 'gr2 * gr2 / ((gr0 * gr0) + (gr1 * gr1) + (gr2 * gr2))'
    execute_on = 'initial timestep_end'
  []
[]

# [BCs]
#   [Periodic]
#     [All]
#       auto_direction = 'x y'
#     []
#   []
# []

[Materials]
  [constants]
    type = GenericConstantMaterial
    prop_names = 'L0             kappa_op   mu      gbe_iso   '
    prop_values = '2.7815e-6    2.100e07    1.305e7    1.125e7   '
    # prop_values = '1.15382e-6    2.07337e7  5.521269e6  4.60748e6'
  []
  [hgb]
    type = SwitchingFunctionGBMaterial
    h_name = hgb
    hgb_threshold = 0
    func_type = COMBINED
    output_properties = 'hgb'
    outputs = 'exodus'
  []
  [new_aniso_mat]
    type = GBCombinedAnisotropyMaterial
    gb_id_method = graintracker
    grain_tracker = grain_tracker
    tj_mode = weighted
    gb_mode = ISO
    bulk_scalar = 0.75
    iso_gbe = 0.75
    alpha_tol = 0
    hgbalpha_tol = 0
    # ifunc_a = 0.1
    kappa_name = kappa_op
    gbe_iso_name = gbe_iso
    enable_misorientation = false
    # euler_angle_provider = euler_file
    w_inc = 1
    w_miso = 1
    output_properties = 'gt_num L gamma_asymm int_width gbe_norm'
    outputs = 'exodus'
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
  # IW and Gamma_asymm
  [avg_iw]
    type = ElementAverageMaterialProperty
    mat_prop = int_width
  []
  [avg_gamma]
    type = ElementAverageMaterialProperty
    mat_prop = gamma_asymm
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
  # num_steps = 0
  end_time = 10
  # dt = 0.01
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01 #0.001
    # cutback_factor = 0.9
    # growth_factor = 1.1
    optimal_iterations = 6
    linear_iteration_ratio = 1e5
  []

  automatic_scaling = true
  compute_scaling_once = true #SET TO TRUE TO SPEED UP
[]

[Outputs]
  # exodus = true # Exodus file will be outputted
  #nemesis = true
  console = true
  # csv = true
  [exodus]
    type = Exodus
    sync_only = true
    sync_times = '1 2 3 4 5 6 7 8 9 10'
  []
  [csv]
    type = CSV
    sync_only = true
    sync_times = '1 2 3 4 5 6 7 8 9 10'
  []
  #
  print_linear_residuals = false
  print_nonlinear_residuals = false
  #
  checkpoint = false
  # [checkpoint]
  #   type = Checkpoint
  #   num_files = 3
  # []
  # [pgraph]
  #   type = PerfGraphOutput
  #   execute_on = 'initial final' # Default is "final"
  #   level = 2 # Default is 1
  #   heaviest_branch = true # Default is false
  #   heaviest_sections = 7 # Default is 0
  # []
  # file_base = 06_aniso_iso_iso${iso}_bulk${bulk}_out
[]
