##############################################################################
# File: 06_full_tex.i
# File Location: /examples/agg/08_combined/15_1000gr_512mesh_8s/06_full_tex
# Created Date: Thursday August 6th 2026
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Thursday August 6th 2026
# Modified By: Brandon Battas
# -----
# Description:
#
#
#
#
##############################################################################

randseed = 20

[Mesh]
  [gmg]
    type = DistributedRectilinearMeshGenerator
    dim = 2
    nx = 768
    ny = 768
    xmin = 0
    xmax = 512
    ymin = 0
    ymax = 512
  []
  parallel_type = DISTRIBUTED
  second_order = false
  uniform_refine = 0
[]

[GlobalParams]
  op_num = 25
  var_name_base = gr
[]

[Variables]
  [PolycrystalVariables]
  []
[]

[ICs]
  [PolycrystalICs]
    [PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
    []
  []
[]

[UserObjects]
  [euler_file]
    type = EulerAngleTxtFileReader
    file_name = '../../00_sub/euler_1k_textured.txt'
  []
  [voronoi]
    type = PolycrystalVoronoi
    coloring_algorithm = bt #jp #bt
    grain_num = 1000
    rand_seed = ${randseed}
    int_width = 3 #6
  []
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
    expression = 'grain_tracker < 50'
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
[]

[Kernels]
  [PolycrystalKernel]
    variable_mobility = false
    # order = SECOND
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
  # [halos]
  #   type = FeatureFloodCountAux
  #   variable = halos
  #   flood_counter = grain_tracker
  #   field_display = HALOS
  #   execute_on = 'initial timestep_end'
  # []
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
    prop_names = 'L0             kappa_op   mu      gbe_iso   '
    prop_values = '2.7815e-6    2.100e07    1.305e7    1.125e7   '
  []
  [hgb]
    type = SwitchingFunctionGBMaterial
    h_name = hgb
    hgb_threshold = 0
    func_type = COMBINED
    # output_properties = 'hgb'
    # outputs = 'exodus'
  []
  [new_aniso_mat]
    type = GBCombinedAnisotropyMaterial
    gb_id_method = graintracker
    grain_tracker = grain_tracker
    gb_mode = FULL
    bulk_scalar = 0.75
    alpha_tol = 8 #10
    hgbalpha_tol = 5
    # ifunc_a = 0.1
    kappa_name = kappa_op
    gbe_iso_name = gbe_iso
    enable_misorientation = true
    euler_angle_provider = euler_file
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
  nl_rel_tol = 1e-6#1e-8 #1e-10 # Relative tolerance for nonlienar solves
  # nl_abs_tol = 1e-10

  start_time = 0.0
  # num_steps = 3
  end_time = 8.04
  # dt = 0.1
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01 #0.001
    # cutback_factor = 0.9
    # growth_factor = 1.1
    optimal_iterations = 6
    linear_iteration_ratio = 1e5
  []

  automatic_scaling = true
  compute_scaling_once = false
  line_search = 'none'
[]

[Outputs]
  # exodus = true # Exodus file will be outputted
  #nemesis = true
  console = true
  # csv = true
  [exodus]
    type = Exodus
    sync_only = true
    sync_times = '0 0.067 0.134 0.201 0.268 0.335 0.402 0.469 0.536 0.603 0.67 0.737 0.804 0.871 0.938 1.005
    1.072 1.139 1.206 1.273 1.34 1.407 1.474 1.541 1.608 1.675 1.742 1.809 1.876 1.943 2.01 2.077
    2.144 2.211 2.278 2.345 2.412 2.479 2.546 2.613 2.68 2.747 2.814 2.881 2.948 3.015 3.082 3.149
    3.216 3.283 3.35 3.417 3.484 3.551 3.618 3.685 3.752 3.819 3.886 3.953 4.02 4.087 4.154 4.221
    4.288 4.355 4.422 4.489 4.556 4.623 4.69 4.757 4.824 4.891 4.958 5.025 5.092 5.159 5.226 5.293
    5.36 5.427 5.494 5.561 5.628 5.695 5.762 5.829 5.896 5.963 6.03 6.097 6.164 6.231 6.298 6.365
    6.432 6.499 6.566 6.633 6.7 6.767 6.834 6.901 6.968 7.035 7.102 7.169 7.236 7.303 7.37 7.437
    7.504 7.571 7.638 7.705 7.772 7.839 7.906 7.973 8.04'
  []
  [csv]
    type = CSV
    sync_only = true
    sync_times = '0 0.067 0.134 0.201 0.268 0.335 0.402 0.469 0.536 0.603 0.67 0.737 0.804 0.871 0.938 1.005
    1.072 1.139 1.206 1.273 1.34 1.407 1.474 1.541 1.608 1.675 1.742 1.809 1.876 1.943 2.01 2.077
    2.144 2.211 2.278 2.345 2.412 2.479 2.546 2.613 2.68 2.747 2.814 2.881 2.948 3.015 3.082 3.149
    3.216 3.283 3.35 3.417 3.484 3.551 3.618 3.685 3.752 3.819 3.886 3.953 4.02 4.087 4.154 4.221
    4.288 4.355 4.422 4.489 4.556 4.623 4.69 4.757 4.824 4.891 4.958 5.025 5.092 5.159 5.226 5.293
    5.36 5.427 5.494 5.561 5.628 5.695 5.762 5.829 5.896 5.963 6.03 6.097 6.164 6.231 6.298 6.365
    6.432 6.499 6.566 6.633 6.7 6.767 6.834 6.901 6.968 7.035 7.102 7.169 7.236 7.303 7.37 7.437
    7.504 7.571 7.638 7.705 7.772 7.839 7.906 7.973 8.04'
  []
  #
  print_linear_residuals = true
  print_nonlinear_residuals = true
  #
  # checkpoint = false
  [checkpoint]
    type = Checkpoint
    num_files = 3
  []
  # [pgraph]
  #   type = PerfGraphOutput
  #   execute_on = 'initial final' # Default is "final"
  #   level = 2 # Default is 1
  #   heaviest_branch = true # Default is false
  #   heaviest_sections = 7 # Default is 0
  # []
  # file_base = 01_base_a${a_tol}_i${i_tol}_out
[]
