##############################################################################
# File: 07_bicr_miso_high_bulkmatch.i
# File Location: /examples/agg/10_fd_comparison/07_bicr_miso/07_bicr_miso_high_bulkmatch
# Created Date: Friday October 16th 2026
# Author: Battas,Brandon Scott (bbattas@ufl.edu)
# -----
# Last Modified: Friday October 16th 2026
# Modified By: Battas,Brandon Scott
# -----
# Description:
#  setting bulk gbe to match the single miso value (0.984)
#
#
#
##############################################################################

[Mesh]
  [gmg]
    type = DistributedRectilinearMeshGenerator
    dim = 2
    nx = 160
    ny = 160
    xmin = 0
    xmax = 160
    ymin = 0
    ymax = 160
  []
  parallel_type = DISTRIBUTED
  second_order = false
  uniform_refine = 0
[]

[GlobalParams]
  op_num = 2
  var_name_base = gr
[]

[Variables]
  [PolycrystalVariables]
  []
[]

[ICs]
  [gr0_IC]
    type = SmoothCircleIC
    invalue = 1
    outvalue = 0
    radius = 60
    variable = gr0
    x1 = 80
    y1 = 80
    int_width = 6
  []
  [gr1_IC]
    type = SmoothCircleIC
    invalue = 0
    outvalue = 1
    radius = 60
    variable = gr1
    x1 = 80
    y1 = 80
    int_width = 6
  []
[]

[UserObjects]
  [euler_file]
    type = EulerAngleTxtFileReader
    file_name = '../00_sub/bicr_miso_high.txt'
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
    gb_mode = MISO
    bulk_scalar = 0.984
    iso_gbe = 0.5
    alpha_tol = 0
    hgbalpha_tol = 0
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
  nl_rel_tol = 1e-8 #1e-10 # Relative tolerance for nonlienar solves
  # nl_abs_tol = 1e-10

  start_time = 0.0
  # num_steps = 0
  end_time = 30
  dt = 0.001
  # [TimeStepper]
  #   type = IterationAdaptiveDT
  #   dt = 0.01 #0.001
  #   # cutback_factor = 0.9
  #   # growth_factor = 1.1
  #   optimal_iterations = 6
  #   linear_iteration_ratio = 1e5
  # []

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
    sync_times = '0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5
      1.6 1.7 1.8 1.9 2 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3 3.1
      3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4 4.1 4.2 4.3 4.4 4.5 4.6 4.7
      4.8 4.9 5 5.1 5.2 5.3 5.4 5.5 5.6 5.7 5.8 5.9 6 6.1 6.2 6.3
      6.4 6.5 6.6 6.7 6.8 6.9 7 7.1 7.2 7.3 7.4 7.5 7.6 7.7 7.8 7.9
      8 8.1 8.2 8.3 8.4 8.5 8.6 8.7 8.8 8.9 9 9.1 9.2 9.3 9.4 9.5
      9.6 9.7 9.8 9.9 10 10.1 10.2 10.3 10.4 10.5 10.6 10.7 10.8 10.9 11 11.1
      11.2 11.3 11.4 11.5 11.6 11.7 11.8 11.9 12 12.1 12.2 12.3 12.4 12.5 12.6 12.7
      12.8 12.9 13 13.1 13.2 13.3 13.4 13.5 13.6 13.7 13.8 13.9 14 14.1 14.2 14.3
      14.4 14.5 14.6 14.7 14.8 14.9 15 15.1 15.2 15.3 15.4 15.5 15.6 15.7 15.8 15.9
      16 16.1 16.2 16.3 16.4 16.5 16.6 16.7 16.8 16.9 17 17.1 17.2 17.3 17.4 17.5
      17.6 17.7 17.8 17.9 18 18.1 18.2 18.3 18.4 18.5 18.6 18.7 18.8 18.9 19 19.1
      19.2 19.3 19.4 19.5 19.6 19.7 19.8 19.9 20 20.1 20.2 20.3 20.4 20.5 20.6 20.7
      20.8 20.9 21 21.1 21.2 21.3 21.4 21.5 21.6 21.7 21.8 21.9 22 22.1 22.2 22.3
      22.4 22.5 22.6 22.7 22.8 22.9 23 23.1 23.2 23.3 23.4 23.5 23.6 23.7 23.8 23.9
      24 24.1 24.2 24.3 24.4 24.5 24.6 24.7 24.8 24.9 25 25.1 25.2 25.3 25.4 25.5
      25.6 25.7 25.8 25.9 26 26.1 26.2 26.3 26.4 26.5 26.6 26.7 26.8 26.9 27 27.1
      27.2 27.3 27.4 27.5 27.6 27.7 27.8 27.9 28 28.1 28.2 28.3 28.4 28.5 28.6 28.7
      28.8 28.9 29 29.1 29.2 29.3 29.4 29.5 29.6 29.7 29.8 29.9 30'
  []
  [csv]
    type = CSV
    sync_only = true
    sync_times = '0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5
      1.6 1.7 1.8 1.9 2 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3 3.1
      3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4 4.1 4.2 4.3 4.4 4.5 4.6 4.7
      4.8 4.9 5 5.1 5.2 5.3 5.4 5.5 5.6 5.7 5.8 5.9 6 6.1 6.2 6.3
      6.4 6.5 6.6 6.7 6.8 6.9 7 7.1 7.2 7.3 7.4 7.5 7.6 7.7 7.8 7.9
      8 8.1 8.2 8.3 8.4 8.5 8.6 8.7 8.8 8.9 9 9.1 9.2 9.3 9.4 9.5
      9.6 9.7 9.8 9.9 10 10.1 10.2 10.3 10.4 10.5 10.6 10.7 10.8 10.9 11 11.1
      11.2 11.3 11.4 11.5 11.6 11.7 11.8 11.9 12 12.1 12.2 12.3 12.4 12.5 12.6 12.7
      12.8 12.9 13 13.1 13.2 13.3 13.4 13.5 13.6 13.7 13.8 13.9 14 14.1 14.2 14.3
      14.4 14.5 14.6 14.7 14.8 14.9 15 15.1 15.2 15.3 15.4 15.5 15.6 15.7 15.8 15.9
      16 16.1 16.2 16.3 16.4 16.5 16.6 16.7 16.8 16.9 17 17.1 17.2 17.3 17.4 17.5
      17.6 17.7 17.8 17.9 18 18.1 18.2 18.3 18.4 18.5 18.6 18.7 18.8 18.9 19 19.1
      19.2 19.3 19.4 19.5 19.6 19.7 19.8 19.9 20 20.1 20.2 20.3 20.4 20.5 20.6 20.7
      20.8 20.9 21 21.1 21.2 21.3 21.4 21.5 21.6 21.7 21.8 21.9 22 22.1 22.2 22.3
      22.4 22.5 22.6 22.7 22.8 22.9 23 23.1 23.2 23.3 23.4 23.5 23.6 23.7 23.8 23.9
      24 24.1 24.2 24.3 24.4 24.5 24.6 24.7 24.8 24.9 25 25.1 25.2 25.3 25.4 25.5
      25.6 25.7 25.8 25.9 26 26.1 26.2 26.3 26.4 26.5 26.6 26.7 26.8 26.9 27 27.1
      27.2 27.3 27.4 27.5 27.6 27.7 27.8 27.9 28 28.1 28.2 28.3 28.4 28.5 28.6 28.7
      28.8 28.9 29 29.1 29.2 29.3 29.4 29.5 29.6 29.7 29.8 29.9 30'
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
  # file_base = 01_base_a${a_tol}_i${i_tol}_out
[]
