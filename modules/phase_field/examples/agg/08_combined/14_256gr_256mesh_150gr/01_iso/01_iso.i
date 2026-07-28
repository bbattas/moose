##############################################################################
# File: 01_iso.i
# File Location: /examples/agg/08_combined/14_256gr_256mesh_150gr/01_iso
# Created Date: Tuesday July 28th 2026
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Tuesday July 28th 2026
# Modified By: Brandon Battas
# -----
# Description:
#  150 grain end case with 120 frames
#
#
#
##############################################################################

randseed = 20

[Mesh]
  [gmg]
    type = DistributedRectilinearMeshGenerator
    dim = 2
    nx = 384
    ny = 384
    xmin = 0
    xmax = 256
    ymin = 0
    ymax = 256
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
  # [euler_file]
  #   type = EulerAngleTxtFileReader
  #   file_name = '../../00_sub/euler_256_mackenzie.txt'
  # []
  [voronoi]
    type = PolycrystalVoronoi
    coloring_algorithm = bt #jp #bt
    grain_num = 256
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
    gb_mode = ISO
    bulk_scalar = 0.75
    iso_gbe = 0.5
    alpha_tol = 10
    hgbalpha_tol = 5
    # ifunc_a = 0.1
    kappa_name = kappa_op
    gbe_iso_name = gbe_iso
    enable_misorientation = false
    # euler_angle_provider = euler_file
    w_inc = 1
    w_miso = 1
    output_properties = 'gt_num L gamma_asymm int_width gbe_norm gbe_gb gb_notj_mask gb_tj_mask'
    outputs = 'exodus'
  []
  [gbe_notj_masked]
    type = ParsedMaterial
    property_name = gbe_notj_masked
    material_property_names = 'gbe_gb gb_notj_mask'
    expression = 'gbe_gb * gb_notj_mask'
    outputs = 'exodus'
  []
  [gbe_tj_masked]
    type = ParsedMaterial
    property_name = gbe_tj_masked
    material_property_names = 'gbe_gb gb_tj_mask'
    expression = 'gbe_gb * gb_tj_mask'
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
    # GBE
  [gbe_notj_total]
    type = ElementIntegralMaterialProperty
    mat_prop = gbe_notj_masked
  []
  [gbe_tj_total]
    type = ElementIntegralMaterialProperty
    mat_prop = gbe_tj_masked
  []
  [gb_notj_area]
    type = ElementIntegralMaterialProperty
    mat_prop = gb_notj_mask
  []
  [gb_tj_area]
    type = ElementIntegralMaterialProperty
    mat_prop = gb_tj_mask
  []
  [gbe_notj_avg]
    type = ParsedPostprocessor
    pp_names = 'gbe_notj_total gb_notj_area'
    expression = 'gbe_notj_total / gb_notj_area'
  []
  [gbe_tj_avg]
    type = ParsedPostprocessor
    pp_names = 'gbe_tj_total gb_tj_area'
    expression = 'gbe_tj_total / gb_tj_area'
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
  end_time = 6.9
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
    sync_times = '0 0.0575 0.115 0.1725 0.23 0.2875 0.345 0.4025 0.46 0.5175 0.575 0.6325 0.69 0.7475 0.805 0.8625
    0.92 0.9775 1.035 1.0925 1.15 1.2075 1.265 1.3225 1.38 1.4375 1.495 1.5525 1.61 1.6675 1.725 1.7825
    1.84 1.8975 1.955 2.0125 2.07 2.1275 2.185 2.2425 2.3 2.3575 2.415 2.4725 2.53 2.5875 2.645 2.7025
    2.76 2.8175 2.875 2.9325 2.99 3.0475 3.105 3.1625 3.22 3.2775 3.335 3.3925 3.45 3.5075 3.565 3.6225
    3.68 3.7375 3.795 3.8525 3.91 3.9675 4.025 4.0825 4.14 4.1975 4.255 4.3125 4.37 4.4275 4.485 4.5425
    4.6 4.6575 4.715 4.7725 4.83 4.8875 4.945 5.0025 5.06 5.1175 5.175 5.2325 5.29 5.3475 5.405 5.4625
    5.52 5.5775 5.635 5.6925 5.75 5.8075 5.865 5.9225 5.98 6.0375 6.095 6.1525 6.21 6.2675 6.325 6.3825
    6.44 6.4975 6.555 6.6125 6.67 6.7275 6.785 6.8425 6.9'
  []
  [csv]
    type = CSV
    sync_times = '0 0.0575 0.115 0.1725 0.23 0.2875 0.345 0.4025 0.46 0.5175 0.575 0.6325 0.69 0.7475 0.805 0.8625
    0.92 0.9775 1.035 1.0925 1.15 1.2075 1.265 1.3225 1.38 1.4375 1.495 1.5525 1.61 1.6675 1.725 1.7825
    1.84 1.8975 1.955 2.0125 2.07 2.1275 2.185 2.2425 2.3 2.3575 2.415 2.4725 2.53 2.5875 2.645 2.7025
    2.76 2.8175 2.875 2.9325 2.99 3.0475 3.105 3.1625 3.22 3.2775 3.335 3.3925 3.45 3.5075 3.565 3.6225
    3.68 3.7375 3.795 3.8525 3.91 3.9675 4.025 4.0825 4.14 4.1975 4.255 4.3125 4.37 4.4275 4.485 4.5425
    4.6 4.6575 4.715 4.7725 4.83 4.8875 4.945 5.0025 5.06 5.1175 5.175 5.2325 5.29 5.3475 5.405 5.4625
    5.52 5.5775 5.635 5.6925 5.75 5.8075 5.865 5.9225 5.98 6.0375 6.095 6.1525 6.21 6.2675 6.325 6.3825
    6.44 6.4975 6.555 6.6125 6.67 6.7275 6.785 6.8425 6.9'
  []
  #
  print_linear_residuals = true
  print_nonlinear_residuals = true
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
