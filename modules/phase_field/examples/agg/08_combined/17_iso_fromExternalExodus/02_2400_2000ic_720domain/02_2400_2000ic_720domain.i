##############################################################################
# File: 02_2400_2000ic_720domain.i
# File Location: /examples/agg/08_combined/17_iso_fromExternalExodus/02_2400_2000ic_720domain
# Created Date: Thursday August 27th 2026
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Thursday August 27th 2026
# Modified By: Brandon Battas
# -----
# Description:
#  Same style but a 2000 op problem in the larger mesh (relaxed from 2400 for the ic)
#
#
#
##############################################################################

# randseed = 20

[Mesh]
  [gmg]
    type = DistributedRectilinearMeshGenerator
    dim = 2
    nx = 1080
    ny = 1080
    xmin = 0
    xmax = 720
    ymin = 0
    ymax = 720
  []
  parallel_type = DISTRIBUTED
  second_order = false
  uniform_refine = 0
[]

[GlobalParams]
  op_num = 20 #25
  var_name_base = gr
[]

[Variables]
  [PolycrystalVariables]
  []
[]

[ICs]
  [ic_gr0]
    type = SolutionIC
    from_variable = gr0
    variable = gr0
    solution_uo = solution_userobject
  []
  [ic_gr1]
    type = SolutionIC
    from_variable = gr1
    variable = gr1
    solution_uo = solution_userobject
  []
  [ic_gr2]
    type = SolutionIC
    from_variable = gr2
    variable = gr2
    solution_uo = solution_userobject
  []
  [ic_gr3]
    type = SolutionIC
    from_variable = gr3
    variable = gr3
    solution_uo = solution_userobject
  []
  [ic_gr4]
    type = SolutionIC
    from_variable = gr4
    variable = gr4
    solution_uo = solution_userobject
  []
  [ic_gr5]
    type = SolutionIC
    from_variable = gr5
    variable = gr5
    solution_uo = solution_userobject
  []
  [ic_gr6]
    type = SolutionIC
    from_variable = gr6
    variable = gr6
    solution_uo = solution_userobject
  []
  [ic_gr7]
    type = SolutionIC
    from_variable = gr7
    variable = gr7
    solution_uo = solution_userobject
  []
  [ic_gr8]
    type = SolutionIC
    from_variable = gr8
    variable = gr8
    solution_uo = solution_userobject
  []
  [ic_gr9]
    type = SolutionIC
    from_variable = gr9
    variable = gr9
    solution_uo = solution_userobject
  []
  [ic_gr10]
    type = SolutionIC
    from_variable = gr10
    variable = gr10
    solution_uo = solution_userobject
  []
  [ic_gr11]
    type = SolutionIC
    from_variable = gr11
    variable = gr11
    solution_uo = solution_userobject
  []
  [ic_gr12]
    type = SolutionIC
    from_variable = gr12
    variable = gr12
    solution_uo = solution_userobject
  []
  [ic_gr13]
    type = SolutionIC
    from_variable = gr13
    variable = gr13
    solution_uo = solution_userobject
  []
  [ic_gr14]
    type = SolutionIC
    from_variable = gr14
    variable = gr14
    solution_uo = solution_userobject
  []
  [ic_gr15]
    type = ConstantIC
    variable = gr15
    value = 0.0
  []
  [ic_gr16]
    type = ConstantIC
    variable = gr16
    value = 0.0
  []
  [ic_gr17]
    type = ConstantIC
    variable = gr17
    value = 0.0
  []
  [ic_gr18]
    type = ConstantIC
    variable = gr18
    value = 0.0
  []
  [ic_gr19]
    type = ConstantIC
    variable = gr19
    value = 0.0
  []
[]

[UserObjects]
  # [euler_file]
  #   type = EulerAngleTxtFileReader
  #   file_name = '../../00_sub/euler_256_mackenzie.txt'
  # []
  # [voronoi]
  #   type = PolycrystalVoronoi
  #   coloring_algorithm = bt #jp #bt
  #   grain_num = 1000
  #   rand_seed = ${randseed}
  #   int_width = 3 #6
  # []
  [solution_userobject]
    type = SolutionUserObject
    mesh = 'relaxed_2000.e'
    timestep = 'LATEST'
    system_variables = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 gr8 gr9 gr10 gr11 gr12 gr13 gr14'
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
    expression = 'grain_tracker < 500'
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
[]

[Outputs]
  # exodus = true # Exodus file will be outputted
  #nemesis = true
  console = true
  # csv = true
  [exodus]
    type = Exodus
    sync_only = true
    sync_times = '0 0.0335 0.067 0.1005 0.134 0.1675 0.201 0.2345 0.268 0.3015 0.335 0.3685 0.402 0.4355 0.469 0.5025
    0.536 0.5695 0.603 0.6365 0.67 0.7035 0.737 0.7705 0.804 0.8375 0.871 0.9045 0.938 0.9715 1.005 1.0385
    1.072 1.1055 1.139 1.1725 1.206 1.2395 1.273 1.3065 1.34 1.3735 1.407 1.4405 1.474 1.5075 1.541 1.5745
    1.608 1.6415 1.675 1.7085 1.742 1.7755 1.809 1.8425 1.876 1.9095 1.943 1.9765 2.01 2.0435 2.077 2.1105
    2.144 2.1775 2.211 2.2445 2.278 2.3115 2.345 2.3785 2.412 2.4455 2.479 2.5125 2.546 2.5795 2.613 2.6465
    2.68 2.7135 2.747 2.7805 2.814 2.8475 2.881 2.9145 2.948 2.9815 3.015 3.0485 3.082 3.1155 3.149 3.1825
    3.216 3.2495 3.283 3.3165 3.35 3.3835 3.417 3.4505 3.484 3.5175 3.551 3.5845 3.618 3.6515 3.685 3.7185
    3.752 3.7855 3.819 3.8525 3.886 3.9195 3.953 3.9865 4.02 4.0535 4.087 4.1205 4.154 4.1875 4.221 4.2545
    4.288 4.3215 4.355 4.3885 4.422 4.4555 4.489 4.5225 4.556 4.5895 4.623 4.6565 4.69 4.7235 4.757 4.7905
    4.824 4.8575 4.891 4.9245 4.958 4.9915 5.025 5.0585 5.092 5.1255 5.159 5.1925 5.226 5.2595 5.293 5.3265
    5.36 5.3935 5.427 5.4605 5.494 5.5275 5.561 5.5945 5.628 5.6615 5.695 5.7285 5.762 5.7955 5.829 5.8625
    5.896 5.9295 5.963 5.9965 6.03 6.0635 6.097 6.1305 6.164 6.1975 6.231 6.2645 6.298 6.3315 6.365 6.3985
    6.432 6.4655 6.499 6.5325 6.566 6.5995 6.633 6.6665 6.7 6.7335 6.767 6.8005 6.834 6.8675 6.901 6.9345
    6.968 7.0015 7.035 7.0685 7.102 7.1355 7.169 7.2025 7.236 7.2695 7.303 7.3365 7.37 7.4035 7.437 7.4705
    7.504 7.5375 7.571 7.6045 7.638 7.6715 7.705 7.7385 7.772 7.8055 7.839 7.8725 7.906 7.9395 7.973 8.0065
    8.04 8.0735 8.107 8.1405 8.174 8.2075 8.241 8.2745 8.308 8.3415 8.375 8.4085 8.442 8.4755 8.509 8.5425
    8.576 8.6095 8.643 8.6765 8.71 8.7435 8.777 8.8105 8.844 8.8775 8.911 8.9445 8.978 9.0115 9.045 9.0785
    9.112 9.1455 9.179 9.2125 9.246 9.2795 9.313 9.3465 9.38 9.4135 9.447 9.4805 9.514 9.5475 9.581 9.6145
    9.648 9.6815 9.715 9.7485 9.782 9.8155 9.849 9.8825 9.916 9.9495 9.983 10.0165 10.05 10.0835 10.117 10.1505
    10.184 10.2175 10.251 10.2845 10.318 10.3515 10.385 10.4185 10.452 10.4855 10.519 10.5525 10.586 10.6195 10.653 10.6865
    10.72 10.7535 10.787 10.8205 10.854 10.8875 10.921 10.9545 10.988 11.0215 11.055 11.0885 11.122 11.1555 11.189 11.2225
    11.256 11.2895 11.323 11.3565 11.39 11.4235 11.457 11.4905 11.524 11.5575 11.591 11.6245 11.658 11.6915 11.725 11.7585
    11.792 11.8255 11.859 11.8925 11.926 11.9595 11.993 12.0265 12.06'
  []
  [csv]
    type = CSV
    sync_only = true
    sync_times = '0 0.0335 0.067 0.1005 0.134 0.1675 0.201 0.2345 0.268 0.3015 0.335 0.3685 0.402 0.4355 0.469 0.5025
    0.536 0.5695 0.603 0.6365 0.67 0.7035 0.737 0.7705 0.804 0.8375 0.871 0.9045 0.938 0.9715 1.005 1.0385
    1.072 1.1055 1.139 1.1725 1.206 1.2395 1.273 1.3065 1.34 1.3735 1.407 1.4405 1.474 1.5075 1.541 1.5745
    1.608 1.6415 1.675 1.7085 1.742 1.7755 1.809 1.8425 1.876 1.9095 1.943 1.9765 2.01 2.0435 2.077 2.1105
    2.144 2.1775 2.211 2.2445 2.278 2.3115 2.345 2.3785 2.412 2.4455 2.479 2.5125 2.546 2.5795 2.613 2.6465
    2.68 2.7135 2.747 2.7805 2.814 2.8475 2.881 2.9145 2.948 2.9815 3.015 3.0485 3.082 3.1155 3.149 3.1825
    3.216 3.2495 3.283 3.3165 3.35 3.3835 3.417 3.4505 3.484 3.5175 3.551 3.5845 3.618 3.6515 3.685 3.7185
    3.752 3.7855 3.819 3.8525 3.886 3.9195 3.953 3.9865 4.02 4.0535 4.087 4.1205 4.154 4.1875 4.221 4.2545
    4.288 4.3215 4.355 4.3885 4.422 4.4555 4.489 4.5225 4.556 4.5895 4.623 4.6565 4.69 4.7235 4.757 4.7905
    4.824 4.8575 4.891 4.9245 4.958 4.9915 5.025 5.0585 5.092 5.1255 5.159 5.1925 5.226 5.2595 5.293 5.3265
    5.36 5.3935 5.427 5.4605 5.494 5.5275 5.561 5.5945 5.628 5.6615 5.695 5.7285 5.762 5.7955 5.829 5.8625
    5.896 5.9295 5.963 5.9965 6.03 6.0635 6.097 6.1305 6.164 6.1975 6.231 6.2645 6.298 6.3315 6.365 6.3985
    6.432 6.4655 6.499 6.5325 6.566 6.5995 6.633 6.6665 6.7 6.7335 6.767 6.8005 6.834 6.8675 6.901 6.9345
    6.968 7.0015 7.035 7.0685 7.102 7.1355 7.169 7.2025 7.236 7.2695 7.303 7.3365 7.37 7.4035 7.437 7.4705
    7.504 7.5375 7.571 7.6045 7.638 7.6715 7.705 7.7385 7.772 7.8055 7.839 7.8725 7.906 7.9395 7.973 8.0065
    8.04 8.0735 8.107 8.1405 8.174 8.2075 8.241 8.2745 8.308 8.3415 8.375 8.4085 8.442 8.4755 8.509 8.5425
    8.576 8.6095 8.643 8.6765 8.71 8.7435 8.777 8.8105 8.844 8.8775 8.911 8.9445 8.978 9.0115 9.045 9.0785
    9.112 9.1455 9.179 9.2125 9.246 9.2795 9.313 9.3465 9.38 9.4135 9.447 9.4805 9.514 9.5475 9.581 9.6145
    9.648 9.6815 9.715 9.7485 9.782 9.8155 9.849 9.8825 9.916 9.9495 9.983 10.0165 10.05 10.0835 10.117 10.1505
    10.184 10.2175 10.251 10.2845 10.318 10.3515 10.385 10.4185 10.452 10.4855 10.519 10.5525 10.586 10.6195 10.653 10.6865
    10.72 10.7535 10.787 10.8205 10.854 10.8875 10.921 10.9545 10.988 11.0215 11.055 11.0885 11.122 11.1555 11.189 11.2225
    11.256 11.2895 11.323 11.3565 11.39 11.4235 11.457 11.4905 11.524 11.5575 11.591 11.6245 11.658 11.6915 11.725 11.7585
    11.792 11.8255 11.859 11.8925 11.926 11.9595 11.993 12.0265 12.06'
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
  [pgraph]
    type = PerfGraphOutput
    execute_on = 'initial final' # Default is "final"
    level = 2 # Default is 1
    heaviest_branch = true # Default is false
    heaviest_sections = 7 # Default is 0
  []
  # file_base = 01_base_a${a_tol}_i${i_tol}_out
[]
