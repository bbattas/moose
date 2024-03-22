##############################################################################
# File: 01_IC_check_1step.i
# File Location: /examples/sintering/paper3/14_large_2D_noIrradiation/01_IC_check_1step
# Created Date: Friday March 22nd 2024
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Friday March 22nd 2024
# Modified By: Brandon Battas
# -----
# Description:
#  0 or 1 timestep input to check what the generated ICs look like
#  2.3m DoFs
#
#
##############################################################################

[Mesh]
  [ebsd_mesh]
    type = EBSDMeshGenerator
    filename = ../00_d3d_txt/2D_100x100um_8umavg_config2_8pore.txt
  []
  [subdomain_external]
    type = ParsedSubdomainMeshGenerator
    input = ebsd_mesh
    combinatorial_geometry = 'x > 95000'
    block_id = 1
  []
  parallel_type = DISTRIBUTED
  uniform_refine = 0
[]

[GlobalParams]
  op_num = 10
  var_name_base = gr
  int_width = 2000 #1000 #min radius is like 2250, element size of 250
  # profile = TANH # not used at the moment? only in circleic?
[]

[Variables]
  [w]
  []
  [phi]
  []
  [PolycrystalVariables]
  []
  # [rhoi_scalar]
  #   family = SCALAR
  #   order = FIRST
  #   initial_condition = 0.0001 #0.05 #0.01
  # []
[]

[AuxVariables]
  [bnds]
  []
  [T]
    order = CONSTANT
    family = MONOMIAL
  []
  [OU]
    order = CONSTANT
    family = MONOMIAL
  []
  [voids]
    order = CONSTANT
    family = MONOMIAL
  []
  [unique_grains]
    order = CONSTANT
    family = MONOMIAL
  []
  [var_indices]
    order = CONSTANT
    family = MONOMIAL
  []
  [ebsd_grains]
    family = MONOMIAL
    order = CONSTANT
  []
  # [rhoi_aux]
  # []
[]

[ICs]
  [PolycrystalICs]
    [PolycrystalColoringIC]
      polycrystal_ic_uo = ebsd
    []
  []
  [VoidIC]
    type = BoundingBoxIC
    variable = phi
    block = 1
    inside = 1
    outside = 0.01
    x1 = 100000
    x2 = 124000
    y1 = -4000
    y2 = 104000
  []
  [voidIC2]
    type = SpecifiedSmoothCircleIC
    variable = phi
    invalue = 1
    outvalue = 0.01
    radii = '4548.11 5312.36 4711.77 5240.34 5026.54 5430.06 4389.2  4289.61'
    x_positions = '58015.35 61941.06 74432.15 45472.52 24616.72 18372.46 29933.71 82823.36'
    y_positions = '53348.45 28238.81 75444.65 74557.06 39339.06 78867.79 17592.7  52979.64'
    z_positions = '0 0 0 0 0'
    block = 0
  []
[]

[BCs]
  [phi]
    type = NeumannBC
    variable = phi
    value = 0
    boundary = 'left right top bottom'
  []
  [gr0]
    type = NeumannBC
    variable = gr0
    value = 0
    boundary = 'left right top bottom'
  []
  [gr1]
    type = NeumannBC
    variable = gr1
    value = 0
    boundary = 'left right top bottom'
  []
  [gr2]
    type = NeumannBC
    variable = gr2
    value = 0
    boundary = 'left right top bottom'
  []
  [gr3]
    type = NeumannBC
    variable = gr3
    value = 0
    boundary = 'left right top bottom'
  []
  [gr4]
    type = NeumannBC
    variable = gr4
    value = 0
    boundary = 'left right top bottom'
  []
  [gr5]
    type = NeumannBC
    variable = gr5
    value = 0
    boundary = 'left right top bottom'
  []
  [gr6]
    type = NeumannBC
    variable = gr6
    value = 0
    boundary = 'left right top bottom'
  []
  [gr7]
    type = NeumannBC
    variable = gr7
    value = 0
    boundary = 'left right top bottom'
  []
  [gr8]
    type = NeumannBC
    variable = gr8
    value = 0
    boundary = 'left right top bottom'
  []
  [gr9]
    type = NeumannBC
    variable = gr9
    value = 0
    boundary = 'left right top bottom'
  []
[]

[Functions]
  [f_T]
    type = ConstantFunction
    value = 1600
  []
  [f_OU]
    type = ConstantFunction
    value = 2.0
  []
  # [rhoi_func]
  #   type = ParsedFunction
  #   symbol_names = 'rhoi_scalar_reporter'
  #   symbol_values = 'rhoi_scalar_reporter'
  #   expression = 'rhoi_scalar_reporter'
  #   execute_on = 'INITIAL LINEAR'
  # []
[]

[Materials]
  # Free energy coefficients for parabolic curves
  [ks]
    type = ParsedMaterial
    property_name = ks
    coupled_variables = 'T'
    constant_names = 'a b'
    constant_expressions = '-0.0025 157.16'
    expression = 'a*T + b'
  []
  [kv]
    type = ParsedMaterial
    property_name = kv
    material_property_names = 'ks'
    expression = '10*ks'
  []
  # Diffusivity and mobilities
  [chiD]
    type = GrandPotentialIsoMaterial
    f_name = chiD
    solid_mobility = L #CHANGED FROM L
    void_mobility = Lv
    chi = chi
    c = phi
    T = T
    D0 = 8.33e9
    GBmob0 = 1.4759e9
    Q = 2.77
    Em = 3.608
    bulkindex = 1
    gbindex = -1 # -1 sets the GB D to the LANL MD Value in GPIsoMat
    surfindex = 1e11
    GBwidth = 1.0
    surf_thickness = 0.5
    iw_scaling = true
  []
  [cv_eq]
    type = UO2CvMaterial
    property_name = cv_eq
    T = T
    c = phi
    OU = OU
    gb_se_csv = '../../../gb_segregation_csv/sigma9_se.csv
                 ../../../gb_segregation_csv/sigma11_se.csv'
    outputs = 'none'
  []
  [sintering]
    type = GrandPotentialSinteringMaterial
    chemical_potential = w
    void_op = phi
    Temperature = T
    surface_energy = 19.7
    grainboundary_energy = 9.86
    void_energy_coefficient = kv
    solid_energy_coefficient = ks
    solid_energy_model = PARABOLIC
    equilibrium_vacancy_concentration = cv_eq
  []
  # Concentration is only meant for output
  [c]
    type = ParsedMaterial
    property_name = c
    material_property_names = 'hs rhos hv rhov'
    constant_names = 'Va'
    constant_expressions = '0.04092'
    expression = 'Va*(hs*rhos + hv*rhov)'
    outputs = none #exodus
  []
  [L_mat]
    type = DerivativeParsedMaterial
    property_name = L_mat
    material_property_names = 'L Lv'
    coupled_variables = 'phi'
    constant_names = 'p0'
    constant_expressions = '0.3'
    expression = 'hv:=if(phi<=0.0,0.0,if(phi>=p0,1.0,6*(phi/p0)^5 - 15*(phi/p0)^4 + 10*(phi/p0)^3));
                hv*Lv + (1-hv)*L'
    outputs = none #'nemesis'
  []
  [combined_rho_vac]
    type = DerivativeParsedMaterial
    property_name = combined_rho_vac
    coupled_variables = 'w phi'
    material_property_names = 'rhov rhos hv(phi)'
    expression = 'hv*rhov + (1-hv)*rhos' #'(1-hv)*rhos' #
    outputs = none #nemesis #'nemesis'
  []
[]

[Modules]
  [PhaseField]
    [GrandPotential]
      switching_function_names = 'hv hs'
      anisotropic = 'false'

      chemical_potentials = 'w'
      mobilities = 'chiD'
      susceptibilities = 'chi'
      free_energies_w = 'rhov rhos'

      gamma_gr = gamma
      mobility_name_gr = L_mat
      kappa_gr = kappa
      free_energies_gr = 'omegav omegas'

      additional_ops = 'phi'
      gamma_grxop = gamma
      mobility_name_op = L_mat
      kappa_op = kappa
      free_energies_op = 'omegav omegas'
    []
  []
[]

[Kernels]
  [barrier_phi]
    type = ACBarrierFunction
    variable = phi
    v = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 gr8 gr9' # gr10 gr11 gr12 gr13 gr14 gr15'
    gamma = gamma
    mob_name = L_mat
  []
  [kappa_phi]
    type = ACKappaFunction
    variable = phi
    mob_name = L_mat
    kappa_name = kappa
  []
[]

[AuxKernels]
  [bnds_aux]
    type = BndsCalcAux
    variable = bnds
    execute_on = 'initial timestep_end'
  []
  [T_aux]
    type = FunctionAux
    variable = T
    function = f_T
  []
  [OU_aux]
    type = FunctionAux
    variable = OU
    function = f_OU
  []
  [voids_aux]
    type = FeatureFloodCountAux
    variable = voids
    flood_counter = void_tracker
    field_display = UNIQUE_REGION
    execute_on = 'INITIAL TIMESTEP_END'
  []
  #NEW Grain_Tracker related stuff
  [unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
    execute_on = 'initial timestep_end'
  []
  [var_indices]
    type = FeatureFloodCountAux
    variable = var_indices
    flood_counter = grain_tracker
    field_display = VARIABLE_COLORING
    execute_on = 'initial timestep_end'
  []
  [grain_aux]
    type = EBSDReaderPointDataAux
    variable = ebsd_grains
    ebsd_reader = ebsd_reader
    data_name = 'feature_id'
    execute_on = 'initial timestep_end'
  []
[]

[Postprocessors]
  # [memoryAll]
  #   type = MemoryUsage
  #   mem_units = megabytes
  #   outputs = csv
  #   execute_on = 'NONLINEAR LINEAR TIMESTEP_END'
  #   report_peak_value = false
  # []
  # [memoryPeak]
  #   type = MemoryUsage
  #   mem_units = megabytes
  #   outputs = csv
  #   execute_on = 'NONLINEAR LINEAR TIMESTEP_END'
  #   report_peak_value = true
  # []
  [memory1CPU]
    type = MemoryUsage
    mem_units = megabytes
    outputs = csv
    execute_on = 'NONLINEAR LINEAR TIMESTEP_END'
    value_type = max_process
  []
  [n_DOFs]
    type = NumDOFs
    outputs = csv
  []
  [c_total]
    type = ElementIntegralMaterialProperty
    mat_prop = c
    outputs = csv
  []
  [nonlinear]
    type = NumNonlinearIterations
    outputs = csv
  []
  [linear]
    type = NumLinearIterations
    outputs = csv
  []
  [residuals]
    type = NumResidualEvaluations
    outputs = csv
  []
  [runtime]
    type = PerfGraphData
    section_name = "Root"
    data_type = TOTAL
  []
  [void_tracker]
    type = FeatureFloodCount
    variable = phi
    threshold = 0.5 # 0.1 #0.2
    connecting_threshold = 0.5 #0.09 #0.08
    compute_var_to_feature_map = true
    execute_on = 'initial timestep_end'
  []
  # [void_tracker_05]
  #   type = FeatureFloodCount
  #   variable = phi
  #   threshold = 0.5 #0.2
  #   connecting_threshold = 0.5 #0.08
  #   compute_var_to_feature_map = true
  #   execute_on = 'initial timestep_end'
  # []
  [grain_tracker_fc]
    type = FeatureFloodCount
    variable = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 gr8 gr9' #'gr0 gr1 gr2'
    threshold = 0.5 #0.2
    connecting_threshold = 0.5 #0.08
    compute_var_to_feature_map = true
    execute_on = 'initial timestep_end'
  []
  [timestep]
    type = TimestepSize
    outputs = csv
  []
  # # Irradiation PPs Used in Materials
  # [hs_average]
  #   type = ElementAverageMaterialProperty
  #   mat_prop = hs
  #   outputs = csv
  # []
  # [average_rho_vac]
  #   type = ElementAverageMaterialProperty
  #   mat_prop = combined_rho_vac
  #   outputs = csv
  # []
  # Other irradiation based PPs
  [total_phi]
    type = ElementIntegralVariablePostprocessor
    variable = phi
    outputs = csv
  []
  # [total_rhoi]
  #   # type = ElementIntegralMaterialProperty
  #   # mat_prop = rho_i_dpm
  #   type = ElementIntegralVariablePostprocessor
  #   variable = rhoi_aux
  #   outputs = csv
  # []
  # [total_rhoi_scalar]
  #   # type = ElementIntegralMaterialProperty
  #   # mat_prop = rho_i_dpm
  #   type = ElementIntegralVariablePostprocessor
  #   variable = rhoi_scalar
  #   outputs = csv
  # []
  [total_rhov]
    type = ElementIntegralMaterialProperty
    mat_prop = combined_rho_vac
    outputs = csv
  []
  # # Scalar Kernel
  # [rhoi_scalar_reporter]
  #   type = ScalarVariable
  #   variable = rhoi_scalar
  #   execute_on = 'INITIAL LINEAR TIMESTEP_END'
  # []
  # [a_r_pp]
  #   type = ElementAverageMaterialProperty
  #   mat_prop = a_r
  #   outputs = csv
  # []
  # [fdot_x10]
  #   type = ElementExtremeMaterialProperty
  #   mat_prop = rho_gen_vac
  #   value_type = max
  #   outputs = csv
  # []
[]

# [VectorPostprocessors]
#   [voids]
#     type = FeatureVolumeVectorPostprocessor
#     flood_counter = void_tracker
#     execute_on = 'initial timestep_end final'
#     output_centroids = true #false #was true
#     outputs = csv
#   []
#   # [alt_voids]
#   #   type = FeatureVolumeVectorPostprocessor
#   #   flood_counter = void_tracker_05
#   #   execute_on = 'initial timestep_end final'
#   #   output_centroids = false #was true
#   #   outputs = csv
#   # []
#   [grain_sizes]
#     type = FeatureVolumeVectorPostprocessor
#     flood_counter = grain_tracker_fc
#     execute_on = 'initial timestep_end final'
#     output_centroids = false #was true
#     outputs = csv
#   []
#   # [vectorMemory]
#   #   type = VectorMemoryUsage
#   #   mem_units = gigabytes
#   #   outputs = csv
#   #   execute_on = 'NONLINEAR LINEAR TIMESTEP_END'
#   # []
# []

[UserObjects]
  [ebsd_reader]
    type = EBSDReader
  []
  [ebsd]
    type = PolycrystalEBSD
    # coloring_algorithm = bt
    ebsd_reader = ebsd_reader
    enable_var_coloring = true
    output_adjacency_matrix = true
    phase = 1
  []
  [grain_tracker]
    type = GrainTracker
    threshold = 0.1 #0.2
    connecting_threshold = 0.09 #0.08
    compute_halo_maps = false #true#false
    verbosity_level = 1
  []
  # [terminator_2void]
  #   type = Terminator
  #   expression = 'void_tracker < 2'
  # []
[]

[Preconditioning]
  [SMP] #slow but good, very slow for 3D (might be another option then)
    type = SMP
    # full = true
    coupled_groups = 'w,phi'
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -sub_pc_type -pc_factor_levels'
  petsc_options_value = 'asm ilu 2'
  # petsc_options_iname = '-pc_type -pc_hypre_type' # -snes_type'
  # petsc_options_value = 'hypre boomeramg' # vinewtonrsls'
  # petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap'
  # petsc_options_value = ' asm      lu           2'
  nl_max_its = 12 #20 #40 too large- optimal_iterations is 6
  l_max_its = 200 #30 #if it seems like its using a lot it might still be fine
  l_tol = 1e-04
  nl_rel_tol = 1e-6 #default is 1e-8
  nl_abs_tol = 1e-6 #only needed when near equilibrium or veeeery small timesteps and things changing FAST
  start_time = 0
  # end_time = 1e8 #5e6 #0.006
  num_steps = 0
  steady_state_detection = true
  # # From tonks ode input
  # automatic_scaling = true
  # compute_scaling_once = false
  # line_search = none
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 6
    dt = 1 #2.5
    linear_iteration_ratio = 1e5 #needed with large linear number for asmilu
    # growth_factor = 1.2
    # cutback_factor = 0.8
    # cutback_factor_at_failure = 0.5 #might be different from the curback_factor
  []
  # [Adaptivity]
  #   refine_fraction = 0.8
  #   coarsen_fraction = 0.05 #minimize this- adds error
  #   max_h_level = 2 #test a short simulation with 1,2,3,4 for this to see where it stops helping
  #   initial_adaptivity = 2
  # []
[]

[Outputs]
  perf_graph = false
  csv = true
  exodus = false
  # nemesis = false
  # fr_1.00e-10_csv/fr_1.00e-10
  # [csv]
  #   type = CSV
  #   file_base = fr_${f_dot}_csv/fr_${f_dot}
  # []
  [nemesis]
    type = Nemesis
    # file_base = fr_${f_dot}_nemesis/fr_${f_dot}_nemesis
    # interval = 3 # this ExodusII will only output every third time step
  []
  print_linear_residuals = false
  # [checkpoint]
  #   type = Checkpoint
  #   file_base = fr_${f_dot}_checkpoint
  #   num_files = 3
  # []
[]
