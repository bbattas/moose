##############################################################################
# File: 08_gbE_equal_surfE.i
# File Location: /examples/sintering/paper3/11_2D_square_pore_debugging/08_gbE_equal_surfE
# Created Date: Thursday February 8th 2024
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Thursday February 8th 2024
# Modified By: Brandon Battas
# -----
# Description:
#  The no irradiaton takes several hours to run even with the small input
#   to the point of squaring the pore, so using HPG to run several at the same time
#  This is the case with the surface_energy = gb_energy = 9.86
#   and no irradiation kernels
##############################################################################

[Mesh]
  [ebsd_mesh]
    type = EBSDMeshGenerator
    filename = ../00_d3d_txt/2D_20x20um_8umavg_allVoids.txt
  []
  [subdomain_external]
    type = ParsedSubdomainMeshGenerator
    input = ebsd_mesh
    combinatorial_geometry = 'x > 19000'
    block_id = 1
  []
  parallel_type = DISTRIBUTED
  uniform_refine = 0
[]

[GlobalParams]
  op_num = 3 #10
  var_name_base = gr
  int_width = 1000 #min radius is like 2250, element size of 250
  profile = TANH # not used at the moment? only in circleic?
[]

[Variables]
  [w]
  []
  [phi]
  []
  [PolycrystalVariables]
  []
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
[]

[ICs]
  [PolycrystalICs]
    [PolycrystalColoringIC]
      polycrystal_ic_uo = ebsd
    []
  []
  # [VoidIC]
  #   type = ReconPhaseVarIC
  #   ebsd_reader = ebsd_reader
  #   variable = phi
  #   phase = 2
  #   block = 1
  # []
  [VoidIC]
    type = BoundingBoxIC
    variable = phi
    block = 1
    inside = 1
    outside = 0.01
    x1 = 20000
    x2 = 25000
    y1 = -2000
    y2 = 22000
  []
  [voidIC2]
    type = SmoothCircleIC
    variable = phi
    invalue = 1
    outvalue = 0.01
    radius = 5046.26504
    x1 = 8893.09397
    y1 = 8144.4245
    z1 = 0
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
  # [gr3]
  #   type = NeumannBC
  #   variable = gr3
  #   value = 0
  #   boundary = 'left right top bottom'
  # []
  # [gr4]
  #   type = NeumannBC
  #   variable = gr4
  #   value = 0
  #   boundary = 'left right top bottom'
  # []
  # [gr5]
  #   type = NeumannBC
  #   variable = gr5
  #   value = 0
  #   boundary = 'left right top bottom'
  # []
  # [gr6]
  #   type = NeumannBC
  #   variable = gr6
  #   value = 0
  #   boundary = 'left right top bottom'
  # []
  # [gr7]
  #   type = NeumannBC
  #   variable = gr7
  #   value = 0
  #   boundary = 'left right top bottom'
  # []
  # [gr8]
  #   type = NeumannBC
  #   variable = gr8
  #   value = 0
  #   boundary = 'left right top bottom'
  # []
  # [gr9]
  #   type = NeumannBC
  #   variable = gr9
  #   value = 0
  #   boundary = 'left right top bottom'
  # []
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
    GBwidth = 1
    surf_thickness = 0.5
  []
  [cv_eq]
    type = UO2CvMaterial
    f_name = cv_eq
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
    surface_energy = 9.86 #19.7
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
    outputs = 'nemesis'
  []
  # Testing Outputs
  [D_mat]
    type = ParsedMaterial
    property_name = D_mat
    material_property_names = 'diffusivity'
    expression = 'diffusivity'
    outputs = 'nemesis'
  []
  # Irradiation and Interstitials
  # Interstitials
  [rho_gen_int]
    type = DerivativeParsedMaterial
    property_name = rho_gen_int
    derivative_order = 1
    constant_names = 'Nc Nd f_dot noise'
    constant_expressions = '2 5 1e-8 1'
    material_property_names = 'hs'
    postprocessor_names = 'hs_average'
    expression = 'f_dot * noise * Nc * Nd * hs_average' # * hs
    outputs = nemesis #'nemesis'
  []
  [a_r]
    type = ParsedMaterial
    property_name = a_r
    constant_names = 'Va Z a_0 Di_0 Ei_B kB'
    constant_expressions = '0.04092 250 0.25 1e13 2 8.617343e-5'
    material_property_names = 'hs'
    coupled_variables = 'T'
    expression = 'Di:=Di_0*exp(-Ei_B/(kB*T));
                  Va * Z * Di / (a_0^2)' #hs *
    outputs = none #'nemesis'
  []
  [combined_rho_vac]
    type = DerivativeParsedMaterial
    property_name = combined_rho_vac
    material_property_names = 'rhov rhos hv(phi)'
    expression = 'hv*rhov + (1-hv)*rhos'
    outputs = nemesis #'nemesis'
  []
  [rho_i_dpm]
    type = DerivativeParsedMaterial
    property_name = rho_i_dpm
    # derivative_order = 0 #2
    # coupled_variables = ''
    material_property_names = 'rho_gen_int a_r'
    postprocessor_names = 'average_rho_vac'
    expression = 'if(average_rho_vac>0.0,(rho_gen_int / (a_r * average_rho_vac)),0.0)'
    outputs = nemesis #'nemesis'
  []
  # Vacancies
  [rho_gen_vac]
    type = DerivativeParsedMaterial
    property_name = rho_gen_vac
    derivative_order = 1
    constant_names = 'Nc Nd f_dot noise'
    constant_expressions = '2 5 1e-8 1'
    material_property_names = 'hs'
    # postprocessor_names = 'hs_average'
    expression = 'f_dot * noise * Nc * Nd * hs' # * hs
    outputs = 'nemesis'
  []
  [rho_v_recombRate]
    type = DerivativeParsedMaterial
    property_name = rho_v_recombRate
    coupled_variables = 'w'
    # additional_derivative_symbols = w
    material_property_names = 'a_r rho_i_dpm combined_rho_vac'
    postprocessor_names = 'total_rhoi' # average_rho_vac'
    expression = 'a_r * rho_i_dpm * combined_rho_vac'
    outputs = 'nemesis'
  []
  [rho_v_mixing]
    type = DerivativeParsedMaterial
    property_name = rho_v_mixing
    coupled_variables = 'w'
    derivative_order = 1
    constant_names = 'Nc Vc f_dot noise tc Dc'
    constant_expressions = '2 268 1e-8 1 1e-11 1e12'
    material_property_names = 'chi'
    expression = 'f_dot * noise * Nc * tc * Vc * Dc * chi' # * hs
    outputs = 'nemesis'
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
    v = 'gr0 gr1 gr2' # gr3 gr4 gr5 gr6 gr7 gr8 gr9' # gr10 gr11 gr12 gr13 gr14 gr15'
    gamma = gamma
    mob_name = L_mat
  []
  [kappa_phi]
    type = ACKappaFunction
    variable = phi
    mob_name = L_mat
    kappa_name = kappa
  []
  # # Irradiation
  # # Source/Generation
  # [source_w]
  #   type = MaskedBodyForce
  #   variable = w
  #   mask = rho_gen_vac
  # []
  # # Sink/Recombination
  # [recombination_w]
  #   type = MatReaction
  #   variable = w
  #   mob_name = rho_v_recombRate
  #   # args = rhoi #but its a constant and material not a variable
  # []
  # # Damage/Mixing
  # [ballistic_mix_w]
  #   type = MatDiffusion
  #   variable = w
  #   diffusivity = rho_v_mixing
  # []
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
  # [./ghosted_entities]
  #   type = FeatureFloodCountAux
  #   variable = ghost_regions
  #   flood_counter = grain_tracker
  #   field_display = GHOSTED_ENTITIES
  #   execute_on = 'initial timestep_end'
  # [../]
  # [./halos]
  #   type = FeatureFloodCountAux
  #   variable = halos
  #   flood_counter = grain_tracker
  #   field_display = HALOS
  #   execute_on = 'initial timestep_end'
  # [../]
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
    threshold = 0.1 #0.2
    connecting_threshold = 0.09 #0.08
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
    variable = 'gr0 gr1 gr2'
    threshold = 0.5 #0.2
    connecting_threshold = 0.5 #0.08
    compute_var_to_feature_map = true
    execute_on = 'initial timestep_end'
  []
  [timestep]
    type = TimestepSize
    outputs = csv
  []
  # Irradiation PPs Used in Materials
  [hs_average]
    type = ElementAverageMaterialProperty
    mat_prop = hs
    outputs = csv
  []
  [average_rho_vac]
    type = ElementAverageMaterialProperty
    mat_prop = combined_rho_vac
    outputs = csv
  []
  # Other irradiation based PPs
  [total_phi]
    type = ElementIntegralVariablePostprocessor
    variable = phi
    outputs = csv
  []
  [total_rhoi]
    type = ElementIntegralMaterialProperty
    mat_prop = rho_i_dpm
    outputs = csv
  []
[]

[VectorPostprocessors]
  [voids]
    type = FeatureVolumeVectorPostprocessor
    flood_counter = void_tracker
    execute_on = 'initial timestep_end final'
    output_centroids = false #was true
    outputs = csv
  []
  # [alt_voids]
  #   type = FeatureVolumeVectorPostprocessor
  #   flood_counter = void_tracker_05
  #   execute_on = 'initial timestep_end final'
  #   output_centroids = false #was true
  #   outputs = csv
  # []
  [grain_sizes]
    type = FeatureVolumeVectorPostprocessor
    flood_counter = grain_tracker_fc
    execute_on = 'initial timestep_end final'
    output_centroids = false #was true
    outputs = csv
  []
  # [vectorMemory]
  #   type = VectorMemoryUsage
  #   mem_units = gigabytes
  #   outputs = csv
  #   execute_on = 'NONLINEAR LINEAR TIMESTEP_END'
  # []
[]

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
  end_time = 5e6 #0.006
  steady_state_detection = true
  # num_steps = 30
  # automatic_scaling = true
  # dt = 0.00002
  # dtmax = 1000
  # dt = 0.0001
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 6
    dt = 100 #2.5
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
  [nemesis]
    type = Nemesis
    # interval = 3 # this ExodusII will only output every third time step
  []
  print_linear_residuals = false
  [checkpoint]
    type = Checkpoint
    num_files = 3
  []
[]
