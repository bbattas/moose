##############################################################################
# File: 03_4grain_gamma1_D1.i
# File Location: /examples/sintering/paper2/05_234_gamma1_D1/03_4grain_gamma1_D1/03_4grain_gamma1_D1.i
# Created Date: Wednesday August 16th 2023
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Monday August 28th 2023
# Modified By: Brandon Battas
# -----
# Description:
#  Full 3 grain structure, del*Dgb/Ds = 1, sigma_gb/sigma_s = 1
#  delta=iw but on scaled D internal value, so backcalc Dgb = 1e6, Ds = 2e7
#  sigma_gb=sigma_s=sig5tiltgb=9.86 eV/nm2
#  Coarsened the mesh to 4 elements horizontal across iw (Tonks advice)
#  100x100x88 had 5.4m DoFs-> 450(12k/ea), 300(18k/ea)
##############################################################################

[Mesh]
  [gmg]
    type = DistributedRectilinearMeshGenerator
    dim = 3
    nx = 100
    ny = 100
    nz = 88
    xmin = 0
    xmax = 500
    ymin = 0
    ymax = 500
    zmin = 0
    zmax = 440
  []
  # parallel_type = DISTRIBUTED
  # uniform_refine = 1
  # second_order = false
[]

[GlobalParams]
  op_num = 4
  var_name_base = gr
  int_width = 20 #particle radius is 100
  profile = TANH
[]

[Variables]
  [w]
  []
  [phi]
  []
  [gr0]
  []
  [gr1]
  []
  [gr2]
  []
  [gr3]
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
  # [./ghost_regions]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./halos]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
[]

[ICs]
  [phi_IC]
    type = SpecifiedSmoothCircleIC
    variable = phi
    x_positions = '150 350 250 250'
    y_positions = '250 250 150 350'
    z_positions = '150 150 290 290'
    radii = '100 100 100 100'
    invalue = 0
    outvalue = 1
  []
  [gr0_IC]
    type = SmoothCircleIC
    variable = gr0
    x1 = 150
    y1 = 250
    z1 = 150
    radius = 100
    invalue = 1
    outvalue = 0
  []
  [gr1_IC]
    type = SmoothCircleIC
    variable = gr1
    x1 = 350
    y1 = 250
    z1 = 150
    radius = 100
    invalue = 1
    outvalue = 0
  []
  [gr2_IC]
    type = SmoothCircleIC
    variable = gr2
    x1 = 250
    y1 = 150
    z1 = 290
    radius = 100
    invalue = 1
    outvalue = 0
  []
  [gr3_IC]
    type = SmoothCircleIC
    variable = gr3
    x1 = 250
    y1 = 350
    z1 = 290
    radius = 100
    invalue = 1
    outvalue = 0
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
  [f_cond]
    type = ParsedFunction
    expression = 'void_tracker = 2'
    symbol_names = void_tracker
    symbol_values = void_tracker
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
    gbindex = 1e6
    surfindex = 2e7
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
    outputs = none
  []
  [Diff] # Diffusivity output for debugging
    type = ParsedMaterial
    property_name = Diff
    material_property_names = 'diffusivity'
    expression = 'diffusivity'
    outputs = nemesis
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
    v = 'gr0 gr1 gr2 gr3' #gr2 gr3'# gr4 gr5'# gr6 gr7 gr8 gr9 gr10 gr11 gr12 gr13 gr14 gr15'
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
  #NEW
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
  # [./memoryAll]
  #   type = MemoryUsage
  #   mem_units = megabytes
  #   outputs = csv
  #   execute_on = 'NONLINEAR LINEAR TIMESTEP_END'
  #   report_peak_value = false
  # [../]
  # [./memoryPeak]
  #   type = MemoryUsage
  #   mem_units = megabytes
  #   outputs = csv
  #   execute_on = 'NONLINEAR LINEAR TIMESTEP_END'
  #   report_peak_value = true
  # [../]
  # [./memory1CPU]
  #   type = MemoryUsage
  #   mem_units = megabytes
  #   outputs = csv
  #   execute_on = 'NONLINEAR LINEAR TIMESTEP_END'
  #   value_type = max_process
  # [../]
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
    threshold = 0.6
    connecting_threshold = 0.5 #was 0.2 and worked fine for iso not tensor
    compute_var_to_feature_map = true
    execute_on = 'initial timestep_end'
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
  #  [./vectorMemory]
  #    type = VectorMemoryUsage
  #    mem_units = gigabytes
  #    outputs = csv
  #    execute_on = 'NONLINEAR LINEAR TIMESTEP_END'
  #  [../]
[]

[Controls]
  # [diff]
  #   type = TimePeriod
  #   disable_objects = 'UserObjects::terminator_void'
  #   start_time = '0'
  #   end_time = '100'
  # []
  [conditional]
    type = ConditionalFunctionEnableControl
    conditional_function = f_cond
    enable_objects = 'UserObjects::terminator_void'
    reverse_on_false = false
  []
[]

[UserObjects]
  [terminator_void] #do i have to specify that this is off so that the control can turn it on?
    type = Terminator
    expression = 'void_tracker = 1'
    execute_on = TIMESTEP_END
    enable = false
    # enable = true
  []
  [terminator] #do i have to specify that this is off so that the control can turn it on?
    type = Terminator
    expression = 'grain_tracker < 4'
    execute_on = TIMESTEP_END
    enable = true
  []
  [grain_tracker]
    type = GrainTracker
    threshold = 0.1 #0.2
    connecting_threshold = 0.09 #0.08
    compute_halo_maps = false #true#false
  []
  # [./circle_IC]
  #   type = PolycrystalCircles
  #   # file_name = '../gt4gr.txt'
  #   read_from_file = true
  #   execute_on = 'initial'
  #   threshold = 0.2
  #   connecting_threshold = 0.08
  #   coloring_algorithm = jp
  # [../]
[]

# [Preconditioning]
#   [./SMP] #slow but good, very slow for 3D (might be another option then)
#     type = SMP
#     coupled_groups = 'w,phi'
#   [../]
# []

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
  # petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap'
  # petsc_options_value = ' asm      lu           2'
  nl_max_its = 20 #40 too large- optimal_iterations is 6
  l_max_its = 30 #if it seems like its using a lot it might still be fine
  l_tol = 1e-04
  nl_rel_tol = 1e-6 #default is 1e-8
  nl_abs_tol = 1e-6 #only needed when near equilibrium or veeeery small timesteps and things changing FAST
  start_time = 0
  end_time = 150
  steady_state_detection = true
  # num_steps = 2
  # dt = 0.0001
  dtmax = 500
  # dt = 0.0001
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 8 #WAS 6
    dt = 0.01
    growth_factor = 1.2
    cutback_factor = 0.8
    cutback_factor_at_failure = 0.5 #might be different from the curback_factor
  []
  #[Adaptivity]
  #  refine_fraction = 0.8
  #  coarsen_fraction = 0.05 #minimize this- adds error
  #  max_h_level = 2 #test a short simulation with 1,2,3,4 for this to see where it stops helping
  #  initial_adaptivity = 2
  #[]
[]

[Outputs]
  perf_graph = false
  csv = true
  exodus = false
  # nemesis = true
  [nemesis]
    type = Nemesis
    # interval = 5              # this ExodusII will only output every third time step
  []
  print_linear_residuals = false
  [checkpoint]
    type = Checkpoint
    num_files = 3
  []
[]
