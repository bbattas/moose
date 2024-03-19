##############################################################################
# File: 01_voronoivoidIC.i
# File Location: /examples/sintering/paper3/01_test_inputs/14_devel_moose_examples_Tanvir
# Created Date: Tuesday March 19th 2024
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Tuesday March 19th 2024
# Modified By: Brandon Battas
# -----
# Description:
#  An example input using voronoivoid IC to generate a grain structure with
#   bubbles inside as well as an external pore/void region on 1 boundary
#  The cv_eq is set for UO_2.0 and T=1600K
#
##############################################################################



[Mesh]
  [gmg]
    type = DistributedRectilinearMeshGenerator
    dim = 2
    nx = 100
    ny = 80
    xmin = 0
    xmax = 25000
    ymin = 0
    ymax = 20000
  []
  [subdomain_external]
    type = ParsedSubdomainMeshGenerator
    input = gmg
    combinatorial_geometry = 'x > 20000'
    block_id = 1
  []
  # parallel_type = DISTRIBUTED
  # uniform_refine = 0
[]

[GlobalParams]
  op_num = 3 #10
  var_name_base = gr
  int_width = 1000 #min radius is like 2250, element size of 250
  profile = TANH # not used at the moment? only in circleic?
  # Voronoi Values
  radius = 5046
  bubspac = 10
  numbub = 4
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
  # [ebsd_grains]
  #   family = MONOMIAL
  #   order = CONSTANT
  # []
[]

[ICs]
  [PolycrystalICs]
    [PolycrystalVoronoiVoidIC]
      invalue = 1.0
      outvalue = 0.0
      polycrystal_ic_uo = voronoi
      rand_seed = 10
      block = 0
    []
  []
  [bubble_IC]
    type = PolycrystalVoronoiVoidIC
    variable = phi
    structure_type = voids
    invalue = 1.0
    outvalue = 0.0
    polycrystal_ic_uo = voronoi
    rand_seed = 10
    block = 0
  []
  [VoidIC]
    type = BoundingBoxIC
    variable = phi
    block = 1
    inside = 1.0
    outside = 0.0
    x1 = 20000
    x2 = 30000
    y1 = -2000
    y2 = 22000
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
    solid_mobility = L
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
  []
  [cv_eq]
    type = ParsedMaterial
    property_name = cv_eq
    coupled_variables = 'T'
    constant_names = 'a b'
    constant_expressions = '-0.0025 157.16'
    expression = 'a*T + b'
  []
  # [cv_eq]
  #   type = ParsedMaterial
  #   # type = UO2CvMaterial
  #   # f_name = cv_eq
  #   # T = T
  #   # c = phi
  #   # OU = OU
  #   # gb_se_csv = '../../../gb_segregation_csv/sigma9_se.csv
  #   #              ../../../gb_segregation_csv/sigma11_se.csv'
  #   # outputs = 'none'
  # []
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
    outputs = none
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
  # [voids_aux]
  #   type = FeatureFloodCountAux
  #   variable = voids
  #   flood_counter = void_tracker
  #   field_display = UNIQUE_REGION
  #   execute_on = 'INITIAL TIMESTEP_END'
  # []
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
  # [grain_aux]
  #   type = EBSDReaderPointDataAux
  #   variable = ebsd_grains
  #   ebsd_reader = ebsd_reader
  #   data_name = 'feature_id'
  #   execute_on = 'initial timestep_end'
  # []
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
  [timestep]
    type = TimestepSize
    outputs = csv
  []
[]


[UserObjects]
  [grain_tracker]
    type = GrainTracker
    threshold = 0.1 #0.2
    connecting_threshold = 0.09 #0.08
    compute_halo_maps = false #true#false
    verbosity_level = 1
  []
  [voronoi]
    type = PolycrystalVoronoi
    grain_num = 3 # Number of grains
    rand_seed = 10
    # int_width = 7 # global param
  []
[]

[Preconditioning]
  [SMP] #slow but good, very slow for 3D (might be another option then)
    type = SMP
    full = true
    # coupled_groups = 'w,phi'
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
  nl_max_its = 20 #40 too large- optimal_iterations is 6
  l_max_its = 200 #30 #if it seems like its using a lot it might still be fine
  l_tol = 1e-04
  nl_rel_tol = 1e-6 #default is 1e-8
  nl_abs_tol = 1e-6 #only needed when near equilibrium or veeeery small timesteps and things changing FAST
  start_time = 0
  # end_time = 50000 #0.006
  steady_state_detection = true
  num_steps = 1
  # automatic_scaling = true
  # dt = 0.00002
  # dtmax = 500
  # dt = 0.0001
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 6
    dt = 50 #5#2.5
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
  # nemesis = true
  [nemesis]
    type = Nemesis
    # interval = 3 # this ExodusII will only output every third time step
  []
  print_linear_residuals = true
  # [checkpoint]
  #   type = Checkpoint
  #   num_files = 3
  # []
[]
