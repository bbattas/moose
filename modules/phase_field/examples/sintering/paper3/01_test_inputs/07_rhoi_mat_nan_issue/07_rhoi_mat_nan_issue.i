##############################################################################
# File: 07_rhoi_mat_nan_issue.i
# File Location: /examples/sintering/paper3/01_test_inputs/07_rhoi_mat_nan_issue
# Created Date: Thursday November 2nd 2023
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Thursday November 2nd 2023
# Modified By: Brandon Battas
# -----
# Description:
#  Debugging why suddenly the rho_i_mat is displaying only nan in paraview
#  while the pp average or integral of the value is a normal number???
#
#
##############################################################################

[Mesh]
  [gmg]
    type = DistributedRectilinearMeshGenerator
    dim = 2
    nx = 90
    ny = 80
    xmin = 0
    xmax = 90
    ymin = 0
    ymax = 80
  []
  [subdomains]
    type = ParsedSubdomainMeshGenerator
    input = gmg
    combinatorial_geometry = 'x > 80'
    block_id = 1
  []
[]

[GlobalParams]
  op_num = 8
  var_name_base = gr
  int_width = 6 #particle radius is 100
  # profile = TANH # not used at the moment? only in circleic?
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
  [gr4]
  []
  [gr5]
  []
  [gr6]
  []
  [gr7]
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
  [unique_grains]
    order = CONSTANT
    family = MONOMIAL
  []
  [var_indices]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[ICs]
  [PolycrystalICs]
    [PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
      block = 0
    []
  []
  [VoidIC]
    # type = ConstantIC
    # block = 1
    # variable = phi
    # value = 1
    type = BoundingBoxIC
    variable = phi
    block = 1
    inside = 1
    outside = 0
    x1 = 80
    x2 = 100
    y1 = -10
    y2 = 90
  []
  # [VoidIC2]
  #   type = ConstantIC
  #   block = 0
  #   variable = phi
  #   value = 1e-4 #0.02
  # []
  [BubbleIC]
    type = MultiSmoothCircleIC
    variable = phi
    invalue = 1
    outvalue = 0.01
    numbub = 5
    radius = 8 #4
    bubspac = 12
    block = 0
    numtries = 10000
  []
  # [Rho_i_IC]
  #   type = ConstantIC
  #   block = 0
  #   variable = rho_i
  #   value = 0 #1e-4 #0.02
  # []
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
  [rhov_avg_func]
    type = ParsedFunction
    symbol_names = 'rhov_pp_avg'
    symbol_values = 'rhov_pp_avg'
    expression = 'rhov_pp_avg'
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
    outputs = nemesis #exodus
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
  # Irradiation and Interstitials
  [rho_testing]
    type = ParsedMaterial
    property_name = rho_testing
    material_property_names = 'rhos rhov'
    expression = 'rhov'
    outputs = 'nemesis'
  []
  [sum_grs]
    type = ParsedMaterial
    property_name = sum_grs
    coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'
    expression = 'gr0 + gr1 + gr2 + gr3 + gr4 + gr5 + gr6 + gr7'
    outputs = 'nemesis'
  []
  # [f_dot]
  #   type = RandomMaterial
  #   f_name = f_dot
  #   noise = noise
  #   outputs = exodus
  #   min_value = 0.8e-8
  #   max_value = 1.2e-8
  # []
  # Interstitials
  [rho_gen]
    type = DerivativeParsedMaterial
    property_name = rho_gen
    derivative_order = 1
    constant_names = 'Nc Nd f_dot noise'
    constant_expressions = '2 5 1e-8 1'
    material_property_names = 'hs'
    expression = 'f_dot * noise * Nc * Nd' # * hs
    outputs = 'nemesis'
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
    outputs = 'nemesis'
  []
  [rho_i_dpm]
    type = DerivativeParsedMaterial
    property_name = rho_i_dpm
    derivative_order = 2
    # coupled_variables = ''
    material_property_names = 'rho_gen a_r'
    postprocessor_names = 'rhov_pp_avg'
    expression = 'rho_gen / a_r' # (a_r * rhov_pp_avg)'
    outputs = 'nemesis'
  []
  # Vacancies
  [rho_v_recombRate]
    type = DerivativeParsedMaterial
    property_name = rho_v_recombRate
    coupled_variables = 'w'
    # additional_derivative_symbols = w
    material_property_names = 'rho_i_dpm a_r rhov'
    # postprocessor_names = 'rhov_pp_avg'
    expression = 'a_r * rhov * rho_i_dpm'
    outputs = 'nemesis'
  []
  # [source]
  #   type = DerivativeParsedMaterial
  #   f_name = source
  #   args = 'phi'
  #   material_property_names = 'hm(phi) f_dot'
  #   constant_names = 'X' #defects generated per cascade
  #   constant_expressions = '3'
  #   function = '2*hm*f_dot*X'
  #   derivative_order = 1
  #   outputs = nemesis
  # []
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
    v = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7' # gr8 gr9 gr10 gr11 gr12 gr13 gr14 gr15'
    gamma = gamma
    mob_name = L_mat
  []
  [kappa_phi]
    type = ACKappaFunction
    variable = phi
    mob_name = L_mat
    kappa_name = kappa
  []
  # Irradiation
  # Source/Generation
  # [source_w]
  #   type = MaskedBodyForce
  #   variable = w
  #   mask = rho_gen
  # []
  # # Recombination/sink
  # [recombination_w]
  #   type = MatReaction
  #   variable = w
  #   mob_name = rho_v_recombRate
  #   # args = rhoi #but its a constant and material not a variable
  # []
  # Damage
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
[]

[Postprocessors]
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
  [total_phi]
    type = ElementIntegralVariablePostprocessor
    variable = phi
    outputs = csv
  []
  [total_rhov]
    type = ElementIntegralMaterialProperty
    mat_prop = rhov
    outputs = csv
  []
  [total_rhos]
    type = ElementIntegralMaterialProperty
    mat_prop = rhos
    outputs = csv
  []
  [total_rhoi]
    type = ElementIntegralMaterialProperty
    mat_prop = rho_i_dpm
    outputs = csv
  []
  [total_grs]
    type = ElementIntegralMaterialProperty
    mat_prop = sum_grs
    outputs = csv
  []
  [rhov_pp_avg]
    # type = ElementAverageValue
    type = ElementAverageMaterialProperty
    mat_prop = rhov
    execute_on = 'initial TIMESTEP_BEGIN'
  []
[]

# [VectorPostprocessors]
# #  [./voids]
# #    type = FeatureVolumeVectorPostprocessor
# #    flood_counter = void_tracker
# #    execute_on = 'initial timestep_end final'
# #    output_centroids = false  #was true
# #    outputs = csv
# #  [../]
#  [./vectorMemory]
#    type = VectorMemoryUsage
#    mem_units = gigabytes
#    outputs = csv
#    execute_on = 'NONLINEAR LINEAR TIMESTEP_END'
#  [../]
# []

[UserObjects]
  [grain_tracker]
    type = GrainTracker
    threshold = 0.1 #0.2
    connecting_threshold = 0.09 #0.08
    compute_halo_maps = false #true#false
  []
  [voronoi]
    type = PolycrystalVoronoi
    grain_num = 8 # Number of grains
    rand_seed = 10
    # int_width = 7 # global param
  []
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
  petsc_options_iname = '-pc_type -pc_hypre_type' # -snes_type'
  petsc_options_value = 'hypre boomeramg' # vinewtonrsls'
  # petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap'
  # petsc_options_value = ' asm      lu           2'
  nl_max_its = 20 #40 too large- optimal_iterations is 6
  l_max_its = 30 #if it seems like its using a lot it might still be fine
  l_tol = 1e-04
  nl_rel_tol = 1e-6 #default is 1e-8
  nl_abs_tol = 1e-6 #only needed when near equilibrium or veeeery small timesteps and things changing FAST
  start_time = 0
  # end_time = 10 #0.006
  steady_state_detection = true
  num_steps = 4
  # dt = 0.00002
  # dtmax = 500
  # dt = 0.0001
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 6
    dt = 0.00002
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
  # [checkpoint]
  #   type = Checkpoint
  #   num_files = 3
  # []
[]
