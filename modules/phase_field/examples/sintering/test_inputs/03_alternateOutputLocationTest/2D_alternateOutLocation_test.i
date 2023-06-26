##############################################################################
# File: 2D_alternateOutLocation_test.i
# File Location: /examples/sintering/test_inputs/03_alternateOutputLocationTest/2D_alternateOutLocation_test.i
# Created Date: Monday June 26th 2023
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Monday June 26th 2023
# Modified By: Brandon Battas
# -----
# Description:
#  Testing the ability to write output to a different location and
#  the naming specifics needed for that (goal is to use this with MOOSE on
#  HPG to keep inputs in ~/moose but output results to /blue)
#  https://github.com/idaholab/moose/discussions/22986
#  WORKs- to do outside base moose directory you have to cd to the output directory
#  and then run it from there!
##############################################################################

[Mesh]
  [gmg]
    type = DistributedRectilinearMeshGenerator
    dim = 2
    nx = 60
    ny = 60
    xmin = 0
    xmax = 600
    ymin = 0
    ymax = 600
  []
  parallel_type = DISTRIBUTED
  uniform_refine = 0
  second_order = false
[]

[GlobalParams]
  op_num = 2
  var_name_base = gr
  int_width = 20
  profile = TANH
[]

[Variables]
  [w]
    # order = SECOND
  []
  [phi]
    # order = SECOND
  []
  [PolycrystalVariables]
    # order = SECOND
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
  #  [./voids]
  #    order = CONSTANT
  #    family = MONOMIAL
  #  [../]
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
    x_positions = '320 600'
    y_positions = '600 320'
    z_positions = '  0   0'
    radii = '200 200'
    invalue = 0
    outvalue = 1
  []
  [gr0_IC]
    type = SmoothCircleIC
    variable = gr0
    x1 = 320
    y1 = 600
    z1 = 0
    radius = 200
    invalue = 1
    outvalue = 0
  []
  [gr1_IC]
    type = SmoothCircleIC
    variable = gr1
    x1 = 600
    y1 = 320
    z1 = 0
    radius = 200
    invalue = 1
    outvalue = 0
  []
[]

[BCs]
  [phi_bc]
    type = NeumannBC
    variable = 'phi'
    boundary = 'right top'
    value = 0
  []
  [gr0_bc]
    type = NeumannBC
    variable = 'gr0'
    boundary = 'right top'
    value = 0
  []
  [gr1_bc]
    type = NeumannBC
    variable = 'gr1'
    boundary = 'right top'
    value = 0
  []
[]

[Functions]
  [f_T]
    type = ConstantFunction
    value = 1800
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
    f_name = ks
    args = 'T'
    constant_names = 'a b'
    constant_expressions = '-0.0025 157.16'
    function = 'a*T + b'
  []
  [kv]
    type = ParsedMaterial
    f_name = kv
    material_property_names = 'ks'
    function = '10*ks'
  []
  # Diffusivity and mobilities
  [chiD]
    type = GrandPotentialIsoMaterial
    f_name = chiD
    solid_mobility = L #CHANGED FROM L
    void_mobility = Lv
    chi = chi
    surface_energy = 19.7
    c = phi
    T = T
    D0 = 8.33e9
    GBmob0 = 1.4759e9
    Q = 2.77
    Em = 3.608
    bulkindex = 1
    gbindex = -1
    surfindex = 1e11
    # iw_scaling = TRUE
  []
  # Should work without gr# in args?
  # [./cv_eq_old]
  #   type = DerivativeParsedMaterial
  #   f_name = cv_eq_old
  #   args = 'T phi bnds'
  #   constant_names = 'Ef_b Ef_gb kB'
  #   constant_expressions = '2.6551964327999986 0.5296941200355877 8.617343e-5'
  #   derivative_order = 2
  #   function = 'c_B:=exp(-Ef_b/kB/T);
  #               c_GB:=exp(-Ef_gb/kB/T);
  #               bounds:=bnds + phi^2;
  #               c_B + 4.0 * (c_GB - c_B) * (1.0 - bounds)^2'
  #   outputs = 'none'
  # [../]
  [cv_eq]
    type = UO2CvMaterial
    f_name = cv_eq
    T = T
    c = phi
    OU = OU
    gb_se_csv = '/home/bbattas/projects/moose/modules/phase_field/examples/sintering/gb_segregation_csv/sigma9_se.csv
                 /home/bbattas/projects/moose/modules/phase_field/examples/sintering/gb_segregation_csv/sigma11_se.csv'
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
    f_name = c
    material_property_names = 'hs rhos hv rhov'
    constant_names = 'Va'
    constant_expressions = '0.04092'
    function = 'Va*(hs*rhos + hv*rhov)'
    outputs = none #exodus
  []
  # [./L_mat1]
  #   type = ParsedMaterial
  #   # type = DerivativeParsedMaterial
  #   f_name = L_mat1
  #   material_property_names = 'L Lv hv'
  #   args = 'phi'
  #   constant_names = 'p0'
  #   constant_expressions = '1.0'
  #   function = 'hv*Lv + (1-hv)*L'
  #   outputs = nemesis #none
  # [../]
  [L_mat]
    # type = ParsedMaterial
    type = DerivativeParsedMaterial
    f_name = L_mat
    material_property_names = 'L Lv'
    args = 'phi'
    constant_names = 'p0'
    constant_expressions = '0.3'
    function = 'hv:=if(phi<=0.0,0.0,if(phi>=p0,1.0,6*(phi/p0)^5 - 15*(phi/p0)^4 + 10*(phi/p0)^3));
                hv*Lv + (1-hv)*L'
    # function = '0.1*Lv'
    outputs = none #none
  []
  # [./L_mat3]
  #   # type = ParsedMaterial
  #   type = DerivativeParsedMaterial
  #   f_name = L_mat3
  #   material_property_names = 'L Lv'
  #   args = 'phi'
  #   constant_names = 'p0'
  #   constant_expressions = '0.3'
  #   function = 'hv:=6*(phi/p0)^5 - 15*(phi/p0)^4 + 10*(phi/p0)^3;
  #               hv*Lv + (1-hv)*L'
  #   outputs = nemesis #none
  # [../]
  # [Diff]
  #   type = ParsedMaterial
  #   # type = DerivativeParsedMaterial
  #   f_name = Diff
  #   material_property_names = 'diffusivity'
  #   function = 'diffusivity'
  #   outputs = nemesis #none
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
      mobility_name_gr = L_mat #_mat2 #CHANGED FROM L
      kappa_gr = kappa
      free_energies_gr = 'omegav omegas'

      additional_ops = 'phi'
      gamma_grxop = gamma
      mobility_name_op = L_mat #_mat2
      kappa_op = kappa
      free_energies_op = 'omegav omegas'
    []
  []
[]

[Kernels]
  [barrier_phi]
    type = ACBarrierFunction
    variable = phi
    v = 'gr0 gr1' #gr2 gr3'# gr4 gr5'# gr6 gr7 gr8 gr9 gr10 gr11 gr12 gr13 gr14 gr15'
    gamma = gamma
    mob_name = L_mat #_mat2 #????
  []
  [kappa_phi]
    type = ACKappaFunction
    variable = phi
    mob_name = L_mat #_mat2
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
  #  [./voids_aux]
  #    type = FeatureFloodCountAux
  #    variable = voids
  #    flood_counter = void_tracker
  #    field_display = UNIQUE_REGION
  #    execute_on = 'INITIAL TIMESTEP_END'
  #  [../]
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
  #  [./void_tracker]
  #    type = FeatureFloodCount
  #    variable = phi
  #    threshold = 0.6
  #    connecting_threshold = 0.5 #was 0.2 and worked fine for iso not tensor
  #    compute_var_to_feature_map = true
  #    execute_on = 'initial timestep_end'
  #  [../]
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
  # [./terminator]  #do i have to specify that this is off so that the control can turn it on?
  #   type = Terminator
  #   expression = 'void_tracker = 1'
  #   execute_on = TIMESTEP_END
  #   enable = true
  # [../]
  [terminator] #do i have to specify that this is off so that the control can turn it on?
    type = Terminator
    expression = 'grain_tracker < 2'
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
  l_tol = 1e-4
  nl_rel_tol = 1e-6 #default is 1e-8
  nl_abs_tol = 1e-6 #only needed when near equilibrium or veeeery small timesteps and things changing FAST
  start_time = 0
  # end_time = 5000
  steady_state_detection = true
  num_steps = 1
  # dt = 0.0001
  # dtmax = 500
  # dt = 0.001
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 8 #WAS 6
    dt = 0.0001
    growth_factor = 1.2
    cutback_factor = 0.8
    cutback_factor_at_failure = 0.5 #might be different from the curback_factor
  []
  [Adaptivity]
    refine_fraction = 0.8
    coarsen_fraction = 0.05 #minimize this- adds error
    max_h_level = 2 #test a short simulation with 1,2,3,4 for this to see where it stops helping
    initial_adaptivity = 2
  []
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
  # print_linear_residuals = false
  [checkpoint]
    type = Checkpoint
    num_files = 3
  []
  # file_base = /home/bbattas/projects/moose/modules/phase_field/examples/sintering/test_inputs/03.1_alternateOutputLocationTest/output_name
  file_base = /home/bbattas/Documents/moose_test/alternateOutputLocation/outputname
[]
