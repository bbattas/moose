##############################################################################
# File: 10_ebsd_conc.i
# File Location: /examples/sintering/test_inputs/10_ebsd_conc
# Created Date: Monday June 17th 2024
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Tuesday June 18th 2024
# Modified By: Brandon Battas
# -----
# Description:
#  Simple test for using EBSD reader to read in concentrations for mass cons
#  version of GP Sintering
#
#
##############################################################################

[Problem]
  type = FEProblem
  solve = false
  kernel_coverage_check = false
[]

[Mesh]
  [ebsd_mesh]
    type = EBSDMeshGenerator
    filename = 00_d3d_txt/2D_20x20nm_cvci_irr.txt
  []
  [subdomain_external]
    type = ParsedSubdomainMeshGenerator
    input = ebsd_mesh
    combinatorial_geometry = 'x > 18'
    block_id = 1
  []
  parallel_type = DISTRIBUTED
  uniform_refine = 0
[]

[GlobalParams]
  op_num = 3
  var_name_base = gr
  int_width = 2 #1000 #min radius is like 2250, element size of 250
  # profile = TANH # not used at the moment? only in circleic?
[]

[Variables]
  [wvac]
  []
  [cvac]
  []
  [cint]
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
    initial_condition = 1600
  []
  [voids]
    order = CONSTANT
    family = MONOMIAL
  []
  [unique_grains]
    order = CONSTANT
    family = MONOMIAL
  []
  [c_aux]
    order = CONSTANT
    family = MONOMIAL
  []
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
    x1 = 20
    x2 = 30
    y1 = -10
    y2 = 30
  []
  [voidIC2]
    type = SpecifiedSmoothCircleIC
    variable = phi
    invalue = 1
    outvalue = 0.01
    radii = '1.56 1.98'
    x_positions = '12.2  9.17'
    y_positions = '16.05 11.54'
    z_positions = '0     0'
    block = 0
  []
  [c_ic]
    type = EBSDConsVarIC
    variable = cvac
    ebsd_reader = ebsd_reader
    column = 0
  []
  [ci_ic]
    type = EBSDConsVarIC
    variable = cint
    ebsd_reader = ebsd_reader
    column = 1
  []
[]

# [BCs]
#   [phi]
#     type = NeumannBC
#     variable = phi
#     value = 0
#     boundary = 'left right top bottom'
#   []
#   [gr0]
#     type = NeumannBC
#     variable = gr0
#     value = 0
#     boundary = 'left right top bottom'
#   []
#   [gr1]
#     type = NeumannBC
#     variable = gr1
#     value = 0
#     boundary = 'left right top bottom'
#   []
#   [gr2]
#     type = NeumannBC
#     variable = gr2
#     value = 0
#     boundary = 'left right top bottom'
#   []
# []

[AuxKernels]
  [bnds_aux]
    type = BndsCalcAux
    variable = bnds
    execute_on = 'initial timestep_end'
  []
  #NEW Grain_Tracker related stuff
  [unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
    execute_on = 'initial timestep_end'
  []
  # [test_aux]
  #   type = EBSDReaderPointDataAux
  #   variable = c_aux
  #   ebsd_reader = ebsd_reader
  #   data_name = "CUSTOM0"
  #   execute_on = 'initial'
  # []
[]

[UserObjects]
  [ebsd_reader]
    type = EBSDReader
    custom_columns = 2
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
    compute_halo_maps = true #true#false
    verbosity_level = 1
  []
[]

# [Preconditioning]
#   [SMP] #slow but good, very slow for 3D (might be another option then)
#     type = SMP
#     full = true
#     # coupled_groups = 'wvac,wint,phi'
#   []
# []

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = PJFNK
  # petsc_options_iname = '-pc_type -sub_pc_type -pc_factor_levels'
  # petsc_options_value = 'asm ilu 2'
  # petsc_options_iname = '-pc_type -pc_hypre_type' # -snes_type'
  # petsc_options_value = 'hypre boomeramg' # vinewtonrsls'
  petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -sub_pc_factor_shift_type'
  petsc_options_value = ' asm      lu           2                nonzero'
  nl_max_its = 30 #12 #20 #40 too large- optimal_iterations is 6
  l_max_its = 60 #200 #30 #200 #30 #if it seems like its using a lot it might still be fine
  l_tol = 1e-06 #4
  nl_rel_tol = 1e-8 #6 #default is 1e-8
  # nl_abs_tol = 1e-14 #only needed when near equilibrium or veeeery small timesteps and things changing FAST
  start_time = 0
  # end_time = 1e6 #1e10 #5e6 #0.006
  num_steps = 1 #00
  # steady_state_detection = true
  # # From tonks ode input
  automatic_scaling = true
  compute_scaling_once = false
  line_search = none
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 6
    dt = 0.001 #2.5
    # linear_iteration_ratio = 1e5 #needed with large linear number for asmilu
  []
[]

[Outputs]
  perf_graph = false
  csv = true
  exodus = true
  checkpoint = false
  print_linear_residuals = false
[]

# [Debug]
#   show_var_residual_norms = true
# []
