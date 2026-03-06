##############################################################################
# File: 01_sub_poly10gr.i
# File Location: /examples/agg/05_large_bicr_L/00_sub
# Created Date: Friday March 6th 2026
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Friday March 6th 2026
# Modified By: Brandon Battas
# -----
# Description:
#  Sub file to pass second order polycrystal variables for main IC
#
#
#
##############################################################################

# i_tol = 0
# a_tol = 5
# amag = 0.0

[Mesh]
  [ebsd_mesh]
    type = EBSDMeshGenerator
    filename = '2D_100x100_10gr.txt'
    # filename = '../../00_sub/2D_100x100_10gr.txt'
  []
  parallel_type = DISTRIBUTED
  second_order = true
  # uniform_refine = 1
[]

[GlobalParams]
  op_num = 10
  var_name_base = gr
[]

[Variables]
  [PolycrystalVariables]
    # order = SECOND
  []
[]

[ICs]
  [PolycrystalICs]
    [PolycrystalColoringIC]
      polycrystal_ic_uo = ebsd
    []
  []
[]

[UserObjects]
  [ebsd_reader]
    type = EBSDReader
  []
  [ebsd]
    type = PolycrystalEBSD
    coloring_algorithm = bt #jp
    ebsd_reader = ebsd_reader
    enable_var_coloring = true
  []
  [grain_tracker]
    type = GrainTracker
    # variable = 'gr0 gr1 gr2 gr3 gr4'
    threshold = 0.01
    connecting_threshold = 0.01
    compute_var_to_feature_map = false
    compute_halo_maps = false # Only necessary for displaying HALOS
    execute_on = 'initial timestep_end'
    halo_level = 6
    # use_less_than_threshold_comparison = true
    # polycrystal_ic_uo = ebsd
  []
[]

[AuxVariables]
  [gr0a]
    order = SECOND
  []
  [gr1a]
    order = SECOND
  []
  [gr2a]
    order = SECOND
  []
  [gr3a]
    order = SECOND
  []
  [gr4a]
    order = SECOND
  []
  [gr5a]
    order = SECOND
  []
  [gr6a]
    order = SECOND
  []
  [gr7a]
    order = SECOND
  []
  [gr8a]
    order = SECOND
  []
  [gr9a]
    order = SECOND
  []
[]

[Kernels]
  [PolycrystalKernel]
    variable_mobility = false
  []
[]

[AuxKernels]
  [auxgr0]
    type = ProjectionAux
    variable = gr0a
    v = gr0
  []
  [auxgr1]
    type = ProjectionAux
    variable = gr1a
    v = gr1
  []
  [auxgr2]
    type = ProjectionAux
    variable = gr2a
    v = gr2
  []
  [auxgr3]
    type = ProjectionAux
    variable = gr3a
    v = gr3
  []
  [auxgr4]
    type = ProjectionAux
    variable = gr4a
    v = gr4
  []
  [auxgr5]
    type = ProjectionAux
    variable = gr5a
    v = gr5
  []
  [auxgr6]
    type = ProjectionAux
    variable = gr6a
    v = gr6
  []
  [auxgr7]
    type = ProjectionAux
    variable = gr7a
    v = gr7
  []
  [auxgr8]
    type = ProjectionAux
    variable = gr8a
    v = gr8
  []
  [auxgr9]
    type = ProjectionAux
    variable = gr9a
    v = gr9
  []
[]

# [BCs]
#   [Periodic]
#     [All]
#       auto_direction = 'x y'
#     []
#   []
# []

[Materials]
  [constants]
    type = GenericConstantMaterial
    prop_names = 'L0             kappa_op  int_width_iso   const_m      sigma0' #'gamma_asymm'
    prop_values = '1.15382e-6  2.07337e7        6      5.521269e6  4.60748e6' #'1.5    '
  []
  [constants2]
    type = GenericConstantMaterial
    prop_names = 'int_width  gamma_asymm  sigma      L    mu'
    prop_values = '6             1.5    4.6075e6 1.153820e-6  5.521269e6'
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
  num_steps = 1
  dt = 1e-6
  automatic_scaling = true
  compute_scaling_once = false
[]

[Outputs]
  # file_base = MLPaperActualSimulations/Case3
  exodus = false
  #nemesis = true
  console = true
  csv = false
  checkpoint = false
[]
