##############################################################################
# File: 04_4cpu_main.i
# File Location: /examples/agg/01_initial_testing/37_suo_parallelization_ic/04_4cpu_main
# Created Date: Monday August 10th 2026
# Author: Battas,Brandon Scott (bbattas@ufl.edu)
# -----
# Last Modified: Tuesday August 25th 2026
# Modified By: Battas,Brandon Scott
# -----
# Description:
#  Using the 4 cpu base for 4 here
#
#
# INFO: Finished. Runtime: 0h 0m 13.32s | Peak memory delta: 208208.0 MB
##############################################################################

[Mesh]
  [gmg]
    type = DistributedRectilinearMeshGenerator
    dim = 2 # Problem dimension
    nx = 12 # Number of elements in the x-direction
    ny = 12 # Number of elements in the y-direction
    nz = 0 # Number of elements in the z-direction
    xmin = 0 # minimum x-coordinate of the mesh
    xmax = 1000 # maximum x-coordinate of the mesh
    ymin = 0 # minimum y-coordinate of the mesh
    ymax = 1000 # maximum y-coordinate of the mesh
    zmin = 0
    zmax = 0
    elem_type = QUAD4 # Type of elements used in the mesh
  []
  uniform_refine = 3 # Initial uniform refinement of the mesh
  parallel_type = DISTRIBUTED # Periodic BCs
[]

[GlobalParams]
  # Parameters used by several kernels that are defined globally to simplify input file
  op_num = 15 # Number of grains
  var_name_base = gr # Base name of grains
[]

[UserObjects]
  # [voronoi]
  #   type = PolycrystalVoronoi
  #   grain_num = 15
  #   rand_seed = 42
  #   coloring_algorithm = bt # We must use bt to force the UserObject to assign one grain to each op
  # []
  [solution_userobject]
    type = SolutionUserObject
    mesh = '../03_4cpu_base/03_4cpu_base_out.e'
    timestep = 'LATEST'
    system_variables = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 gr8 gr9 gr10 gr11 gr12 gr13 gr14'
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
[]

[Variables]
  # Variable block, where all variables in the simulation are declared
  [PolycrystalVariables]
    # Custom action that created all of the grain variables
    order = FIRST # element type used by each grain variable
    family = LAGRANGE
  []
[]

[AuxVariables]
  #active = ''
  # Dependent variables
  [bnds]
    # Variable used to visualize the grain boundaries in the simulation
    order = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  # Kernel block, where the kernels defining the residual equations are set up.
  [PolycrystalKernel]
    # Custom action creating all necessary kernels for grain growth.  All input parameters are up in GlobalParams
  []
[]

[AuxKernels]
  #active = ''
  # AuxKernel block, defining the equations used to calculate the auxvars
  [bnds_aux]
    # AuxKernel that calculates the GB term
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  []
[]

[BCs]
  # Boundary Condition block
  [Periodic]
    [top_bottom]
      auto_direction = 'x y' # Makes problem periodic in the x and y directions
    []
  []
[]

[Materials]
  [CuGrGr]
    # Material properties
    type = GBEvolution # Quantitative material properties for copper grain growth.  Dimensions are nm and ns
    GBmob0 = 2.5e-6 #Mobility prefactor for Cu from schonfelder1997molecular bibtex entry
    GBenergy = 0.708 #GB energy for Cu from schonfelder1997molecular bibtex entry
    Q = 0.23 #Activation energy for grain growth from Schonfelder 1997
    T = 450 # K   #Constant temperature of the simulation (for mobility calculation)
    wGB = 14 # nm      #Width of the diffuse GB
  []
[]

[Postprocessors]
  active = 'dt '
  # Scalar postprocessors
  [dt]
    # Outputs the current time step
    type = TimestepSize
  []
[]

[Executioner]
  type = Transient # Type of executioner, here it is transient with an adaptive time step
  scheme = bdf2 # Type of time integration (2nd order backward euler), defaults to 1st order backward euler

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -mat_mffd_type'
  petsc_options_value = 'hypre    boomeramg      101                ds'

  l_max_its = 30 # Max number of linear iterations
  l_tol = 1e-4 # Relative tolerance for linear solves
  nl_max_its = 40 # Max number of nonlinear iterations
  nl_abs_tol = 1e-11 # Relative tolerance for nonlienar solves
  nl_rel_tol = 1e-8 # Absolute tolerance for nonlienar solves

  start_time = 0.0
  end_time = 100 #4000
  dt = 25

  # [TimeStepper]
  #   type = IterationAdaptiveDT
  #   dt = 25 # Initial time step.  In this simulation it changes.
  #   optimal_iterations = 6 #Time step will adapt to maintain this number of nonlinear iterations
  # []

  # [Adaptivity]
  #   # Block that turns on mesh adaptivity. Note that mesh will never coarsen beyond initial mesh (before uniform refinement)
  #   initial_adaptivity = 2 # Number of times mesh is adapted to initial condition
  #   refine_fraction = 0.7 # Fraction of high error that will be refined
  #   coarsen_fraction = 0.1 # Fraction of low error that will coarsened
  #   max_h_level = 4 # Max number of refinements used, starting from initial mesh (before uniform refinement)
  # []
[]

[Outputs]
  exodus = true
  csv = true
  [console]
    type = Console
    max_rows = 20
  []
  [pgraph]
    type = PerfGraphOutput
    execute_on = 'initial final' # Default is "final"
    level = 2 # Default is 1
    heaviest_branch = true # Default is false
    heaviest_sections = 7 # Default is 0
  []
[]
