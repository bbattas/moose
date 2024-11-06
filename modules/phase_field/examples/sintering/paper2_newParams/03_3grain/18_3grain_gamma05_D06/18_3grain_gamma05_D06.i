##############################################################################
# File: 18_3grain_gamma05_D06.i
# File Location: /examples/sintering/paper2_newParams/03_3grain/18_3grain_gamma05_D06
# Created Date: Monday October 14th 2024
# Author: Battas,Brandon Scott (bbattas@ufl.edu)
# -----
# Last Modified: Wednesday November 6th 2024
# Modified By: Brandon Battas
# -----
# Description:
#  NEW: finer mesh 3 grain input with intermediary D
#   Dgb = 5.8422e7, Ds = iw*Dgb/0.6 = 1.94740e9, sigma_s = 2*sigma_gb
#
#  OLD: Additional 3 grain input with intermediary delDgb/Ds
#  Full 3 grain structure, del*Dgb/Ds = 0.6, sigma_gb/sigma_s = 0.5
#  delta=iw but on scaled D internal value, so backcalc Dgb = 1e6, Ds = 3.33e7
#  sigma_gb=sig5tiltgb=9.86 eV/nm2, sig_s = 2*sig_gb = 19.7
#  02 but with different gamma/D
#  cut down the end time and refined the mesh a bit compared to gamma=1
#  now elements are 4 in xy and 4.0254... (475/118) in z instead of 5
#  125x75x118 had 5.7m DoFs-> 300(19k/ea), 450(12.7k/ea)
##############################################################################

[Mesh]
  [gmg]
    type = DistributedRectilinearMeshGenerator
    dim = 3
    nx = 125
    ny = 75
    nz = 118
    xmin = 0
    xmax = 500
    ymin = 0
    ymax = 300
    zmin = 0
    zmax = 475
  []
  # parallel_type = DISTRIBUTED
  # uniform_refine = 1
  # second_order = false
[]

[GlobalParams]
  op_num = 3
  var_name_base = gr
  int_width = 20 #particle radius is 100
  # profile = TANH
[]

[MultiApps]
  [NMC_1step]
    type = FullSolveMultiApp
    execute_on = initial
    positions = '0 0 0'
    input_files = ../../00_sub/04_sub_3grain_fine.i
  []
[]

[Transfers]
  [from_NMC]
    type = MultiAppCopyTransfer
    from_multi_app = NMC_1step
    source_variable = 'gr0 gr1 gr2 phi cvac_aux_th'
    variable = 'gr0 gr1 gr2 phi cvac_var'
  []
[]

[Variables]
  [w]
    initial_condition = 0.0
  []
  [cvac_var]
  []
  [phi]
  []
  [gr0]
  []
  [gr1]
  []
  [gr2]
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

[Materials]
  # Free energy coefficients for parabolic curves
  [k_constants]
    type = GenericConstantMaterial
    prop_names = 'ks kv' # Using the GB based values (lowest of mine)
    # prop_values = '7.751e2 7.751e2' # Irradiation
    prop_values = '6.569e2 6.569e2' # No Irradiation
  []
  [gb_e_mat] # eV/nm^2
    type = ParsedMaterial
    property_name = gb_e_mat
    coupled_variables = 'T'
    constant_names = 'a b c'
    constant_expressions = '-5.87e-4 1.56 6.24151'
    expression = 'c*(a*T + b)'
    outputs = none
  []
  [surf_e_mat] # eV/nm^2
    type = ParsedMaterial
    property_name = surf_e_mat
    material_property_names = gb_e_mat
    expression = '2 * gb_e_mat'
    outputs = none
  []
  [hgb]
    type = SwitchingFunctionGBMaterial
    h_name = hgb
    grain_ops = 'gr0 gr1 gr2'
    hgb_threshold = 0.0 #0.0001
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
    GBmob0 = 3.42828e10 # nm4/eVs #1.4759e9 # new value from Tonks/PC/Jake GG Paper
    Q = 3.01 #2.77 # new value from Tonks/PC/Jake GG Paper
    D0 = 4.2488e11 #8.33e9
    Em = 4.23317 #3.608
    bulkindex = 1
    gbindex = 5.8422e7 # GB avg from 9/11 at 1600 # -1 sets the GB D to the LANL MD Value in GPIsoMat
    surfindex = 1.94740e9 #iw*DGB/0.6  #5.8422e10
    GBwidth = 3.5 # based on avg of two lanl values
    surf_thickness = 3.5 # keeping equal to gb for simplicity
    iw_scaling = true
  []
  [cv_eq]
    type = UO2CeqMaterial
    ceq_name = cv_eq
    hgb = hgb
    vi = VAC_TH
  []
  [sintering]
    type = GrandPotentialSinteringMaterial
    chemical_potential = w
    void_op = phi
    Temperature = T
    surface_energy = surf_e_mat #9.86 #19.7
    grainboundary_energy = gb_e_mat #9.86
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
    [GrandPotentialAlt]
      switching_function_names = 'hv hs'
      # Chempot
      chemical_potentials = 'w'
      mobilities = 'chiD' #cons_mob
      anisotropic = 'false'
      susceptibilities = 'chi'
      free_energies_w = 'rhov rhos '
      # Grains
      mobility_name_gr = L_mat
      kappa_gr = kappa
      gamma_gr = gamma
      free_energies_gr = 'omegav omegas'
      # Other OPs (Phi)
      additional_ops = 'phi'
      mobility_name_op = L_mat
      kappa_op = kappa
      gamma_grxop = gamma
      free_energies_op = 'omegav omegas' #empty when no phi'omegaa omegab'
      # Mass Conservation
      mass_conservation = true
      concentrations = 'cvac_var'
      hj_over_kVa = 'hv_over_kVa hs_over_kVa'
      hj_c_min = 'hv_c_min hs_c_min'
    []
  []
[]

[Kernels]
  [barrier_phi]
    type = ACBarrierFunction
    variable = phi
    v = 'gr0 gr1 gr2' #gr2 gr3'# gr4 gr5'# gr6 gr7 gr8 gr9 gr10 gr11 gr12 gr13 gr14 gr15'
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

[VectorPostprocessors]
  [ctrline]
    type = LineValueSampler
    variable = phi
    start_point = '250 0 207.5'
    end_point = '250 300 207.5'
    sort_by = y
    num_points = 31
    outputs = csv
  []
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
[]

[UserObjects]
  # [./terminator]  #do i have to specify that this is off so that the control can turn it on?
  #   type = Terminator
  #   expression = 'void_tracker = 1'
  #   execute_on = TIMESTEP_END
  #   enable = true
  # [../]
  [terminator] #do i have to specify that this is off so that the control can turn it on?
    type = Terminator
    expression = 'grain_tracker < 3'
    execute_on = TIMESTEP_END
    # enable = true
  []
  [grain_tracker]
    type = GrainTracker
    threshold = 0.1 #0.2
    connecting_threshold = 0.09 #0.08
    # compute_halo_maps = false #true#false
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

[Preconditioning]
  [./SMP] #slow but good, very slow for 3D (might be another option then)
    type = SMP
    coupled_groups = 'w,phi'
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = PJFNK
  # petsc_options_iname = '-pc_type -pc_hypre_type'
  # petsc_options_value = 'hypre boomeramg'
  petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -sub_pc_factor_shift_type'
  petsc_options_value = ' asm      lu           2                nonzero'
  nl_max_its = 20 #40 too large- optimal_iterations is 6
  l_max_its = 60 #if it seems like its using a lot it might still be fine
  l_tol = 1e-06#4
  nl_rel_tol = 1e-6 #default is 1e-8
  # nl_abs_tol = 1e-6 #only needed when near equilibrium or veeeery small timesteps and things changing FAST
  start_time = 0
  end_time = 7500
  steady_state_detection = true
  automatic_scaling = true
  compute_scaling_once = false
  # num_steps = 2
  # dt = 0.0001
  dtmax = 50
  # dt = 0.0001
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 6 #WAS 6
    dt = 0.01
    growth_factor = 1.2
    cutback_factor = 0.8
    cutback_factor_at_failure = 0.5 #might be different from the curback_factor
    linear_iteration_ratio = 1e5
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

