##############################################################################
# File: 04_sub_3grain_fine.i
# File Location: /examples/sintering/paper2_newParams/00_sub
# Created Date: Monday October 14th 2024
# Author: Battas,Brandon Scott (bbattas@ufl.edu)
# -----
# Last Modified: Monday October 14th 2024
# Modified By: Battas,Brandon Scott
# -----
# Description:
#  Finer mesh i used on the gamma 0.5 inputs for 5 elements across interface
#   instead of 4
#
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
[]

[AuxVariables]
  [T]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 1600
  []
  # THEQ
  # Vacancy mat to variable
  [cvac_aux_th]
    family = LAGRANGE
    order = FIRST
  []
  [cvac_aux_elem_th]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[ICs]
  [phi_IC]
    type = SpecifiedSmoothCircleIC
    variable = phi
    x_positions = '150 350 250'
    y_positions = '150 150 150'
    z_positions = '150 150 322.5'
    radii = '100 100 100'
    invalue = 0
    outvalue = 1
  []
  [gr0_IC]
    type = SmoothCircleIC
    variable = gr0
    x1 = 150
    y1 = 150
    z1 = 150
    radius = 100
    invalue = 1
    outvalue = 0
  []
  [gr1_IC]
    type = SmoothCircleIC
    variable = gr1
    x1 = 350
    y1 = 150
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
    z1 = 322.5
    radius = 100
    invalue = 1
    outvalue = 0
  []
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
    # D0 = 8.33e9
    # GBmob0 = 1.4759e9
    # Q = 2.77
    # Em = 3.608
    D0 = 4.2488e11 #8.33e9
    Em = 4.23317 #3.608
    GBmob0 = 3.42828e10 # nm4/eVs #1.4759e9 # new value from Tonks/PC/Jake GG Paper
    Q = 3.01 #2.77 # new value from Tonks/PC/Jake GG Paper
    bulkindex = 1
    gbindex = -1 # -1 sets the GB D to the LANL MD Value in GPIsoMat
    surfindex = 1.16844e9 #iw*DGB  #5.8422e10
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
    surface_energy = gb_e_mat # 9.86 #19.7
    grainboundary_energy = gb_e_mat # 9.86
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
  # Cvac th is for auxvar for IC
  [cvac_th]
    type = ParsedMaterial
    property_name = cvac_th
    material_property_names = 'hs cv_eq hv'
    expression = 'cv:=(hs*cv_eq ) + hv*1.0;
                  if(cv<0.0, 0.0, cv)'
    # outputs = nemesis
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
    [GrandPotentialAlt]
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
      # Mass Conservation
      mass_conservation = false
      # concentrations = 'cvac_var cint_var'
      # hj_over_kVa = 'hv_over_kvVa hs_over_kvVa hv_over_kiVa hs_over_kiVa' #'hoverk_vu hoverk_su hoverk_vi hoverk_si'
      # hj_c_min = 'hv_cv_min hs_cv_min hv_ci_min hs_ci_min' #'cvueq_mask csueq_mask cvieq_mask csieq_mask'
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
  # THEQ
  # Vacancy
  [cvac_aux_elem_th]
    type = MaterialRealAux
    variable = 'cvac_aux_elem_th'
    property = 'cvac_th'
  []
  [cvac_aux_th]
    type = ProjectionAux
    variable = cvac_aux_th
    v = cvac_aux_elem_th
  []
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
    # outputs = csv
  []
  [c_total]
    type = ElementIntegralMaterialProperty
    mat_prop = c
    # outputs = csv
  []
  [nonlinear]
    type = NumNonlinearIterations
    # outputs = csv
  []
  [linear]
    type = NumLinearIterations
    # outputs = csv
  []
  [residuals]
    type = NumResidualEvaluations
    # outputs = csv
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

[Preconditioning]
  [SMP] #slow but good, very slow for 3D (might be another option then)
    type = SMP
    coupled_groups = 'w,phi'
  []
[]

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
  # end_time = 25000
  # steady_state_detection = true
  automatic_scaling = true
  compute_scaling_once = false
  num_steps = 1
  dt = 0.0001
  # dtmax = 200
  # dt = 0.0001
  # [TimeStepper]
  #   type = IterationAdaptiveDT
  #   optimal_iterations = 6 #WAS 6
  #   dt = 0.01
  #   growth_factor = 1.2
  #   cutback_factor = 0.8
  #   cutback_factor_at_failure = 0.5 #might be different from the curback_factor
  # []
  #[Adaptivity]
  #  refine_fraction = 0.8
  #  coarsen_fraction = 0.05 #minimize this- adds error
  #  max_h_level = 2 #test a short simulation with 1,2,3,4 for this to see where it stops helping
  #  initial_adaptivity = 2
  #[]
[]

[Outputs]
  perf_graph = false
  csv = false
  exodus = false
  nemesis = false
  # [nemesis]
  #   type = Nemesis
  #   # interval = 5              # this ExodusII will only output every third time step
  # []
  print_linear_residuals = false
  # [checkpoint]
  #   type = Checkpoint
  #   num_files = 3
  # []
[]
