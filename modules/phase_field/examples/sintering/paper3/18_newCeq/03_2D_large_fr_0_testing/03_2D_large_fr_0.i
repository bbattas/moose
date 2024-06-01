##############################################################################
# File: 01_2D_1pore_fr_0.i
# File Location: /examples/sintering/paper3/18_newCeq/01_2D_1pore_fr_0
# Created Date: Thursday May 30th 2024
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Friday May 31st 2024
# Modified By: Brandon Battas
# -----
# Description:
#  Testing stuff with the no irradiation for new large 2D IC (20 pores 6%vol)
#
#  Using the 5% GB volume version of parabolic coeffs for 1600K
#  1,735,929 DoFs
##############################################################################

f_dot = 0
# Thermal 5%GB Volume 1600K
ks_vac = 1.302e4
ks_int = 1.092e9
# # Irradiation 5%GB Volume 1600K
# ks_vac = 5.753e3
# ks_int = 1.116e7

[Mesh]
  [ebsd_mesh]
    type = EBSDMeshGenerator
    filename = ../00_d3d_txt/2D_100x100um_8umavg_20pore.txt
  []
  [subdomain_external]
    type = ParsedSubdomainMeshGenerator
    input = ebsd_mesh
    combinatorial_geometry = 'x > 99000'
    block_id = 1
  []
  parallel_type = DISTRIBUTED
  uniform_refine = 0
[]

[GlobalParams]
  op_num = 6
  var_name_base = gr
  int_width = 2000 #1000 #min radius is like 2250, element size of 250
  # profile = TANH # not used at the moment? only in circleic?
[]

[Variables]
  [wvac]
  []
  [wint]
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
  [VoidIC] #External
    type = BoundingBoxIC
    variable = phi
    block = 1
    inside = 1
    outside = 0.01
    x1 = 100000
    x2 = 150000
    y1 = -5000
    y2 = 125000
  []
  [voidIC2] #Internal
    type = SpecifiedSmoothCircleIC
    variable = phi
    invalue = 1
    outvalue = 0.01
    radii = '3387.41 3204.63 3097.33 3157.62 2608.26 2573.66 3162.96 2580.68 2922.15 3499.
    3438.79 3442.07 3025.42 2697.36 3107.91 3175.05 3355.8  2983.18 2935.   3191.81'
    x_positions = '65392.28 87521.94 74200.75 33056.45 38578.91 73717.72 79407.04 61057.73 65502.19 34105.4
    15495.84 39069.52 47270.14 81182.91 19481.53 53804.82 16891.68 80337.59 52473.89 51584.05'
    y_positions = '57040.12 66912.04 21248.38 73556.73 87021.67 91228.2  35742.18 89239.61 9333.9
    11173.35 81164.89 38512.26 52276.59 55357.99 23911.9  24603.55 55058.71 80676.15 65521.33 78256.87'
    z_positions = '0.       0.       0.       0.       0.       0.       0.       0. 0.       0.
    0.       0.       0.       0.       0.       0. 0.       0.       0.       0.'
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
[]

[Functions]
  [f_T]
    type = ConstantFunction
    value = 1600
  []
[]

[Materials]
  # Free energy coefficients for parabolic curves
  [k_constants]
    type = GenericConstantMaterial
    prop_names = 'ksu kvu ksi kvi'
    prop_values = '${ks_vac} ${ks_vac} ${ks_int} ${ks_int}' #'26.8 26.8 3.6e21 3.6e21'#154.16 154.16
  []
  # New GB and Surface energy as a function of T
  [gb_e_mat] # eV/nm^2
    type = ParsedMaterial
    property_name = gb_e_mat
    coupled_variables = 'T'
    constant_names = 'a b c'
    constant_expressions = '-5.87e-4 1.56 6.24151'
    expression = 'c*(a*T + b)'
    outputs = none
  []
  # [surf_e_mat]
  #   type = ParsedMaterial
  #   property_name = surf_e_mat
  #   material_property_names = 'gb_e_mat'
  #   expression = '2 * gb_e_mat' #2*
  #   outputs = nemesis
  # []
  # Diffusivity and mobilities
  [chiuD]
    type = GrandPotentialIsoMaterial
    f_name = chiuD
    solid_mobility = L #CHANGED FROM L
    void_mobility = Lv
    chi = chiu
    c = phi
    T = T
    D0 = 4.2488e11 #8.33e9
    Em = 4.23317 #3.608
    GBmob0 = 3.42828e10 # nm4/eVs #1.4759e9 # new value from Tonks/PC/Jake GG Paper
    Q = 3.01 #2.77 # new value from Tonks/PC/Jake GG Paper
    bulkindex = 1
    gbindex = -1 # -1 sets the GB D to the LANL MD Value in GPIsoMat
    surfindex = -1 #1e11
    GBwidth = 3.5 #LANL Creep DFT JNM
    surf_thickness = 3.5 #LANL Creep DFT JNM
    iw_scaling = true
    D_out_name = vac_diffus
  []
  [chiiD]
    type = GrandPotentialIsoIntMaterial
    f_name = chiiD
    chi = chii
    c = phi
    T = T
    D0 = 4.0767e11 #1e13 #Ian's irradiation paper (Matzke 1987)
    Em = 4.08453089 #2 #Ian's irradiation paper (Matzke 1987)
    bulkindex = 1
    gbindex = -1 #10 # -1 sets the GB D to the LANL MD Value in GPIsoMat
    surfindex = -1 #100 #1e11
    GBwidth = 3.5 #LANL Creep DFT JNM
    surf_thickness = 3.5 #LANL Creep DFT JNM
    iw_scaling = true
    D_out_name = int_diffus
  []
  [cv_eq]
    type = DerivativeParsedMaterial
    property_name = cv_eq
    coupled_variables = 'gr0 gr1 gr2 phi T wvac'
    material_property_names = 'hgb(phi,gr0,gr1,gr2,gr3,gr4,gr5)' #'rhovi(vac) rhosi hv(phi)'
    constant_names = 'cb cgb'
    # constant_expressions = '3.877e-04 4.347e-03' #Irradiation
    constant_expressions = '2.424e-06 5.130e-03' #No Irradiation- LANL
    expression = 'cgb * hgb + (1 - hgb)*cb'
    outputs = none #'nemesis' # + phi^2
  []
  [ci_eq]
    type = DerivativeParsedMaterial
    property_name = ci_eq
    coupled_variables = 'gr0 gr1 gr2 phi T wint'
    material_property_names = 'hgb(phi,gr0,gr1,gr2,gr3,gr4,gr5)' # 'rhovi(wint) rhosi(wint) hv(phi)'
    constant_names = 'cb cgb'
    # constant_expressions = '7.258e-09 5.900e-06' #Irradiation
    constant_expressions = '1.667e-32 6.170e-08' #No Irradiation- LANL
    expression = 'cgb * hgb + (1 - hgb)*cb'
    outputs = none #'nemesis' #+ phi^2
  []
  [densificaiton]
    type = GrandPotentialMultiSinteringMaterial
    chemical_potential_vac = wvac
    chemical_potential_int = wint
    void_op = phi
    Temperature = T
    surface_energy = gb_e_mat #surf_e_mat #19.7
    grainboundary_energy = gb_e_mat #9.86
    vac_solid_energy_coefficient = ksu
    int_solid_energy_coefficient = ksi
    vac_void_energy_coefficient = kvu
    int_void_energy_coefficient = kvi
    equilibrium_vacancy_concentration = cv_eq
    equilibrium_interstitial_concentration = ci_eq
    solid_energy_model = PARABOLIC
  []
  # Concentration is only meant for output
  [cvac]
    type = ParsedMaterial
    property_name = cvac
    material_property_names = 'hs rhosu hv rhovu'
    constant_names = 'Va'
    constant_expressions = '0.04092'
    expression = 'Va*(hs*rhosu + hv*rhovu)'
    outputs = nemesis #none #exodus
  []
  [cint]
    type = ParsedMaterial
    property_name = cint
    material_property_names = 'hs rhosi hv rhovi'
    constant_names = 'Va'
    constant_expressions = '0.04092'
    expression = 'Va*(hs*rhosi + hv*rhovi)'
    outputs = nemesis #none #exodus
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
    outputs = none #'nemesis'
  []
  [combined_rho_vac]
    type = DerivativeParsedMaterial
    property_name = combined_rho_vac
    coupled_variables = 'wvac phi'
    derivative_order = 2
    material_property_names = 'rhovu(wvac) rhosu(wvac) hv(phi)'
    expression = 'hv*rhovu + (1-hv)*rhosu' #'(1-hv)*rhos' #
    outputs = none #'nemesis'
  []
  [combined_rho_int]
    type = DerivativeParsedMaterial
    property_name = combined_rho_int
    coupled_variables = 'wint phi'
    derivative_order = 2
    material_property_names = 'rhovi(wint) rhosi(wint) hv(phi)'
    expression = 'hv*rhovi + (1-hv)*rhosi' #'(1-hv)*rhos' #
    outputs = none #'nemesis'
  []
  # Irradiation materials
  [a_r]
    type = ParsedMaterial
    property_name = a_r
    constant_names = 'Va Z a_0 kB Di_0 Ei_B' # Di_0 Ei_B'
    constant_expressions = '0.04092 250 0.25 8.617343e-5 4.0767e11 4.08453089' #1e13 2
    material_property_names = 'hs(phi)'
    coupled_variables = 'phi T gr0 gr1 gr2'
    expression = 'dint:=Di_0 * exp(-Ei_B / (kB * T));
                  hs * Va * Z * dint / (a_0^2)'
    outputs = none #nemesis #nemesis #'nemesis'
  []
  [rho_gen]
    type = DerivativeParsedMaterial
    property_name = rho_gen
    coupled_variables = 'phi'
    derivative_order = 1
    constant_names = 'Nc Nd noise f_dot'
    constant_expressions = '2 5 1 ${f_dot}'
    material_property_names = 'hs(phi)'
    expression = 'f_dot * noise * Nc * Nd * hs'
    outputs = none #'nemesis'
  []
  [rho_recomb] #This one is off on GB?
    type = DerivativeParsedMaterial
    property_name = rho_recomb
    coupled_variables = 'wvac wint phi' # gr0 gr1' #gr0 gr1 gr2
    derivative_order = 2
    # additional_derivative_symbols = w # combined_rho_vac combined_rho_int
    material_property_names = 'a_r(phi) combined_rho_vac(wvac,phi) combined_rho_int(wint,phi)'
    expression = 'out:=a_r * combined_rho_vac * combined_rho_int;
                  if(out>0.0,0.0-out,0.0)'
    outputs = none #'nemesis'
  []
  [rho_mixing_vac]
    type = DerivativeParsedMaterial
    property_name = rho_mixing_vac
    coupled_variables = 'wvac'
    derivative_order = 2
    constant_names = 'Nc Vc noise tc Dc f_dot'
    constant_expressions = '2 268 1 1e-11 1e12 ${f_dot}'
    material_property_names = 'chiu(phi,wvac)'
    expression = 'f_dot * noise * Nc * tc * Vc * Dc * chiu' # * hs
    outputs = none #nemesis #'nemesis'
  []
  [rho_mixing_int]
    type = DerivativeParsedMaterial
    property_name = rho_mixing_int
    coupled_variables = 'wint'
    derivative_order = 2
    constant_names = 'Nc Vc noise tc Dc f_dot'
    constant_expressions = '2 268 1 1e-11 1e12 ${f_dot}'
    material_property_names = 'chii(phi,wint)'
    expression = 'f_dot * noise * Nc * tc * Vc * Dc * chii' # * hs
    outputs = none #nemesis #'nemesis'
  []
  [hgb]
    type = DerivativeParsedMaterial
    property_name = hgb
    coupled_variables = 'phi gr0 gr1 gr2 gr3 gr4 gr5'
    expression = '16 * ( (gr0 * gr1)^2 + (gr0 * gr2)^2 + (gr0 * gr3)^2 + (gr0 * gr4)^2 + (gr0 * gr5)^2 +
    (gr1 * gr2)^2 + (gr1 * gr3)^2 + (gr1 * gr4)^2 + (gr1 * gr5)^2 +
    (gr2 * gr3)^2 + (gr2 * gr4)^2 + (gr2 * gr5)^2 +
    (gr3 * gr4)^2 + (gr3 * gr5)^2 + (gr4 * gr5)^2 )'
    outputs = none
  []
[]

[Modules]
  [PhaseField]
    [GrandPotential]
      switching_function_names = 'hv hs'
      anisotropic = 'false false'

      chemical_potentials = 'wvac wint'
      mobilities = 'chiuD chiiD'
      susceptibilities = 'chiu chii'
      free_energies_w = 'rhovu rhosu rhovi rhosi'

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
  # Dont need these two if surface and gb energy are equal
  # [barrier_phi]
  #   type = ACBarrierFunction
  #   variable = phi
  #   v = 'gr0 gr1 gr2' # gr3 gr4 gr5 gr6 gr7 gr8 gr9' # gr10 gr11 gr12 gr13 gr14 gr15'
  #   gamma = gamma
  #   mob_name = L_mat
  # []
  # [kappa_phi]
  #   type = ACKappaFunction
  #   variable = phi
  #   mob_name = L_mat
  #   kappa_name = kappa
  # []
  # Irradiation
  # Source/Generation
  [source_vac]
    type = MaskedBodyForce
    variable = wvac
    mask = rho_gen #change_vac #rho_gen
    coupled_variables = 'phi' # wint'
  []
  [source_int]
    type = MaskedBodyForce
    variable = wint
    mask = rho_gen #change_vac #rho_gen
    coupled_variables = 'phi' # wvac'
  []
  # Sink/Recombination
  [recombination_vac]
    type = MatReaction
    variable = wvac
    mob_name = rho_recomb
    args = 'wint phi' # gr0 gr1 gr2'
    #coupled_variables = 'phi gr0 gr1 gr2'
  []
  [recombination_int]
    type = MatReaction
    variable = wint
    mob_name = rho_recomb
    args = 'wvac phi' # gr0 gr1 gr2'
    #coupled_variables = 'wvac phi gr0 gr1 gr2'
  []
  # Damage/Mixing
  [ballistic_mix_vac]
    type = MatDiffusion
    variable = wvac
    diffusivity = rho_mixing_vac
  []
  [ballistic_mix_int]
    type = MatDiffusion
    variable = wint
    diffusivity = rho_mixing_int
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
  # [var_indices]
  #   type = FeatureFloodCountAux
  #   variable = var_indices
  #   flood_counter = grain_tracker
  #   field_display = VARIABLE_COLORING
  #   execute_on = 'initial timestep_end'
  # []
  # [grain_aux]
  #   type = EBSDReaderPointDataAux
  #   variable = ebsd_grains
  #   ebsd_reader = ebsd_reader
  #   data_name = 'feature_id'
  #   execute_on = 'initial timestep_end'
  # []
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
  [c_total_vac]
    type = ElementIntegralMaterialProperty
    mat_prop = cvac
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
  # FFC
  [void_tracker]
    type = FeatureFloodCount
    variable = phi
    threshold = 0.5 # 0.1 #0.2
    connecting_threshold = 0.5 #0.09 #0.08
    compute_var_to_feature_map = true
    execute_on = 'initial timestep_end'
  []
  [grain_tracker_fc]
    type = FeatureFloodCount
    variable = 'gr0 gr1 gr2 gr3 gr4 gr5' # gr3 gr4 gr5 gr6 gr7 gr8 gr9' #'gr0 gr1 gr2'
    threshold = 0.5 #0.2
    connecting_threshold = 0.5 #0.08
    compute_var_to_feature_map = true
    execute_on = 'initial timestep_end'
  []
  # # Irradiation PPs Used in Materials
  # [hs_average]
  #   type = ElementAverageMaterialProperty
  #   mat_prop = hs
  #   outputs = csv
  # []
  # [average_rho_vac]
  #   type = ElementAverageMaterialProperty
  #   mat_prop = combined_rho_vac
  #   outputs = csv
  # []
  # Other irradiation based PPs
  [total_phi]
    type = ElementIntegralVariablePostprocessor
    variable = phi
    outputs = csv
  []
  [avg_rhov]
    type = ElementAverageMaterialProperty
    mat_prop = combined_rho_vac
    outputs = csv
  []
  [avg_rhoi]
    type = ElementAverageMaterialProperty
    mat_prop = combined_rho_int
    outputs = csv
  []
  [total_rhov]
    type = ElementIntegralMaterialProperty
    mat_prop = combined_rho_vac
    outputs = csv
  []
  [total_rhoi]
    type = ElementIntegralMaterialProperty
    mat_prop = combined_rho_int
    outputs = csv
  []
[]

[VectorPostprocessors]
  [voids]
    type = FeatureVolumeVectorPostprocessor
    flood_counter = void_tracker
    execute_on = 'initial timestep_end final'
    output_centroids = true #false #was true
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

# [Controls]
#   # [extraGen]
#   #   type = TimePeriod
#   #   disable_objects = 'Kernels::source_vac Kernels::source_int
#   #   Kernels::recombination_vac Kernels::recombination_int'
#   #   # Kernels::ballistic_mix_vac Kernels::ballistic_mix_int'
#   #   enable_objects = 'Kernels::extra_source_int Kernels::extra_source_vac'
#   #   start_time = 0
#   #   end_time = 2000
#   #   execute_on = 'INITIAL TIMESTEP_END'
#   # []
#   # # [irr_kernels]
#   # #   type = TimePeriod
#   # #   # disable_objects = 'Kernels::source_vac Kernels::source_int
#   # #   # Kernels::recombination_vac Kernels::recombination_int
#   # #   # Kernels::ballistic_mix_vac Kernels::ballistic_mix_int'
#   # #   enable_objects = 'Kernels::source_vac Kernels::source_int'
#   # #   disable_objects = 'Kernels::recombination_vac Kernels::recombination_int'
#   # #   start_time = 0
#   # #   end_time = 1e5 #1e6
#   # #   execute_on = 'INITIAL TIMESTEP_END'
#   # # []
#   # [irr_newdt]
#   #   type = TimePeriod
#   #   enable_objects = 'TimeStepper::tstepper1'
#   #   disable_objects = 'TimeStepper::tstepper2'
#   #   start_time = 0
#   #   end_time = 1e5
#   #   execute_on = 'INITIAL TIMESTEP_END'
#   # []
# []

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
  [terminator_2void]
    type = Terminator
    expression = 'void_tracker < 2'
  []
[]

[Preconditioning]
  [SMP] #slow but good, very slow for 3D (might be another option then)
    type = SMP
    full = true
    # coupled_groups = 'wvac,wint,phi'
  []
[]

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
  end_time = 1e6 #1e10 #5e6 #0.006
  # num_steps = 20 #00
  # steady_state_detection = true
  # # From tonks ode input
  automatic_scaling = true
  compute_scaling_once = false
  line_search = none
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 6
    dt = 10 #2.5
    # linear_iteration_ratio = 1e5 #needed with large linear number for asmilu
    # growth_factor = 1.8
    # cutback_factor = 0.5
    # cutback_factor_at_failure = 0.5 #might be different from the curback_factor
  []
  # [TimeSteppers]
  #   [tstepper1]
  #     type = IterationAdaptiveDT
  #     optimal_iterations = 6
  #     dt = 10 #2.5
  #     # linear_iteration_ratio = 1e5 #needed with large linear number for asmilu
  #   []
  #   [tstepper2]
  #     type = IterationAdaptiveDT
  #     optimal_iterations = 6
  #     reset_dt = true
  #     dt = 10 #2.5
  #     # linear_iteration_ratio = 1e5 #needed with large linear number for asmilu
  #   []
  # []
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
  checkpoint = false
  # file_base = fr_${f_dot}_kIrrGB_1000R/fr_${f_dot}_kIrrGB_1000R
  # nemesis = false
  # fr_1.00e-10_csv/fr_1.00e-10
  # [csv]
  #   type = CSV
  #   # file_base = 02_2D_8pore_config1_csv/02_2D_8pore_config1
  #   time_step_interval = 3
  # []
  [nemesis]
    type = Nemesis
    # file_base = 02_2D_8pore_config1_nemesis/02_2D_8pore_config1_nemesis
    # interval = 3 # this ExodusII will only output every third time step
    # time_step_interval = 3
  []
  print_linear_residuals = false
  # [checkpoint]
  #   type = Checkpoint
  #   # file_base = 02_2D_8pore_config1_checkpoint
  #   num_files = 5
  # []
[]

# [Debug]
#   show_var_residual_norms = true
# []
