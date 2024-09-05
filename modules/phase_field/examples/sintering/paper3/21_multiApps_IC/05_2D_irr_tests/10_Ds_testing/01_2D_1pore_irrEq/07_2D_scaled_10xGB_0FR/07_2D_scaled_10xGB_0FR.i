##############################################################################
# File: 07_2D_scaled_10xGB_0FR.i
# File Location: /examples/sintering/paper3/21_multiApps_IC/05_2D_irr_tests/10_Ds_testing/01_2D_1pore_irrEq/07_2D_scaled_10xGB_0FR
# Created Date: Thursday September 5th 2024
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Thursday September 5th 2024
# Modified By: Brandon Battas
# -----
# Description:
#  Testing the base iw scaled Ds=10xDgb with irr eq c adn k, but with a fission
#  rate of 0 so its basically recombination only no irradiation
#
#
##############################################################################

f_dot = 0

[Mesh]
  [gmg]
    type = DistributedRectilinearMeshGenerator
    dim = 2
    nx = 120
    ny = 80
    xmin = 0
    xmax = 30000
    ymin = 0
    ymax = 20000
  []
  uniform_refine = 0 #2
  # second_order = true
[]

[GlobalParams]
  # profile = TANH
  int_width = 2000
  op_num = 2
  var_name_base = gr
[]

[MultiApps]
  [full_NMC]
    type = FullSolveMultiApp
    execute_on = initial
    positions = '0 0 0'
    input_files = ../../00_sub/01_sub_2D_1e-8FR.i
  []
[]

[Transfers]
  [from_NMC]
    type = MultiAppCopyTransfer
    from_multi_app = full_NMC
    source_variable = 'gr0 gr1 phi cvac_aux cint_aux'
    variable = 'gr0 gr1 phi cvac_var cint_var'
  []
[]

[Variables]
  [wvac]
    initial_condition = 0
  []
  [wint]
    initial_condition = 0
  []
  [cvac_var]
  []
  [cint_var]
  []
  [gr0]
  []
  [gr1]
  []
  # [gr2]
  # []
  # [gr3]
  # []
  # [gr4]
  # []
  # [gr5]
  # []
  # [gr6]
  # []
  # [gr7]
  # []
  [phi]
  []
[]

[AuxVariables]
  [T]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 1600
  []
[]

[BCs]
  [phi_bc]
    type = NeumannBC
    variable = 'phi'
    boundary = 'left right top bottom'
    value = 0
  []
  [gr0_bc]
    type = NeumannBC
    variable = 'gr0'
    boundary = 'left right top bottom'
    value = 0
  []
  # [gr1_bc]
  #   type = NeumannBC
  #   variable = 'gr1'
  #   boundary = 'left right top bottom'
  #   value = 0
  # []
  # [gr2_bc]
  #   type = NeumannBC
  #   variable = 'gr2'
  #   boundary = 'left right top bottom'
  #   value = 0
  # []
  # [gr3_bc]
  #   type = NeumannBC
  #   variable = 'gr3'
  #   boundary = 'left right top bottom'
  #   value = 0
  # []
  # [gr4_bc]
  #   type = NeumannBC
  #   variable = 'gr4'
  #   boundary = 'left right top bottom'
  #   value = 0
  # []
  # [gr5_bc]
  #   type = NeumannBC
  #   variable = 'gr5'
  #   boundary = 'left right top bottom'
  #   value = 0
  # []
  # [gr6_bc]
  #   type = NeumannBC
  #   variable = 'gr6'
  #   boundary = 'left right top bottom'
  #   value = 0
  # []
  # [gr7_bc]
  #   type = NeumannBC
  #   variable = 'gr7'
  #   boundary = 'left right top bottom'
  #   value = 0
  # []
  [wv_bc]
    type = NeumannBC
    variable = 'wvac'
    boundary = 'left right top bottom'
    value = 0
  []
  [wi_bc]
    type = NeumannBC
    variable = 'wint'
    boundary = 'left right top bottom'
    value = 0
  []
  [cv_bc]
    type = NeumannBC
    variable = 'cvac_var'
    boundary = 'left right top bottom'
    value = 0
  []
  [ci_bc]
    type = NeumannBC
    variable = 'cint_var'
    boundary = 'left right top bottom'
    value = 0
  []
[]

[Modules]
  [PhaseField]
    [GrandPotentialAlt]
      switching_function_names = 'hv hs'
      # Chempot
      chemical_potentials = 'wvac wint'
      mobilities = 'chiuD chiiD' #cons_mob
      anisotropic = 'false false'
      susceptibilities = 'chiu chii'
      free_energies_w = 'rhovu rhosu rhovi rhosi'
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
      concentrations = 'cvac_var cint_var'
      hj_over_kVa = 'hv_over_kvVa hs_over_kvVa hv_over_kiVa hs_over_kiVa' #'hoverk_vu hoverk_su hoverk_vi hoverk_si'
      hj_c_min = 'hv_cv_min hs_cv_min hv_ci_min hs_ci_min' #'cvueq_mask csueq_mask cvieq_mask csieq_mask'
    []
  []
[]

[Kernels]
  #   # Sintering Terms for gb/surface energy
  #   [barrier_phi]
  #     type = ACBarrierFunction
  #     variable = phi
  #     v = 'gr0 gr1 gr2' #gr2 gr3'# gr4 gr5'# gr6 gr7 gr8 gr9 gr10 gr11 gr12 gr13 gr14 gr15'
  #     gamma = gamma
  #     mob_name = L_mat
  #   []
  #   [kappa_phi]
  #     type = ACKappaFunction
  #     variable = phi
  #     mob_name = L_mat
  #     kappa_name = kappa
  #   []
  # Irradiation
  # Source/Generation
  [source_vac]
    type = MaskedBodyForce
    variable = cvac_var
    mask = rho_gen #change_vac #rho_gen
    coupled_variables = 'phi wvac' # wint'
  []
  [source_int]
    type = MaskedBodyForce
    variable = cint_var
    mask = rho_gen #change_vac #rho_gen
    coupled_variables = 'phi wint' # wvac'
  []
  # Sink/Recombination
  [recombination_vac]
    type = MatReaction
    variable = cvac_var
    mob_name = rho_recomb
    args = 'wvac wint phi gr0 gr1' # gr2 gr3 gr4 gr5 gr6 gr7' # grain ops mean i should include them in combined_rho
    #coupled_variables = 'phi gr0 gr1 gr2'
  []
  [recombination_int]
    type = MatReaction
    variable = cint_var
    mob_name = rho_recomb
    args = 'wvac wint phi gr0 gr1' # gr2 gr3 gr4 gr5 gr6 gr7' # grain ops mean i should include them in combined_rho
    #coupled_variables = 'wvac phi gr0 gr1 gr2'
  []
  # Damage/Mixing
  [ballistic_mix_vac]
    type = MatDiffusion
    variable = cvac_var
    diffusivity = rho_mixing_vac
    args = 'wvac'
  []
  [ballistic_mix_int]
    type = MatDiffusion
    variable = cint_var
    diffusivity = rho_mixing_int
    args = 'wint'
  []
[]

[Materials]
  [consts]
    type = GenericConstantMaterial
    prop_names = 'Va negOverVa cvieq_mask'
    prop_values = '0.04092 -24.4379 0.0'
  []
  [k_constants]
    type = GenericConstantMaterial
    prop_names = 'ksu kvu ksi kvi' # Using the GB based values (lowest of mine)
    prop_values = '7.751e2 7.751e2 5.711e5 5.711e5' # Irradiation
    # prop_values = '6.569e2 6.569e2  5.461e7 5.461e7' # No Irradiation
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
  # [surf_e_mat] # eV/nm^2
  #   type = ParsedMaterial
  #   property_name = surf_e_mat
  #   material_property_names = gb_e_mat
  #   expression = '2 * gb_e_mat'
  #   outputs = none
  # []
  [hgb]
    type = DerivativeParsedMaterial
    property_name = hgb
    derivative_order = 2
    coupled_variables = 'gr0 gr1'
    # expression = '16 * ( (gr0 * gr1)^2 )'
    expression = 'hg:=16 * ( (gr0 * gr1)^2 );
                  if(hg>1e-4,hg,0.0)' #+ + (gr0 * gr2)^2 + (gr1 * gr2)^2
    # expression = '4*(1 - (gr0^2 + gr1^2 + etad0^2))^2'
    # outputs = exodus
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
  [chiuD]
    type = GrandPotentialIsoMaterial
    f_name = chiuD #chiuD
    solid_mobility = L #CHANGED FROM L
    void_mobility = Lv
    chi = chiu
    c = phi
    T = T
    D0 = 4.2488e11 #8.33e9
    Em = 4.23317 #3.608
    GBmob0 = 3.42828e10 # nm4/eVs #1.4759e9 # new value from Tonks/PC/Jake GG Paper
    Q = 3.01 #2.77 # new value from Tonks/PC/Jake GG Paper
    vaporindex = 1
    bulkindex = 1
    gbindex = -1 # -1 sets the GB D to the LANL MD Value in GPIsoMat
    surfindex = -1 #1e11
    GBwidth = 3.5 # based on avg of two lanl values
    surf_thickness = 3.5 # keeping equal to gb for simplicity
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
    vaporindex = 1
    bulkindex = 1
    gbindex = -1 #10 # -1 sets the GB D to the LANL MD Value in GPIsoMat
    surfindex = -1 #100 #1e11
    GBwidth = 3.5 # based on avg of two lanl values
    surf_thickness = 3.5 # keeping equal to gb for simplicity
    iw_scaling = true
    D_out_name = int_diffus
  []
  [densificaiton]
    type = GrandPotentialMultiSinteringMaterial
    chemical_potential_vac = wvac
    chemical_potential_int = wint
    void_op = phi
    Temperature = T
    surface_energy = gb_e_mat #gb_e_mat #surf_e_mat #19.7
    grainboundary_energy = gb_e_mat #9.86
    vac_solid_energy_coefficient = ksu
    int_solid_energy_coefficient = ksi
    vac_void_energy_coefficient = kvu
    int_void_energy_coefficient = kvi
    equilibrium_vacancy_concentration = cv_eq
    equilibrium_interstitial_concentration = ci_eq
    solid_energy_model = PARABOLIC
    mass_conservation = true
  []
  [cv_eq]
    type = DerivativeParsedMaterial
    property_name = cv_eq
    derivative_order = 2
    coupled_variables = 'gr0 gr1 phi wvac'
    material_property_names = 'hgb(phi,gr0,gr1)' #'rhovi(vac) rhosi hv(phi)'
    constant_names = 'cb cgb'
    constant_expressions = '3.877e-04 4.347e-03' #Irradiation
    # constant_expressions = '2.424e-06 5.130e-03' #No Irradiation- LANL
    expression = 'cgb * hgb + (1 - hgb)*cb'
    outputs = none # + phi^2
  []
  [ci_eq]
    type = DerivativeParsedMaterial
    property_name = ci_eq
    derivative_order = 2
    coupled_variables = 'gr0 gr1 phi wint'
    material_property_names = 'hgb(phi,gr0,gr1)' # 'rhovi(wint) rhosi(wint) hv(phi)'
    constant_names = 'cb cgb'
    constant_expressions = '7.258e-09 5.900e-06' #Irradiation
    # constant_expressions = '1.667e-32 6.170e-08' #'1.667e-32 6.170e-08' #No Irradiation- LANL
    # constant_expressions = '2.424e-06 5.130e-03' #No Irradiation VACANCY- LANL
    expression = 'cgb * hgb + (1 - hgb)*cb'
    outputs = none #+ phi^2
  []
  # Outputs for visualization
  [cvac]
    type = ParsedMaterial
    property_name = cvac
    material_property_names = 'hs rhosu hv rhovu Va'
    expression = 'Va*(hs*rhosu + hv*rhovu)'
    outputs = nemesis
  []
  [cint]
    type = ParsedMaterial
    property_name = cint
    material_property_names = 'hs rhosi hv rhovi Va'
    expression = 'Va*(hs*rhosi + hv*rhovi)'
    outputs = nemesis
  []
  [hgb_out]
    type = ParsedMaterial
    property_name = hgb_out
    material_property_names = 'hs hgb hv'
    expression = 'hgb'
    outputs = nemesis
  []
  [net_defect]
    type = ParsedMaterial
    property_name = net_defect
    material_property_names = 'rho_gen rho_recomb'
    expression = 'rho_gen + rho_recomb'
    outputs = nemesis
  []
  # IRRADIATION
  [combined_rho_vac]
    type = DerivativeParsedMaterial
    property_name = combined_rho_vac
    coupled_variables = 'wvac phi'
    derivative_order = 2
    material_property_names = 'rhovu(wvac) rhosu(wvac) hv(phi) hs(phi)'
    expression = 'hv*rhovu + hs*rhosu' #'(1-hv)*rhos' #
    outputs = none #'nemesis'
  []
  [combined_rho_int]
    type = DerivativeParsedMaterial
    property_name = combined_rho_int
    coupled_variables = 'wint phi'
    derivative_order = 2
    material_property_names = 'rhovi(wint) rhosi(wint) hv(phi) hs(phi)'
    expression = 'hv*rhovi + hs*rhosi' #'(1-hv)*rhos' #
    outputs = none #'nemesis'
  []
  [a_r]
    type = ParsedMaterial
    property_name = a_r
    constant_names = 'Va Z a_0 kB Di_0 Ei_B' # Di_0 Ei_B'
    constant_expressions = '0.04092 250 0.25 8.617343e-5 4.0767e11 4.08453089' #1e13 2
    material_property_names = 'hs(phi)'
    coupled_variables = 'phi T gr0 gr1'
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
    material_property_names = 'hs(phi) Va hv(phi) htj'
    expression = 'rg:=f_dot * noise * Nc * Nd * hs * Va;
                  if(hv<=1e-6,rg,0)'
    # expression = 'f_dot * noise * Nc * Nd * hs * Va * (1-htj)'
    outputs = none #emesis #'nemesis'
  []
  [rho_recomb] #This one is off on GB?
    type = DerivativeParsedMaterial
    property_name = rho_recomb
    coupled_variables = 'wvac wint phi' # gr0 gr1' #gr0 gr1 gr2
    derivative_order = 2
    # additional_derivative_symbols = w # combined_rho_vac combined_rho_int
    material_property_names = 'Va a_r(phi) combined_rho_vac(wvac,phi) combined_rho_int(wint,phi) hv'
    expression = 'out:=a_r * combined_rho_vac * combined_rho_int * Va;
                  if((out>0.0 & hv<=1e-6),0.0-out,0.0)'
    outputs = none #nemesis #'nemesis'
  []
  [rho_mixing_vac]
    type = DerivativeParsedMaterial
    property_name = rho_mixing_vac
    coupled_variables = 'wvac'
    derivative_order = 2
    constant_names = 'Nc Vc noise tc Dc f_dot'
    constant_expressions = '2 268 1 1e-11 1e12 ${f_dot}'
    material_property_names = 'chiu(phi,wvac) Va'
    expression = 'f_dot * noise * Nc * tc * Vc * Dc * chiu * Va' # * hs
    outputs = none #nemesis #'nemesis'
  []
  [rho_mixing_int]
    type = DerivativeParsedMaterial
    property_name = rho_mixing_int
    coupled_variables = 'wint'
    derivative_order = 2
    constant_names = 'Nc Vc noise tc Dc f_dot'
    constant_expressions = '2 268 1 1e-11 1e12 ${f_dot}'
    material_property_names = 'chii(phi,wint) Va'
    expression = 'f_dot * noise * Nc * tc * Vc * Dc * chii * Va' # * hs
    outputs = none #nemesis #'nemesis'
  []
[]

[Postprocessors]
  [max_mpi_memory]
    type = MemoryUsage
    value_type = max_process
    report_peak_value = True
    mem_units = gigabytes
    execute_on = 'NONLINEAR TIMESTEP_END'
  []
  [cv_var_total]
    type = ElementIntegralVariablePostprocessor
    variable = cvac_var
  []
  [wvac_total]
    type = ElementIntegralVariablePostprocessor
    variable = wvac
  []
  [ci_var_total]
    type = ElementIntegralVariablePostprocessor
    variable = cint_var
  []
  [wint_total]
    type = ElementIntegralVariablePostprocessor
    variable = wint
  []
  [phi_total]
    type = ElementIntegralVariablePostprocessor
    variable = phi
  []
  [gr0_total]
    type = ElementIntegralVariablePostprocessor
    variable = gr0
  []
  [gr1_total]
    type = ElementIntegralVariablePostprocessor
    variable = gr1
  []
  [total_defect_change]
    type = ElementIntegralMaterialProperty
    mat_prop = net_defect
  []
  # [gr2_total]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = gr2
  # []
  # [gr3_total]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = gr3
  # []
  # [gr4_total]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = gr4
  # []
  # [gr5_total]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = gr5
  # []
  # [gr6_total]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = gr6
  # []
  # [gr7_total]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = gr7
  # []
  [hgb_total]
    type = ElementIntegralMaterialProperty
    mat_prop = hgb
  []
  [cvac_total]
    type = ElementIntegralMaterialProperty
    mat_prop = cvac
  []
  [cint_total]
    type = ElementIntegralMaterialProperty
    mat_prop = cint
  []
  [runtime]
    type = PerfGraphData
    section_name = "Root"
    data_type = TOTAL
  []
  [void_tracker]
    type = FeatureFloodCount
    variable = phi
    threshold = 0.2 # 0.1 #0.2
    connecting_threshold = 0.08 #0.09 #0.08
    compute_var_to_feature_map = true
    execute_on = 'initial timestep_end'
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
[]

[UserObjects]
  [terminator]
    type = Terminator
    expression = 'void_tracker < 2'
    execute_on = TIMESTEP_END
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -sub_pc_factor_shift_type'
  petsc_options_value = ' asm      lu           2                nonzero'
  nl_max_its = 30
  l_max_its = 60
  l_tol = 1e-06 #4
  nl_rel_tol = 1e-6 #6 #default is 1e-8
  # nl_abs_tol = 1e-14 #only needed when near equilibrium or veeeery small dt
  # start_time = 0
  end_time = 1e7 #1e8
  # num_steps = 5
  # steady_state_detection = true
  # # From tonks ode input
  automatic_scaling = true
  compute_scaling_once = false
  # line_search = none
  # dt = 1.0
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 6
    dt = 0.1
  []
[]

[Outputs]
  csv = true
  exodus = false
  nemesis = true
  checkpoint = false
  # file_base = MC_GifRif_highDs
[]

# [Debug]
#   show_var_residual_norms = true
# []
