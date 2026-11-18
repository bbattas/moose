##############################################################################
# File: 01_sub_3D_25_5pore.i
# File Location: /examples/sintering/paper3/21_multiApps_IC/12_3D_25/00_sub
# Created Date: Wednesday November 18th 2026
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Wednesday November 18th 2026
# Modified By: Brandon Battas
# -----
# Description:
#  Subfile for new larger 3D tests
#  77 grains, 5% porosity by volume with 12 pores
#
#
##############################################################################

[Mesh]
  [ebsd_mesh]
    type = EBSDMeshGenerator
    filename = 3D_25x25x25_8avg_1_plusVoid.txt
  []
  [subdomain_external]
    type = ParsedSubdomainMeshGenerator
    input = ebsd_mesh
    combinatorial_geometry = 'x > 25000'
    block_id = 1
  []
  parallel_type = DISTRIBUTED
  uniform_refine = 0
[]

[GlobalParams]
  # profile = TANH
  int_width = 2000
  op_num = 20
  var_name_base = gr
[]

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
    # use_kdtree = true
  []
  [grain_tracker]
    type = GrainTracker
    threshold = 0.1 #0.2
    connecting_threshold = 0.09 #0.08
    # compute_halo_maps = false #true#false
    verbosity_level = 1
  []
  #   [terminator]
  #     type = Terminator
  #     expression = 'void_tracker < 2'
  #     execute_on = TIMESTEP_END
  #   []
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
    outside = 0.0 #0.01
    x1 = 25000
    x2 = 40000
    y1 = -5000
    y2 = 35000
  []
  [voidIC2] #Internal
    type = SpecifiedSmoothCircleIC
    variable = phi
    invalue = 1
    outvalue = 0 #0.01
    radii = '2738.0 2279.0 2589.0 2741.0 2192.0 2430.0 2349.0 2878.0 2365.0 2424.0 2187.0 2550.0'
    x_positions = '13087.0 12124.0 17745.0 9255.0 18815.0 19631.0 4934.0 13681.0 10555.0 7674.0 4619.0 5099.0'
    y_positions = '10255.0 5401.0 6644.0 17961.0 13025.0 17038.0 18670.0 18100.0 9188.0 12930.0 12850.0 7544.0'
    z_positions = '13419.0 20208.0 8490.0 9978.0 19866.0 11035.0 15942.0 16096.0 6392.0 19127.0 4645.0 10489.0'
    block = 0
  []
[]

[Variables]
  [wvac]
    initial_condition = 0
  []
  [wint]
    initial_condition = 0
  []
  [PolycrystalVariables]
  []
  [phi]
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
  # Interstitial mat to variable
  [cint_aux_th]
    family = LAGRANGE
    order = FIRST
  []
  [cint_aux_elem_th]
    family = MONOMIAL
    order = CONSTANT
  []
  # IRR EQ
  # Vacancy mat to variable
  [cvac_aux_ir]
    family = LAGRANGE
    order = FIRST
  []
  [cvac_aux_elem_ir]
    family = MONOMIAL
    order = CONSTANT
  []
  # Interstitial mat to variable
  [cint_aux_ir]
    family = LAGRANGE
    order = FIRST
  []
  [cint_aux_elem_ir]
    family = MONOMIAL
    order = CONSTANT
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
  # Interstitial
  [cint_aux_elem_th]
    type = MaterialRealAux
    variable = 'cint_aux_elem_th'
    property = 'cint_th'
  []
  [cint_aux_th]
    type = ProjectionAux
    variable = cint_aux_th
    v = cint_aux_elem_th
  []
  # THEQ
  # Vacancy
  [cvac_aux_elem_ir]
    type = MaterialRealAux
    variable = 'cvac_aux_elem_ir'
    property = 'cvac_ir'
  []
  [cvac_aux_ir]
    type = ProjectionAux
    variable = cvac_aux_ir
    v = cvac_aux_elem_ir
  []
  # Interstitial
  [cint_aux_elem_ir]
    type = MaterialRealAux
    variable = 'cint_aux_elem_ir'
    property = 'cint_ir'
  []
  [cint_aux_ir]
    type = ProjectionAux
    variable = cint_aux_ir
    v = cint_aux_elem_ir
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
      mass_conservation = false
      # concentrations = 'cvac_var cint_var'
      # hj_over_kVa = 'hv_over_kvVa hs_over_kvVa hv_over_kiVa hs_over_kiVa' #'hoverk_vu hoverk_su hoverk_vi hoverk_si'
      # hj_c_min = 'hv_cv_min hs_cv_min hv_ci_min hs_ci_min' #'cvueq_mask csueq_mask cvieq_mask csieq_mask'
    []
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
    # prop_values = '7.751e2 7.751e2 5.711e5 5.711e5' # Irradiation
    prop_values = '6.569e2 6.569e2  5.461e7 5.461e7' # No Irradiation
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
    type = SwitchingFunctionGBMaterial
    h_name = hgb
    grain_ops = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'
    hgb_threshold = 0.0001
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
    surfindex = 5.8422e10 #-1 #1e11
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
    surfindex = 1.3989e11 #-1 #100 #1e11
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
    mass_conservation = false
  []
  [cv_eq]
    type = UO2CeqMaterial
    ceq_name = cv_eq
    hgb = hgb
    vi = VAC_TH
  []
  [ci_eq]
    type = UO2CeqMaterial
    ceq_name = ci_eq
    hgb = hgb
    vi = INT_TH
  []
  [cv_eq_ir]
    type = UO2CeqMaterial
    ceq_name = cv_eq_ir
    hgb = hgb
    vi = VAC_IRR
  []
  [ci_eq_ir]
    type = UO2CeqMaterial
    ceq_name = ci_eq_ir
    hgb = hgb
    vi = INT_IRR
  []
  # Outputs for visualization
  [cvac_th]
    type = ParsedMaterial
    property_name = cvac_th
    material_property_names = 'hs cv_eq hv Va'
    expression = 'cv:=(hs*cv_eq ) + hv*1.0;
                  if(cv<0.0, 0.0, cv)'
    # outputs = nemesis
  []
  [cint_th]
    type = ParsedMaterial
    property_name = cint_th
    material_property_names = 'hs ci_eq hv Va'
    expression = 'ci:=(hs*ci_eq) + hv*0.0;
                  if(ci<0.0, 0.0, ci)'
    # outputs = nemesis
  []
  [cvac_ir]
    type = ParsedMaterial
    property_name = cvac_ir
    material_property_names = 'hs cv_eq_ir hv Va'
    expression = 'cv:=(hs*cv_eq_ir ) + hv*1.0;
                  if(cv<0.0, 0.0, cv)'
    # outputs = nemesis
  []
  [cint_ir]
    type = ParsedMaterial
    property_name = cint_ir
    material_property_names = 'hs ci_eq_ir hv Va'
    expression = 'ci:=(hs*ci_eq_ir) + hv*0.0;
                  if(ci<0.0, 0.0, ci)'
    # outputs = nemesis
  []
  # # hgb check for interface width
  # [hgb_out]
  #   type = ParsedMaterial
  #   property_name = hgb_out
  #   material_property_names = 'hgb'
  #   expression = 'hgb'
  #   outputs = nemesis
  # []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = false
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = PJFNK
  # petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -sub_pc_factor_shift_type'
  # petsc_options_value = ' asm      lu           2                nonzero'
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -mat_mffd_type'
  petsc_options_value = 'hypre boomeramg 101 ds'
  nl_max_its = 30
  l_max_its = 60
  l_tol = 1e-04 #4
  nl_rel_tol = 1e-6 #6 #default is 1e-8
  # nl_abs_tol = 1e-14 #only needed when near equilibrium or veeeery small dt
  # start_time = 0
  end_time = 10000 #1e8
  # num_steps = 0
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
    linear_iteration_ratio = 1e5
  []
[]

[Outputs]
  csv = false
  exodus = false
  nemesis = false
  checkpoint = false
  # file_base = MC_GifRif_highDs
[]

# [Debug]
#   show_var_residual_norms = true
# []
