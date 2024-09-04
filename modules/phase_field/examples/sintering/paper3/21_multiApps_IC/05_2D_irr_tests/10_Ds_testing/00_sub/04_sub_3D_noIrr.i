##############################################################################
# File: 04_sub_3D_noIrr.i
# File Location: /examples/sintering/paper3/21_multiApps_IC/05_2D_irr_tests/10_Ds_testing/00_sub
# Created Date: Wednesday September 4th 2024
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Wednesday September 4th 2024
# Modified By: Brandon Battas
# -----
# Description:
#  Manually defined 3D version of the simple 1 pore input, still 2 grains with
#   a GB plane with spherical pore in the center
#  Using the thermal eq c values
#  Has about 1.8m DoFs so like 150-180 cpus should be fine (12k-10k)
##############################################################################

[Mesh]
  [gmg]
    type = DistributedRectilinearMeshGenerator
    dim = 3
    nx = 100#120
    ny = 60#80
    nz = 60
    xmin = 0
    xmax = 25000#30000
    ymin = 0
    ymax = 15000#20000
    zmin = 0
    zmax = 15000
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

[Variables]
  [wvac]
    initial_condition = 0
  []
  [wint]
    initial_condition = 0
  []
  # [cvac_var]
  # []
  # [cint_var]
  #   # order = SECOND
  # []
  [gr0]
  []
  [gr1]
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
  # Vacancy mat to variable
  [cvac_aux]
    family = LAGRANGE
    order = FIRST
  []
  [cvac_aux_elem]
    family = MONOMIAL
    order = CONSTANT
  []
  # Interstitial mat to variable
  [cint_aux]
    family = LAGRANGE
    order = FIRST
  []
  [cint_aux_elem]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  # Vacancy
  [cvac_aux_elem]
    type = MaterialRealAux
    variable = 'cvac_aux_elem'
    property = 'cvac'
  []
  [cvac_aux]
    type = ProjectionAux
    variable = cvac_aux
    v = cvac_aux_elem
  []
  # Interstitial
  [cint_aux_elem]
    type = MaterialRealAux
    variable = 'cint_aux_elem'
    property = 'cint'
  []
  [cint_aux]
    type = ProjectionAux
    variable = cint_aux
    v = cint_aux_elem
  []
[]

[ICs]
  # Grains
  [gr0_IC]
    type = FunctionIC
    variable = gr0
    function = ic_func_gr0
    # block = '0 1' #2
  []
  [gr1_IC]
    type = FunctionIC
    variable = gr1
    function = ic_func_gr1
    # block = '0 1' #2
  []
  [phi_IC]
    type = FunctionIC
    variable = phi
    function = ic_func_pore
    # block = '0 1'
  []
[]

[Functions]
  [ic_func_pore] # Pore and exterior right
    type = ParsedFunction
    symbol_names = 'iw x0 y0 z0 r xr'
    symbol_values = '1000 10000 7500 7500 2200 20000'
    expression = 'd:=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2);
    1-0.5*(1.0-tanh((r-d)/iw)) + 0.5*(1.0+tanh((x-xr)/iw))'
  []
  [ic_func_gr0] # Top grain
    type = ParsedFunction
    symbol_names = 'iw x0 y0 z0 r xr'
    symbol_values = '1000 10000 7500 7500 2200 20000'
    expression = 'd:=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2);
    0.5*(1.0-tanh((r-d)/iw)) * 0.5*(1.0-tanh((x-xr)/iw)) * 0.5*(1.0-tanh((y-y0)/iw))'
  []
  [ic_func_gr1] # Top grain
    type = ParsedFunction
    symbol_names = 'iw x0 y0 z0 r xr'
    symbol_values = '1000 10000 7500 7500 2200 20000'
    expression = 'd:=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2);
    0.5*(1.0-tanh((r-d)/iw)) * 0.5*(1.0-tanh((x-xr)/iw)) * 0.5*(1.0+tanh((y-y0)/iw))'
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
      # hj_over_kVa = 'hoverk_vu hoverk_su hoverk_vi hoverk_si' #'hv_over_kVa hs_over_kVa' #
      # hj_c_min = 'cvueq_mask csueq_mask cvieq_mask csieq_mask' #cvueq_mask=hv*1 'hv_c_min hs_c_min' #
    []
  []
[]

[Materials]
  [consts]
    type = GenericConstantMaterial
    prop_names = 'Va negOverVa' #cvieq_mask cvueq_mask csieq_mask csueq_mask'
    prop_values = '0.04092 -24.4379' #0.0 0.0 0.0 0.0'
  []
  [k_constants]
    type = GenericConstantMaterial
    prop_names = 'ksu kvu ksi kvi' # Using the GB based values (lowest of mine)
    # prop_values = '7.751e2 7.751e2 5.711e5 5.711e5' # Irradiation
    prop_values = '6.569e2 6.569e2 5.461e7 5.461e7' # No Irradiation
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
  [hgb]
    type = DerivativeParsedMaterial
    property_name = hgb
    # derivative_order = 2
    coupled_variables = 'gr0 gr1'
    # expression = '16 * ( (gr0 * gr1)^2 )'
    expression = 'hg:=16 * ( (gr0 * gr1)^2 );
                  if(hg>1e-4,hg,0.0)' # + (gr0 * gr2)^2 + (gr1 * gr2)^2
    # expression = '4*(1 - (gr0^2 + gr1^2 + etad0^2))^2' 0.0001
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
    GBwidth = 1.0
    surf_thickness = 1.0 #0.5
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
    GBwidth = 1.0
    surf_thickness = 1.0 #0.5
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
    type = DerivativeParsedMaterial
    property_name = cv_eq
    # derivative_order = 2
    coupled_variables = 'gr0 gr1 phi wvac'
    material_property_names = 'hgb(phi,gr0,gr1)' #'rhovi(vac) rhosi hv(phi)'
    constant_names = 'cb cgb'
    # constant_expressions = '3.877e-04 4.347e-03' #Irradiation
    constant_expressions = '2.424e-06 5.130e-03' #No Irradiation- LANL
    expression = 'cgb * hgb + (1 - hgb)*cb'
    outputs = none #'nemesis' # + phi^2
  []
  [ci_eq]
    type = DerivativeParsedMaterial
    property_name = ci_eq
    # derivative_order = 2
    coupled_variables = 'gr0 gr1 phi wint'
    material_property_names = 'hgb(phi,gr0,gr1)' # 'rhovi(wint) rhosi(wint) hv(phi)'
    constant_names = 'cb cgb'
    # constant_expressions = '7.258e-09 5.900e-06' #Irradiation
    constant_expressions = '1.667e-32 6.170e-08' #'1.667e-32 6.170e-08' #No Irradiation- LANL
    # constant_expressions = '2.424e-06 5.130e-03' #No Irradiation VACANCY- LANL
    expression = 'cgb * hgb + (1 - hgb)*cb'
    outputs = none #'nemesis' #+ phi^2
  []
  # [cvac_nmc]
  #   type = ParsedMaterial
  #   property_name = cvac_nmc
  #   material_property_names = 'hs rhosu hv rhovu Va'
  #   expression = 'Va*(hs*rhosu + hv*rhovu)'
  #   # outputs = exodus
  # []
  # [cint_nmc]
  #   type = ParsedMaterial
  #   property_name = cint_nmc
  #   material_property_names = 'hs rhosi hv rhovi Va'
  #   expression = 'Va*(hs*rhosi + hv*rhovi)'
  #   # outputs = exodus
  # []
  [cvac]
    type = ParsedMaterial
    property_name = cvac
    material_property_names = 'hs cv_eq hv Va'
    expression = 'cv:=(hs*cv_eq ) + hv*1.0;
                  if(cv<0.0, 0.0, cv)'
    # outputs = exodus
  []
  [cint]
    type = ParsedMaterial
    property_name = cint
    material_property_names = 'hs ci_eq hv Va'
    expression = 'ci:=(hs*ci_eq) + hv*0.0;
                  if(ci<0.0, 0.0, ci)'
    # outputs = exodus
  []
  # # CONSERVATION
  # # h / (V*k)
  # [hoverk_vu]
  #   type = DerivativeParsedMaterial
  #   property_name = hoverk_vu
  #   coupled_variables = 'phi gr0 gr1'
  #   derivative_order = 2
  #   material_property_names = 'hv(phi) Va kvu'
  #   expression = 'hv / (Va * kvu)'
  #   # outputs = exodus
  # []
  # [hoverk_su]
  #   type = DerivativeParsedMaterial
  #   property_name = hoverk_su
  #   coupled_variables = 'phi gr0 gr1'
  #   derivative_order = 2
  #   material_property_names = 'hs(phi) Va ksu'
  #   expression = 'hs / (Va * ksu)'
  #   # outputs = exodus
  # []
  # # h*ceq Masks
  # [cvueq_mask] # cvueq_mask = hv*1
  #   type = DerivativeParsedMaterial
  #   property_name = cvueq_mask
  #   coupled_variables = 'phi gr0 gr1'
  #   derivative_order = 2
  #   material_property_names = 'hv(phi)'
  #   expression = 'hv'
  #   # outputs = exodus
  # []
  # [csueq_mask]
  #   type = DerivativeParsedMaterial
  #   property_name = csueq_mask
  #   coupled_variables = 'phi gr0 gr1'
  #   derivative_order = 2
  #   material_property_names = 'hs(phi) cv_eq(phi,gr0,gr1)'
  #   expression = 'hs * cv_eq'
  #   # outputs = exodus
  # []
[]

[Postprocessors]
  # [cv_var_total]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = cvac_var
  # []
  [wvac_total]
    type = ElementIntegralVariablePostprocessor
    variable = wvac
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
  # [gr2_total]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = gr2
  # []
  [hgb_total]
    type = ElementIntegralMaterialProperty
    mat_prop = hgb
  []
  [cvac_total]
    type = ElementIntegralMaterialProperty
    mat_prop = cvac
  []
  [runtime]
    type = PerfGraphData
    section_name = "Root"
    data_type = TOTAL
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
  start_time = 0
  # end_time = 100
  num_steps = 1
  dt = 0.001
  # steady_state_detection = true
  # # From tonks ode input
  automatic_scaling = true
  compute_scaling_once = false
  # line_search = none
  # dt = 1.0
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 6
    dt = 0.1 #0.001
  []
[]

[Outputs]
  csv = false
  exodus = false
  checkpoint = false
  # file_base = fullList
[]

# [Debug]
#   show_var_residual_norms = true
# []
