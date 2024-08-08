##############################################################################
# File: 01_normal.i
# File Location: /examples/sintering/paper3/01_test_inputs/15_IC_IW/01_normal
# Created Date: Tuesday August 6th 2024
# Author: Battas,Brandon Scott (bbattas@ufl.edu)
# -----
# Last Modified: Thursday August 8th 2024
# Modified By: Brandon Battas
# -----
# Description:
#  Testing variations of the int_width vs IC int_width to see if i can make it
#   not have to relax or narrow from the IC
#  Both at 1000 here
#
##############################################################################

[Mesh]
  [gmg]
    type = DistributedRectilinearMeshGenerator
    dim = 2
    nx = 80
    ny = 30
    xmin = 0
    xmax = 40000
    ymin = 0
    ymax = 15000
  []
  [subdomain_right]
    type = ParsedSubdomainMeshGenerator
    input = gmg
    combinatorial_geometry = 'x > 20000'
    block_id = 1
  []
  [subdomain_botright]
    type = ParsedSubdomainMeshGenerator
    input = subdomain_right
    combinatorial_geometry = 'y < 7500'
    block_id = 2
    excluded_subdomains = 0
  []
  uniform_refine = 1 #2
  # second_order = true
[]

[GlobalParams]
  profile = TANH
  int_width = 2000
  op_num = 2
  var_name_base = gr
[]

[Variables]
  [wvac]
    initial_condition = 0
  []
  # [wint]
  #   initial_condition = 0
  #   # order = SECOND
  # []
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
  [cvac_aux]
    family = LAGRANGE
    order = FIRST
  []
  [cvac_aux_elem]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [cvac_aux_elem] #NEEDS WORK
    type = MaterialRealAux
    variable = 'cvac_aux_elem'
    property = 'cvac'
  []
  [cvac_aux]
    type = ProjectionAux
    variable = cvac_aux
    v = cvac_aux_elem
  []
[]

[ICs]
  # Grains
  [gr0_IC]
    type = FunctionIC
    variable = gr0
    function = ic_func_gr0
    block = '0 1 2'
  []
  [gr1_IC_L]
    type = FunctionIC
    variable = gr1
    function = ic_func_gr1
    block = '0'
  []
  [gr1_IC_R]
    type = SmoothCircleIC
    variable = gr1
    x1 = 30000
    y1 = 7500
    radius = 5000
    invalue = 0
    outvalue = 1
    block = '1 2'
    int_width = 2000
  []
  [phi_IC]
    type = SmoothCircleIC
    variable = phi
    x1 = 30000
    y1 = 7500
    radius = 5000
    invalue = 1
    outvalue = 0.0
    block = '0 1 2'
    int_width = 2000
  []
[]

[Functions]
  [ic_func_gr0] # Left grain
    type = ParsedFunction
    symbol_names = 'iw x0 y0 r'
    symbol_values = '2000 -5000 7500 15000'
    expression = 'd:=sqrt((x-x0)^2+(y-y0)^2);1-0.5*(1.0-tanh((r-d)/iw))'
  []
  [ic_func_gr1] # Right grain
    type = ParsedFunction
    symbol_names = 'iw x0 y0 r'
    symbol_values = '2000 -5000 7500 15000'
    expression = 'd:=sqrt((x-x0)^2+(y-y0)^2);0.5*(1.0-tanh((r-d)/iw))'
  []
[]

[Modules]
  [PhaseField]
    [GrandPotentialAlt]
      switching_function_names = 'hv hs'
      # Chempot
      chemical_potentials = 'wvac'
      mobilities = 'chiD' #cons_mob
      anisotropic = 'false'
      susceptibilities = 'chi'
      free_energies_w = 'rhov rhos'
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
    prop_values = '6.569e2 6.569e2 6.569e2 6.569e2' # 5.461e7 5.461e7' # No Irradiation
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
    expression = '16 * ( (gr0 * gr1)^2 )'
    # expression = 'hg:=16 * ( (gr0 * gr1)^2 );
    #               if(hg>1e-8,hg,0.0)' #+ (gr0 * etad0)^2 + (gr1 * etad0)^2)
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
  [chiD]
    type = GrandPotentialIsoMaterial
    f_name = chiD #chiuD
    solid_mobility = L #CHANGED FROM L
    void_mobility = Lv
    chi = chi
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
  [sintering]
    type = GrandPotentialSinteringMaterial
    chemical_potential = wvac
    void_op = phi
    Temperature = T
    surface_energy = gb_e_mat
    grainboundary_energy = gb_e_mat
    void_energy_coefficient = kvu
    solid_energy_coefficient = ksu
    solid_energy_model = PARABOLIC
    equilibrium_vacancy_concentration = cv_eq
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
  [cvac]
    type = ParsedMaterial
    property_name = cvac
    material_property_names = 'hs rhos hv rhov Va'
    expression = 'Va*(hs*rhos + hv*rhov)'
    outputs = exodus
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
  end_time = 1e5
  # num_steps = 1
  # dt = 0.001
  # steady_state_detection = true
  # # From tonks ode input
  automatic_scaling = true
  compute_scaling_once = false
  # line_search = none
  # dt = 1.0
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 6
    dt = 1 #0.001
  []
[]

[Outputs]
  csv = false
  exodus = true
  checkpoint = false
  # file_base = fullList
[]

# [Debug]
#   show_var_residual_norms = true
# []
