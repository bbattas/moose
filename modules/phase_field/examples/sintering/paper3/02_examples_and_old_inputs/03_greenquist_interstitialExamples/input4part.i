[Mesh]
  type = GeneratedMesh
  dim = 2
  xmin = -560
  xmax =  560
  nx = 14
  ymin = -560
  ymax =  560
  ny = 14
  uniform_refine = 1
[]

[GlobalParams]
  op_num = 4
  var_name_base = gr
  int_width = 80
  profile = TANH
[]

[Variables]
  [./wv]
  [../]
  [./wi]
  [../]
  [./PolycrystalVariables]
  [../]
  [./phi]
  [../]
[]

[Modules]
  [./PhaseField]
    [./GrandPotential]
      switching_function_names = 'hm hb'

      chemical_potentials = 'wi    wv'
      mobilities =          'chiDi chiDv'
      anisotropy =          '0     1'
      susceptibilities =    'chii  chiv'
      free_energies_w = 'rhomi rhobi rhomv rhobv'

      free_energies_gr = 'omegam omegab'
      mobility_name_gr = Lm
      gamma_gr = gamma
      kappa_gr = kappa

      additional_ops = 'phi'
      free_energies_op = 'omegam omegab'
      mobility_name_op = Lb
      gamma_grxop = gamma
      kappa_op = kappa
    [../]
  [../]
[]

[Kernels]
  [./barrier_phi]
    type = ACBarrierFunction
    variable = phi
    v = 'gr0 gr1 gr2 gr3'
    gamma = gamma
    mob_name = Lb
  [../]
  [./kappa_phi]
    type = ACKappaFunction
    variable = phi
    mob_name = Lb
    kappa_name = kappa
  [../]
  [./recom_wv]
    type = GrandPotentialRecombination
    variable = wv
    rho = rhomv
    rho_r = rhomi
    D = Di
    omega = 0.04092
    a0 = 0.25
    Z = 50
    hm = hm
    args = 'wv wi gr0 gr1 gr2 gr3'
  [../]
  [./recom_wi]
    type = GrandPotentialRecombination
    variable = wi
    rho = rhomi
    rho_r = rhomv
    D = Di
    omega = 0.04092
    a0 = 0.25
    Z = 50
    hm = hm
    args = 'wv wi gr0 gr1 gr2 gr3'
  [../]
[]

[AuxVariables]
  [./bnds]
  [../]
  [./T]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./box]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./bnds_aux]
    type = BndsCalcAux
    variable = bnds
    execute_on = 'Initial Timestep_end'
  [../]
  [./T_aux]
    type = FunctionAux
    variable = T
    function = f_T
    execute_on = 'Timestep_begin'
  [../]
[]

[ICs]
  #VARIABLES
  [./gr0_IC]
    type = SmoothCircleIC
    variable = gr0
    x1 = -200
    y1 = -200
    radius = 200
    invalue = 1
    outvalue = 0
  [../]
  [./gr1_IC]
    type = SmoothCircleIC
    variable = gr1
    x1 =  200
    y1 = -200
    radius = 200
    invalue = 1
    outvalue = 0
  [../]
  [./gr2_IC]
    type = SmoothCircleIC
    variable = gr2
    x1 = -200
    y1 = 200
    radius = 200
    invalue = 1
    outvalue = 0
  [../]
  [./gr3_IC]
    type = SmoothCircleIC
    variable = gr3
    x1 = 200
    y1 = 200
    radius = 200
    invalue = 1
    outvalue = 0
  [../]
  [./phi_IC]
    type = SpecifiedSmoothCircleIC
    variable = phi
    x_positions = '-200 200 -200 200'
    y_positions = '-200 -200 200 200'
    z_positions = '   0   0    0   0'
    radii = '200 200 200 200'
    invalue = 0
    outvalue = 1
  [../]
  #AUXVARIABLES
  [./T_IC]
    type = FunctionIC
    variable = T
    function = f_T
  [../]
  [./box_IC]
    type = BoundingBoxIC
    variable = box
    x1 = -200
    x2 = 200
    y1 = -200
    y2 = 200
    inside = 1
    outside = 0
  [../]
[]

[Functions]
  [./f_T]
    type = PiecewiseLinear
    x = '0 50 200'
    y = '1000 1700 1700'
  [../]
[]

[Materials]
  [./constants]
    type = GenericConstantMaterial
    prop_names =  'kB'
    prop_values = '8.617343e-5'
  [../]
  [./kmv]
    type = ParsedMaterial
    f_name = kmv
    args = 'T'
    constant_names = 'a b'
    constant_expressions = '-0.0025 157.16'
    function = 'a * T + b'
  [../]
  [./kmi]
    type = ParsedMaterial
    f_name = kmi
    args = 'T'
    constant_names = 'a b'
    constant_expressions = '-0.0062 504.98'
    function = 'a * T + b'
  [../]
  [./kbv]
    type = ParsedMaterial
    f_name = kbv
    material_property_names = 'kmv'
    function = '10 * kmv'
  [../]
  [./kbi]
    type = ParsedMaterial
    f_name = kbi
    material_property_names = 'kmi'
    function = '10 * kmi'
  [../]
  [./cmv_eq]
    type = GBSegregation
    c_name = cmv_eq
    formation_energy = -2.9714
    c_prefactor = 1.1722e-15
    GB_energy = 1.5
    Temperature = T
  [../]
  [./cmi_eq]
    type = GBSegregation
    c_name = cmi_eq
    formation_energy = 6.9
    GB_energy = 1.5
    Temperature = T
  [../]
  [./sintering]
    type = GrandPotentialRadiationDamage
    chemical_potentials = 'wi wv'
    void_op = phi
    Temperature = T
    void_energy_coefficient = 'kbi kbv'
    solid_energy_coefficient = 'kmi kmv'
    solid_energy_model = PARABOLIC
    surface_energy = 19.7
    GB_energy = 9.86
    atomic_volume = 0.04092
    solid_equilibrium_concentrations = 'cmi_eq cmv_eq'
    void_equilibrium_concentrations = '0.0 1.0'
    interface_switch_value = 0.3
    switching_base = h
    potential_base = omega
    void_subscript = b
    solid_subscript = m
    potential_subscripts = 'i v'
    susceptibility_base = chi
    kappa_name = kappa
    gamma_name = gamma
  [../]
  [./chiDv]
    type = GrandPotentialTensorMaterial
    f_name = chiDv
    solid_mobility = Lm
    void_mobility = Lb
    chi = chiv
    surface_energy = 19.7
    c = phi
    T = T
    D0 = 8.33e9
    GBmob0 = 1.4759e9
    Q = 2.77
    Em = 3.608
    bulkindex = 1
    gbindex = 1e6
    surfindex = 1e9
  [../]
  [./chiDi]
    type = GrandPotentialScalarMaterial
    D_name = Di
    M_name = chiDi
    chemical_potential = wi
    void_variable = phi
    Temperature = T
    susceptibility = chii
    D0 = 1e13
    migration_energy = 2.0 #Matzke 1987
    bulk_index = 1
    GB_index = 1e6
    surface_index = 1e9
  [../]
  # Output Only Materials
  [./cv]
    type = ParsedMaterial
    f_name = cv
    material_property_names = 'hm rhomv hb rhobv'
    constant_names = 'Va'
    constant_expressions = '0.04092'
    function = 'Va * (hm * rhomv + hb * rhobv)'
    outputs = exodus
  [../]
  [./ci]
    type = ParsedMaterial
    f_name = ci
    material_property_names = 'hm rhomi hb rhobi'
    constant_names = 'Va'
    constant_expressions = '0.04092'
    function = 'Va * (hm * rhomi + hb * rhobi)'
    outputs = exodus
  [../]
  [./cv_var]
    type = ParsedMaterial
    f_name = cv_var
    material_property_names = 'cv cmv_eq hm hb'
    function = 'cv - (hm * cmv_eq + hb)'
    outputs = exodus
  [../]
  [./ci_var]
    type = ParsedMaterial
    f_name = ci_var
    material_property_names = 'ci cmi_eq hm hb'
    function = 'ci - (hm * cmi_eq)'
    outputs = exodus
  [../]
  [./boxphi]
    type = ParsedMaterial
    f_name = boxphi
    args = 'phi box'
    function = 'phi * box'
  [../]
[]

[Postprocessors]
  [./memory]
    type = MemoryUsage
    mem_type = physical_memory
  [../]
  [./DOFs]
    type = NumDOFs
  [../]
  [./cv_tot]
    type = ElementIntegralMaterialProperty
    mat_prop = cv
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./ci_tot]
    type = ElementIntegralMaterialProperty
    mat_prop = ci
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./pore_vol]
    type = ElementIntegralMaterialProperty
    mat_prop = boxphi
    execute_on = 'initial timestep_end'
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -ksp_gmres_restart -sub_ksp_type'
  petsc_options_value = ' asm      lu           1               31                 preonly'
  nl_max_its = 25
  l_max_its = 30
  l_tol = 1e-4
  nl_rel_tol = 1e-7
  nl_abs_tol = 1e-5
  start_time = 0
  end_time = 200
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-3
    optimal_iterations = 10
    iteration_window = 2
    growth_factor = 1.5
    cutback_factor = 0.5
  [../]
  [./Adaptivity]
    max_h_level = 2
    initial_adaptivity = 1
    coarsen_fraction = 0.05
    refine_fraction = 0.95
    weight_names = 'wv wi phi gr0 gr1 gr2 gr3'
    weight_values = '0  0  1   1   1   1   1'
  [../]
[]

[Outputs]
  exodus = true
  csv = true
[]
