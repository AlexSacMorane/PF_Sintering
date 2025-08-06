# from https://github.com/idaholab/moose/discussions/20805
## The phase field model to simulate sintering of two particles in 2D

[GlobalParams]        # Automatically define order parameters 
  var_name_base = eta
  op_num = 2
[]

[Mesh]                # Generate mesh
  type = GeneratedMesh
  dim = 2
  nx = 80             # n°s element in x. 
  ny = 40             # n°s element in y.
  xmin = 0  
  xmax = 40      # µm
  ymin = 0
  ymax = 20      # µm
  elem_type = QUAD4
  uniform_refine = 2
[]

[Variables]
  [./c]               # particle concentration (conserved variable)
    order = FIRST
    family = LAGRANGE
  [../]
  [./w]               # functional derivative of the total free energy  
    order = FIRST
    family = LAGRANGE
  [../]
  [./PolycrystalVariables] # Automatically creates order parameter variables
  [../]
[]

[Kernels]
  # variable c 
  [./c_res]
  # compute the residual therm of total free energy as function of c.
    type = SplitCHParsed
    variable = c
    f_name = F
    kappa_name = kappa_c
    args = 'eta0 eta1'
    w = w
  [../]
  [./w_res]
  # compute the residual therm of functional free energy as function of w and M.
    type = SplitCHWRes
    variable = w
    mob_name = M
  [../]
  [./time]
  # compute the time derivative term of w and c.
    type = CoupledTimeDerivative
    variable = w
    v = c
  [../]
  # variable eta0
  [./eta0_dot]
    type = TimeDerivative
    variable = eta0
  [../]
  [./acint_eta0]
    type = ACInterface
    variable = eta0
    mob_name = L
    args = c
    kappa_name = kappa_eta
  [../]
  [./acbulk_eta0]
    type = AllenCahn
    variable = eta0
    mob_name = L
    f_name = F
    args = 'c eta1'
  [../]
  # variable eta1
  [./eta1_dot]
    type = TimeDerivative
    variable = eta1
  [../]
  [./acint_eta1]
    type = ACInterface
    variable = eta1
    mob_name = L
    args = c
    kappa_name = kappa_eta
  [../]
  [./acbulk_eta1]
    type = AllenCahn
    variable = eta1
    mob_name = L
    f_name = F
    args = 'c eta0'
  [../]
[]

[Materials]  
  [./constants]
    type = GenericFunctionMaterial
    prop_names = 'kappa_c  kappa_eta  A     B    L
                    Dvol   Dvap       Dsuf  Dgb  M'
    
    prop_values = '1.0     0.5        16.0  1.0  1.0
                   0.01    0.001      4.0   0.4  1.0'
  [../]
  [./free_energy]            # Total free energy function of c and eta.
    type = DerivativeParsedMaterial
    f_name = F
    args = 'c eta0 eta1'          # Must be changed as op_num changes.
    function = 'A*c^2*(1-c)^2 + B*(c^2+6*(1-c)*(eta0^2+eta1^2)-4*(2-c)*(eta0^3+eta1^3)
                + 3*(eta0^2+eta1^2)^2)'
    material_property_names = 'A  B'
    derivative_order = 2
  [../]
  #[./CHmobility]    # Cahn and Hilliard Mobility incluing all diffusin path (m^5/J/s)
  #  type = DerivativeParsedMaterial
  #  f_name = M
  #  args = 'c eta0 eta1'
  #  function = 'Dvol*c^3*(10-15*c+6*c^2) + Dvap*(1-(c^3*(10-15*c+6*c^2))) + Dsuf*(3*c^2-2*c^3)*(1-3*c^2-2*c^3) + Dgb*(eta0*eta1)'
  #  material_property_names = 'Dvol Dvap Dsuf Dgb'
  #  derivative_order = 1
  #[../]
[]

[ICs]
  [./ic_eta0]
  # initial condition for each order parameter 
    type = SmoothCircleIC
    variable = eta0
    x1 = 10.5
    y1 = 10.0
    radius = 7.5
    invalue = 1.0
    outvalue = 0.0
    int_width = 2.0
  [../]
  [./IC_eta1]
    type = SmoothCircleIC
    variable = eta1
    x1 = 25.5
    y1 = 10.0
    radius = 7.5
    invalue = 1.0
    outvalue = 0.0
    int_width = 2.0
  [../]
  [./ic_c]
  # initial condition for conserved parameter 
    type = SpecifiedSmoothCircleIC
    invalue = 1.0
    outvalue = 0.0
    int_width = 2.0
    x_positions = '10.5 25.5'
    y_positions = '10.0 10.0'
    z_positions = '0.0 0.0'
    radii       = '7.5 7.5'
    3D_spheres = false
    variable = c
    block = 0
  [../]
[]

[Preconditioning]
  [./coupled]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = NEWTON
  automatic_scaling = true
  petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_ksp_type
                         -sub_pc_type -pc_asm_overlap'
  petsc_options_value = 'asm      31                  preonly
                         lu          1'
  l_tol = 1e-05
  nl_max_its = 30
  l_max_its = 30
  nl_rel_tol = 1e-07
  nl_abs_tol = 1e-09
  start_time = 0.0
  end_time = 100
  #dt = 0.05
  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.5
    dt = 0.01
    growth_factor = 2
    optimal_iterations = 8
  [../]
  [./Adaptivity]
    initial_adaptivity = 1
    refine_fraction = 0.6
    coarsen_fraction = 0.01
    max_h_level = 2
    print_changed_info = true
  [../]
[]

[Outputs]
  exodus = true
  #[./debug]
  #  type = VariableResidualNormsDebugOutput
  #[../]
  [./other]
    type = VTK
    execute_on = 'TIMESTEP_END'
  [../]
[]