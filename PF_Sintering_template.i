[Mesh]
    type = GeneratedMesh
    dim = 2
    nx =
    ny =
    nz = 0
    xmin =
    xmax =
    ymin =
    ymax =
    elem_type = QUAD4
[]

[Variables]

    [./c]
        outputs = exodus
        [./InitialCondition]
            type = FunctionIC
            function = c_txt
        [../]
    [../]
    [./w] # functional derivative of the total free energy  
        order = FIRST
        family = LAGRANGE
    [../]
[]

[Kernels]

    # Order parameter c
    [./c_res]
        type = SplitCHParsed
        variable = c
        f_name = f_tot
        kappa_name = kappa_c
        coupled_variables = 
        w = w
    [../]
    # variable w (derivative of the total free energy)
    [./w_res]
        type = SplitCHWRes
        variable = w
        mob_name = M
    [../]
    [./time]
        type = CoupledTimeDerivative
        variable = w
        v = c
    [../]
[]

[BCs]

[]

[Materials]
    [./consts]
        type = GenericConstantMaterial
        prop_names  = 'L kappa_eta M kappa_c'
        prop_values =
    [../]
    [./free_energy]
        type = DerivativeParsedMaterial
        block = 0
        property_name = f_tot
        coupled_variables =
        constant_names = 'A B'
        constant_expressions = 
        expression = 
        enable_jit = true
        derivative_order = 2
        #outputs = exodus
    [../]
[]

[Functions]

    [c_txt]
        type = PiecewiseMultilinear
        data_file = data/c.txt
    []
[]

[Preconditioning]
    # This preconditioner makes sure the Jacobian Matrix is fully populated. Our
    # kernels compute all Jacobian matrix entries.
    # This allows us to use the Newton solver below.
    [./SMP]
        type = SMP
        full = true
    [../]
[]

[Executioner]
    type = Transient
    scheme = 'bdf2'

    # Automatic differentiation provides a _full_ Jacobian in this example
    # so we can safely use NEWTON for a fast solve
    solve_type = 'NEWTON'

    l_max_its = 20
    l_tol =
    l_abs_tol =

    nl_max_its = 10
    nl_rel_tol =
    nl_abs_tol =

    start_time = 0.0
    num_steps =

    [./TimeStepper]
        type = SolutionTimeAdaptiveDT
        dt =
    [../]
[]

[Postprocessors]

    [c_pp]
      type = ElementAverageValue
      variable = c
    []
[]

[UserObjects]
    [./grain_tracker]
        type = GrainTracker
        variable = 
        outputs = none
        compute_var_to_feature_map = true
        execute_on = 'initial timestep_begin'
    [../]
[]

[Outputs]
    execute_on = 'initial timestep_end'
    exodus = true
    [./other]
        type = VTK
        execute_on = 'TIMESTEP_END'
    [../]
    [console]
        type = Console
        execute_on = 'nonlinear'
        all_variable_norms = true
        max_rows = 3
        show = 
    []
    [./csv]
        type = CSV
        execute_on = 'TIMESTEP_END'
    [../]
[]
