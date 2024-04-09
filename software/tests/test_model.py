import diaux.model
import numpy as np 

def assert_FPA_try_block(suballocation, exception_type, **kwargs):
    """
    Helper function that applies a try/catch block of an FPA, seeking a specific
    exception
    """
    try:
        diaux.model.FluxParityAllocator(suballocation, **kwargs)
        assert False
    except exception_type:
        assert True

def test_initialize_FPA():
    """
    Tests that the flux-parity allocator appropriately initiates with supplied 
    constants.
    """
    suballocation = {'nu_max':[0, 1]}

    # Ensure initialization fails when not specifying the type of allocation
    assert_FPA_try_block(suballocation, RuntimeError)

    ### STATIC STRATEGY TESTS
    # Ensure that static definition fails early if alpha is not provided.
    suballocation['strategy'] = 'static'
    assert_FPA_try_block(suballocation, RuntimeError)

    # Ensure failure if number of alphas is not equal to the number of nutrients
    alphas = [[0, 0, 1], [1]]
    for a in alphas:
        suballocation['alpha'] = a
        assert_FPA_try_block(suballocation, ValueError)

    # Ensure failure if suballocation does not sum to 1
    alphas = [1, 2]
    suballocation['alpha'] = a
    assert_FPA_try_block(suballocation, ValueError)

    ## HIERARCHICAL STRATEGY TESTS
    suballocation['strategy'] = 'hierarchical'
    suballocation['K'] = [1, 1]
    assert_FPA_try_block(suballocation, KeyError)
    del suballocation['K'] 
    suballocation['n'] = [1,1]
    assert_FPA_try_block(suballocation, KeyError)

    ## PROPORTIONAL AND METABOLIC HIERARCHY TESTS
    # If metabolic hierarchy is passed, ensure checks are in place to calculate
    # metabolic hierarchy 
    suballocation['strategy'] = 'proportional'
    alphas = [0, 1]
    suballocation['K'] = [1, 1]
    del suballocation['n']
    assert_FPA_try_block(suballocation, KeyError, metabolic_hierarchy=True)
    del suballocation['K'] 
    suballocation['n']=[1, 1]
    assert_FPA_try_block(suballocation, KeyError, metabolic_hierarchy=True)

    ## CONSTANT OVERRIDES AND BOUND ENFORCEMENT
    # Assert that the useful fraction must be bounded between 0 and o
    suballocation['frac_useful'] = [-1, 1]
    assert_FPA_try_block(suballocation, ValueError)
    suballocation['frac_useful'] = [1, 10]
    assert_FPA_try_block(suballocation, ValueError)
    del suballocation['frac_useful'] 

    # Ensure that constants are being overriden if provided
    constants = {'gamma_max': 0,
                 'phi_O': 0,
                 'tau': 0,
                 'kappa_max': 0,
                 'Km_u':0,
                 'Km_c':0,
                 'Km': [0, 0],
                 'Y': [0, 0]}
    FPA = diaux.model.FluxParityAllocator(suballocation, constants=constants)
    for k, v in constants.items():
        assert getattr(FPA, k) == v

def test_FPA_property_calculation():
    """
    Asserts that calculation of the properties in the flux-parity framework 
    yield expected values. 
    """

    # Define the suballocation for the self replicator
    suballocation = {'strategy': 'static',
                     'alpha': [0.0, 1.0],
                     'nu_max': [1.0, 2.0],
                     'hierarchy': [1, 0]}
    # Define simple constants to ensure values are easily calculated
    cst = {'gamma_max': 1,
           'phi_O': 0.5,
           'tau': 2,
           'kappa_max': 3,
           'Km_u': 4,
           'Km_c': 5,
           'Km': [6, 7],
           'Y': [8, 9]}
    allocator = diaux.model.FluxParityAllocator(suballocation, constants=cst)

    # Ensure that tRNA concentration bound limits result in approprate corrections
    # to allocation
    tRNA_u = 0
    tRNA_c = 1
    nutrients = [10, 11]
    allocator.compute_properties(tRNA_c, tRNA_u, nutrients)
    assert allocator.phi_Rb == 1 - cst['phi_O']
    assert allocator.ratio == np.inf
    tRNA_c = 0
    allocator.compute_properties(tRNA_c, tRNA_u, nutrients)
    assert allocator.phi_Rb == 0
    assert allocator.ratio == 0

    # Assert that all properties are calculated correctly
    tRNA_c = 1.4
    tRNA_u = 1.9
    allocator.compute_properties(tRNA_c, tRNA_u, nutrients)
    assert allocator.ratio == tRNA_c / tRNA_u
    assert allocator.phi_Rb == (1 - cst['phi_O']) * (allocator.ratio / (cst['tau'] + allocator.ratio))
    assert allocator.kappa == cst['kappa_max'] * (allocator.ratio / (cst['tau'] + allocator.ratio))
    assert allocator.gamma == cst['gamma_max'] * (tRNA_c / (tRNA_c + cst['Km_c']))

    # Test that the inverted hierarchy is respected and metabolic rates are appropriately tuned 
    metab_factor = tRNA_u / (tRNA_u + cst['Km_u'])
    env_factor = np.array([nutrients[1] / (nutrients[1] + cst['Km'][1]),
                  nutrients[0] / (nutrients[0] + cst['Km'][0])])
    nu = np.array([env_factor[i] * metab_factor * suballocation['nu_max'][i] for i in range(2)])
    for i in range(2):
        assert allocator.nu[i] == nu[i]

    # Confirm that running this with the metabolic hierarchy in place results in 
    # a hierarchically controlled metabolic rate
    suballocation['K'] = [12.0,13.0]
    suballocation['n'] = [1,1]
    mh_allocator = diaux.model.FluxParityAllocator(suballocation,
                                                   constants=cst,
                                                   metabolic_hierarchy=True)
    mh_allocator.compute_properties(tRNA_c, tRNA_u, nutrients)
    numer = (nutrients[1] / suballocation['K'][1])**suballocation['n'][1]
    hierarchy_factor = 1 - (numer / (1 + numer))
    assert mh_allocator.nu[0] == hierarchy_factor * nu[0]
    del suballocation['K']
    del suballocation['n']

    ## METABOLIC ALLOCATION CHECKS
    # First define fixed things true for all strategies
    ratio = tRNA_c / tRNA_u
    phi_Rb = (1 - cst['phi_O']) * (ratio / (ratio + cst['tau']))
    phi_Mb = 1 - cst['phi_O'] - phi_Rb

    # Check that metabolic suballocation follows prescribed alphas.
    suballocation['strategy'] = 'static'
    suballocation['alpha'] = np.array([0.1, 0.9])
    allocator = diaux.model.FluxParityAllocator(suballocation, constants=cst)
    allocator.compute_properties(tRNA_c, tRNA_u, nutrients) 
    phi_Mbs = suballocation['alpha'] * phi_Mb 
    for phi_true, phi_test in zip(phi_Mbs, allocator.phi_Mb): 
        assert phi_test == phi_true
    
    # Set the proportional allocator case with skewed alphas
    suballocation['strategy'] = 'proportional' 
    allocator = diaux.model.FluxParityAllocator(suballocation, constants=cst)
    nutrients = [1.0, 9.0]
    alpha_true = nutrients / np.sum(nutrients)
    allocator.compute_properties(tRNA_c, tRNA_u, nutrients) 
    phi_Mbs = alpha_true * phi_Mb
    for phi_true, phi_test in zip(phi_Mbs, allocator.phi_Mb):
        assert phi_test == phi_true

    # Test that the proportional case works with zero nutrients.
    allocator.compute_properties(tRNA_c, tRNA_u, [0.0 ,0.0])
    assert allocator.alpha[0] == 1
    assert allocator.alpha[1] == 0
    
    # Set the hierarchical allocator case
    suballocation['strategy'] = 'hierarchical'
    suballocation['K'] = [12, 13]
    suballocation['n'] = [1, 1]
    allocator = diaux.model.FluxParityAllocator(suballocation, constants=cst)
    nutrients = [3.0, 4.0] 
    allocator.compute_properties(tRNA_c, tRNA_u, nutrients)

    numer = (nutrients[1]/suballocation['K'][1])**suballocation['n'][1]
    alpha = numer / (1 + numer)
    alpha_true = np.array([1 - alpha, alpha])
    phi_Mbs = alpha_true * phi_Mb
    for phi_true, phi_test in zip(phi_Mbs, allocator.phi_Mb):
        assert phi_test == phi_true


def test_allocation_shuffled_hierarchy():
    """
    Will test if the metabolic suballocation works as advertised for a long and 
    shuffled hierarchical consumption nutrients.
    """
    hierarchy = [4, 2, 0, 3, 1]
    nu_max = [1, 2, 3, 4, 5]

    suballocation = {'strategy': 'hierarchical', 
                     'nu_max': nu_max,
                     'hierarchy': hierarchy,
                     'n': [1, 1, 2, 2, 2],
                     'K': [5, 4, 3, 2, 1]}
    cst = {'gamma_max': 1,
           'phi_O': 0.5,
           'tau': 2,
           'kappa_max': 3,
           'Km_u': 4,
           'Km_c': 5,
           'Km': [6, 7, 8, 9, 10],
           'Y': [10, 9, 7, 8, 6]}
    nutrients = np.array([1, 2, 3, 4, 5]).astype(float)
    tRNA_u = 0.013
    tRNA_c = 0.013
    allocator = diaux.model.FluxParityAllocator(suballocation, constants=cst)
    allocator.compute_properties(tRNA_c, tRNA_u, nutrients) 

    # Determine the suballocation functions for each 
    numers = [(c/suballocation['K'][i])**suballocation['n'][i] for i, c in enumerate(nutrients)]
    alloc_fn = [(numer / (1 + numer)) for numer in numers]
    allocated_idx = []

    # USing defined heirarchy, explictly compute the alpha
    hier_0 = alloc_fn[2]
    hier_1 = alloc_fn[4] * (1 - hier_0)
    hier_2 = alloc_fn[1] * (1 - hier_0 - hier_1)
    hier_3 = alloc_fn[3] * (1 - hier_0 - hier_1 - hier_2)
    hier_4 = 1 - hier_0 - hier_1 - hier_2 - hier_3  
    alpha = [hier_4, hier_2, hier_0, hier_3,  hier_1]

    # FIXME: There appears to be numeric precision errors here, maybe due to the
    # way that the multiplication and constraints are being enforced. I don't think 
    # that this really functionally matters, since they are the same to within 8 
    # decimals, but I would like to understand why.
    for alpha_true, alpha_test in zip(alpha, allocator.alpha):
        assert np.round(alpha_true, decimals=8) == np.round(alpha_test, decimals=8)
    
    # Check the total metabolic allocation is correct
    ratio = tRNA_c/tRNA_u
    phi_Rb = (1 - cst['phi_O']) * (ratio / (cst['tau'] + ratio))
    phi_Mb = (1 - cst['phi_O'] - phi_Rb) * np.array(alpha)
    for phi_true, phi_test in zip(phi_Mb, allocator.phi_Mb):
        assert np.round(phi_true, decimals=8) == np.round(phi_test, decimals=8) 


def test_FPA_derivative_calc():
    """
    Tests that the derivatives are correctly calculated and returned
    """
    suballocation = {'strategy': 'hierarchical', 
                     'nu_max':  [1, 2],
                     'hierarchy': [1, 0],
                     'K': [5, 7],
                     'n': [1, 2],
                     'death_rate': 0.1}

    # Define simple constants to ensure values are easily calculated
    cst = {'gamma_max': 1,
           'phi_O': 0.5,
           'tau': 2,
           'kappa_max': 3,
           'Km_u': 4,
           'Km_c': 5,
           'Km': [6, 7],
           'Y': [8, 9]}
    nutrients = [3, 4]
    # tRNA_c = 0.01
    # tRNA_u = 0.003
    # Instantiate the system
    allocator = diaux.model.FluxParityAllocator(suballocation, constants=cst)
    
    ## TEST POSITIVITY BOUNDING
    nutrients = [-3, -4]
    masses = np.array([1, -1, 1, -1, 1, -1])
    _ = allocator.compute_derivatives(masses, nutrients)
    assert np.all(allocator.nutrients >= 0)
    assert np.all(allocator.M_Mb == [1, 0])
    assert allocator.M_Rb == 1
    assert allocator.M_O == 0
    assert allocator.M == 2 

    ## TEST CALCULATION OF DERIVATIVES
    nutrients = np.array([4, 8]) 
    masses = np.array([1, 2, 3, 4, 0.5, 0.6])

    # Define tRNA concentrations from provided masses 
    tRNA_u = masses[-2]/np.sum(masses[:-2])
    tRNA_c = masses[-1]/np.sum(masses[:-2])

    # Compute the properites, which has already been tested to be accurate
    allocator.compute_properties(tRNA_c, tRNA_u, nutrients) 
    M = np.sum(masses[:-2])
    dM_dt = allocator.gamma * masses[2]
    dM_Rb = allocator.phi_Rb * dM_dt - allocator.death_rate * masses[2]
    dM_Mb = [allocator.phi_Mb[i] * dM_dt  - allocator.death_rate * masses[i] for i in range(2)]
    dM_O = cst['phi_O'] * dM_dt - allocator.death_rate * masses[3]
    dm_u = allocator.kappa * M - np.sum(allocator.nu * masses[:2]) + dM_dt - allocator.death_rate * masses[-2]
    dm_c = np.sum(allocator.nu * masses[:2]) - dM_dt - allocator.death_rate * masses[-1]
    mass_derivs = [dM_Mb[0], dM_Mb[1], dM_Rb, dM_O, dm_u, dm_c]
    nut_derivs = (-allocator.nu * masses[:2])/allocator.Y

    dM, dc = allocator.compute_derivatives(masses, nutrients)
    assert np.all(dM == mass_derivs)
    assert np.all(nut_derivs == dc)