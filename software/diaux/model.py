import numpy as np
import scipy.integrate
OD_CONV = 1.15E17 # Amino acids per OD600 unit (approximately)

def steady_state_precursors(gamma_max, 
                            phi_Rb,
                            nu_max,
                            Kd_cpc,
                            phi_O):
    """
    Computes the steady state value of the charged-tRNA abundance.

    Parameters
    ----------
    gamma_max: positive float
        The maximum translational efficiency in units of inverse time.
    phi_Rb: float [0, 1]
        The fraction of the proteome occupied by ribosomal proteins.
    nu_max : positive float 
        The maximum nutritional capacity in units of inverse time. 
    Kd_cpc : positive float
        The effective dissociation constant of the precursors to the ribosome. 
    phi_O : float [0, 1]
        Allocation towards other proteins.
    Returns
    -------
    c_pc : float
        The charged tRNA abunadance given the model parameters. This is defined
        as the mass of the charged-tRNA relative to the total biomass.

    Notes
    -----
    This function assumes that in steady state, the nutrients are in such abundance 
    that the nutritional capacy is equal to its maximal value. 

    """
    ss_lam = steady_state_growth_rate(
        gamma_max, phi_Rb, nu_max, Kd_cpc, phi_O=phi_O)
    cpc = (nu_max * (1 - phi_Rb - phi_O) / ss_lam) - 1
    return cpc


def steady_state_growth_rate(gamma_max,
                             phi_Rb,
                             nu_max,
                             Kd_cpc,
                             phi_O):
    """
    Computes the steady-state growth rate of the self-replicator model. 

    Parameters
    ----------
    gamma_max : positive float 
        The maximum translational capacity in units of inverse time.
    phi_Rb : float [0, 1]
        The fraction of the proteome occupied by ribosomal protein mass
    nu_max : positive float 
        The maximum nutritional capacity in units of inverse time. 
    Kd_cpc :  positive float
        The effective dissociation constant of charged tRNA to the elongating
        ribosome.
    phi_O : float [0, 1]
        Allocation towards other proteins.

    Returns
    -------
    lam : float 
        The physically meaningful growth rate (in units of inverse time) given 
        the provided model parameeters.

    Notes
    -----
    This function assumes that in steady state, the nutrients are in such abundance 
    that the nutritional capacy is equal to its maximal value. 
    """
    Nu = nu_max * (1 - phi_Rb - phi_O)
    Gamma = gamma_max * phi_Rb
    numer = Nu + Gamma - \
        np.sqrt((Nu + Gamma)**2 - 4 * (1 - Kd_cpc) * Nu * Gamma)
    denom = 2 * (1 - Kd_cpc)
    lam = numer / denom
    return lam


def steady_state_gamma(gamma_max,
                       phi_Rb,
                       nu_max,
                       Kd_cpc,
                       phi_O=0):
    """
    Computes the steady-state translational efficiency, gamma.

    Parameters 
    -----------
    gamma_max : positive float
        The maximum translational capacity in units of inverse time.
    phi_Rb : float [0, 1]
        The fraction of the proteome occupied by ribosomal protein mass.
    nu_max : positive float 
        The maximum nutritional capacity in units of inverse time.
    Kd_cpc : positive float 
        The effective dissociation constant of charged tRNA to the elongating
        ribosome.
    phi_O : float [0, 1]
        Allocation towards other proteins

    Returns
    -------
    gamma : positive float
        The translational efficiency in units of inverse time
    """

    c_pc = steady_state_precursors(
        gamma_max, phi_Rb, nu_max, Kd_cpc, phi_O=phi_O)
    return gamma_max * (c_pc / (c_pc + Kd_cpc))


def phiRb_optimal_allocation(gamma_max,
                             nu_max,
                             Kd_cpc,
                             phi_O):
    """
    Computes the optimal fraction of proteome that is occupied by ribosomal 
    proteins which maximizes the growth rate. 

    Parameters
    ----------
    gamma_max : positive float 
        The maximum translational efficiency in units of inverse time.
    nu_max : positive float
        The maximum nutritional capacity in units of inverse time.
    Kd_cpc: positive float 
        The effective dissociation constant of charged tRNA to the elongating 
        ribosome.
    phi_O : float [0, 1]
        Allocation towards other proteins.
    Returns
    -------
    phi_Rb_opt : positive float [0, 1]
        The optimal allocation to ribosomes.
    """
    numer = nu_max * (-2 * Kd_cpc * gamma_max + gamma_max + nu_max) +\
        np.sqrt(Kd_cpc * gamma_max * nu_max) * (gamma_max - nu_max)
    denom = -4 * Kd_cpc * gamma_max * nu_max + \
        gamma_max**2 + 2 * gamma_max * nu_max + nu_max**2
    phi_Rb_opt = (1 - phi_O) * numer / denom
    return phi_Rb_opt

class FluxParityAllocator():
    """Base class for a self replicator obeying flux-parity allocation."""
    def __init__(self, 
                 suballocation, 
                 constants={}, 
                 label=0):
        """
        Instantiates a self replicating organism that undergoes flux-parity regulation 
        of its resource allocation. 

        Parameters
        ----------
        suballocation : dict 
            A dictionary defining the metabolic suballocation strategy of the 
            self replicator. Must have the following keys:
                strategy: str, `'dynamic'` or `'static'`
                    The type of suballocation. If `'dynamic'`, the suballocation
                    of each metabolic sector follows a Monod function with a 
                    Monod constant `K` and sensitivty `n`. If `'static'`,
                    suballocation is deemed to be at a fixed value, supplied
                    with key `alpha.`    
                metabolic_rates: numpy.ndarray of floats [0, inf)
                    The metabolic rate for consumption of each nutrient in units 
                    of inverse hours.
                frac_useful: numpy.ndarray
                    The fraction of each metabolic class that is deemed to be 
                    useful and contributes to the net metabolic flux.
                    Default is assumed to be 1 (all is useful.)
                hierarchy : numpy.ndarray of ints
                    The hierarchy order of the nutrients beginning at 0. This is only used if 
                    `strategy: 'dynamic'` is supplied. If not supplied, hierarchy 
                    is assumed to be in the order of the provided nutrients.
                alpha : `numpy.ndarray` of float, [0, 1]
                    The fixed suballocation of the metabolic strategy. Values 
                    must sum to `1.0`. This is only used if `strategy: 'static'` 
                    is supplied.
                K : numpy.ndarray of float [0, inf)
                    The Monod constant for the suballocation, in units of concentration.
                    This is only used if `strategy: 'dynamic'` is supplied.
                n : numpy.ndarray of int [1, inf) 
                    The sensitivity parameter of the Monod function. This is only 
                    used if `strategy: 'dynamic'` is supplied.

               
        constants : dict, optional
            A dictionary of the constant model parameters. Default values
            (described below) can be overridden by providing a key,value pair.
                gamma_max : float
                    The maximum translation rate in inverse hours. Default 
                    value is 9.65 inv. hr. 
                phi_O : float [0, 1]
                    The allocation towards all non-ribosomal and non-metabolic 
                    proteins. Default value is 0.55. 
               tau : float [0, inf)                               
                    The sensitivity parameter of the flux-parity ribosomal 
                    allocation function. Default is 1.0
               kappa_max : float [0, inf)
                    The maximum transcription rate of uncharged tRNA in 
                    mass abundance units per hour. Default is 1.15E-3 inverse hr.
               Km_u : float [0, inf)                
                    The Michaelis-Menten constant for binding of uncharged tRNA 
                    to the metabolic components in units of relative mass abundance. 
                    Default value is 3E-5.
               Km_c : float [0, inf)                
                    The Michaelis-Menten constant for binding of charged tRNA 
                    to the ribosomal components in units of relative mass abundance. 
                    Default value is 3E-5.
               Km : numpy.ndarray [0, inf)
                    The Monod constant for metabolic activity as a function of
                    the external nutrient concentration in units of concentration. 
                    Default is the same for each nutrient with 5E-6 M.
               Y : numpy.ndarray [0, inf)
                    The yield coefficient for growth on each nutrient. Default 
                    is 2.95E19 for each nutrient supplied.              
               death_rate : float
                    The death rate of the species in units of mass per unit time 
                    (typically hours). Default is 0.

        """
        self.num_metab = len(suballocation['nu_max'])
        self.label = label
        self.extinct = False
        self._properties = False
        # Set the properties of the self replicator.
        _constants = {'gamma_max': 9.65,
                     'phi_O': 0.55,
                     'tau': 1,
                     'kappa_max':1.15E-3,
                     'Km_u': 3E-5,
                     'Km_c': 3E-5,
                     'Km': [5E-6 for _ in range(self.num_metab)],
                     'Y': [2.95E19 for _ in range(self.num_metab)],
                     'nu_max': suballocation['nu_max'],
                     'death_rate': 0,
                     'frac_useful':[1 for _ in range(self.num_metab)]}
        self._overridden_pars = {}
        for k, v in constants.items():
            _constants[k] = v
            self._overridden_pars[k] = v

        # Set attributes from the constant dictionary.
        for k, v in _constants.items():
            setattr(self, k, v)

        # Set the suballocation details
        self.strategy = suballocation['strategy']
        if self.strategy == 'static':
            if np.sum(suballocation['alpha']) != 1:
                raise ValueError(f"Suballocation parameters must sum to 1!")
            self.alpha = suballocation['alpha'] 
        elif self.strategy == 'dynamic':
            self.K = suballocation['K']
            self.n = suballocation['n']
            if 'hierarchy' not in suballocation.keys():
                self.hierarchy = np.arange(self.num_metab).astype(int)
            else:
                self.hierarchy = np.array([int(k) for k in suballocation['hierarchy']]).astype(int)
        else:
            raise ValueError("Supplied metabolic strategy must be either `static` or `dynamic`.") 

    def __repr__(self):
        rep = f"""
=============== Self Replicator #{self.label} ================
extinct?                    : {self.extinct}
metabolic_strategy          : {self.strategy}
number of metabolic classes : {self.num_metab}
metabolic rates             : {self.nu_max}
hierarchy                   : {self.hierarchy}
"""
        if self.strategy == 'dynamic':
            rep += f"\tK's [M] : {self.K}"
            rep += f"\n\tn's : {self.n}"
        elif self.strategy == 'static':
            rep += f"\talpha's: {self.alpha}"
        if len(self._overridden_pars) > 0:
            rep += f"\nOverridden Default Parameters:" 
            for k, v in self._overridden_pars.items():
                rep += f"\n\t{k}: {v}"
        if self._properties:
            rep+=f"""\n
---------------------Allocation------------------------------
tRNA ratio (charged/uncharged) : {self.ratio}
gamma (translation rate)       : {self.gamma} [inv. hr.]
nu (metabolic rate)            : {self.nu} [inv. hr.]
kappa (transcription rate)     : {self.kappa} [inv. hr.]
phi_Rb (ribosomal allocation)  : {self.phi_Rb}   
phi_Mb (metabolic allocation)  : {self.phi_Mb} 
            """ 
        return rep 

    def compute_properties(self, tRNA_c, tRNA_u, nutrients):
        """
        Computes the self replicator properties, including the allocation parameters 
        and the corresponding rates. 

        # Parameters
        # ----------
        # tRNA_c : float [0, inf)
        #     The charged tRNA concentration in relative mass abundance units.
        # tRNA_u : float[0, inf) 
        #     The uncharged tRNA concentration in relative mass abundance units.
        # nutrients : numpy.ndarray float 
        #     The concentration of the nutrients in the environment for calculation
        #     of rates. 
        """

        # Determine if the uncharged tRNA concentration is 0. If so, set the 
        # ribosomal allocation to 1 and do not compute the ratio
        if tRNA_u <= 0: #Include less than zero for numerical underflow
            self.phi_Rb = (1 - self.phi_O)
            self.ratio = np.inf
        else:
            self.ratio = tRNA_c / tRNA_u
            self.phi_Rb = (1 - self.phi_O) * (self.ratio / (self.ratio + self.tau))

        # Adjust the flux parity regulation on transcription
        self.kappa = self.kappa_max * self.phi_Rb / (1 - self.phi_O)

        # Set the translation rate
        self.gamma = self.gamma_max * tRNA_c / (tRNA_c + self.Km_c)
        
        # Set the metabolic rates
        self.nu  = np.zeros_like(self.nu_max)
        metab_factor = tRNA_u / (tRNA_u + self.Km_u)
        for i in range(self.num_metab):
            env_factor = nutrients[i] / (nutrients[i] + self.Km[i])
            self.nu[i] = self.nu_max[i] * env_factor * metab_factor

        # Compute the metabolic allocation
        self.phi_Mb = np.zeros(self.num_metab)
        if self.strategy == 'static':
            # Set the fixed suballocations as prescribed by the alpha parameters.
            for i in range(self.num_metab):
                self.phi_Mb[i] = (1 - self.phi_O - self.phi_Rb) * self.alpha[i] 

        elif self.strategy == 'dynamic':
            # set the suballocation based on the nutrient concentrations
            occupied_phi_Mb = 0 
            idx_sort = self.hierarchy.argsort()
            for i, v in enumerate(idx_sort):  
                if self.hierarchy[v] == (self.num_metab - 1):
                    # Evaluate the remaining suballocation to the final nutrient in the hierarchy.
                    self.phi_Mb[i] = (1 - self.phi_O - self.phi_Rb - occupied_phi_Mb) 
                else:
                    # Set the suballocation given the supplied monod function properties
                    numer = (nutrients[v] / self.K[i])**self.n[v]
                    factor = numer / (numer + 1)
                    self.phi_Mb[v] =  (1 - self.phi_O - self.phi_Rb - occupied_phi_Mb) * factor
                    occupied_phi_Mb += self.phi_Mb[v]
        self._properties = True

    def compute_derivatives(self, 
                            masses, 
                            nutrients):
        """
        Computes the mass and concentration derivatives given the internal flux-parity 
        allocation state.

        Parameters
        ----------
        masses : numpy.ndarray float
            The masses of the protein sectors and the concentrations of the charged 
            and uncharged tRNAs. This should be in the order of the metabolic masses,
            ribosomal mass, "other" protein masses, uncharged tRNAs, and charged tRNAs

        nutrients: numpy.ndarray float 
            The concentrations of the external nutrients.

        Returns
        -------
        masses_dot : numpy.ndarray float 
            The time derivative of the masses given the self replicator dynamics.
        nutrients_dt : numpy.ndarray float
            The time derivative the nutrient concentrations consumed by a the 
            self replicator given the dynamics.    
        """

        # Unpack the variables 
        self.M_Mb = masses[:self.num_metab]
        self.M_Rb = masses[self.num_metab]
        self.M_O = masses[self.num_metab + 1]
        self.M = np.sum(masses[:self.num_metab + 2])
        self.tRNA_u = masses[self.num_metab + 2]
        self.tRNA_c = masses[self.num_metab + 3]
        self.nutrients = np.array(nutrients) * (np.array(nutrients) >= 0)

        # Evaluate the properties
        self.compute_properties(self.tRNA_c, self.tRNA_u, self.nutrients)

        # Evalutate the derivatives
        dM_Rb = self.phi_Rb * self.gamma * self.M_Rb - self.M_Rb * self.death_rate
        dM_O = self.phi_O * self.gamma * self.M_Rb  -  self.M_O * self.death_rate
        dM_Mb = [phi_Mb * self.gamma * self.M_Rb - self.M_Mb[i] * self.death_rate for i, phi_Mb in enumerate(self.phi_Mb)]
        dtRNA_u = self.kappa + self.gamma * self.M_Rb * (1 - self.tRNA_u) / self.M - np.sum(self.nu * self.M_Mb * self.frac_useful) / self.M - self.tRNA_u * self.death_rate
        dtRNA_c = np.sum(self.nu * self.M_Mb * self.frac_useful) / self.M - self.gamma * self.M_Rb * (1 + self.tRNA_c) / self.M - self.tRNA_c * self.death_rate
        dc_nt = -self.nu * self.M_Mb / self.Y
        masses_dt = []
        for d in dM_Mb:
            masses_dt.append(d)
        masses_dt.append(dM_Rb)
        masses_dt.append(dM_O)
        masses_dt.append(dtRNA_u)
        masses_dt.append(dtRNA_c)
        return [np.array(masses_dt), dc_nt]


class Ecosystem():
    """Base class for an ecosystem of self-replicators"""
    def __init__(self, 
                 species, 
                 nutrients):
        """
        Set up an ecosystem of J self replicators growing on K nutrients

        Parameters
        ----------
        species : list of `FluxParityAllocator`
            A list of the consitituent species of the ecosystem. Each species
            should be specified as a FluxParityAllocator object.

        nutrients: list of dict 
            A list of nutrients with dictionary specifying the details about the 
            nutrient loads with the following keys.
                init_concs : numpy.ndarray
                    The initial concentrations of the nutrients in the environment.
                feed_concs : numpy.ndarray
                    The concentration of the nutrients in the feed stock. If not 
                    provided, it will be assumed to be the same as the initial 
                    concentration.
                inflow_rates : numpy.ndarray 
                    The inflow rates of each nutrient in units of concentration 
                    per unit time.  If not supplied, it will be assumed to be 0.
                degradation_rates : numpy.ndarray
                    The degradation rates of each nutrient in units of concentraiton 
                    per unit time. If not supplied, it will be assumed to be 0.
 
        """
        self.species = species
        self.num_species =  len(species)
        self.num_nutrients = len(nutrients['init_concs'])

        # Set the defaults
        if 'feed_concs' not in nutrients.keys():
            nutrients['feed_concs'] = nutrients['init_concs']
        for k in ['inflow_rates', 'degradation_rates']:
            if k not in nutrients.keys():
                nutrients[k] = np.zeros(self.num_nutrients) 
        for k, v in nutrients.items():
            setattr(self, k, v)

    # TODO: Allow for approximate steady state initiation. 
    def initialize(self, 
                   freqs=None, 
                   od_init=0.04):
        """
        Initialize the ecosystem by setting the starting masses and tRNA 
        concentrations with a prescribed starting frequency. 

        Parameters
        ----------
        freqs : numpy.ndarray or None
            The desired starting frequencies for each species. If `None`, all 
            species are assumed to start with equal frequency.
        od_init : float
            The approximate initial optical density of the culture. This is 
            converted to mass of amino acids using a conversion factor of 
            1.15E17 AA per OD600. Default value is 0.04.
        """
        if freqs is None:
            freqs = np.ones(self.num_species) / self.num_species
        M0 = od_init * OD_CONV * freqs
        params = []
        for i, _species in enumerate(self.species):
            # Set the concentrations of the charged and uncharged tRNA and 
            # compute the initial properties
            _species.compute_properties(1E-5, 1E-5, self.init_concs)

            #TODO: Here is where we could set up an integration to find the 
            # steady-state solution and use that to initialize the physiology.

            # Set the initial parameters in the order of metabolic proteins,
            # ribosomal proteins, other proteins, uncharged tRNA, and charged tRNA
            for j, phi_Mb in enumerate(_species.phi_Mb):
                params.append(phi_Mb * M0[i])
            params.append(_species.phi_Rb * M0[i])
            params.append(_species.phi_O * M0[i])
            params.append(1E-5)
            params.append(1E-5)
        for c in self.init_concs:
            params.append(c) 
        self._p0 = params 

    def _dynamical_system(self, t, 
                          params, 
                          args):
        """
        Computes the derivatives of the masses of the community members.
        """ 

        # Re-access dimensional information
        num_params = self.num_species * (4 + self.num_nutrients) 
        nutrients = params[-self.num_nutrients:]  

        # Re-access ecosystem parameters
        inflow = self.inflow_rates
        deg = self.degradation_rates
        feed = self.feed_concs

        # Unpack the parameters for easy iteration
        if self.num_species > 1:
            masses = [params[i*num_params:i*num_params + 1] for i in range(self.num_species)]
        else:
            masses = [params[:-self.num_nutrients]]

        # Initialize storage vectors for the derivatives
        derivs = []
        nutrient_derivs = np.zeros_like(nutrients)
        for i, m in enumerate(masses):
            species = args['species'][i]
            _masses, _nutrients = species.compute_derivatives(m, nutrients)
            for _m in _masses:
                derivs.append(_m) 
            nutrient_derivs += _nutrients 
    
        for i, n in enumerate(nutrient_derivs):
            # Account for chemostatic change
            _n = inflow[i] * feed[i] + n - deg[i] * nutrients[i] 
            derivs.append(_n)
        return derivs 

    def integrate(self, 
                  time_range=[0, 10], 
                  **solver_kwargs):
        """
        Integrate the temporal dynamics of the ecosystem. 
        """
        args = {'species': self.species}
        print('Integrating dynamics...', end='')
        self.sol = scipy.integrate.solve_ivp(self._dynamical_system, time_range, self._p0,
                                            args=(args,))
        print('done!')
        return self.sol